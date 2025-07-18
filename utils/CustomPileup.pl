#!/usr/bin/env perl
use warnings;
use strict;
use Bio::DB::HTS;
use Class::Struct;
use File::Basename;
use Getopt::Std;

#Local Libs
use FindBin;
use lib "$FindBin::Bin/../PerlLib";
use MVPipe::Util::File qw(OpenFileHandle LoadMapFile);
use MVPipe::Util::Array qw(Unique);

struct (PaddedVariantSet => {ref => '$', alleles => '@', len => '$'});
#alt_idx_by_sample is a sample keyed hashref of comma sep strings of alt indexes
#   alt numbering goes from 1 - n_allele
struct (CallSite => {   chrom => '$', pos => '$', ref => '$', s_call => '@',
                        ref_len => '$', n_allele => '$', indel => '$',
                        max_pad_len => '$', each_var_set => '@', 
                        alt_idx_str_by_smpl => '%'}); 
struct (SiteData => { DP => '$', AD => '@', ACO => '@', AF => '@', CCA => '$'});
#incomplete:
#   -1: missing left part of allele; 0: complete, 1: missing right part of allele
struct (AlignedAllele => {subject => '$', query => '$', incomplete => '$'});

sub LoadCommonCalls($@);
sub GetHTSObjects($$);
sub CallSite2Str($);
sub CallSite2VCFStr($);
sub LocateAlignedSite($$$);
sub GetAlignedAllele($$);
sub OutputVCFHeader($@);
sub ProcessSampleCallSite($$\%);
sub SiteInfo2Str($);

my %CigarOpsConsumingReference = map {($_ => 1)} qw(I S H P);

my %opts;
getopts('D',\%opts);
if(exists $opts{D}){
    exit 0;
}

if(@ARGV < 3){
    print STDERR "Usage: ".basename($0). " refFasta commonCalls BamMap\n".
                "\trefFasta - Path to a fasta file containing the reference sequence\n".
                "\tcommonCalls - tab delim, headers of CHROM,POS,REF,ALT,Callers,and Sample\n".
                "\t\tall calls are alleleic primatives\n".
                "\tBamMap\tTab delim mapping between sample IDs and paths to bams\n";
    exit 0;
}

sub main {
    my ($refFile, $commonCallsFile, $bamMapFile, $maxReadLen) = @_;
    my %htsMap = GetHTSObjects($bamMapFile,$refFile);
    my @sampleNames = sort keys %htsMap;
    my %sampleIdxMap;
    @sampleIdxMap{@sampleNames} = (0 .. $#sampleNames);
    my @callSites = LoadCommonCalls($commonCallsFile,@sampleNames);
    #foreach my $site (@callSites){
    #    print CallSite2Str($site),"\n";
    #}; return;
    OutputVCFHeader($refFile,@sampleNames);
    #Iterate over call Sites
    foreach my $callSite (sort {$a->pos <=> $b->pos} @callSites) {
        #Iterate over hts Objects in order of name
        my $formatStr = "DP:AD:AF:ACO:CCA";
        foreach my $smplName (@sampleNames){
            my $siteInfo = ProcessSampleCallSite($callSite,$smplName,%htsMap); 
            $formatStr .= "\t".SiteInfo2Str($siteInfo);
        }
        my @smplIdx = sort {$a <=> $b} @sampleIdxMap{@{$callSite->s_call}};
        my $infoStr = "NCALL=".scalar(@smplIdx);
        if(@smplIdx){
            $infoStr .= ";SCALL=".join(",",@smplIdx);
        }
        $infoStr .= ";INDEL" if($callSite->indel);
        print  join("\t",(  CallSite2VCFStr($callSite),('.') x 2,
                            $infoStr,$formatStr)),"\n";
    }
} main (@ARGV);

#Load a common calls file and return the uniq sites to process
#Note: Common Calls is expected to be in the form of alleleic primatives
#       With no multi-allelic entries
#      All calls at one site are combined into a single call site
#Inputs - A path to a common Calls File:
#           tab delim, with header of CHROM POS, etc
#Output - An array of uniq call sites in the calls file
sub LoadCommonCalls($@){
    my ($file,@sampleNames) = @_;
    my $fh = OpenFileHandle($file,"CommonCalls");
    my %uniqSites;
    #label keyed hash of ref allele keyed hash ref of alt alleles keyed array of samples
    my %refAltSetDict;
    my $header = <$fh>;
    while(my $line = <$fh>){
        chomp($line);
        my ($chrom,$pos,$ref,$alt,$vCaller,$sample) = split(/\t/,$line);
        $sample=undef if($sample eq "");
        my $label = "$chrom.$pos";
        my $bIndel = (length($ref) != length($alt)) ? 1 : 0;
        unless(exists $uniqSites{$label}){
            $uniqSites{$label} = CallSite->new( chrom => $chrom, pos => $pos,
                                                ref => undef, each_var_set => [],
                                                indel => $bIndel, s_call => []);
        }
        unless(exists $refAltSetDict{$label}){
            $refAltSetDict{$label} = {};
        }
        unless(exists $refAltSetDict{$label}->{$ref}){
            $refAltSetDict{$label}->{$ref} = {}
        }
        $refAltSetDict{$label}->{$ref}->{$alt} = [] unless(exists $refAltSetDict{$label}->{$ref}->{$alt});
        if(defined $sample){
            push(@{$refAltSetDict{$label}->{$ref}->{$alt}},$sample);
        }
    }
    #Ensure that all reference alleles are of the same form across samples
    #   padding is added at this point which can cause multiple ref forms
    #   at one site
    while (my ($label, $rRefAltSet) = each %refAltSetDict){
        my $callSite = $uniqSites{$label};
        my ($commonRef) = sort { length($b) <=> length($a) } (keys %{$rRefAltSet});
        my $commonRefLen = length($commonRef);
        $callSite->ref($commonRef);
        $callSite->ref_len($commonRefLen);
        push(   @{$callSite->each_var_set},
                PaddedVariantSet->new(  ref => $commonRef,
                                        len => $commonRefLen,
                                        alleles => [$commonRef]));
        #Used to Count the number of unique samples in which this site was called
        # sample keyed hash of padded allele setRef
        my %sampleSet;
        #Used to map unique alt alleles to samples they in which they were called
        #paddedRef Keyed hash of paddedAlt setRef
        my %padSetDict;
        while (my ($ref,$rAltSet) = each %{$rRefAltSet}){
            #Determine the suffix needed to make all calls at the site of a
            #   consistent form
            my $refLen = length($ref);
            my $suffix = "";
            if($refLen < $commonRefLen){
                $suffix = substr($commonRef,$refLen);
            }
            #Pad each ref-alt pair prior to adding the suffix
            #Then record the unique pad ref forms observed and their associated
            # padded alts   
            while( my ($alt,$rSampleList) = each %{$rAltSet}){
                my $lenDiff = length($alt) - $refLen;
                my $padRef = $ref;
                $padRef .= '-' x abs($lenDiff)  if($lenDiff > 0);
                $padRef .= $suffix;
                my $padAlt = $alt;
                $padAlt .= '-' x abs($lenDiff)  if($lenDiff < 0);
                $padAlt .= $suffix;
                $padSetDict{$padRef} = {} unless(exists $padSetDict{$padRef});
                $padSetDict{$padRef}->{$padAlt} = 1;
                #For every sample in which this alt appeared,
                # add this padded alt to the list of padded Alts for this sample
                foreach my $sample (@{$rSampleList}){
                    $sampleSet{$sample} = {} unless(exists $sampleSet{$sample});
                    $sampleSet{$sample}->{$padAlt} = 1;
                }
            }
        }
        $callSite->s_call([keys %sampleSet]);
        my $nAllele = 1;
        my %alleleIdxBySample = map {($_ => [])} @sampleNames;
        my $altIdx = 1;
        while(my ($padRef,$rPadAltSet) = each %padSetDict){
            #The ref allele will always be included
            delete $rPadAltSet->{$padRef} if(exists $rPadAltSet->{$padRef});
            my $varSet = PaddedVariantSet->new( ref => $padRef,
                            len => length($padRef),
                            alleles => [sort keys %{$rPadAltSet}]);
            push(@{$callSite->each_var_set}, $varSet);
            foreach my $padAlt (@{$varSet->alleles}){
                while (my ($sample, $rPadSet) = each %sampleSet){
                    if(exists $rPadSet->{$padAlt}) {
                        push(@{$alleleIdxBySample{$sample}},$altIdx);
                    }
                }
            } continue {$altIdx++}
            $nAllele += scalar(keys %{$rPadAltSet});
        }
        my @varSetList = @{$callSite->each_var_set};
        my ($maxObj) = sort {$b->len <=> $a->len} @varSetList;
        $callSite->n_allele($nAllele);
        $callSite->max_pad_len($maxObj->len);
        #Construct altIdx Strings
        my %altIdxStrDict = map {($_ => '.')} @sampleNames;
        foreach my $smpl (@sampleNames){
            my $rAltIdxList = $alleleIdxBySample{$smpl};
            if(scalar @{$rAltIdxList}) {
                $altIdxStrDict{$smpl} = join(',',
                    sort {$a <=> $b} @{$rAltIdxList});
            }
        }
        $callSite->alt_idx_str_by_smpl(\%altIdxStrDict);
    }
    return values %uniqSites;
}


sub GetHTSObjects($$) {
    my $bamMapFile=shift;
    my $refFile=shift;
    my %BamMap = LoadMapFile($bamMapFile,"BamMap","\t");
    my %htsMap;
    while (my ($name,$path) = each %BamMap){
        $htsMap{$name} = Bio::DB::HTS->new(-fasta => $refFile, -bam => $path);
    }
    return %htsMap;
}


sub VarSet2Str($){
    my $obj = shift;
    my $alleleStr = join(",",@{$obj->alleles});
    return $obj->len.':'.$obj->ref.' -> '.$alleleStr;
}

sub CallSite2Str($){
    my $obj = shift;
    my $varSetStr = join(";",map {VarSet2Str($_)} @{$obj->each_var_set});
    return join("\t",($obj->chrom,$obj->pos,$obj->ref,$obj->n_allele,
                      $obj->max_pad_len,$varSetStr));
}

sub CallSite2VCFStr($){
    my $obj = shift;
    my %alleleSet;
    my @varSetList = @{$obj->each_var_set};
    #Strip of the first var set (reference)
    shift @varSetList;
    my @altAlleles;
    foreach my $varSet (@varSetList){
        my @alleles = @{$varSet->alleles};
        for(my $i = 0; $i < @alleles; $i++){
            $alleles[$i] =~ tr/-//d;
        }
        push(@altAlleles,@alleles);
    }
    my $altStr = join(",",@altAlleles);
    return join("\t",(  $obj->chrom,$obj->pos,".",$obj->ref,$altStr));
}

#Given a cigar array and an offset within an alignment, determine the position within an alignment where a variant position is
#Inputs - An offset in the unaligned sequence
#       - An array ref of two elelement array refs,
#           the elements are CIGAR op and oplen
#Output - an offset within the aligned sequence
sub LocateAlignedSite($$$) {
    my ($off,$rCigarArray,$bLeft) = @_;
    my $alignOff = $off;
    my $i = 0;
    while($i < $#{$rCigarArray} and $off >= 0) {
        my $opLen = $rCigarArray->[$i]->[1];
        if(exists $CigarOpsConsumingReference{$rCigarArray->[$i]->[0]}){
            $alignOff += $opLen;
        } else {
            $off -= $opLen;
        }
    } continue {$i++}
    return $alignOff;
}

#Given a call Site and an alignment retreives the portion of the read which 
#covers the allele, without gaps
#First finds the in alignment position of the site
#Then pulls a subsequence of the read from that position, up to the length of the longest allele, without any gaps
#Inputs - a CallSite Object
#       - a Bio::DB::HTS::AlignmentWrapper Object
#Output - an AlignedAlleleObject
sub GetAlignedAllele($$){
    my ($callSite, $align) = @_;
    my $alnStart = $align->start;
    my $siteOff = LocateAlignedSite($callSite->pos - $align->start,
                                    $align->cigar_array,0);
    my $lengthOff = 0;
    my $incomplete = 0;
    if($siteOff < 0) {
        $lengthOff = abs($siteOff);
        $siteOff = 0;
        $incomplete = -1;
    }
    my $callEnd = $callSite->pos + $callSite->ref_len - 1;
    #Find the first position in the alignment after the end of the call site
    # then return the position just before that, which allows for gaps
    my $endOff = LocateAlignedSite($callEnd - $align->start + 1,
                                    $align->cigar_array,1);
    my ($subject, $matches, $query) = $align->padded_alignment;
    my $alnLen = length($subject);
    #take the read suffix from the start of the allele
    my $alleleLen = $endOff-$siteOff-$lengthOff;
    $query = substr($query,$siteOff,$alleleLen);
    $subject = substr($subject,$siteOff,$alleleLen);

    if(!$incomplete and (   length($query) < $alleleLen or
                            $endOff >= $alnLen))
    {
        $incomplete = 1;
    }
    return AlignedAllele->new(  subject => $subject, query => $query,
                                incomplete => $incomplete);
}


sub OutputVCFHeader($@) {
    my ($ref,@sampleNames) = @_;
    print   "##fileformat=VCFv4.2\n".
            "##reference=file://@{[basename($ref)]}\n";
    my %contigLenMap = LoadMapFile("$ref.fai","faidx","\t");
    foreach my $contig (sort keys %contigLenMap){
        print "##contig=<ID=$contig,length=$contigLenMap{$contig}>\n";
    }
    print   "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates variant is an INDEL\">\n".
            "##INFO=<ID=NCALL,Number=1,Type=Integer,Description=\"Number of samples in which this site was called\">\n".
            "##INFO=<ID=SCALL,Number=.,Type=Integer,Description=\"The (zero-indexed) sample numbers in which this site was called\">\n".
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of reads overlapping a position\">\n".
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n".
            "##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Allelic Frequency\">\n".
            "##FORMAT=<ID=ACO,Number=R,Type=Integer,Description=\"Number of reads consistent with an allele\">\n".
            "##FORMAT=<ID=CCA,Number=.,Type=Integer,Description=\"The indexes of the alt allele(s), if any, called in this sample\">\n".
            "#".join("\t",(qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT),@sampleNames))."\n"
            ;
}


sub ProcessSampleCallSite($$\%) {
    my $callSite = shift;
    my $smplName = shift;
    my $rHTSMap = shift;
    my $hts = $rHTSMap->{$smplName};
    #Contruct default site info with zero coverage
    my $altIdxStr = $callSite->alt_idx_str_by_smpl($smplName);
    my $nAllele = $callSite->n_allele;
    my $siteInfo = SiteData->new(   DP => 0, AD => [(0) x $nAllele],
                                    ACO => [(0) x $nAllele], CCA => $altIdxStr,
                                    AF => [('.') x $nAllele]);
    #If a call site is outside of the reference genome -
    # a user provided a set of variants only some of which apply to this genome -
    # then segment will be undefined
    my $segment = $hts->segment($callSite->chrom,$callSite->pos,
                                $callSite->pos + $callSite->ref_len - 1);
    return $siteInfo unless(defined $segment);
    #Find all reads overlapping this position
    my $iterator = $segment->features(-iterator => 1, -type => 'match');
    #Iterate over the alignments, extracting the alt sequences
    while(my $align = $iterator->next_seq){
        $siteInfo->DP($siteInfo->DP + 1);
        my $siteOff = LocateAlignedSite($callSite->pos - $align->start,
                                        $align->cigar_array,0);
        my $alnAlleleObj = GetAlignedAllele($callSite,$align);
        #Iterate over variant sets with the same reference form 
        my $adIdx = undef;
        my $alleleIdx = 0;
        foreach my $varSetObj (@{$callSite->each_var_set}){
            my $refAln = $alnAlleleObj->subject;
            my $queryAln = $alnAlleleObj->query;
            #Iterate over alleles
            foreach my $allele (@{$varSetObj->alleles}){
                my $alleleQ = $allele;
                my $alleleRef = $varSetObj->ref;
                my $diff = length($alleleQ) - length($refAln);
                if($diff > 0){
                    next if(! $alnAlleleObj->incomplete);
                    if($alnAlleleObj->incomplete < 0){
                        $alleleRef = substr($alleleRef,$diff);
                        $alleleQ = substr($allele,$diff);
                    } else {
                        $alleleRef = substr($alleleRef,0,-$diff);
                        $alleleQ = substr($allele,0,-$diff);
                    }
                }
                #Check if the read is consistent with the allele
                if($alleleRef eq $refAln and $alleleQ eq $queryAln){
                    #Set the Ad to this allele, if multiple alleles
                    # are consistent, mark as no AD
                    $adIdx = (defined $adIdx) ? -1 : $alleleIdx;
                    #Increment consistency counter
                    $siteInfo->ACO->[$alleleIdx]++;
                }
            } continue {$alleleIdx++;}
        }
        #If the read is consistent with only one allele increment the depth
        if(defined $adIdx and $adIdx >= 0){
            $siteInfo->AD->[$adIdx]++;
        }
    }
    #Calculate Allele Frequencies
    for(my $i = 0; $i < $nAllele;$i++){
        my $ad = $siteInfo->AD->[$i];
        if($siteInfo->DP > 0){
            $siteInfo->AF->[$i] = sprintf("%0.03f",$ad/$siteInfo->DP);
        }
    }
    return $siteInfo;
}


sub SiteInfo2Str($){
    my $self = shift;
    return join(":",(   $self->DP,
                        join(",",@{$self->AD}),
                        join(",",@{$self->AF}),
                        join(",",@{$self->ACO}),
                        $self->CCA));
}
