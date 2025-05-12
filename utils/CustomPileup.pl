#!/usr/bin/env perl
use warnings;
use strict;
use Bio::DB::HTS;
use Class::Struct;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../PerlLib";
use MVPipe::Util::File qw(OpenFileHandle LoadMapFile);
use MVPipe::Util::Array qw(Unique);

struct (PaddedVariantSet => {ref => '$', alleles => '@', len => '$'});
struct (CallSite => {   chrom => '$', pos => '$', ref => '$',
                        ref_len => '$', n_allele => '$',
                        max_pad_len => '$', each_var_set => '@'}); 
struct (AlignedAllele => {subject => '$', query => '$'});

sub LoadCommonCalls($);
sub GetHTSObjects($$);
sub CallSite2Str($);
sub LocateAlignedSite($$$);
sub GetAlignedAllele($$);
#sub pileupCallback($$$$);

my %CigarOpsConsumingReference = map {($_ => 1)} qw(I S H P);

if(@ARGV < 3){
    die "Usage: ".basename($0). " refFasta commonCalls BamMap\n";
}

sub main {
    my ($refFile, $commonCallsFile, $bamMapFile, $maxReadLen) = @_;
    my %htsMap = GetHTSObjects($bamMapFile,$refFile);
    my @callSites = LoadCommonCalls($commonCallsFile);
    #foreach my $site (@callSites){
    #    print CallSite2Str($site),"\n";
    #}; return;
    my @sampleNames = sort keys %htsMap;
    #Iterate over call Sites
    foreach my $callSite (sort {$a->pos <=> $b->pos} @callSites) {
        #Iterate over hts Objects in order of name
        foreach my $smplName (@sampleNames){
            my $hts = $htsMap{$smplName};
            my $maxPaddedLen = $callSite->max_pad_len;
            my $segment = $hts->segment($callSite->chrom,$callSite->pos,
                                        $callSite->pos + $callSite->ref_len - 1);
            #Find all reads overlapping this position
            my $iterator = $segment->features(-iterator => 1, -type => 'match');
            #Iterate over the alignments, extracting the alt sequences
            my @allelicDepth = (0) x $callSite->n_allele;
            my @allelicConsistency = (0) x $callSite->n_allele;
            while(my $align = $iterator->next_seq){
                #TODO: This still doesn't work because C->C is always consistent with C--- -> CXXX, but C--- wont ever be there alone
                my $siteOff = LocateAlignedSite($callSite->pos - $align->start,
                                                $align->cigar_array,0);
                my $alnAlleleObj = GetAlignedAllele($callSite,$align);
                #print "aa| ", join("\t",($alnAlleleObj->subject,$alnAlleleObj->query)),"\n";
                #TODO: extract the cigar adjusted Ref and Alt strings
                my $lengthOff = ($siteOff < 0) ? abs($siteOff) : 0;
                $siteOff = 0 if($siteOff < 0);
                my ($refseq, $matches, $queryseq) = $align->padded_alignment;
                my $adIdx = undef;
                my $alleleIdx = 0;
                foreach my $varSetObj (@{$callSite->each_var_set}){
                    #my $len = $varSetObj->len - $lengthOff;
                    ##If the read is entirely after the alleles it can't be used
                    #if($len < 1){
                    #    $alleleIdx += scalar @{$varSetObj->alleles};
                    #    next;
                    #}
                    #my $refAln = substr($refseq,$siteOff,$len);
                    #print "$refAln $siteOff,$len, (@{[VarSet2Str($varSetObj)]})\n";
                    #my $queryAln = substr($queryseq,$siteOff,$len);
                    my $refAln = $alnAlleleObj->subject;
                    my $queryAln = $alnAlleleObj->query;
                    foreach my $allele (@{$varSetObj->alleles}){
                        my $alleleQ = $allele;
                        my $alleleRef = $varSetObj->ref;
                        my $alleleLen = length($alleleQ);
                        my $alnLen = length($refAln);
                        if($alleleLen > $alnLen){
                            my $diff = $alleleLen - $alnLen;
                            if($lengthOff) { #allele starts before read
                                $alleleRef = substr($alleleRef,$diff);
                                $alleleQ = substr($allele,$diff);
                            } else { #allele ends after read, leave off diff
                                $alleleRef = substr($alleleRef,0,-$diff);
                                $alleleQ = substr($allele,0,-$diff);
                            }
                            
                        }
                        print "$alleleRef vs $refAln and $alleleQ vs $queryAln\n";
                        #Check if the read is consistent with the allele
                        if($alleleRef eq $refAln and $alleleQ eq $queryAln){
                            #Set the Ad to this allele, if multiple alleles
                            # are consistent, mark as no AD
                            $adIdx = (defined $adIdx) ? -1 : $alleleIdx;
                            #Increment consistency counter
                            $allelicConsistency[$alleleIdx]++;
                        }
                    } continue {$alleleIdx++;}
                }
                if(defined $adIdx and $adIdx >= 0){
                    $allelicDepth[$adIdx]++;
                }
                print "$smplName\n";
                print "CS) ", CallSite2Str($callSite),"\n";
                print $align->query->name," @ $siteOff, $lengthOff\n";
                print join("\n",($align->padded_alignment)),"\n";
                print "@allelicDepth\n";
                print "@allelicConsistency\n";

                #TODO: Figure out how best to parse the alignments to get the
                #   expanded data we want
                #Idea: Iterate over all uniq alleles called at the site
                #   Count the number of reads which are CONSISTENT with that allele
                #   For each allele, Note the number of reads which are 
                #       consistent only with it
                #   The former could be the support for the read
                #   The latter would be the allelic depth
                #   At the end any allele which has no allelic depth can be
                #       dropped (except the reference)
                #TODO: Probably should try to combine things together at this
                #   point into haplotypes where possible
                #   Perhaps a first pass to identify reads which span multiple 
                #   sites, those can be combined into a single call site entry and
                #       treated the same as above
                #TODO: Work would then need to be done to mark the samples
                #   which had actual alleles called

                #my $alnStart = $align->start;

                #if($smplName eq "JC3" && $callSite->pos == 2930) {
                #    print join(" ",$align->start,$align->end,$align->strand),"\n";
                #}
            }
        }
    }
} main (@ARGV);

#Load a common calls file and return the uniq sites to process
#Note: Common Calls is expected to be in the form of alleleic primatives
#       With no multi-allelic entries
#      All calls at one site are combined into a single call site
#Inputs - A path to a common Calls File:
#           tab delim, with header of CHROM POS, etc
#Output - An array of uniq call sites in the calls file
sub LoadCommonCalls($){
    my $file = shift;
    my $fh = OpenFileHandle($file,"CommonCalls");
    my %uniqSites;
    #label keyed hash of ref allele keyed hash ref of alt alleles set ref
    my %refAltSetDict;
    my $header = <$fh>;
    while(my $line = <$fh>){
        chomp($line);
        my ($chrom,$pos,$ref,$alt,@other) = split(/\t/,$line);
        my $label = "$chrom.$pos";
        unless(exists $uniqSites{$label}){
            $uniqSites{$label} = CallSite->new( chrom => $chrom, pos => $pos,
                                                ref => undef, each_var_set => []);
        }
        unless(exists $refAltSetDict{$label}){
            $refAltSetDict{$label} = {};
        }
        unless(exists $refAltSetDict{$label}->{$ref}){
            $refAltSetDict{$label}->{$ref} = {}
        }
        $refAltSetDict{$label}->{$ref}->{$alt} = 1
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
        my %alleleSet;
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
            foreach my $alt (keys %{$rAltSet}){
                my $lenDiff = length($alt) - $refLen;
                my $padRef = $ref;
                $padRef .= '-' x abs($lenDiff)  if($lenDiff > 0);
                $padRef .= $suffix;
                my $padAlt = $alt;
                $padAlt .= '-' x abs($lenDiff)  if($lenDiff < 0);
                $padAlt .= $suffix;
                $padSetDict{$padRef} = {} unless(exists $padSetDict{$padRef});
                $padSetDict{$padRef}->{$padAlt} = 1;
            }
        }
        my $nAllele = 1;
        while(my ($padRef,$rPadAltSet) = each %padSetDict){
            #The ref allele will always be included
            delete $rPadAltSet->{$padRef} if(exists $rPadAltSet->{$padRef});
            push(   @{$callSite->each_var_set},
                    PaddedVariantSet->new(ref => $padRef, len => length($padRef),
                                          alleles => [sort keys %{$rPadAltSet}]));
            $nAllele += scalar(keys %{$rPadAltSet});
        }
        my ($maxObj) = sort {$b->len <=> $a->len} @{$callSite->each_var_set};
        $callSite->n_allele($nAllele);
        $callSite->max_pad_len($maxObj->len);
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
    if($siteOff < 0) {
        $lengthOff = abs($siteOff);
        $siteOff = 0;
    }
    my $callEnd = $callSite->pos + $callSite->ref_len - 1;
    #Find the first position in the alignment after the end of the call site
    # then return the position just before that, which allows for gaps
    my $endOff = LocateAlignedSite($callEnd - $align->start + 1,
                                    $align->cigar_array,1);
    print "I) $alnStart; @{[$callSite->pos]} -> $siteOff  $callEnd -> $endOff\n";
    my ($subject, $matches, $query) = $align->padded_alignment;
    #take the read suffix from the start of the allele
    $query = substr($query,$siteOff,$endOff-$siteOff-$lengthOff);
    $subject = substr($subject,$siteOff,$endOff-$siteOff-$lengthOff);
    return AlignedAllele->new(subject => $subject, query => $query);
}

##
#sub pileupCallback($$$$){
#    my ($seqid,$pos,$pileup,$hts) = @_;
#
#}
