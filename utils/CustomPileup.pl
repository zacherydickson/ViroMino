#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::HTS;
use Class::Struct;
use File::Basename;
use Zach::Util::File qw(OpenFileHandle LoadMapFile);
use Zach::Util::Array qw(Unique);

struct (CallSite => {chrom => '$', pos => '$', ref => '$', alleles => '@'}); 

sub LoadCommonCalls($);
sub GetHTSObjects($$);
#sub pileupCallback($$$$);

if(@ARGV < 3){
    die "Usage: ".basename($0). " refFasta commonCalls BamMap \n";
}

sub main {
    my ($refFile, $commonCallsFile, $bamMapFile) = @_;
    my @callSites = LoadCommonCalls($commonCallsFile);
    my %htsMap = GetHTSObjects($bamMapFile,$refFile);
    my @sampleNames = sort keys %htsMap;
    #Iterate over call Sites
    foreach my $callSite (@callSites) {
        #Iterate over hts Objects in order of name
        foreach my $smplName (@sampleNames){
            my $hts = $htsMap{$smplName};
            my $segment = $hts->segment($callSite->chrom,$callSite->pos,$callSite->pos);
            my $ref = $segment->dna;
            #Find all reads overlapping this position
            my $iterator = $segment->features(-iterator => 1, -type => 'match');
            #Iterate over the alignments, extracting the alt sequences
            while(my $align = $iterator->next_seq){
                my $alnStart = $align->start;
                my $siteOff = $callSite->pos - $align->start;
                my ($ref, $matches, $query) = $align->padded_alignment;
                #TODO: Figure out how best to parse the alignments to get the expanded data we want
                #Idea: Iterate over all uniq alleles called at the site
                #   Count the number of reads which are CONSISTENT with that allele
                #   For each allele, Note the number of reads which are consistent only with it
                #   The former could be the support for the read
                #   The latter would be the allelic depth
                #   At the end any allele which has no allelic depth can be dropped (except the reference)
                #TODO: Probably should try to combine things together at this point into haplotypes where possible
                #   Perhaps a first pass to identify reads which span multiple sites,
                #   those can be combined into a single call site entry and treated the same as above
                #TODO: Work would then need to be done to mark the samples which had actual alleles called
                if($smplName eq "JC3" && $callSite->pos == 2147) {
                    print join(" ",$align->start,$align->end,$align->strand),"\n";
                    print $align->query->name,"@ $siteOff\n";
                    print join("\n",($align->padded_alignment)),"\n";
                }
            }
        }
    }
} main (@ARGV);

#Load a common calls file and return the uniq sites to process
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
        my ($chrom,$pos,$ref,$altstr,@other) = split(/\t/,$line);
        my $label = "$chrom.$pos";
        unless(exists $uniqSites{$label}){
            $uniqSites{$label} = CallSite->new(chrom => $chrom, pos => $pos, ref => undef, alleles => []);
        }
        unless(exists $RefAltSetsDict{$label}){
            $RefAllSetsDict{$label} = {};
        }
        unless(exists $RefAllSetsDict{$label}->{$ref}){
            $RefAllSetsDict{$label}->{$ref} = {}
        }
        foreach my $alt (split(/,/,$altstr)){
            $RefAllSetsDict{$label}->{$ref}->{$alt} = 1;
        }
    }
    while (my ($label, $rRefAltSet)){
    #TODO:
    #   Ensure that all reference alleles are of the same form across samples
        $uniqSites{$label}->ref($commonRef);
        $uniqSites{$label}->alleles(\@alleles);
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

##
#sub pileupCallback($$$$){
#    my ($seqid,$pos,$pileup,$hts) = @_;
#
#}
