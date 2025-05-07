#!/usr/bin/env perl
use warnings;
use strict;
use File::Basename;
use Class::Struct;
use POSIX qw(ceil);
use Zach::Util::File qw(OpenFileHandle);

struct (Interval => {start => '$', end => '$'});
struct (SuppTag => {id => '$', type => '$', value => '$'});
struct (SAMEntry => {	qname => '$', flag => '$', rname => '$', pos => '$',
			mqual => '$', cigar => '$', rnext => '$', pnext => '$',
			tlen => '$', seq => '$', qual => '$', 
			tag_by_id => '%'});
struct (Variant => {chrom => '$', pos => '$', alleleMap => '%'});

my $BDebug=0;

sub BuildVariant;
sub BuildSAMEntry;
sub BuildRef2QueryMap($);
sub OpenAlnMap($$);
sub LocateVariant($$);
sub IsValidInterval($){my $iv = shift; return ($iv->end >= $iv->start);}
sub VariantToStr($);
sub EntryToStr($);

if(@ARGV < 2){
    die "Usage: ".basename($0). " BAM Chrom Pos Ref Alt1 ... > out\n".
	"\tAlignmentMap\tIndexed BAM file\n".
	"\tAll other arguments specify the variant to count\n".
	"\tOutput is tab delim with headers:\n".
	"\t Chrom, Pos, AlleleID, Allele, Positions\n".
	"\t\tAlleleID is Ref, or AltN, defined by order\n".
	"\t\tAllele is the actual allele provided\n".
	"\t\tPositions is a comma separated list of read positions\n".
	"\t\t For variants which span outside of a read, the position is the\n".
	"\t\t leftmost position in the read\n".
	"\t\tPositions are normalized to be a proportion of the read length\n";
}

sub main {
    my $alnMapFile = shift;
    my $varObj = BuildVariant(@_);
    my $fh = OpenAlnMap($alnMapFile,$varObj);
    my %allelePos;
    while(my $line = <$fh>){
	chomp($line);
	my $entry = BuildSAMEntry($line);
	print STDERR VariantToStr($varObj),"\n" if($BDebug);
	print STDERR EntryToStr($entry),"\n" if($BDebug);
	my $iv = LocateVariant($entry,$varObj);
	print STDERR $iv->start."-".$iv->end."\n" if($BDebug);
	next unless(IsValidInterval($iv));
	my $allele = substr($entry->seq,$iv->start-1,$iv->end-$iv->start+1);
	print STDERR $allele,"\n" if($BDebug);
	$allelePos{$allele} = [] unless(exists $allelePos{$allele});
	my $normPos = ($iv->start-1)/length($entry->seq);
	push(@{$allelePos{$allele}}, $normPos);
    }
    print "Chrom\tPos\tAlleleID\tAllele\tPositions\n";
    foreach my $id (sort keys %{$varObj->alleleMap}){
        my $allele = $varObj->alleleMap->{$id};
        my $posStr ="";
        if(exists $allelePos{$allele}) {
	    $posStr = join(",",@{$allelePos{$allele}});
        }
	print join("\t",($varObj->chrom,$varObj->pos,$id,$allele,$posStr)),"\n";
    }
    close($fh);
} main (@ARGV);

#Given an array of inputs constructs a variant object
#Inputs - an array with elements in this order:
#	    chrom, pos, ref, alt1, ..., altN
#Ouput - a Variant Object
sub BuildVariant{
    my ($chrom,$pos,$ref,@alts) = @_;
    my $varObj = Variant->new(chrom => $chrom, pos => $pos);
    $varObj->alleleMap->{'Ref'} = $ref;
    my $digits = ceil(sprintf("%0.1f",log(@alts) / log(10)));
    for(my $i = 0; $i < @alts; $i++){
	$varObj->alleleMap->{sprintf("Alt%*d",$digits,$i+1)} = $alts[$i];
    }
    return $varObj;
}

#Given a SAM Entry constructs a mapping for reference positions to all corresponding query positions
#Inputs a SAMEntry objects
#Output a refposition keyed hash of array refs of query positions
sub BuildRef2QueryMap($){
    my $entry = shift;
    my $refPos = $entry->pos - 1;
    my $queryPos = 0;
    my $cigar = $entry->cigar;
    my %map;
    print STDERR $cigar,"\n" if($BDebug);
    while($cigar =~ s/^([0-9]+)(.)//){
	my $opLen = $1;
	my $op = $2;
	print STDERR "$opLen,$op,$cigar\n" if($BDebug);
	my $rUpdate = ($op =~ m/[MDN=X]/) ? 1 : 0;
	my $qUpdate = ($op =~ m/[MIS=X]/) ? 1 : 0;
	while($opLen-- > 0){
	    if($refPos and $queryPos){
		$map{$refPos} = [] unless(exists $map{$refPos});
		push(@{$map{$refPos}},$queryPos);
	    }
	    $refPos += $rUpdate;
	    $queryPos += $qUpdate;
	}
    }
    return %map;
}

#Given a string corresponding to a SAM entry, constructs a SAMEntry object
#Inputs - a string
#Output - a SAMEntry object
sub BuildSAMEntry {
    my $line = shift;
    my ($qname,$flag,$rname,$pos,$mqual,$cigar,$rnext,$pnext,$tlen,$seq,$qual,@tags) = split(/\t/,$line);
    my $entry = SAMEntry->new(	qname => $qname, flag => $flag, rname => $rname,
				pos => $pos, mqual => $mqual, cigar => $cigar,
				rnext => $rnext, pnext => $pnext, tlen => $tlen,
				seq => $seq, qual => $qual);
    foreach my $tagStr (@tags){
	my @tagFields = split(/:/,$tagStr);
	my $tag = SuppTag->new(	id => $tagFields[0],type => $tagFields[1],
				value => $tagFields[2]);
	$entry->tag_by_id->{$tag->id} = $tag;
    }
    return $entry;
}

#Opens a file handle to the entries in an alignment map file relevant to
# a variant
#Inputs - a path to an alignment map file
#	- a Variant object to search for
#Outputs - an open File handle;
sub OpenAlnMap($$) {
    my ($file,$obj) = @_;
    my $refLen = length($obj->alleleMap->{"Ref"});
    my $regEnd = $obj->pos + $refLen - 1;
    my $reg = $obj->chrom.":".$obj->pos."-$regEnd";
    return OpenFileHandle("samtools view -F0xF04 $file $reg |","samtools Cmd");
}

#Finds the interval of position in the entry which correspond to the positions of the reference allele
#Inputs - a SAMEntry object
#	- a Variant object
#Output - an interval object (1 indexed)
sub LocateVariant($$){
    my ($entry,$var) = @_;
    my $refLen = length($var->alleleMap->{"Ref"});
    my $nullIV = Interval->new(start => 1, end => 0);
    if(	$entry->rname ne $var->chrom or
	$entry->pos > $var->pos+$refLen - 1){
	return $nullIV;
    }
    my %map = BuildRef2QueryMap($entry);
    foreach my $refPos (sort {$a <=> $b} keys %map){
        print STDERR "$refPos => [@{$map{$refPos}}]\n" if($BDebug);
    }
    my @posList;
    for(my $offset = 0; $offset < $refLen; $offset++){
	my $pos = $var->pos + $offset;
	push(@posList,@{$map{$pos}}) if(exists $map{$pos});
    }
    return $nullIV if(!@posList);
    @posList = sort {$a <=> $b} @posList;
    print STDERR "@posList\n" if($BDebug);
    my $iv = Interval->new(start => $posList[0], end => $posList[$#posList]);
    return $iv;
}

#Collapse a variant object to a string 
#Inputs - a Variant Object
#Output - a string representing the objects's data
sub VariantToStr($) {
    my $var = shift;
    my $str = join("\t",($var->chrom,$var->pos));
    my @alleleStrs;
    foreach my $id (sort keys %{$var->alleleMap}){
	push(@alleleStrs,"$id:".$var->alleleMap->{$id});
    }
    return join("\t",($str,join("|",@alleleStrs)));
}


#Collapse a SAMEntry object to a string 
#Inputs - a SAMEntry Object
#Output - a string representing the objects's data
sub EntryToStr($){
    my $entry = shift;
    my $str = join("\t",(   $entry->qname,$entry->flag,$entry->rname,
			    $entry->pos,$entry->mqual,$entry->cigar,
			    $entry->rnext,$entry->pnext,$entry->tlen,
			    $entry->seq,$entry->qual));
    my @tagStrs;
    foreach my $id (sort keys %{$entry->tag_by_id}){
	my $tag = $entry->tag_by_id->{$id};
	push(@tagStrs,join(":",($tag->id,$tag->type,$tag->value)));
    }
    return join("\t",($str,@tagStrs));
}
