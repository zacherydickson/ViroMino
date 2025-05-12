package MVPipe::Util::Array;

use strict;
use warnings;
require Exporter;

our $VERSION = 1.00;
our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(BinarySearch Unique Zip Unzip Split);
our %EXPORT_TAGS = (All => [qw(&BinarySearch &Unique &Zip &Unzip &Split)]);

#Given a value and a sorted array of values which have defined behaviour for equality and
#inequality operators, finds all matches for the value if any
#   By default it compares numerically
#Input  - a scalar
#       - a reference to an array
#       - a code reference which takes two values and return -1 for less than, 0 for equal,
#       and +1 for greater than; the first value is the query, the second is an element
#       form the array
#Output - a potentially empty array of indexes in the input array which match the input
#value
sub BinarySearch($\@;&){
    my ($x, $rArr, $method) = @_;
    $method = sub { $_[0] <=> $_[1] } unless(defined $method);
    my ($l, $r) = (0, $#{$rArr});
    my @idx;
    while($l <= $r){
        my $m = int(($l + $r) / 2);
        my $cmp = $method->($x,$rArr->[$m]);
        if(!$cmp){
            my $i = $m-1;
            $i-- while($i >= 0 and !$method->($x,$rArr->[$i]));
            $i++;
            push(@idx,$i++) while($i < @{$rArr} and !$method->($x,$rArr->[$i]));
            last;
        } elsif($cmp < 0){
            $r = $m - 1;
        } else {
            $l = $m + 1;
        }
    }
    return @idx;
}


#Given an array of values, returns an array containing only unique values
sub Unique {
    my %hash = map {$_ => undef} @_;
    return keys %hash;
}

#Given an array, which may be two concatenated array zips the first half and the second
#half together
#Output - A single array which alternates between elements from the the first and second input half
sub Zip {
    my @arr = @_;
    return map {($arr[$_],$arr[$_ + @arr / 2])} 0 .. $#arr;
}

#Given an array reorders such that all even elements come before all odd elements
#   with order maintained
#Indexing starts at 0 :: 0,1,2,3,4 -> 0,2,4,1,3 
sub Unzip {
    my @arr = @_;
    my @even = @arr[ grep {$_ % 2 == 0} (0 .. $#arr)];
    my @odd = @arr[ grep {$_ % 2 == 1} (0 .. $#arr)];
    return (@even,@odd);
}

#Given two arrays of equal size returns a hash with keys of uniq 2nd array values, and values of
#corresponding first array values
sub Split {
    my @arr = @_;
    my %hash;
    for(my $i = 0; $i < @arr / 2; $i++){
        $hash{$arr[$i + @arr/2]} = [] unless(exists $hash{$arr[$i + @arr/2]});
        push(@{$hash{$arr[$i + @arr/2]}},$arr[$i]);
    }
}

1;

