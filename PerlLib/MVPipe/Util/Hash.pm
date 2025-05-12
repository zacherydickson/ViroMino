package MVPipe::Util::Hash;

use strict;
use warnings;
require Exporter;

our $VERSION = 1.00;
our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(FilterHash);
our %EXPORT_TAGS = (All => [qw(&FilterHash)]);

#Given an hash and a hash containing white and black lists, removes key-value pairs from
#the hash whose keys do not appear in the whitelist, or do appear in the blacklist
#Input  - A hashref
#       - A hash with optional keys of white and black, both of which have values of array ref
#Output - None, the hash is modified directly
sub FilterHash($%){
    my ($hashRef,%filter) = @_;
    $filter{black} = [] unless(exists $filter{black});
    if(exists $filter{white}){
        my %whiteSet = map {($_ => undef)} @{$filter{white}};
        push(@{$filter{black}},grep {!exists $whiteSet{$_}} keys %{$hashRef});
    }
    delete @{$hashRef}{@{$filter{black}}};
}

1;
