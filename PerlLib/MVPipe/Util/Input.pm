package MVPipe::Util::Input;

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
require Exporter;

our $VERSION = 1.00;
our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(GetProcessorCount ProcessNumericOption ValidateThreadCount);
our %EXPORT_TAGS = (All => [qw(&GetProcessorCount &ProcessNumericOption &ValidateThreadCount)]);

#Returns the maximum number of processors on a linux based system
#Output: -1 if /proc/cpuinfo not found, the number of processors otherwise
sub GetProcessorCount(){
    my $cpu_count = -1;
    if(open my $handle, "/proc/cpuinfo"){
        $cpu_count = scalar(map /^processor/, <$handle>);
        close($handle);
    }
    return $cpu_count;
}

sub ValidateThreadCount($){
    my $max_proc = shift;
    my $cpu_count = GetProcessorCount();
    my $message = undef;
    if($cpu_count == -1){
        $message = "Could not determine Sys_max threads: $!\n\t Can't fully validate max_threads\n";
        if($max_proc == 0){
            $message = "Sys_max threads unknown: proceeding with 1 thread";
            $max_proc = 1;
        }
    } elsif($max_proc == 0 or $max_proc > $cpu_count){
        if($max_proc > $cpu_count){
            $message = "Max threads is greater than sys_max: proceeding with $cpu_count threads";
        }
        $max_proc = $cpu_count;
    }
    if(defined $message){
        my ($sec,$min,$hour) = localtime;
        printf STDERR "%02d:%02d:%02d - [%s] %s\n", ($hour,$min,$sec,"WARNING",$message);
    }
    return $max_proc;
}

sub ProcessNumericOption($$$$$$){
    my($val,$default,$min,$max,$bInt,$varName) = @_;
    return $default unless(defined $val);
    if(looks_like_number($val)){
        $val = int($val) if($bInt);
        if(!defined $min or $val >= $min){
            if(!defined $max or $val <= $max){
                return $val;
            }
        }
    }
    my $message = sprintf("%s must be a%s between %s and %s",$varName,
       ($bInt ? "n integer" : " value"),(defined $min ? $min : "-∞"),(defined $max ? $max : "∞"));
    my ($sec,$minute,$hour) = localtime;
    printf STDERR "%02d:%02d:%02d - [%s] %s\n", ($hour,$minute,$sec,"ERROR",$message);
    exit 1;
}

1;
