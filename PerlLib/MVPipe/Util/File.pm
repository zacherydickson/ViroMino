package MVPipe::Util::File;

use strict;
use warnings;
require Exporter;

our $VERSION = 1.00;
our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(OpenFileHandle LoadConfig LoadListFile LoadMapFile LoadMultiMapFile LoadBiMapFile TestFile);
our %EXPORT_TAGS = (All => [qw(&OpenFileHandle &LoadConfig &LoadListFile &LoadMapFile &LoadMultiMapFile &LoadBiMapFile &TestFile)]);

#Given a file path, opens it for reading and fails with the provided exit function if unsuccessful
# The default exit function is die
#Input	File - Path to the file to open
#	Type - A descriptor of the type of file being opened
# (opt)	level - the logging level for failure to open the file
#Output	A file handle if successful, 0 otherwise.
#Note: Calling process must close the file handle
sub OpenFileHandle($$;$){
    my ($file,$type,$level) = @_;
    $level = "ERROR" unless(defined $level);
    my $cmd = ($file =~ m/\.gz$/) ? "zcat $file |" : $file;
    if(open(my $fh, $cmd) ){
        return $fh;
    } else {
        my $message = "Could not open $type file ($file): $!";
        my ($sec,$min,$hour) = localtime;
        printf STDERR "%02d:%02d:%02d - [%s] %s\n", ($hour,$min,$sec,$level,$message);
        exit 1 if($level eq "ERROR");
        return 0;
    }
}

sub LoadConfig($;$){
    my $file = shift;
    my $level = shift;
    $level = "ERROR" unless(defined $level);
    my $fh = OpenFileHandle($file,"Config",$level);
    my %Config;
    while(my $line = <$fh>){
        next if($line =~ /^#/);
        chomp($line);
        next if($line eq "");
        my($key,$value, @other) = split(/\t/,$line);
        next unless(defined $key and defined $value and scalar(@other) == 0);
        if(exists $Config{$key}){
            my $message = "Duplicate config entry for $key on line $. ignored\n";
            my ($sec,$minute,$hour) = localtime;
            printf STDERR "%02d:%02d:%02d - [%s] %s\n", ($hour,$minute,$sec,"WARNING",$message);
        }
        $Config{$key} = $value;
    }
    close($fh);
    return %Config;
}

sub LoadListFile($$;$){
    my $file = shift;
    my $type = shift;
    my $level = shift || "ERROR";
    my $fh = OpenFileHandle($file,$type,$level);
    my @lines = <$fh>;
    chomp(@lines);
    @lines = grep {$_ ne ""} @lines;
    close($fh);
    return @lines;
}

sub LoadMapFile($$$;$){
    my $file = shift;
    my $type = shift;
    my $delim = shift;
    my $level = shift || "ERROR";
    my $fh = OpenFileHandle($file,$type,$level);
    return () unless($fh);
    my %map;
    while(my $line = <$fh>){
        chomp($line);
        next if($line eq "");
        my ($key,$value) = split(/$delim/,$line);
        $map{$key} = $value;        
    }
    close($fh);
    return %map;
}

sub LoadMultiMapFile($$$;$){
    my $file = shift;
    my $type = shift;
    my $delim = shift;
    my $level = shift || "ERROR";
    my $fh = OpenFileHandle($file,$type,$level);
    my %map;
    while(my $line = <$fh>){
        chomp($line);
        next if($line eq "");
        my ($key,$value) = split(/$delim/,$line);
        $map{$key} = [] unless(exists $map{$key});
        push(@{$map{$key}},$value);
    }
    close($fh);
    return %map;
}

sub LoadBiMapFile($$$;$){
    my $file = shift;
    my $type = shift;
    my $delim = shift;
    my $level = shift || "ERROR";
    my $fh = OpenFileHandle($file,$type,$level);
    my %map;
    while(my $line = <$fh>){
        chomp($line);
        next if($line eq "");
        my @vals = split(/$delim/,$line);
        for(my $dir = 0; $dir < 2; $dir++){
            $map{$vals[$dir]} = [] unless(exists $map{$vals[$dir]});
            push(@{$map{$vals[$dir]}},$vals[!$dir]);
        }
    }
    close($fh);
    return %map;
}

#Checks if a file meets requirements 
#Inputs - a string of characters corresponding to tests
#	    see perldoc -X
#Output - a string of characters corresponding to the tests the file failed
sub TestFile($;$){
    my $file = shift;
    my $testStr = shift || "es";
    my @tests = split(//,$testStr);
    my %tests = map { $_ => eval("sub { -$_ shift }")} @tests;
    my $failed = "";
    foreach my $test (@tests){
	$failed .= $test unless($tests{$test}->($file));
    }
    return $failed;
}

1;
