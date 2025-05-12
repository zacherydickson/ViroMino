package MVPipe::Logger;
use strict;

use Carp;
use Scalar::Util qw(openhandle);

our %_LOG_LEVEL = (ERROR => 0, WARNING => 1, INFO => 2, DEBUG => 3);
our %LEVEL_LOG = reverse(%_LOG_LEVEL);
our $_DEFAULT_LEVEL = "WARNING";

sub new {
    my ($package, %Args) = @_;
    my $loglevel = (exists $Args{level}) ? uc($Args{level}) : $_DEFAULT_LEVEL;
    my $handle = (exists $Args{handle}) ? $Args{handle} : \*STDERR;
    my $title = (exists $Args{title}) ? $Args{title} : undef;
    unless(exists $_LOG_LEVEL{$loglevel}){
        my @keys = sort keys %_LOG_LEVEL;
        croak("Logger level must be in @keys");
    }
    delete $Args{level};
    delete $Args{handle};
    delete $Args{title};
    carp("Unrecognized arguments to Logger::new (@{[keys %Args]})") if(scalar(keys %Args));

    my $self = {'Logger::level' => $_LOG_LEVEL{$loglevel}, 'Logger::handle' => $handle, 'Logger::title' => $title};
    bless $self => $package;

    return $self;
}

#Gets/sets the level of logs which will be written
sub level (){
    my $self = shift;
    my $level;
    if(@_){
        $level = uc(shift);
        carp("Too many arguments to Logger::level") if(@_);
    }
    if(defined $level){
        unless(exists $_LOG_LEVEL{$level}){
            my @keys = sort keys %_LOG_LEVEL;
            croak("Logger level must be in @keys");
        }
        $self->{'Logger::level'} = $_LOG_LEVEL{$level};
    }
    return {reverse %_LOG_LEVEL}->{$self->{'Logger::level'}};
}

#Gets/sets the handle to which logs should be written
sub handle(){
    my $self = shift;
    my $handle;
    if(@_){
        $handle = openhandle(shift);
        unless(defined $handle){
            croak("Logger::handle must be provided an open file handle");
        }
        carp("Too many arguments to Logger::handle") if(@_);
    }
    if(defined $handle){
        $self->{'Logger::handle'} = $handle;
    }
    return $self->{'Logger::level'};
}

#Gets/sets the handle to which logs should be written
sub title(){
    my $self = shift;
    my $title;
    if(@_){
        $title = shift;
        carp("Too many arguments to Logger::title") if(@_);
    }
    if(defined $title){
        $self->{'Logger::title'} = $title;
    }
    return $self->{'Logger::title'};
}

#Writes a message to STDERR with a timestamp and the level of severity, exits if necessary
sub Log(){
    my ($self, $message, $level) = @_;
    croak("Logger::Log requires a message") unless defined $message;
    $level = (defined $level) ? uc($level) : $_DEFAULT_LEVEL;
    carp("Too many arguments to Logger:Log") if(@_ > 3);
    return if($_LOG_LEVEL{$level} > $self->{'Logger::level'});
    my ($sec,$min,$hour) = localtime;
    my $prevHandle = select($self->{'Logger::handle'});
    my $titleStr = (defined $self->{'Logger::title'}) ? "(".$self->{'Logger::title'}.") " : "";
    printf("%02d:%02d:%02d - %s[%s] %s\n",$hour,$min,$sec,$titleStr,$level,$message);
    exit 1 if($_LOG_LEVEL{$level} == 0);
    select($prevHandle);
}

#Given a hash of parameters, outputs a summary
#Summary is printed in the order in which arguments are provided
#Allowing parameters to be provided as an array with keys in even elements and values in odd elements
sub LogParameters(){
    my ($self,%params) = @_;
    my @order;
    for(my $i = 1; $i < @_; $i += 2){
        push(@order,@_[$i]);
    }
    return if($self->{'Logger::level'} < $_LOG_LEVEL{INFO});
    my $prevHandle = select($self->{'Logger::handle'});
    print join("",(('=') x 60)),"\n";
    print "Initialized on ".(localtime)."\n";
    print "Parameters:\n";
    my $nchar = 10;
    foreach (keys %params){
        $nchar = length($_) if(length($_) > $nchar);
    }
    foreach my $param (@order){
        printf("%*s : %s\n",$nchar + 1,$param,$params{$param});
    }
    print join("",(('=') x 60)),"\n";
    select($prevHandle);
}

1; #Exits properly
