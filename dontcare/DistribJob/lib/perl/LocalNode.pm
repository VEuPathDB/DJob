package DistribJob::LocalNode;

use DistribJob::Node;
use Common::Utils;


@ISA = qw(DistribJob::Node);

use strict;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotcount) = @_;
    my $self = &DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotcount);
    $self->{nodeDir} = "$nodeDir/$nodeNum";
    $self->_initNodeDir();  
    return $self;
}

sub _initNodeDir {
    my($self) = @_;
    if (-e $self->{nodeDir}) {
	$self->runCmd("/bin/rm -r $self->{nodeDir}");
    }
    $self->runCmd("mkdir -p $self->{nodeDir}");
}

sub runCmd {
    my ($self, $cmd) = @_;
    return &Common::Utils::runCmd($cmd);
}
    
sub execSubTask {
    my ($self, $nodeRunDir, $serverResultDir, $cmd) = @_;
    
    $cmd = "subtaskInvoker $nodeRunDir $serverResultDir " . $cmd;
    my $rc = system("$cmd &");
    my $status = $rc >> 8;
    die "Failed with status $status running '$cmd'\n" if $status;
}

1;



