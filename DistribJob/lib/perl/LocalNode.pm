package DJob::DistribJob::LocalNode;

use DJob::DistribJob::Node;
use CBIL::Util::Utils;


@ISA = qw(DJob::DistribJob::Node);

use strict;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotcount) = @_;
    my $self = &DJob::DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotcount);
    $self->{nodeDir} = "$nodeDir/$nodeNum";
    $self->_initNodeDir();  
    return $self;
}


sub runCmd {
    my ($self, $cmd) = @_;
    return &CBIL::Util::Utils::runCmd($cmd);
}
    
sub execSubTask {
    my ($self, $nodeRunDir, $serverSubtaskDir, $cmd) = @_;
    
    $cmd = "subtaskInvoker $nodeRunDir $serverSubtaskDir/result " . $cmd;
    my $rc = system("$cmd &");
    my $status = $rc >> 8;
    die "Failed with status $status running '$cmd'\n" if $status;
}

1;



