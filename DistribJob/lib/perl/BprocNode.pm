package DJob::DistribJob::BprocNode;

use DJob::DistribJob::Node;
use CBIL::Util::Utils;

@ISA = qw(DJob::DistribJob::Node);

use strict;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotcount) = @_;
    my $self = &DJob::DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotcount);
    $self->{nodeDir} = $nodeDir;
    $self->_initNodeDir();  
    return $self;
}

sub _initNodeDir {
    my ($self) = @_;

    if ($self->_fileExists($self->{nodeDir})) {
	$self->runCmd("/bin/rm -r $self->{nodeDir}");
    }

    my $try = 0;
    until($self->_fileExists($self->{nodeDir})){
	die "Can't create $self->{nodeDir} on node $self->{nodeNum}" 
	    if ($try++ > 3);
	$self->runCmd("mkdir -p $self->{nodeDir}");
    }
}

sub runCmd {
    my ($self, $cmd) = @_;
    
    return &CBIL::Util::Utils::runCmd("bpsh -n $self->{nodeNum} $cmd");
}

sub execSubTask {
    my ($self, $nodeRunDir, $serverResultDir, $cmd) = @_;
    
    $cmd = "subtaskInvoker $nodeRunDir $serverResultDir " . $cmd;

    my $rc = system("bpsh -n $self->{nodeNum} $cmd &");
    my $status = $rc >> 8;
    die "Failed with status $status running '$cmd'\n" if $status;
}

sub _fileExists {
  my($self, $file) = @_;

  for(my $a = 0; $a < 2; $a++){
    my $test = `bpsh -n $self->{nodeNum} find $file 2> /dev/null`;
    return 1 if $test =~ /$file/;
  }
}

1;
