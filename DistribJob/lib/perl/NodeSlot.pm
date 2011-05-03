package DJob::DistribJob::NodeSlot;

use strict;

sub new {
    my ($class, $node, $slotNum) = @_;

    my $self = {};
    bless $self;
    $self->{node} = $node;
    $self->{slotNum} = $slotNum;
    $self->{nodeSlotDir} = $node->getDir() . "/slot_$slotNum";
    $self->_initDir();

    return $self;
}

sub assignNewTask {
    my ($self, $task) = @_;
    $self->{task} = $task;
}

sub getTask {
  my $self = shift;
  return $self->{task};
}

sub taskComplete {
    my ($self) = @_;

    my $complete = 1;
    if ($self->{task}) {
	if ($self->{task}->complete()) {
	    $self->{task} = undef;
	    $complete = 1;
	} else {
	    $complete = 0;
	}
    }
    return $complete;
}

sub isRunning {
    my ($self) = @_;
    return $self->{task};
}

sub getDir {
    my ($self) = @_;
    return $self->{nodeSlotDir};
}

sub getNodeNum {
    my ($self) = @_;
    return $self->{node}->getNum();
}

sub getNode {
    my ($self) = @_;
    return $self->{node};
}

sub getNum {
    my ($self) = @_;
    return $self->{slotNum};
}

sub _initDir {
    my ($self) = @_;
    
    $self->{node}->runCmd("mkdir -p $self->{nodeSlotDir}");
}

sub isFinished {
    my $self = shift;
    return $self->{finished};
}

sub cleanUp {
    my ($self) = @_;
    $self->{finished} = 1;
    $self->{node}->cleanUp(); 
}

1;
