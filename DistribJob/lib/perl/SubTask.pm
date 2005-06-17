package DJob::DistribJob::SubTask;

use strict;

sub new {
    my ($class, $subTaskNum, $subTaskDir, $nodeSlot, $task) = @_;

    my $self = {};
    bless $self;
    $self->{subTaskNum} = $subTaskNum;
    $self->{subTaskDir} = $subTaskDir;
    $self->{task} = $task;
    $self->{subTaskResultDir} = "$subTaskDir/result";
    $self->{nodeSlot} = $nodeSlot;
    $self->{startTime} = time;
    return $self;
}

sub setNodeSlot {
  my($self,$ns) = @_;
  $self->{nodeSlot} = $ns;
}

sub getNodeSlot {
  my($self) = @_;
  return $self->{nodeSlot};
}

sub getDir {
    my ($self) = @_;
    return $self->{subTaskDir};
}

sub getResultDir {
    my ($self) = @_;
    return $self->{subTaskResultDir};
}

sub getNum {
    my ($self) = @_;
    return $self->{subTaskNum};
}

sub complete {
    my ($self) = @_;

    my $complete = $self->_isDone();
    if ($complete) {
      $self->{task}->addSubtaskTime($self->getRunningTime());
	if ($self->_isFailure()) {
	    $self->{task}->failSubTask($self);
	} else {
	    $self->{task}->passSubTask($self);
	}
    }
    return $complete;
}

sub getRunningTime {
  my $self = shift;
  return time - $self->{startTime};
}

sub _isDone {
    my ($self) = @_;
    return (-e "$self->{subTaskResultDir}/done");
}

sub _isFailure {
    my ($self) = @_;
    return (-e "$self->{subTaskResultDir}/failed");
}

1;
