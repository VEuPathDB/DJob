package DJob::DistribJob::Node;

use DJob::DistribJob::NodeSlot;

use strict;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotCount) = @_;

    my $self = {};
    bless($self, $class);
    $self->{nodeNum} = $nodeNum;
    $self->{nodeDir} = $nodeDir;
    $self->{slotCount} = $slotCount;

    return $self;
}

sub getDir {
    my ($self) = @_;
    return $self->{nodeDir};
}

sub getNum {
    my ($self) = @_;
    return $self->{nodeNum};
}

sub getSlots {
    my ($self) = @_;
    unless ($self->{slots}) {
	$self->{slots} = [];
	for(my $i=0; $i<$self->{slotCount}; $i++) {
	    push(@{$self->{slots}}, DJob::DistribJob::NodeSlot->new($self, $i+1));
	}
    }

    return $self->{slots};
}

sub cleanUp {
    my ($self) = @_;

    unless ($self->{clean}) {
	$self->runCmd("rm -r $self->{nodeDir}");
	$self->{clean} = 1;
    }
}

1;
