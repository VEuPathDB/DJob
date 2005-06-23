package DJob::DistribJob::Node;

use DJob::DistribJob::NodeSlot;
use POSIX ":sys_wait_h";
use strict;
require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'states' => [ qw( $NOCONNECTION $QUEUED $READYTOINITTASK $INITIALIZINGTASK $RUNNINGTASK $COMPLETE $FAILEDNODE ) ] ); 

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'states'} } );

## Node states
our $NOCONNECTION = 0;
our $QUEUED = 1;
our $READYTOINITTASK = 2;
our $INITIALIZINGTASK = 3;
our $RUNNINGTASK = 4;
our $COMPLETE = 5;
our $FAILEDNODE = 6;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotCount, $runTime, $fileName) = @_;

    my $self = {};
    bless($self, $class);
    $self->{nodeNum} = $nodeNum;
    $self->{nodeDir} = $nodeDir;
    $self->{slotCount} = $slotCount;
    $self->{runTime} = $runTime;
    $self->{fileName} = $fileName;

    $self->setState($NOCONNECTION);

    return $self;
}

sub initialize {
  my $self = shift;
  if($self->{initPid}){  ##am already forked so  check to see if complete
    if(waitpid($self->{initPid},1)){
      $self->setState($READYTOINITTASK);
    }
    return; 
  }
  $self->setState($QUEUED);
  my $pid;
 FORK: {
    if($pid = fork){
      $self->{initPid} = $pid;
    }elsif(defined $pid) {  
      $self->_init();
      exit;
    }elsif($! =~ /No more process/){
      print STDERR "Forking failure: $!\n";
      sleep 1;
      redo FORK;
    } else {
      die "Unable to fork: $!\n";
    }
  } 

}

# must be over ridden in node objects if specific initialization is necessary
sub _init {
  my $self = shift;
  ##put node specific code here

  $self->_initNodeDir();
}

sub _initNodeDir {
    my($self) = @_;
    if (-e $self->{nodeDir}) {
	$self->runCmd("/bin/rm -r $self->{nodeDir}");
    }
    $self->runCmd("mkdir -p $self->{nodeDir}");
}

sub getDir {
    my ($self) = @_;
    return $self->{nodeDir};
}

sub getNum {
    my ($self) = @_;
    $self->initialize() unless defined $self->{nodeNum};
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

sub getSlotCount {
    my($self) = @_;
    return $self->{slotCount};
}

##want to track state here...
#possible states are:
#  $NOCONNECTIOIN || 0 => have no connection to a compute node
#  $QUEUED || 1 => have requested a connection and am waiting for response from scheduler
#  $READYTOINITTASK || 2 => have connection but task is not inititalized
#  $INITIALIZINGTASK || 3 => task is initializing
#  $RUNNINGTASK || 4 => task is initialized and node is ready to run (or running)
#  $COMPLETE || 5 => node has completed work and released link to compute node (cleaned up)
#  $FAILEDNODE || 6 => node is non functional, skipping

sub getState {
  my $self = shift;
  ##need to deal with the in between states..check if the thread is still running
  if($self->{state} == $QUEUED){
    ##call initialize here...check to see state there
    $self->initialize();
  }elsif($self->{state} == $INITIALIZINGTASK){
    $self->setState($RUNNINGTASK) if waitpid($self->{taskPid},1);
  }
  return $self->{state}
}

sub setState {
  my($self,$s) = @_;
  $self->{state} = $s;
}


# also need to have node ask the task to initNode so can do it in a thread
sub initializeTask {
  my($self,$task,$inputDir) = @_;
  $self->setState($INITIALIZINGTASK);
  my $pid;
 FORK: {
    if($pid = fork){
      $self->{taskPid} = $pid;
    }elsif(defined $pid) {  
      $self->_initTask($task,$inputDir);
      exit;
    }elsif($! =~ /No more process/){
      print STDERR "Forking failure: $!\n";
      sleep 1;
      redo FORK;
    } else {
      die "Unable to fork: $!\n";
    }
  } 
}

sub _initTask {
  my($self,$task,$inputDir) = @_;
  $task->initNode($self,$inputDir);
  print "Task initialized on node ".$self->getNum()."\n";
}

sub DESTROY {
  my($self) = @_;
  kill(1, $self->{initPid}) unless waitpid($self->{initPid},1);
  kill(1, $self->{taskPid}) unless waitpid($self->{taskPid},1);
}

##want to only clean up if both slots are finished
sub cleanUp {
    my ($self,$force, $state) = @_;  ##note that if $force is true then will not check if slots are finished
    
    return if $self->getState() >= $COMPLETE;  ##already cleaned up
    
    if(!$force){
      foreach my $slot (@{$self->getSlots()}){
        return unless $slot->isFinished();
      }
    }

    ##want to kill any child processes still running to quit cleanly
    if($self->getState() == $QUEUED && $self->{initPid}){
      kill(1, $self->{initPid}) unless waitpid($self->{initPid},1);
    }elsif($self->getState() == $INITIALIZINGTASK && $self->{taskPid}){
      kill(1, $self->{taskPid}) unless waitpid($self->{taskPid},1);
    }

    $self->runCmd("/bin/rm -r $self->{nodeDir}");
    $self->setState($state ? $state : $COMPLETE);  ##complete
}

1;
