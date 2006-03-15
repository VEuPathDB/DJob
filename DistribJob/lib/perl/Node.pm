package DJob::DistribJob::Node;

use DJob::DistribJob::NodeSlot;
use POSIX ":sys_wait_h";
use strict;
require Exporter;
use Cwd;
use IO::Socket;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'states' => [ qw( $NOCONNECTION $QUEUED $READYTORUN $READYTOINITTASK $INITIALIZINGTASK $RUNNINGTASK $COMPLETE $FAILEDNODE ) ] ); 

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'states'} } );

## Node states
our $NOCONNECTION = 0;
our $QUEUED = 1;
our $READYTORUN = 2;
our $READYTOINITTASK = 3;
our $INITIALIZINGTASK = 4;
our $RUNNINGTASK = 5;
our $COMPLETE = 6;
our $FAILEDNODE = 7;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotCount, $runTime, $fileName, $serverHost, $serverPort) = @_;

    my $self = {};
    bless($self, $class);
    $self->{nodeNum} = $nodeNum;
    $self->{nodeDir} = $nodeDir;
    $self->{slotCount} = $slotCount;
    $self->{runTime} = $runTime;
    $self->{fileName} = $fileName;
    $self->{serverHost} = $serverHost;
    $self->{serverPort} = $serverPort;
    $self->{cwd} = getcwd();

    $self->setState($NOCONNECTION);

    return $self;
}

##queue the node here...will need to be node specific...
sub queueNode {
  my $self = shift;
  if(!$self->{jobid}){
    ##submit node here to the queue..over-ride in subclasses
  }
  $self->setState($QUEUED);
}

sub initialize {
  my($self) = @_;
  if(!$self->getJobid()){
    $self->queueNode();
  }
  return unless ($self->getState() == $READYTORUN); 
  if($self->_init()){
    print "Node $self->{nodeNum} initialized\n";
    $self->setState($READYTOINITTASK);
  }else{
    print STDERR  "\nERROR: unable to initialize node $self->{nodeNum} ... marking FAILEDNODE\n";
  }
}

# must be over ridden in node objects if specific initialization is necessary
sub _init {
  my $self = shift;
  my $ct = 0;
  while (1) {
    last if $self->getNodeAddress();
    sleep $ct < 8 ? 15 : 120;
    $ct++;
  }
  if (!$self->checkNode()) {
    print "Node $self->{nodeNum} is not responding to commands....skipping\n";
    $self->cleanUp(1, $FAILEDNODE);  
    return;
  }
  return $self->_initNodeDir(); 
}

sub _initNodeDir {
  my ($self) = @_;

  if ($self->_fileExists($self->{nodeDir})) {
    $self->runCmd("/bin/rm -r $self->{nodeDir}",1);
  }

  my $try = 0;
  do {
    if($try++ > 3){
      $self->cleanUp(1, $FAILEDNODE);  
      return 0;
    }
#    die "Can't create $self->{nodeDir} on node $self->{nodeNum}" if ($try++ > 3);
    $self->runCmd("mkdir -p $self->{nodeDir}",1);
  } until ($self->_fileExists($self->{nodeDir})); 
  return 1;
}

my $endMatchString = 'FooCmdEnd';
my $endCmdString = "echo \$?.$endMatchString";
sub runCmd {
  my ($self, $cmd, $ignoreErr) = @_;
  my $sock = $self->getPort();
  print $sock "$cmd\n";
  my $res = "";
  while(<$sock>){
    if(/^(\d+)\.$endMatchString/){
      if($1 && !$ignoreErr){
        print STDERR "Node ".$self->getNodeAddress().": Failed with status $1 running '$cmd' ... Inactivating Node\n" ;
#        sleep 100;  ##uncomment if need to test problems on the nodes
        print $sock "closeAndExit\n";  #exits the script on the node..
        close $sock;
        $self->cleanUp(1, $FAILEDNODE);  
#        exit(1);
      }
      last;
    }
    $res .= $_;
  }
  return $res;
}

sub getPort {
  my $self = shift;
  if(!$self->{portCon}){
    ##note..need to try a couple of times here because the script may be running but the port is not ready on the nodes to receive connections...seems to take some time
#    print STDERR "Creating new port connection\n";
    my $sock;
    my $ct = 0;
    until($sock){
      $sock = new IO::Socket::INET (
                                    PeerAddr => $self->getNodeAddress(),
                                    PeerPort => $self->getLocalPort(),
                                    Proto => 'tcp',
                                   );
      unless($sock){
        die "Could not create socket: $!\n" if $ct++ > 5;
        sleep 3;
        next;
      }
      $self->{portCon} = $sock;
    }
  }
  return $self->{portCon};
}

 sub closePort {
  my $self = shift;
  close $self->{portCon};
  undef $self->{portCon};
}
sub execSubTask {
  my ($self, $nodeRunDir, $serverSubtaskDir, $cmd) = @_;
    
  $cmd = "subtaskInvoker $nodeRunDir $serverSubtaskDir/result $self->{jobid} $self->{serverHost} $self->{serverPort} $cmd &";

  return $self->runCmd($cmd);
}

sub _fileExists {
  my($self, $file) = @_;
  for (my $a = 0; $a < 2; $a++) {
    my $test = $self->runCmd("find $file 2> /dev/null", 1);
    return 1 if $test =~ /$file/;
  }
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

sub setNum {
  my($self,$num) = @_;
  $self->{nodeNum} = $num;
}

sub setLocalPort {
  my($self,$port) = @_;
  $self->{localPort} = $port;
}

sub getLocalPort {
  my $self = shift;
  return $self->{localPort};
}

sub getSlots {
    my ($self) = @_;
    unless ($self->{slots}) {
	$self->{slots} = [];
        $self->{slotHash};
	for(my $i=1; $i <= $self->{slotCount}; $i++) {
          my $slot = DJob::DistribJob::NodeSlot->new($self, $i);
          push(@{$self->{slots}}, $slot);
          $self->{slotHash}->{"slot_$i"} = $slot;
	}
    }

    return $self->{slots};
}

sub getSlot {
  my($self,$s) = @_;
  $self->getSlots() unless $self->{slots};
  return $self->{slotHash}->{$s};
}

sub getSlotCount {
    my($self) = @_;
    return $self->{slotCount};
}

##want to track state here...
#possible states are:
#  $NOCONNECTIOIN || 0 => have no connection or am queued 
#  $READYTORUN || 1 => have been placed into the running queue
#  $READYTOINITTASK || 2 => have connection but task is not inititalized
#  $INITIALIZINGTASK || 3 => task is initializing
#  $RUNNINGTASK || 4 => task is initialized and node is ready to run (or running)
#  $COMPLETE || 5 => node has completed work and released link to compute node (cleaned up)
#  $FAILEDNODE || 6 => node is non functional, skipping

sub getState {
  my $self = shift;
  ##need to deal with the in between states..check if the thread is still running
  if($self->{state} == $INITIALIZINGTASK){
    $self->setState($RUNNINGTASK) if waitpid($self->{taskPid},1);
  }
  return $self->{state}
}

sub setState {
  my($self,$s) = @_;
  $self->{state} = $s;
}

sub setJobid { my ($self,$jid) = @_; $self->{jobid} = $jid; }
sub getJobid { my $self = shift; return $self->{jobid}; }


# also need to have node ask the task to initNode so can do it in a thread
sub initializeTask {
  my($self,$task,$inputDir) = @_;
  $self->setState($INITIALIZINGTASK);
  $self->{task} = $task;
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

##returns stored task object for use later...
sub getTask {
  my $self = shift;
  return $self->{task};
}

sub _initTask {
  my($self,$task,$inputDir) = @_;
  $task->initNode($self,$inputDir);
  print "Task initialized on node ".$self->getNum()."\n";
}

sub failedSoGetSubtasks {
  my $self = shift;
  return if $self->{retrievedFailedSubtasks};
  my @st;
  foreach my $ns (@{$self->getSlots()}){
    push(@st,$ns->getTask()) if $ns->getTask();
  }
  $self->{retrievedFailedSubtasks} = 1;
  return @st;
}

sub DESTROY {
  my($self) = @_;
  kill(1, $self->{taskPid}) unless waitpid($self->{taskPid},1);
}

sub checkNode {
  my($self) = @_;
  return 1;  ##implement if there are problems..
}


##want to only clean up if both slots are finished
sub cleanUp {
  my ($self,$force, $state) = @_;

  return if $self->getState() >= $COMPLETE; ##already cleaned up
    
  if (!$force) {
    foreach my $slot (@{$self->getSlots()}) {
      return unless $slot->isFinished();
    }
  }

  ##want to kill any child processes still running to quit cleanly
  if($self->getState() == $INITIALIZINGTASK && $self->{taskPid}){
    kill(1, $self->{taskPid}) unless waitpid($self->{taskPid},1);
  }

  print "Cleaning up node $self->{nodeNum}...\n";

  if($self->{portCon}){
    ##now call the task->cleanUpNode method to enable  users to stop processes running on node
    $self->getTask()->cleanUpNode($self);
    $self->runCmd("/bin/rm -r $self->{nodeDir}", 1);
    $self->runCmd("closeAndExit");
    $self->closePort();
  }
  $self->setState($state ? $state : $COMPLETE); ##complete
}

1;
