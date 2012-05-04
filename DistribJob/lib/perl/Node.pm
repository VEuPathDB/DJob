package DJob::DistribJob::Node;

use DJob::DistribJob::NodeSlot;
use POSIX ":sys_wait_h";
use strict;
require Exporter;
use Cwd;
use IO::Socket;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'states' => [ qw( $NOCONNECTION $QUEUED $READYTORUN $READYTOINITTASK $INITIALIZINGTASK $RUNNINGTASK $COMPLETE $FAILEDNODE @failedNodes) ] ); 

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'states'} } );

## Node states
our $NOCONNECTION = 0;
our $QUEUED = 1;
our $READYTORUN = 2;
our $READYTOINITTASK = 3;
our $INITIALIZINGTASK = 4;
our $RUNNINGTASK = 5;
our $COMPLETE = 6;
our $FAILEDNODE = -1;
our @failedNodes;  ##put failed nodes here and let the controller clean them up so don't immediately get released back into circulation so get again.

my $endMatchString = 'FooCmdEnd';
my $endCmdString = "echo \$?.$endMatchString";
my %gfServerPort;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotCount, $runTime, $fileName, $serverHost, $serverPort, $procsPerNode, $memPerNode, $queue, $masterDir) = @_;

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
    $self->{procsPerNode} = $procsPerNode;
    $self->{memPerNode} = $memPerNode;
    $self->{queue} = $queue if $queue;
    $self->{masterDir} = $masterDir;

    $self->setState($NOCONNECTION);

    $self->setSaveForCleanup(0);

    ##bit of a hack but want to generate a unique high number to use for gfServer port;
    my $port = int(rand(5000)) + 8000;
    while($gfServerPort{$port}){
      $port = int(rand(5000)) + 8000;
    }
    $gfServerPort{$port} = 1;
    $self->{gfport} = $port;
      

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

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    die "ERROR: getNodeAddress failed ... no nodeNum defined ... check nodeMsgs in Controller.pm and / or over-ride the getNodeAddress method in the specific Node class to retrieve it\n";
  }
  return $self->{nodeNum};
}

# must be over ridden in node objects if specific initialization is necessary
sub _init {
  my $self = shift;
  my $ct = 0;
  while (1) {
    last if $self->getNodeAddress();
    sleep $ct < 8 ? 5 : 60;
    $ct++;
  }
  if (!$self->checkNode()) {
    print "Node $self->{nodeNum} is not responding to commands....skipping\n";
    $self->failNode();
    return 0;
  }
  return $self->_initNodeDir(); 
}

sub failNode {
  my($self) = @_;
  return if $self->getState() == $FAILEDNODE; ##already failed so don't do again!
  $self->setState($FAILEDNODE);
  push(@failedNodes,[time(),$self]);
}

sub _initNodeDir {
  my ($self) = @_;

if ($self->_fileExists($self->{nodeDir})) {
    $self->runCmd("/bin/rm -r $self->{nodeDir}",1);
  }

  my $try = 0;
  do {
    if($try++ > 3){
      $self->failNode();
      return 0;
    }
#    die "Can't create $self->{nodeDir} on node $self->{nodeNum}" if ($try++ > 3);
    $self->runCmd("mkdir -p $self->{nodeDir}",1);
  } until ($self->_fileExists($self->{nodeDir})); 
  return 1;
}

sub runCmd {
  my ($self, $cmd, $ignoreErr) = @_;
  return if $self->getState() >= $COMPLETE || $self->getState() == $FAILEDNODE;
  my $sock = $self->getPort();
  if(!$sock){
    print "Failed to get Sock for $self->{nodeNum}\n";
    $self->failNode();
    $self->setErr("No Socket");
    return undef;
  }
  print $sock "$cmd\n";
  my $res = "";
  while(<$sock>){
    if(/^(\d+)\.$endMatchString/){
      my $err = $1;
      $self->setErr($err);
      if($err && !$ignoreErr){
        print "Node ".$self->getNodeAddress().": Failed with status $err running '$cmd' ... ";
        if($self->checkNode()){
          print "node is OK so not inactivating\n";
          return undef;
        }else{
          print "Inactivating Node\n";
          $self->failNode();
          return undef;
        }
        $self->setErr($err); ##need to set back to previous since checkNode will over-write
      }
      last;
    }
    $res .= $_;
  }
  return $res;
}

sub setErr {
  my($self,$err) = @_;
  $self->{cmdErr} = $err;
}

sub getErr {
  my($self) = @_;
  return $self->{cmdErr};
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
        if($ct++ > 5){
          print STDERR "Could not create socket: $!\nInactivating node".$self->getNum()."\n" ;
          $self->failNode();
          last;
        }
        sleep 2;
        next;
      }
      $self->{portCon} = $sock;
    }
  }
  if($self->{portCon} && $self->{portCon}->connected()){
    return $self->{portCon};
  }else{
    return undef;
  }
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

sub setDir {
  my($self,$dir) = @_;
  $self->{nodeDir} = $dir if $dir;  ## don't overwrite if null
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
    if($ns->getTask()){
      my $subt = $ns->getTask();
      $subt->setRedoSubtask(1);
      push(@st,$subt);
    }
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
  my $res = $self->runCmd("ls $self->{masterDir}",1);
  return $res ? 1 : 0;  
}

## saving node for cleanup
sub setSaveForCleanup {
  my($self,$val) = @_;
  $self->{saveForCleanup} = $val;
}

sub getSaveForCleanup {
  my($self) = @_;
  return $self->{saveForCleanup};
}


##want to only clean up if both slots are finished
sub cleanUp {
  my ($self,$force, $state) = @_;

  if(!$self->getSaveForCleanup() || ($self->getSaveForCleanup() && !$force)) { 
    return if $self->getState() >= $COMPLETE && !$force; ##already cleaned up
  }
    
  if (!$force) {
    foreach my $slot (@{$self->getSlots()}) {
      return unless $slot->isFinished();
    }
  }

  ##want to kill any child processes still running to quit cleanly
  if($self->getState() == $INITIALIZINGTASK && $self->{taskPid}){
    kill(1, $self->{taskPid}) unless waitpid($self->{taskPid},1);
  }

  $self->setState($state ? $state : $COMPLETE); ##complete

  ## if saving this one so don't clean up further and release
  return if($self->getSaveForCleanup() && !$force);  

  if($state != $FAILEDNODE){  ## if the node has failed don't want to run commands on it ...

    print "Cleaning up node $self->{nodeNum}...\n";
    
    ##now call the task->cleanUpNode method to enable  users to stop processes running on node
    my $task = $self->getTask();
    $task->cleanUpNode($self) if $task;

    if($self->{nodeNum} && $self->getPort()){
      $self->runCmd("/bin/rm -r $self->{nodeDir}", 1);
      $self->runCmd("closeAndExit",1);
      $self->closePort();
    }
  }
}

1;
