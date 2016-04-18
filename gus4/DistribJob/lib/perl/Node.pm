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
    my ($class, $nodeNum, $nodeWorkingDirsHome, $slotCount, $runTime, $fileName, $serverHost, $serverPort, $procsPerNode, $memPerNode, $queue, $masterDir) = @_;

    my $self = {};
    bless($self, $class);
    $self->{nodeNum} = $nodeNum;
    $self->{nodeWorkingDirsHome} = $nodeWorkingDirsHome;
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

sub getQueue {
  my $self = shift;
  return $self->{queue};
}

##queue the node here...will need to be node specific...
sub queueNode {
  my $self = shift;
  if(!$self->{jobid}){
    ##submit node here to the queue..over-ride in subclasses
  }
  $self->setState($QUEUED);
}

##can call if node has failed in order to get new active job in cluster
sub reQueueNode {
  my $self = shift;
  print "     reQueueing node ".$self->getJobid()." - new JobId = ";
  $self->getTask()->addRedoSubtasks($self->failedSoGetSubtasks()) if $self->getTask();
  $self->setJobid("");
  undef $self->{portCon};
  $self->queueNode();
  print $self->getJobid()."\n";
}

##put all the logic here for what to do if getQueueState doesn't return true
sub manageNodeBasedOnQueueState {
  my $self = shift;
  return 1 if ($self->getState() == $FAILEDNODE || $self->getState() == $COMPLETE);
  my $ret = $self->getQueueState();
  if($ret == 0){
    print STDERR "ERROR: Node '".$self->getJobid()."' failed (no longer in queue)\n";
    if($self->{countCheckFailures} < 3) {
      $self->reQueueNode();
    }else{
      $self->failNode();
    }
    $self->{countCheckFailures}++;
  }
}

##NOTE: this method simply returns whether the job scheduler still has this node in the queue irrexpective of it's state: could be running, queued, waiting etc but returns 1.
sub getQueueState {
  my $self = shift;

  return 1 if $self->getState() == $FAILEDNODE || $self->getState() == $COMPLETE; 

  ##should not be in queue
  my $jobid = $self->getJobid();
  if(!$jobid){
    print STDERR "SgeNode->getQueueState: unable to checkQueueStatus as can't retrieve JobID\n";
    return 0;
  }
  return $self->runJobStatusCheck($jobid);
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
    print "\nERROR: unable to initialize node $self->{nodeNum} ... marking FAILEDNODE\n";
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
  return $self->_initWorkingDir(); 
}

sub failNode {
  my($self) = @_;
  return if $self->getState() == $FAILEDNODE; ##already failed so don't do again!
  print "ERROR: node ".$self->getJobid()." is no longer functional ... failing node\n";
  $self->getTask()->addRedoSubtasks($self->failedSoGetSubtasks()) if $self->getTask();
  $self->setState($FAILEDNODE);
  push(@failedNodes,[time(),$self]);
}

sub _initWorkingDir {
  my ($self) = @_;

if ($self->_fileExists($self->{workingDir})) {
    $self->runCmd("/bin/rm -r $self->{workingDir}",1);
  }

  my $try = 0;
  do {
    if($try++ > 3){
      $self->failNode();
      return 0;
    }
#    die "Can't create node job dir $self->{workingDir} on node $self->{nodeNum}" if ($try++ > 3);
    $self->runCmd("mkdir -p $self->{workingDir}",1);
  } until ($self->_fileExists($self->{workingDir})); 
  return 1;
}

sub runCmdExitIfFail {
  my ($self, $cmd, $ignoreErr,$checkingNode) = @_;
  $self->runCmd($cmd, $ignoreErr,$checkingNode,1);
}

sub runCmd {
  my ($self, $cmd, $ignoreErr,$checkingNode,$exitIfFail) = @_;
#  print "runCmd: '$cmd',$ignoreErr,$checkingNode\n";
  return undef if $self->getState() >= $COMPLETE || $self->getState() == $FAILEDNODE;
  my $sock = $self->getPort();
  if(!$sock){
    print "Failed to get Sock for $self->{nodeNum} ($self->{jobid}) running command '$cmd'\n";
    $self->failNode();
    $self->setErr(1);
    die "Failure making socket connection\n" if $exitIfFail;
    return undef;
  }
  print $sock "$cmd\n";
  my $res = "";
  while(<$sock>){
    if(/^(\d+)\.$endMatchString/){
      my $err = $1;
      $self->setErr($err) unless $checkingNode;
      if($err && !$ignoreErr){
        print "Node ".$self->getNodeAddress()." (".$self->getJobid()."): Failed with status $err running '$cmd' ... ";
        if($exitIfFail){
          die "Command called in exit if failure mode so exiting\n";
          
        }elsif($self->checkNode()){
          print "node is OK so not inactivating\n";
          return undef;
        }else{
          print "Inactivating Node\n";
          $self->failNode();
          return undef;
        }
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
#    print STDERR $self->getNum().": Creating new port connection PeerAddr => ".$self->getNodeAddress().", PeerPort => ".$self->getLocalPort()."\n";
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
          print "Could not create socket: $!\nInactivating node: ".$self->getNum()."\n" ;
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

# the directory the job runs in.
sub getWorkingDir {
    my ($self) = @_;
    return $self->{workingDir}; 
}

sub setWorkingDir {
    my ($self, $dir) = @_;
    $self->{workingDir} = $dir if $dir;  ## don't overwrite if null
}

# the dir that holds job dirs.  set at initialization, by config, and not changed, unless using tmp dir
sub getNodeWorkingDirsHome {
    my ($self) = @_;
    return $self->{nodeWorkingDirsHome}; 
}

sub setNodeWorkingDirsHome {
  my($self,$dir) = @_;
  $self->{nodeWorkingDirsHome} = $dir if $dir;  ## don't overwrite if null
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
    $self->{slotHash} = {};
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
  $self->_initTask($task,$inputDir);
## FORK: {
##    if($pid = fork){
##      $self->{taskPid} = $pid;
##    }elsif(defined $pid) {  
##      $self->_initTask($task,$inputDir);
##      exit;
##    }elsif($! =~ /No more process/){
##      print "Forking failure: $!\n";
##      sleep 1;
##      redo FORK;
##    } else {
##      die "Unable to fork: $!\n";
##    }
##  } 
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
  my @st;
  my @nums;
  foreach my $ns (@{$self->getSlots()}){
    if($ns->getTask()){
      my $subt = $ns->getTask();
      $subt->setRedoSubtask(1);
      push(@st,$subt);
      push(@nums,$subt->getNum());
    }
  }
  undef $self->{slots};
  print "$self->{jobid}: retrieving ".scalar(@st)." subtasks (".join(", ", @nums).")\n";
  return @st;
}

sub DESTROY {
  my($self) = @_;
  kill(1, $self->{taskPid}) unless waitpid($self->{taskPid},1);
}

sub checkNode {
  my($self) = @_;
  my $res = $self->runCmd("ls $self->{masterDir}",1,1);
  if($res){
    return 1;
  }else{
    $self->failNode();
    return 0;
  }
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

  return if $self->{cleanedUp}; #already cleaned up

  if (!$force) {
    foreach my $slot (@{$self->getSlots()}) {
      return unless $slot->isFinished();
    }
  }

  ##want to kill any child processes still running to quit cleanly
  if($self->getState() == $INITIALIZINGTASK && $self->{taskPid}){
    kill(1, $self->{taskPid}) unless waitpid($self->{taskPid},1);
  }

  ## if saving this one so don't clean up further and release
  if($self->getSaveForCleanup() && !$force){
    $self->setState($COMPLETE);  ##note that controller monitors state and resets to running once all subtasks are finished.
    return;
  }

  $self->{cleanedUp} = 1;  ##indicates that have cleaned up this node already

  print "Cleaning up node $self->{nodeNum} ($self->{jobid})\n";
  if($state != $FAILEDNODE){  ## if the node has failed don't want to run commands on it
    my $task = $self->getTask();
    $task->cleanUpNode($self) if $task;

    if($self->{nodeNum} && $self->getState() > $QUEUED && $self->getPort()){
      $self->runCmd("/bin/rm -rf $self->{workingDir}",1);
      $self->runCmd("closeAndExit",1);
      $self->closePort();
    }
  }

  if($self->getState() == $FAILEDNODE){ ##don't want to change if is failed node
    $state = $FAILEDNODE;
  }else{
    $self->setState($state == $FAILEDNODE ? $state : $COMPLETE); ##complete
  }

  ## if subclass is inclined, give it a chance to report statistics
  $self->reportJobStats();

  # node should be off queue already, because we stopped running the node script
  # but, to be safe, try to remove it from queue.
  if ($self->getQueueState()) { $self->removeFromQueue() }
}

# remove this node's job from the queue
# must be implemented by subclass
sub removeJobFromQueue {
  my ($self) = @_;
  die "removeJobFromQueue must be implemented by Node subclass";
}

# an optional method for subclasses to implement
# called at the end of node->cleanUp
# can query the que to return stats about this run
# print results to stdout
sub reportJobStats {
  my ($self) = @_;
}

# run a command to check the status of the job with this job id.
# return 1 if job id is found and valid, otherwise return 0
sub runJobStatusCheck {
  my ($self, $jobId) = @_;
  die "must be implemented by subclass";
}


# static method
sub getQueueSubmitCommand {
  my ($class, $queue) = @_;
  die "must be overridden by subclass";
}

# static method to extract Job Id from job id file text
# used to get job id for distribjob itself
sub getJobIdFromJobInfoString {
  my ($class, $jobInfoString) = @_;

  die "getJobIdFromJobInfoString() must be overridden by subclass";
}

# static method to provide command to run to get status of a job
# used to get status of distribjob itself
sub getCheckStatusCmd {
  my ($class, $jobId) = @_;

  die "parseJobIdFile() must be overridden by subclass";
}

# static method to extract status from status file
# used to check status of distribjob itself
# return 1 if still running
sub checkJobStatus {
  my ($class, $statusFileString, $jobId) = @_;

  die "checkJobStatus() must be overridden by subclass";
}


sub deleteLogFilesAndTmpDir {
  die "must override the deleteLogFilesAndTmpDir method in the subclass\n";
}

1;
