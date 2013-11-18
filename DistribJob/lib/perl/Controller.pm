package DJob::DistribJob::Controller;

use strict;
use CBIL::Util::PropertySet;
use CBIL::Util::Utils;
use DJob::DistribJob::Node ":states";
use IO::Socket;
use IO::Select;
use File::Basename;

my $KILLNOW = 2;


$| = 1;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["inputdir",     "",  "Directory in which the controller can find the input"],
 ["masterdir",    "",  "The controller's master directory"],
 ["nodedir",      "",  "Directory on the node's disk to use as a working dir"],
 ["slotspernode", "",  "Number of subtasks to run on a node simultaneously"],
 ["subtasksize",  "",  "Number of input elements to include in each subtask"],
 ["taskclass",    "",  "Subclass of DJob::DistribJob::Task that manages the task"],
 ["nodeclass",    "",  "Subclass of DJob::DistribJob::Node that handles the node"],
 ["keepNodeForPostProcessing", "no",  "yes/no: keep a node for cleanup at end"],
 ["useTmpDir", "yes",  "[yes]/no: set the nodeDir to the tmpDir assigned by scheduler"],
 );

my @nodes;

sub new {
  my ($class, $propfile, $nodenumlist, $kill, $runTime, $parInit, $fileName, $hostname, $procsPerNode, $memPerNode, $queue) = @_;

  my $self = {};
  bless $self;

  $self->{runTime} = $runTime;
  $self->{parInit} = $parInit;
  $self->{fileName} = $fileName;
  $self->{hostname} = $hostname;
  $self->{procsPerNode} = $procsPerNode;
  $self->{memPerNode} = $memPerNode;
  $self->{queue} = $queue;
  $self->{propFile} = $propfile;
  $self->{kill} = $kill;

  ($self->{inputDir}, $self->{masterDir}, $self->{nodeDir}, $self->{slotsPerNode}, $self->{subTaskSize}, $self->{taskClass}, $self->{nodeClass}, $self->{keepNodeForPostProcessing}, $self->{useTmpDir}) = $self->readPropFile($propfile, \@properties);

  return if ($self->kill($kill));

  $self->{processIdFile} = "$self->{masterDir}/distribjobProcessId"; 

  my $restart = 0;

  # restart case
  if (-e $self->{masterDir}) {
      
      my $runningProcessId = $self->processAlreadyRunning();
      die "This job is already running in process $runningProcessId (found in file $self->{processIdFile}).  \nExiting.\n" if $runningProcessId;

      # remove existing running/ dir.  whatever was running is not any longer
      &runCmd("rm -r $self->{masterDir}/running") if (-e "$self->{masterDir}/running");
      $restart = 1;

      if (-e "$self->{masterDir}/failures" && scalar(glob("$self->{masterDir}/failures/*"))) {
	  my $restartInstructions = $self->getRestartInstructions();
	  die "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Existing failures found in $self->{masterDir}/failures.

Probably that is because you have run a job, had some failures, and are trying to restart.

If so, then:
$restartInstructions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
";
      } 

  # regular case
  } else {
    &runCmd("mkdir -p $self->{masterDir}");
  }

  # make process id file
  open(F, ">$self->{processIdFile}") || die "Can't open processIdFile '$self->{processIdFile}' for writing\n";
  print F "$$\n";
  close(F);
 
  my $taskPath = $self->{taskClass};
  $taskPath =~ s/::/\//g;       # work around perl 'require' weirdness
  require "$taskPath.pm";

  my $nodePath = $self->{nodeClass};
  $nodePath =~ s/::/\//g;       # work around perl 'require' weirdness
  require "$nodePath.pm";


  ##open socket to listen for active nodes...
  $self->{hostname} = `hostname -s` unless $self->{hostname};
  chomp $self->{hostname};
  
  my $numTries = 0;
  my $sock;
  do {
    die "Unable to create port on server\n" if $numTries++ > 5;
    $self->{localPort} = int(rand(3000)) + 5000;
    $sock = new IO::Socket::INET (
 #                                    LocalHost => $self->{hostname},
                                     LocalHost => '0.0.0.0',
                                     LocalPort => $self->{localPort},
                                     Proto => 'tcp',
                                     Listen => 100,
                                     Reuse => 1,
                                    );
  } until ( $sock );
  my $sel = IO::Select->new($sock);

  my $amRunning = 0;
  my $runpid;
  if(ref($nodenumlist) =~ /ARRAY/){
    foreach my $nodenum (@$nodenumlist) {
      my $node = $self->{nodeClass}->new($nodenum, $self->{nodeDir}, $self->{slotsPerNode}, $self->{runTime}, $self->{fileName}, $self->{hostname}, $self->{localPort}, $self->{procsPerNode}, $self->{memPerNode}, $self->{queue},$self->{masterDir});
      if (!$node) {               ##failed to initialize so is null..
        print "  Unable to create node $nodenum....skipping\n";
        next;
      }
      push(@nodes,$node);
    }
  }else{
    for(my$a=1;$a<=scalar($nodenumlist);$a++){
      my $node = $self->{nodeClass}->new(undef, $self->{nodeDir}, $self->{slotsPerNode}, $self->{runTime}, $self->{fileName}, $self->{hostname}, $self->{localPort}, $self->{procsPerNode}, $self->{memPerNode}, $self->{queue},$self->{masterDir});
      if (!$node) {               ##failed to initialize so is null..
        print "Unable to create new node number $a....skipping\n";
        next;
      }
      push(@nodes,$node);
    }
  }

  ## get an init node before initializing the task so this can be done on a node
  print "Waiting for init node\n";
  my $initNode;
  if(!$nodes[0]->getState()){
    print "Submitting node to scheduler ";
    ##NOTE: need to set the memory of the init node here iif differet than rest.

    $nodes[0]->queueNode();
  }

  my $ctInitNode = 0;
  until($initNode){
    $ctInitNode++;
    if($ctInitNode % 10 == 0){ print "."; }
    die "ERROR Queueing init node ... check error logs\n" if $ctInitNode % 20 == 0 && !$nodes[0]->getQueueState();
    $self->getNodeMsgs($sel,$sock);
    if($nodes[0]->getState() >= $READYTORUN || $nodes[0]->getState() == $FAILEDNODE){
      if($nodes[0]->checkNode()){
        $initNode = $nodes[0];
        print "\n";
      }else{
        my $tmpNode = $self->{nodeClass}->new(undef, $self->{nodeDir}, $self->{slotsPerNode}, $self->{runTime}, $self->{fileName}, $self->{hostname}, $self->{localPort}, $self->{procsPerNode}, $self->{memPerNode}, $self->{queue},$self->{masterDir}); 
        print "\nNew node created to replace failed node (".$nodes[0]->getJobid().")\n";
        $nodes[0] = $tmpNode;
        print "Submitting node to scheduler ";
        $nodes[0]->queueNode();
      }
    }
    sleep 1;
  }

  my $task = $self->{taskClass}->new($self->{inputDir}, $self->{subTaskSize}, $restart, $self->{masterDir},$initNode);

  print "Initializing server...\n\n";
  $task->initServer($self->{inputDir});
  ##add a bit of code here to print an error statement if the size is not properly set.
  if(!$task->{size}){
    die "ERROR: the inputSetSize is not properly set. This must either be computed in the getInputSetSize() method or in the initServer() method.  If  you compute it in the initServer() method, you must set the variable\n\n \$self->{size} = <your computed set size>;\n\nin that same method.\n\n";
  }

  ##now start running....
  $self->run($task, $propfile, $sel, $sock);
} 

sub run {
    my ($self, $task, $propfile, $sel, $sock) = @_;


    my $running = 1;
    my $kill;
    $self->{parInit} = 1 unless $self->{parInit};
    my $ctLoops = 0;

    do {
	
	print  ($kill ? "!" : ".") if $ctLoops % 10 == 0;

	$kill = $self->checkKill($kill);

        my $ctRunning = 0;
        my $ctInitTask = 0;
        my $ctFinished = 0;

        foreach my $node (@nodes){

          $self->getNodeMsgs($sel,$sock);

          if($node->getState() < $RUNNINGTASK && !$task->countRemainingSubtasks()){  ##no more tasks...clean up these  nodes...
            $node->cleanUp(1);
            next;
          }
          
          $ctFinished++ if $node->getState() > $RUNNINGTASK;
          
          if(!$node->getState()){
            print "Submitting node to scheduler ...\n";
            $node->queueNode();
          }elsif($node->getState() == $FAILEDNODE){

            # print "NODE (".$node->getJobid().") FAILED: getting subtasks so can reassign\n";

            $self->{haveCleanupNode} = 0 if $self->{keepNodeForPostProcessing} eq 'yes' && $node->getSaveForCleanup();
            ### remove from the array of @nodes
            $self->removeNode($node->getJobid());
            ### want to get a new node here and add to list of nodes but only do once!!
            my $tmpNode = $self->{nodeClass}->new(undef, $self->{nodeDir}, $self->{slotsPerNode}, $self->{runTime}, $self->{fileName}, $self->{hostname}, $self->{localPort}, $self->{procsPerNode}, $self->{memPerNode}, $self->{queue},$self->{masterDir});
            if (!$tmpNode) {               ##failed to initialize so is null..
              print "Unable to create new node to replace failed node ".$node->getJobid()."\n";
              next;
            }
            print "New node created to replace failed node (".$node->getJobid().")\n";
            push(@nodes,$tmpNode);
            next;
          }elsif($node->getState() == $READYTOINITTASK){  ##has connection but task on node has not been initialized
            next if($ctInitTask >= $self->{parInit});
            $ctInitTask++;
            print "Initializing task on node ".$node->getNum()."...".`date`;
            $node->initializeTask($task,$self->{inputDir});
          }elsif($node->getState() == $INITIALIZINGTASK){  ##still initializing task
            $ctInitTask++;
          }elsif($node->getState() == $RUNNINGTASK){  ##running...
            $ctRunning++;
            ##set node to be saved if not already set to be saved and count < saveforcleanup
            if ($self->{keepNodeForPostProcessing} eq 'yes' && !$self->{haveCleanupNode} && !$node->getSaveForCleanup()){
              $node->setSaveForCleanup(1);
              $self->{haveCleanupNode} = 1;
            }
            foreach my $nodeSlot (@{$node->getSlots()}) {
              if ($nodeSlot->taskComplete() && !$kill) {
                last if $node->getState() == $FAILEDNODE;  ##taskComplete can fail if results can't be integrated so need to stop processing if that happens.
                $nodeSlot->assignNewTask($task->nextSubTask($nodeSlot));
              }
#remove              $complete |= !$nodeSlot->isRunning();  ##sets to non-zero if nodeslot is not running
              ##check to see if the node is still functional 
              if($ctLoops % 20 == 0 && $nodeSlot->isRunning()){
                if(($task->getSubtaskTime() && $nodeSlot->getTask()->getRunningTime() > 2 * $task->getSubtaskTime()) || $ctLoops % 1000 == 0){
                  if(!$node->checkNode()){
                    print "ERROR:  Node ".$node->getNum()." (".$node->getJobid().") is no longer functional\n";
                    last;
                  }
                }
              }
            }
          }
        } 
        $running =  $task->countRemainingSubtasks() || $ctRunning;  ##set to 0 if $complete > 0 and $ctRunning == 0
#        if($running && $ctFinished >= scalar(@nodes)){  ##if running and no nodes not finished then stop
#          $running = 0;
#          print "\nERROR:  No nodes are available but subtasks remain to complete ... exiting\nRestart to run the remaining subtasks.\n\n";
#        }
        $ctLoops++;
        if($ctLoops % 20 == 0){
          $self->manageNodesBasedOnQueueState();
          $self->manageFailedNodes();
        }
        ##new need to check on the nodes in the queue in case have failed silently ...
        sleep(1);

      } while ($running && $kill != $KILLNOW);

    print "Cleaning up nodes...\n";
    print "Average time to complete subtasks = ",(int($task->getSubtaskTime()*10)/10)." seconds\n";
    foreach my $node (@nodes) {
      $node->cleanUp(1) unless $node->getSaveForCleanup();
    }
    $self->manageFailedNodes(1);

    my $failures = $self->reportFailures($propfile);

    my $numRedoRemaining = $task->haveRedoSubtasks();
    $failures += $numRedoRemaining;
    if($numRedoRemaining){  ##there are still tasks from failed nodes that didn't get assigned
      print "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ERROR: the following subtasks were not run\n";
      while(my $st = $task->getNextRedoSubtask()){
        last unless $st;
        print "  subTask ",$st->getNum()," did not get reassigned to a new node\n";
      }
      print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    }

    ##get node that have been saved for cleanup here so can pass to cleanUpServer method
    my $cNode;
    foreach my $node (@nodes){
      if($node->getSaveForCleanup()){
        $node->setState($RUNNINGTASK);
        $cNode = $node;
        last;
      }
    }
    print "Cleaning up server ... ". ($cNode ? "running post-processing steps on node ".$cNode->getNum()."\n" : "\n");
    $task->cleanUpServer($self->{inputDir},"$self->{masterDir}/mainresult",$cNode); ##allows user to clean up at end of run if desired
    
    $cNode->cleanUp(1) if $cNode; ##cleanup this node if have it

    ##delete the script file ...
    print "Cleaning up files on server\n";
    my $delScript = "/bin/rm $nodes[0]->{script} > /dev/null 2>&1";
    system($delScript);
    close($sock);
    ##also delete those error files that are of no use since we capture
    sleep 10;  ##give time for the nodes to all be released
    foreach my $n (@nodes){
      if($n->{script}){
        my $errBase = basename($n->{script});
        my $delCmd = "/bin/rm $errBase.* > /dev/null 2>&1";
#        print "$delCmd\n";
        system("$delCmd"); 
        last;
      }else{
        print "ERROR MSG for basename ... script name = '$n->{script}'\n";
      }
    }
    print "Done\n" unless $failures;
}

sub getNodeMsgs {
  my($self,$sel,$sock) = @_;
  while($sel->can_read(0)) {
    my $fh = $sock->accept();
    if(!$fh){
      print "ERROR: getNodeMsgs: There is no file handle from socket\n";
      next;
    }
    my $s = <$fh>;
    chomp $s;
#    print "getNodeMsgs: ($s)\n"; ## unless $subtask;
    my ($jobid,$slot,$status,$tmpDir) = split(" ",$s);
    close($fh);
    if($slot =~ /slot_/){ ##subtask has completed in this slot...setState
      my $subtask =  $self->{nodes}->{$jobid}->getSlot($slot)->getTask();
      $subtask->setState($status);
    }else{ ##node is ready to run...
      foreach my $n (@nodes){
        if($n->getJobid() eq $jobid){
          $n->setState($READYTORUN);
          $n->setNum($slot) if $slot;
          $n->setDir($tmpDir) if ($tmpDir && $self->{useTmpDir} eq 'yes');
          print "Node $slot nodedir set to ".$n->getDir()."\n";
          $n->setLocalPort($status) if $status;
          $n->initialize();
          $self->{nodes}->{$jobid} = $n;  ##put into  hash for nodes
          last;
        }
      }
    }
  }
}

# return undef if failed
sub parseArgs {

    my ($propfile, @nodelist, $kill);

    return ("-help") if $ARGV[0] eq "-help";

    return () if scalar(@ARGV) < 2;


    $propfile = shift @ARGV;
    if ($ARGV[0] =~ /killslow/) {
	$kill = 1;
    } elsif ($ARGV[0] =~ /kill/) {
	$kill = $KILLNOW;
    } else {
	@nodelist = @ARGV;
	return () if $#nodelist < 0;
    }
    
    return ($propfile, \@nodelist, $kill);
}

sub readPropFile {
    my ($self, $propfile, $propDeclaration) = @_;

    my $props  = CBIL::Util::PropertySet->new($propfile, $propDeclaration);

    die "\nError: property 'nodeDir' in $propfile must be a full path\n"
	unless $props->getProp('nodedir') =~ /^\//;

    die "\nError: property 'slotspernode' in $propfile must be between 1 and 10\n" 
	if ($props->getProp('slotspernode') < 1 
	    || $props->getProp('slotspernode') > 10);
	
    die "\nError: property 'subtasksize' in $propfile must be between 1 and 100000\n" 
	if ($props->getProp('subtasksize') < 1 
	    || $props->getProp('subtasksize') > 10000000);

    die "\nError: property 'taskclass' must be set to the name of a subclass of Task\n"
	unless $props->getProp('taskclass');
	
    die "\nError: property 'noceclass' must be set to the name of a subclass of Node\n"
	unless $props->getProp('nodeclass');
	
    print "Controller Properties\n".$props->toString()."\n";
	
    return ($props->getProp('inputdir'), $props->getProp('masterdir'), 
	    $props->getProp('nodedir'), 
	    $props->getProp('slotspernode'), $props->getProp('subtasksize'), 
	    $props->getProp('taskclass'), 
	    $props->getProp('nodeclass'), $props->getProp('keepNodeForPostProcessing'), 
            $props->getProp('useTmpDir') );
}

sub checkKill {
    my ($self, $kill) = @_;

    if (kill != $KILLNOW) {

	if (-e "$self->{masterDir}/kill") {
	    print  "\nKilled.  Gracefully terminating\n";
	    $kill = $KILLNOW;
	} elsif (!$kill && -e "$self->{masterDir}/killslow") {
	    print  "\nKilled.  Terminating as soon as running subTasks complete\n";
	    $kill = 1;
	}
    }
    return $kill;
}

sub resetKill {
    my ($self) = @_;

    if (-e "$self->{masterDir}/killslow") {
	&runCmd("/bin/rm $self->{masterDir}/killslow");
    }
    
    if (-e "$self->{masterDir}/kill") {
	&runCmd("/bin/rm $self->{masterDir}/kill");
    }    
}

sub kill {
    my ($self, $kill) = @_;

    if (-e $self->{masterDir}) {
	if ($kill == $KILLNOW) {
	    open(F, ">$self->{masterDir}/kill");
	    print F "kill\n";
	    close(F);
	    print  "Killing now: will not wait for running jobs to complete\n";
	} elsif ($kill) {
	    open(F, ">$self->{masterDir}/killslow");
	    print F "killslow\n";
	    close(F);
	    print  "Killing slow: will wait for running jobs to complete\n";
	}
    } elsif ($kill) {
	print  "Can't kill: masterDir $self->{masterDir} doesn't exist\n";
    }
    return $kill;    
}

# return the number of failures
sub reportFailures {
    my ($self, $propfile) = @_;

    my @failures = split("\n", &runCmd("ls -l $self->{masterDir}/failures"));
    my $count = scalar(@failures) - 1;
    if ($count > 0) {
	my $restartInstructions = $self->getRestartInstructions();
	print "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Failure: $count subtasks failed
$restartInstructions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
";
    }
    return $count;
}

sub processAlreadyRunning {
  my $self = shift;
  my $f = "$self->{processIdFile}";
  return 0 unless -e "$f";
  my $pid = `cat $f`;
  chomp $pid;
  my $status = 0;
  if ($pid) {
    system("ps -p $pid > /dev/null");
    $status = $? >> 8;
    return $pid if !$status;
  }
  return 0;
}

sub getRestartInstructions {
    my ($self) = @_;

return "
Please look in $self->{masterDir}/failures/*/result

After analyzing and correcting failures:
  1. mv $self->{masterDir}/failures $self->{masterDir}/failures.save
  2. restart the job

";
}

sub removeNode {
  my($self,$jid) = @_;
  my @tmp;
  foreach my $node (@nodes) {
    push(@tmp,$node)  unless $node->getJobid() eq $jid;
  }
  @nodes = @tmp;

}

sub manageFailedNodes {
  my($self,$force) = @_;
#  return unless $force;
#  print "--- manageFailedNodes($force) ---\n";
  my %tmp;
  my $time = time();
  foreach my $n (@failedNodes){
#    print "  ".$n->[1]->getJobid().": $time -> nodeTime = $n->[0]\n";
    if($n->[0] + 300 < $time || $force){
      print "  Releasing failed node ".$n->[1]->getJobid()."\n";
      $n->[1]->cleanUp(1);
    }else{
      $tmp{$n->[1]->getJobid()} = $n;
    }
  }
  @failedNodes = values(%tmp);
#  print "  Remaining failed nodes: ".scalar(@failedNodes)."\n";
}

##want to exit entirely if more than one node and all are either failednodes or getQueueState == 0
##if don't exit, for any where getQueueState == 0 want to run manageNodeBasedOnQueueState
sub manageNodesBasedOnQueueState {
  my ($self) = @_;
  my @bad;  ##put any here that aren't in queue
  my $ct = 0;
  foreach my $node (@nodes){
    $ct++;
    next if ($node->getState() == $FAILEDNODE || $node->getState() == $COMPLETE);
    push(@bad,$node) unless $node->getQueueState();
  }
  ##want to exit gracefully if number of bad nodes == $ct || $ct == 0
  ##unless scalar(@nodes) == 1 then want to manage nodes on node and allow to requeue.
  if($ct == 0 || (scalar(@bad) == $ct && scalar(@nodes) > 1)){
    ##need to die here
    foreach my $n (@nodes){
      $n->cleanUp(1);
    }
    die "ERROR: there are no queued or running nodes available in QUEUE ... exiting\n  Check all error logs, fix problems before rerunning\n\n";
  }
  foreach my $badNode (@bad){
    $badNode->manageNodeBasedOnQueueState();
  }
}

1;
