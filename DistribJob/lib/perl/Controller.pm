package DJob::DistribJob::Controller;

use strict;
use CBIL::Util::PropertySet;
use CBIL::Util::Utils;
use DJob::DistribJob::Node ":states";
use IO::Socket;
use IO::Select;

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
 ["restart",      "",  "yes/no: restart a task"]
 );

my @nodes;
my @redoSubtasks;  ##make global variable so can add to it in the getNodeMsgs method

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
  $self->{kill} = $kill;
  my $restart;

  ($self->{inputDir}, $self->{masterDir}, $self->{nodeDir}, $self->{slotsPerNode}, $self->{subTaskSize}, $self->{taskClass}, $self->{nodeClass}, $self->{keepNodeForPostProcessing}, $restart) = $self->readPropFile($propfile, \@properties);

  return if ($self->kill($kill));

  if ($restart) {
    die "masterDir $self->{masterDir} must exist to restart.\n" unless -e $self->{masterDir};
    $self->resetKill();
  } else {
    die "masterDir $self->{masterDir} already exists. Delete it or use -restart.\n" 
      if -e $self->{masterDir};
    &runCmd("mkdir -p $self->{masterDir}");
  }

  my $taskPath = $self->{taskClass};
  $taskPath =~ s/::/\//g;       # work around perl 'require' weirdness
  require "$taskPath.pm";

  my $nodePath = $self->{nodeClass};
  $nodePath =~ s/::/\//g;       # work around perl 'require' weirdness
  require "$nodePath.pm";

  my $task = $self->{taskClass}->new($self->{inputDir}, $self->{subTaskSize}, $restart, $self->{masterDir});

  print "Initializing server...\n\n";
  $task->initServer($self->{inputDir});

  ##open socket to listen for active nodes...
  $self->{hostname} = `hostname -s` unless $self->{hostname};
  chomp $self->{hostname};
  
  my $numTries = 0;
  my $sock;
  do {
    die "Unable to create port on server\n" if $numTries++ > 5;
    $self->{localPort} = int(rand(3000)) + 5000;
    $sock = new IO::Socket::INET (
                                     LocalHost => $self->{hostname},
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
      my $node = $self->{nodeClass}->new($nodenum, $self->{nodeDir}, $self->{slotsPerNode}, $self->{runTime}, $self->{fileName}, $self->{hostname}, $self->{localPort}, $self->{procsPerNode}, $self->{memPerNode}, $self->{queue});
      if (!$node) {               ##failed to initialize so is null..
        print "  Unable to create node $nodenum....skipping\n";
        next;
      }
      push(@nodes,$node);
    }
  }else{
    for(my$a=1;$a<=scalar($nodenumlist);$a++){
      my $node = $self->{nodeClass}->new(undef, $self->{nodeDir}, $self->{slotsPerNode}, $self->{runTime}, $self->{fileName}, $self->{hostname}, $self->{localPort}, $self->{procsPerNode}, $self->{memPerNode}, $self->{queue});
      if (!$node) {               ##failed to initialize so is null..
        print "Unable to create new node number $a....skipping\n";
        next;
      }
      push(@nodes,$node);
    }
  }
  ##now start running....
  $self->run($task, $propfile, $sel, $sock);
} 

sub run {
    my ($self, $task, $propfile, $sel, $sock) = @_;

    my $running = 1;
    my $kill;
    $self->{parInit} = 1 unless $self->{parInit};
    my $complete = 0;
    my $ctLoops = 0;

    do {
	
	print  ($kill ? "!" : ".");

	$kill = $self->checkKill($kill);

        my $ctRunning = 0;
        my $ctInitTask = 0;
        my $ctFinished = 0;

        foreach my $node (@nodes){

          $self->getNodeMsgs($sel,$sock);

          if($node->getState() < $RUNNINGTASK && $complete){  ##no more tasks...clean up these  nodes...
            $node->cleanUp(1);
            next;
          }
          
          $ctFinished++ if $node->getState() > $RUNNINGTASK;
          
          if(!$node->getState()){
            print "Submitting node to scheduler ...\n";
            $node->queueNode();
          }elsif($node->getState() == $FAILEDNODE){

            push(@redoSubtasks,$node->failedSoGetSubtasks());
            ### may make sense to store nodes in a hash rather than array so can remove failed nodes.
            ### alternatively, could create a method to remove the node from the array... probably safer.

            ### cleanup this node
            $node->cleanUp(1);
            $self->{haveCleanupNode} = 0 if $self->{keepNodeForPostProcessing} eq 'yes';
            ### remove from the array of @nodes
            $self->removeNode($node->getJobid());
            ### want to get a new node here and add to list of nodes but only do once!!
            my $tmpNode = $self->{nodeClass}->new(undef, $self->{nodeDir}, $self->{slotsPerNode}, $self->{runTime}, $self->{fileName}, $self->{hostname}, $self->{localPort}, $self->{procsPerNode}, $self->{memPerNode}, $self->{queue});
            if (!$tmpNode) {               ##failed to initialize so is null..
              print "Unable to create new node to replace failed node ".$node->getJobid()."\n";
              next;
            }
            print "New node created to replace failed node (".$node->getJobid().")\n";
            push(@nodes,$tmpNode);
            next;
          }elsif($node->getState() == $READYTOINITTASK){  ##has connection but task on node has not been initialized
            next if($ctInitTask > $self->{parInit});
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
                if(scalar(@redoSubtasks) > 0){
                  my $st = shift @redoSubtasks;
                  print "Reassigning subtask_".$st->getNum()." to node ".$node->getNum().".".$nodeSlot->getNum()."\n";
                  $st->resetStartTime();
                  $st->setNodeSlot($nodeSlot);
                  $nodeSlot->assignNewTask($st);
                  $task->runNextSubtask($st);
                }else{
                  $nodeSlot->assignNewTask($task->nextSubTask($nodeSlot));
                }
              }
              $complete |= !$nodeSlot->isRunning();  ##sets to non-zero if nodeslot is not running
              ##check to see if the node is still functional 
              if($ctLoops % 20 == 0 && $nodeSlot->isRunning()){
                if($task->getSubtaskTime() && $nodeSlot->getTask()->getRunningTime() > 2 * $task->getSubtaskTime()){
                  if(!$node->runCmd("ls ".$nodeSlot->getDir(),1)){
                    print "ERROR:  Node ".$node->getNum()." is no longer functional\n";
                    ##need to get all subtasks from this node and assign to another...
                    push(@redoSubtasks,$node->failedSoGetSubtasks());
                    $node->cleanUp(1,$FAILEDNODE);
                    $self->{haveCleanupNode} = 0 if $self->{keepNodeForPostProcessing} eq 'yes';
                    last;
                  }
                }
              }
            }
          }
        } 
        $running =  !$complete || $ctRunning;  ##set to 0 if $complete > 0 and $ctRunning == 0
        if($running && $ctFinished >= scalar(@nodes)){  ##if running and no nodes not finished then stop
          $running = 0;
          print STDERR "\nERROR:  No nodes are available but subtasks remain to complete ... exiting\nSet restart=yes in the controller.prop file and rerun to run the remaining subtasks.\n\n";
        }
        $ctLoops++;
        sleep(1);

      } while ($running && $kill != $KILLNOW);

    print "Cleaning up nodes...\n";
    foreach my $node (@nodes) {
      $node->cleanUp(1) unless $node->getSaveForCleanup();
    }

    my $failures = $self->reportFailures($propfile);

    ##get node that have been saved for cleanup here so can pass to cleanUpServer method
    my $cNode;
    foreach my $node (@nodes){
      if($node->getSaveForCleanup()){
        $node->setState($RUNNINGTASK);
        $cNode = $node;
        last;
      }
    }
    
    print "Cleaning up server...\n";
    $task->cleanUpServer($self->{inputDir},"$self->{masterDir}/mainresult",$cNode); ##allows user to clean up at end of run if desired
    
    $cNode->cleanUp(1) if $cNode; ##cleanup this node if have it


    if(scalar(@redoSubtasks) > 0){
      print "\nNodes running the following subtasks failed\n";
      foreach my $subtask (@redoSubtasks){
        next unless defined $subtask;
        print "subtask_".$subtask->getNum()."\n";
      }
      print "Set restart=yes in the controller.prop file and re run to run these failed subtasks.\n\n";
    }

    print "Done\n" unless $failures;
    ##delete the script file ...
    my $delScript = "/bin/rm $nodes[0]->{script} > /dev/null 2>&1";
    system($delScript);
    close($sock);
}

sub getNodeMsgs {
  my($self,$sel,$sock) = @_;
  while($sel->can_read(0)) {
    my $fh = $sock->accept();
    if(!$fh){
      print STDERR "ERROR: getNodeMsgs: There is no file handle from socket\n";
      next;
    }
    my $s = <$fh>;
    chomp $s;
    my ($jobid,$slot,$status,$tmpDir) = split(" ",$s);
    close($fh);
    if($slot =~ /slot_/){ ##subtask has completed in this slot...setState
      my $subtask =  $self->{nodes}->{$jobid}->getSlot($slot)->getTask();
      ## having problems with cluster nodes missing perl modules ....
      ## if status is failed and the subtask time is very short then could inactivate this node?
#      if($subtask->getRunningTime() < 10 && $status eq 'failed'){
#        my $node =  $self->{nodes}->{$jobid};
#        print "ERROR:  Node ".$node->getNum()." can not run task ... inactivating.\n";
#        ##need to get all subtasks from this node and assign to another...
#        push(@redoSubtasks,$node->failedSoGetSubtasks());
#        $node->cleanUp(1,$FAILEDNODE);
#      }
      $subtask->setState($status);
    }else{ ##node is ready to run...
      foreach my $n (@nodes){
        if($n->getJobid() eq $jobid){
          $n->setState($READYTORUN);
          $n->setNum($slot) if $slot;
          $n->setDir($tmpDir) if $tmpDir;
          print STDERR "NOTE: nodedir for node $slot set to ".$n->getDir()."\n";
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
	
    die "\nError: property 'restart' in $propfile must be 'yes' or 'no'\n" 
	if ($props->getProp('restart') ne 'yes' 
	    && $props->getProp('restart') ne 'no');

    my $restart = $props->getProp('restart') eq "yes";

    print "Controller Properties\n".$props->toString()."\n";
	
    return ($props->getProp('inputdir'), $props->getProp('masterdir'), 
	    $props->getProp('nodedir'), 
	    $props->getProp('slotspernode'), $props->getProp('subtasksize'), 
	    $props->getProp('taskclass'), 
	    $props->getProp('nodeclass'), $props->getProp('keepNodeForPostProcessing'), $restart);
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
	print "
Failure: $count subtasks failed
Please look in $self->{masterDir}/failures/*/result
After analyzing and correcting failures:
  1. mv $self->{masterDir}/failures $self->{masterDir}/failures.save
  2. set restart=yes in $propfile
  3. restart the job

";
    }
    return $count;
}

sub removeNode {
  my($self,$jid) = @_;
  my @tmp;
  foreach my $node (@nodes) {
    push(@tmp,$node)  unless $node->getJobid() eq $jid;
  }
  @nodes = @tmp;

}

1;
