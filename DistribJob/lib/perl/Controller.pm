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
 ["restart",      "",  "yes/no: restart a task"],
 );

my @nodes;

sub new {
  my ($class, $propfile, $nodenumlist, $kill, $runTime, $parInit, $fileName, $hostname) = @_;

  my $self = {};
  bless $self;

  my ($inputDir, $masterDir, $nodeDir, $slotsPerNode, $subTaskSize, $taskClass, $nodeClass, $restart) = $self->readPropFile($propfile, \@properties);

  return if ($self->kill($kill, $masterDir));

  if ($restart) {
    die "masterDir $masterDir must exist to restart." unless -e $masterDir;
    $self->resetKill($masterDir);
  } else {
    die "masterDir $masterDir already exists. Delete it or use -restart." 
      if -e $masterDir;
    &runCmd("mkdir -p $masterDir");
  }

  my $taskPath = $taskClass;
  $taskPath =~ s/::/\//g;       # work around perl 'require' weirdness
  require "$taskPath.pm";

  my $nodePath = $nodeClass;
  $nodePath =~ s/::/\//g;       # work around perl 'require' weirdness
  require "$nodePath.pm";

  my $task = $taskClass->new($inputDir, $subTaskSize, $restart, $masterDir);

  print "Initializing server...\n\n";
  $task->initServer($inputDir);

  ##open socket to listen for active nodes...
  my $localPort;
  $hostname = `hostname -s` unless $hostname;
  chomp $hostname;
  
  my $numTries = 0;
  my $sock;
  do {
    die "Unable to create port on server\n" if $numTries++ > 5;
    $localPort = int(rand(3000)) + 5000;
    print "Name of headnode: $hostname, port: $localPort\n\n";
    $sock = new IO::Socket::INET (
                                     LocalHost => $hostname,
                                     LocalPort => $localPort,
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
      my $node = $nodeClass->new($nodenum, $nodeDir, $slotsPerNode, $runTime);
      if (!$node) {               ##failed to initialize so is null..
        print "  Unable to create node $nodenum....skipping\n";
        next;
      }
      push(@nodes,$node);
    }
  }else{
    for(my$a=1;$a<=scalar($nodenumlist);$a++){
      my $node = $nodeClass->new(undef, $nodeDir, $slotsPerNode, $runTime, $fileName, $hostname, $localPort);
      if (!$node) {               ##failed to initialize so is null..
        print "Unable to create new node number $a....skipping\n";
        next;
      }
      push(@nodes,$node);
    }
  }
  ##now start running....
  $self->run($task, $inputDir, $masterDir, $propfile,$nodeClass,$nodeDir,$slotsPerNode, $parInit, $sel, $sock);
} 

sub run {
    my ($self, $task, $inputDir,$masterDir, $propfile,$nodeClass,$nodeDir,$slotsPerNode, $parInit, $sel, $sock) = @_;

    my $running = 1;
    my $kill;
    $parInit = 1 unless $parInit;
    my $complete = 0;
    my $ctLoops = 0;
    my @redoSubtasks;

    do {
	
	print  ($kill ? "!" : ".");

	$kill = $self->checkKill($kill, $masterDir);

        my $ctRunning = 0;
        my $ctInitTask = 0;

        ##now check to see if any new nodes are ready to run...
        $self->getNodeMsgs($sel,$sock);

        foreach my $node (@nodes){

          if($node->getState() < $RUNNINGTASK && $complete){  ##no more tasks...clean up these  nodes...
            $node->cleanUp(1);
            next;
          }
          
          if(!$node->getState()){
            print "Submitting node to scheduler ...\n";
            $node->queueNode();
          }elsif($node->getState() == $READYTOINITTASK){  ##has connection but task on node has not been initialized
            next if($ctInitTask > $parInit);
            $ctInitTask++;
            print "Initializing task on node ".$node->getNum()."...".`date`;
            $node->initializeTask($task,$inputDir);
          }elsif($node->getState() == $INITIALIZINGTASK){  ##still initializing task
            $ctInitTask++;
          }elsif($node->getState() == $RUNNINGTASK){  ##running...
            $ctRunning++;
            foreach my $nodeSlot (@{$node->getSlots()}) {
              if ($nodeSlot->taskComplete() && !$kill) {
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
              if($ctLoops % 10 == 0 && $nodeSlot->isRunning()){
                if($task->getSubtaskTime() && $nodeSlot->getTask()->getRunningTime() > 2 * $task->getSubtaskTime()){
                  if(!$node->runCmd("ls ".$nodeSlot->getDir(),1)){
                    print "ERROR:  Node ".$node->getNum()." is no longer functional\n";
                    $node->setState($FAILEDNODE);
                    ##need to get all subtasks from this node and assign to another...
                    foreach my $ns (@{$node->getSlots()}){
                      push(@redoSubtasks,$ns->getTask());
                    }
                    last;
                  }
                }
              }
            }
          }
        } 
        $running = !$complete || $ctRunning;  ##set to 0 if $complete > 0 and $ctRunning == 0
        $ctLoops++;
        sleep(1);

      } while ($running && $kill != $KILLNOW);

    print "Cleaning up nodes...\n";
    foreach my $node (@nodes) {
	$node->cleanUp(1);
    }

    close($sock);

    $self->reportFailures($masterDir, $propfile);

    print "Cleaning up server...\n";
    $task->cleanUpServer($inputDir,"$masterDir/mainresult"); ##allows user to clean up at end of run if desired

    if(scalar(@redoSubtasks) > 0){
      print "\nNodes running the following subtasks failed\n";
      foreach my $subtask (@redoSubtasks){
        next unless defined $subtask;
        print "subtask_".$subtask->getNum()."\n";
      }
      print "Set restart=yes in the controller.prop file and re run to run these failed subtasks.\n\n";
    }

    print "Done\n";
}

sub getNodeMsgs {
  my($self,$sel,$sock) = @_;
  while($sel->can_read(0)) {
    my $fh = $sock->accept();
    my $s = <$fh>;
    chomp $s;
    my ($jobid,$slot,$status) = split(" ",$s);
    close($fh);
    if($slot =~ /slot_/){ ##subtask has completed in this slot...setState
      $self->{nodes}->{$jobid}->getSlot($slot)->getTask()->setState($status);
    }else{ ##node is ready to run...
      foreach my $n (@nodes){
        if($n->getJobid() eq $jobid){
          $n->setState($READYTORUN);
          $n->setNum($slot) if $slot;
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
	    || $props->getProp('subtasksize') > 100000);

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
	    $props->getProp('nodeclass'), $restart);
}

sub checkKill {
    my ($self, $kill, $masterDir) = @_;

    if (kill != $KILLNOW) {

	if (-e "$masterDir/kill") {
	    print  "\nKilled.  Gracefully terminating\n";
	    $kill = $KILLNOW;
	} elsif (!$kill && -e "$masterDir/killslow") {
	    print  "\nKilled.  Terminating as soon as running subTasks complete\n";
	    $kill = 1;
	}
    }
    return $kill;
}

sub resetKill {
    my ($self, $masterDir) = @_;

    if (-e "$masterDir/killslow") {
	&runCmd("/bin/rm $masterDir/killslow");
    }
    
    if (-e "$masterDir/kill") {
	&runCmd("/bin/rm $masterDir/kill");
    }    
}

sub kill {
    my ($self, $kill, $masterDir) = @_;

    if (-e $masterDir) {
	if ($kill == $KILLNOW) {
	    open(F, ">$masterDir/kill");
	    print F "kill\n";
	    close(F);
	    print  "Killing now: will not wait for running jobs to complete\n";
	} elsif ($kill) {
	    open(F, ">$masterDir/killslow");
	    print F "killslow\n";
	    close(F);
	    print  "Killing slow: will wait for running jobs to complete\n";
	}
    } elsif ($kill) {
	print  "Can't kill: masterDir $masterDir doesn't exist\n";
    }
    return $kill;    
}

sub reportFailures {
    my ($self, $masterDir, $propfile) = @_;

    my @failures = split("\n", &runCmd("ls -l $masterDir/failures"));
    my $count = scalar(@failures) - 1;
    if ($count > 0) {
	print "
Failure: $count subtasks failed
Please look in $masterDir/failures/*/result
After analyzing and correcting failures:
  1. mv $masterDir/failures $masterDir/failures.save
  2. set restart=yes in $propfile
  3. restart the job

";
    }
}

1;
