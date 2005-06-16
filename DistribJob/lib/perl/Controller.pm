package DJob::DistribJob::Controller;

use strict;
use CBIL::Util::PropertySet;
use CBIL::Util::Utils;
use DJob::DistribJob::Node ":states";

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
  my ($class, $propfile, $nodenumlist, $kill, $runTime, $parInit, $fileName) = @_;

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
  $taskPath =~ s/::/\//g;       # work around perl 'require' wierdness
  require "$taskPath.pm";

  my $nodePath = $nodeClass;
  $nodePath =~ s/::/\//g;       # work around perl 'require' wierdness
  require "$nodePath.pm";

  my $task = $taskClass->new($inputDir, $subTaskSize, $restart, $masterDir);

  print "Initializing server...\n\n";
  $task->initServer($inputDir);

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
      my $node = $nodeClass->new(undef, $nodeDir, $slotsPerNode, $runTime, $fileName);
      if (!$node) {               ##failed to initialize so is null..
        print "Unable to create new node number $a....skipping\n";
        next;
      }
      push(@nodes,$node);
    }
  }
  ##now start running....
  $self->run($task, $inputDir, $masterDir, $propfile,$nodeClass,$nodeDir,$slotsPerNode, $parInit);
} 

sub run {
    my ($self, $task, $inputDir,$masterDir, $propfile,$nodeClass,$nodeDir,$slotsPerNode, $parInit) = @_;

    my $running = 1;
    my $kill;
    $parInit = 1 unless $parInit;
    my $complete = 0;

    do {
	
	print  ($kill ? "!" : ".");

	$kill = $self->checkKill($kill, $masterDir);

        my $ctRunning = 0;
        my $ctInitTask = 0;

        foreach my $node (@nodes){

          if($node->getState() < 3 && $complete){  ##no more tasks...clean up these  nodes...
            $node->cleanUp(1);
            next;
          }
          
          if(!$node->getState()){  ##does not have connection to node..
            print "Initializing node\n";
            $node->initialize();
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
		$nodeSlot->assignNewTask($task->nextSubTask($nodeSlot));
              }
              $complete |= !$nodeSlot->isRunning();  ##sets to non-zero if nodeslot is not running
            }
          }
        } 
        $running = !$complete || $ctRunning;  ##set to 0 if $complete > 0 and $ctRunning == 0
        sleep(1);

      } while ($running && $kill != $KILLNOW);

    print "Cleaning up nodes...\n";
    foreach my $node (@nodes) {
	$node->cleanUp(1);
    }

    $self->reportFailures($masterDir, $propfile);

    print "Cleaning up server...\n";
    $task->cleanUpServer($inputDir,"$masterDir/mainresult"); ##allows user to clean up at end of run if desired

    print "Done\n";
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
