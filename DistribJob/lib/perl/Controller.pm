#!/usr/local/bin/perl

package DJob::DistribJob::Controller;

use strict;
use CBIL::Util::PropertySet;
use CBIL::Util::Utils;

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

sub new {
    my ($class, $propfile, $nodenumlist, $kill) = @_;

    my $self = {};
    bless $self;

    my ($inputDir, $masterDir, $nodeDir, $slotsPerNode, $subTaskSize, 
	$taskClass, $nodeClass, $restart) = 
	    $self->readPropFile($propfile, \@properties);

    return if ($self->kill($kill, $masterDir));

    if ($restart) {
	die "masterDir $masterDir must exist to restart." unless -e $masterDir;
	$self->resetKill($masterDir);
    } else {
	die "masterDir $masterDir already exists. Delete it or use -restart." 
	    if  -e $masterDir;
	&runCmd("mkdir -p $masterDir");
    }

    my $taskPath = $taskClass;
    $taskPath =~ s/::/\//g;  # work around perl 'require' wierdness
    require "$taskPath.pm";

    my $nodePath = $nodeClass;
    $nodePath =~ s/::/\//g;  # work around perl 'require' wierdness
    require "$nodePath.pm";

    my $task = $taskClass->new($inputDir, $subTaskSize, $restart, $masterDir);

    print "Initializing server...\n\n";
    $task->initServer($inputDir);

    my @nodeSlots;
    foreach my $nodenum (@$nodenumlist) {
	print "Initializing node $nodenum...\n";
	my $node = $nodeClass->new($nodenum, $nodeDir, $slotsPerNode);
        ##First check the node to make certain it is functional...
        if(!$node){  ##failed to initialize so is null..
          print "  Unable to initialize node $nodenum....skipping\n";
          next;
        }
	push(@nodeSlots, @{$node->getSlots()});
	$task->initNode($node, $inputDir);
    }
    $self->run($task, \@nodeSlots, $masterDir, $propfile); 
}

sub run {
    my ($self, $task, $nodeSlots, $masterDir, $propfile) = @_;

    my ($running, $kill);
    
    do {
	
	print  ($kill? "!" : ".");

	$kill = $self->checkKill($kill, $masterDir);

	$running = 0;
    
	foreach my $nodeSlot (@$nodeSlots) {

	    if ($nodeSlot->taskComplete() && !$kill) {
		$nodeSlot->assignNewTask($task->nextSubTask($nodeSlot));
	    }

	    $running |= $nodeSlot->isRunning();
	}
	sleep(1);

    } while ($running && $kill != $KILLNOW);

    print "Cleaning up nodes...\n";
    foreach my $nodeSlot (@$nodeSlots) {
	$nodeSlot->cleanUp();
    }

    $self->reportFailures($masterDir, $propfile);

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
    if ($count) {
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
