package DistribJob::Task;

use strict;
use Common::Utils;
use Common::PropertySet;
use DistribJob::SubTask;

my $haveReadJobscript = 0;

sub new {
    my ($class, $inputDir, $subTaskSize, $restart, $masterDir,
	$props) = @_;

    my $self = {};
    bless $self, $class;

    $self->{props} = Common::PropertySet->new("$inputDir/task.prop", $props);

    print "Finding input set size\n";
    $self->{size} = $self->getInputSetSize($inputDir);
    if ($self->{size} == 0) {
	print "Error: Input set size is 0\n";
	exit 1;
    }
    print "Input set size is $self->{size}\n";
    $self->{subTaskSize} = $subTaskSize;
    my $c = int($self->{size} / $self->{subTaskSize});
    $c += 1 if $self->{size} % $self->{subTaskSize};
    print "Subtask set size is $self->{subTaskSize} ($c subtasks)\n";
    $self->{start} = -$self->{subTaskSize};
    $self->{end} = 0;
    $self->{subTaskNum} = 0;

    $self->{inputDir} = $inputDir;
    $self->{masterDir} = $masterDir;
    $self->{runningDir} = "$masterDir/running";
    $self->{failedDir} = "$masterDir/failures";
    $self->{completeLog} = "$masterDir/completedSubtasks.log";
    $self->{mainResultDir} = "$masterDir/mainresult";

    $self->_initMasterDir();
    $self->{completedInPastRun} = $self->_findCompletedFromLog($restart);
    
    return $self;
}

sub nextSubTask {
    my ($self, $nodeSlot) = @_;

    my $nextSubTask;
    my $node = $nodeSlot->getNode();

    if ($self->_subTaskAvailable()) {

	my $subTaskDir = $self->_newSubTaskDir($nodeSlot->getNodeNum(),
					       $nodeSlot->getNum());

	$nextSubTask = DistribJob::SubTask->new($self->{subTaskNum}, $subTaskDir, $self);

	my $nodeNum = $nodeSlot->getNodeNum();
	my $slotNum = $nodeSlot->getNum();
	my $nodeDir = $nodeSlot->getDir();

	my $date = `date`;
	chomp $date;
	print "\n[$date] subTask $self->{subTaskNum} dispatching to node $nodeNum.$slotNum\n";

	$node->runCmd("/bin/rm -rf $nodeDir/*");

	$self->initSubTask($self->{start}, $self->{end}, $node, 
			   $self->{inputDir}, $subTaskDir, $nodeDir);

	$self->runSubTask($node, $self->{inputDir}, $subTaskDir, $nodeDir);
    }
    return $nextSubTask;
}

sub failSubTask {
    my ($self, $subTask) = @_;
    my $subTaskNum = $subTask->getNum();
    my $subTaskDir = $subTask->getDir();

    my $date = `date`;
    chomp $date;
    print "\n[$date] subTask $subTaskNum failed\n";
    &runCmd("mv $subTaskDir $self->{failedDir}");
}

sub passSubTask {
    my ($self, $subTask) = @_;
    my $subTaskNum = $subTask->getNum();
    my $subTaskDir = $subTask->getDir();
    my $subTaskResultDir = $subTask->getResultDir();

    my $date = `date`;
    chomp $date;
    print "\n[$date] subTask $subTaskNum succeeded\n";

    $self->integrateSubTaskResults($subTaskNum, $subTaskResultDir, 
				   $self->{mainResultDir});
    &runCmd("/bin/rm -rf $subTaskDir");
    &log($self->{completeLog}, "$subTaskNum\n");
}

sub _newSubTaskDir {
    my ($self, $nodeNum, $nodeSlotNum) = @_;

    my $subTaskDir = "$self->{runningDir}/subtask_$self->{subTaskNum}";
    &runCmd("/bin/rm -rf $subTaskDir");
    &runCmd("mkdir -p $subTaskDir/result");

    open(FILE, ">$subTaskDir/subTaskInfo.txt");
    print FILE "Node: $nodeNum  Slot: $nodeSlotNum   Start: $self->{start}    End: $self->{end}";
    close(FILE);
    return $subTaskDir;
}

sub _findCompletedFromLog {
    my ($self, $restart) = @_;
    my %completed;
    if ($restart) {
	print "Reading list of previously completed subtasks from $self->{completeLog}\n";
	open(FILE,"$self->{completeLog}") || die "can't open log $self->{completeLog}";
	while (<FILE>) {
	    if (/(\d+)/) {
		$completed{$1} = 1;
	    }
	}
	close(FILE);
    }
    return \%completed;
}

sub _subTaskAvailable {
    my ($self) = @_;

    my $found = 0;
    while ($self->{end} < $self->{size}) { 
	$self->{start} += $self->{subTaskSize};
	$self->{end} += $self->{subTaskSize};
	$self->{end} = $self->{size} if $self->{end} > $self->{size};
	$self->{subTaskNum}++;
	if (-e "$self->{failedDir}/subtask_$self->{subTaskNum}") {
	    print "\nsubTask $self->{subTaskNum} is an uncorrected failure.\n";
	} elsif ($self->{completedInPastRun}->{$self->{subTaskNum}}) {
	    print "\nsubTask $self->{subTaskNum} succeeded in past run.\n";
	} else {
	    $found = 1;
	    last;
	}
    }
    return $found;
}

sub _initMasterDir {
    my ($self) = @_;

    if (-e $self->{runningDir}) {
	&runCmd("/bin/rm -r $self->{runningDir}");
    }
    &runCmd("mkdir -p $self->{runningDir}");
    &runCmd("mkdir -p $self->{failedDir}");
    &runCmd("mkdir -p $self->{mainResultDir}");
}

1;
