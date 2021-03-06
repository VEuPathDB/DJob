package DJob::DistribJob::Task;

use strict;
use CBIL::Util::Utils;
use CBIL::Util::PropertySet;
use DJob::DistribJob::SubTask;

my $haveReadJobscript = 0;
$| = 1;

sub new {
    my ($class, $inputDir, $subTaskSize, $restart, $masterDir, $nodeForInit, $props, $sizeAfterInitServer) = @_;

    my $self = {};
    bless $self, $class;

    $self->{props} = CBIL::Util::PropertySet->new("$inputDir/task.prop", $props);
    print "Task Properties\n".$self->{props}->toString()."\n";

    $self->{subTaskSize} = $subTaskSize;
    $self->{start} = -$self->{subTaskSize};
    $self->{end} = 0;
    $self->{subTaskNum} = 0;
    $self->{nodeForInit} = $nodeForInit;

    if($sizeAfterInitServer){
      print "Will find input set size after initializing server\n";
    }else{
      print "Finding input set size\n";
      $self->{size} = $self->getInputSetSize($inputDir);
      if ($self->{size} == 0) {
	print "Input set size is 0 .... it must either be calculated in the getInputSetSize() method or in the initServer() method\n";
      }else{
        print "Input set size is $self->{size}\n";
        my $c = int($self->{size} / $self->{subTaskSize});
        $c += 1 if $self->{size} % $self->{subTaskSize};
        print "Subtask set size is $self->{subTaskSize} ($c subtasks)\n";
      }
    }

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
  
  if($self->haveRedoSubtasks()){
    $nextSubTask = $self->getNextRedoSubtask();
    my $node = $nodeSlot->getNode();
    print "Reassigning subTask ".$nextSubTask->getNum()." to node ".$node->getNum().".".$nodeSlot->getNum()."(".$node->getJobid().") \n";
    $nextSubTask->resetStartTime();
    $nextSubTask->setNodeSlot($nodeSlot);
    $nodeSlot->assignNewTask($nextSubTask);
    $self->runNextSubtask($nextSubTask);
  }elsif ($self->_subTaskAvailable()) {
    
    my $serverSubTaskDir = $self->_newSubTaskDir($nodeSlot->getNodeNum(),
                                                 $nodeSlot->getNum());
    
    $nextSubTask = DJob::DistribJob::SubTask->new($self->{subTaskNum}, $serverSubTaskDir, $nodeSlot, $self);
    
    $self->runNextSubtask($nextSubTask);
    
  }else{
    $nodeSlot->cleanUp(); ##clean up node here as will release it if using SGE
  }
  return $nextSubTask;
}

sub haveRedoSubtasks {
  my $self = shift;
  return defined($self->{redoSubtasks}) ? scalar(@{$self->{redoSubtasks}}) : 0;
}

sub addRedoSubtasks {
  my($self,@st) = @_;
  foreach my $st (@st){
    push(@{$self->{redoSubtasks}},$st);
  }
}

sub getNextRedoSubtask {
  my($self) = @_;
  my $st = shift @{$self->{redoSubtasks}};
  return $st;  
}

##returns the number of subtasks that haven't been submitted 
sub countRemainingSubtasks {
  my $self = shift;
  my $num = $self->{size} - $self->{end} + (defined $self->{redoSubtasks} ? $self->{subTaskSize} * scalar(@{$self->{redoSubtasks}}) : 0);
##  print "Remaining subtasks=$num\n";
  return $num;
}

sub runCmdOnNode {
  my ($self, $node, $cmd, $ignoreError) = @_;
  my $value = $node->runCmd($cmd, $ignoreError);
  my $err = $node->getErr();
  die "Failed running cmd \n'$cmd'\nError: '$err'\n" if ($err && !$ignoreError);
  return $value;
}

sub runNextSubtask {
  my($self,$nextSubTask) = @_;

  eval {
    my $nodeSlot = $nextSubTask->getNodeSlot();
    my $node = $nodeSlot->getNode();
    my $nodeNum = $node->getNum();
    my $slotNum = $nodeSlot->getNum();
    my $nodeSlotDir = $nodeSlot->getDir();
    my $serverSubTaskDir = $nextSubTask->getDir();
    my $subtaskNumber = $nextSubTask->getNum();
  
    my $date = `date`;
    chomp $date;
    print "\n[$date] node:$nodeNum slot:$slotNum ($node->{jobid})\tsubTask $subtaskNumber dispatching to node\n";
  
    $self->runCmdOnNode($node, "/bin/rm -rf $nodeSlotDir");
    $self->runCmdOnNode($node, "mkdir $nodeSlotDir");
    ##also touch the serverSubtaskDir to refresh nfs mounts to node .. problem on rcluster
    $self->runCmdOnNode($node, "touch $serverSubTaskDir.touch",1);
    $self->runCmdOnNode($node, "/bin/rm $serverSubTaskDir.touch",1);
    $self->initSubTask($nextSubTask->getStart(), $nextSubTask->getEnd(), $node, 
		       $self->{inputDir}, $serverSubTaskDir, $nodeSlotDir, $nextSubTask);
  
    $nextSubTask->setRedoSubtask(0);

    my $cmd = $self->makeSubTaskCommand($node, $self->{inputDir}, $nodeSlotDir,$subtaskNumber,$self->{mainResultDir}, $subtaskNumber);
    ##note: with perl subtaskInvoker there seems to be a problem with quotes
    #  $cmd =~ s/\'/\\\'/g;
    #  $cmd =~ s/\s\'/ \'\\\'/g;
    #  $cmd =~ s/\'\s/\\\'\' /g;
    #  $cmd =~ s/\"/\"\\\"/g;
    print "Task.pm command: $cmd\n" if $subtaskNumber == 1;
    $node->execSubTask($nodeSlotDir, $serverSubTaskDir, $cmd);
  };
  if ($@) {
    print "Failing subtask because of: $@";
    $self->failSubTask($nextSubTask);
  }
}

sub failSubTask {
    my ($self, $subTask) = @_;
    my $subTaskNum = $subTask->getNum();
    my $subTaskDir = $subTask->getDir();
    my $node = $subTask->getNodeSlot()->getNode();
    my $slotNum = $subTask->getNodeSlot()->getNum();
    my $nodeNum = $node->getNum();
    my $date = `date`;
    chomp $date;
    print "\n[$date] node:$nodeNum slot:$slotNum ($node->{jobid})\t subTask $subTaskNum FAILED.\n";
    &runCmd("mv $subTaskDir $self->{failedDir}");
}

sub passSubTask {
    my ($self, $subTask) = @_;
    my $subTaskNum = $subTask->getNum();
    my $subTaskDir = $subTask->getDir();
    my $subTaskResultDir = $subTask->getResultDir();
    my $nodeSlot = $subTask->getNodeSlot();
    my $node = $nodeSlot->getNode();
    my $nodeSlotDir = $nodeSlot->getDir();
    my $nodeNum = $node->getNum();
    my $slotNum = $nodeSlot->getNum();
    
    if(!$node->checkNode()){
#      $self->failSubTask($subTask);  ## don't fail here ... rather the node is bad so will get assigned to another node
#      print "passSubTask (".$node->getJobid()."): Failed node so will reassign\n";
      return undef;
    }

    eval {
      $self->integrateSubTaskResults($subTaskNum, $node, $nodeSlotDir,
				     $self->{mainResultDir});
    };
    ##default is for $integRes to be null if command succeeded ... if returns 1 then want to fail subtask
    if($@){
      print "Method integrateSubTask failed because of:\n $@";
      $self->failSubTask($subTask);
      return;
    }

    my $date = `date`;
    chomp $date;
    print "\n[$date] node:$nodeNum slot:$slotNum ($node->{jobid})\tsubTask $subTaskNum succeeded in ".$subTask->getRunningTime()." seconds\n";

    &runCmd("/bin/rm -rf $subTaskDir");
    &appendToLogFile($self->{completeLog}, "$subTaskNum\n");
}

sub getProperty {
  my($self,$prop) = @_;
  return $self->{props}->getProp($prop);
}

sub setProperty {
  my($self,$prop,$value) = @_;
  $self->{props}->setProp($prop,$value);
}

sub cleanUpServer {
  my($self, $inputDir, $mainResultDir) = @_;
  return 1;
}

sub cleanUpNode {
  my($self,$node) = @_;
  return 1;
}

sub addSubtaskTime {
  my($self,$time) = @_;
  $self->{ctCompleteSubtasks}++;
  $self->{totSubtaskTime} += $time;
} 

sub getSubtaskTime {
  my $self = shift;
  return $self->{ctCompleteSubtasks} ? $self->{totSubtaskTime} / $self->{ctCompleteSubtasks} : 0
}

sub _newSubTaskDir {
    my ($self, $nodeNum, $nodeSlotNum) = @_;

    my $subTaskDir = "$self->{runningDir}/subtask_$self->{subTaskNum}";
    ##if subtaskDir already exists then don't delete and make again
    return $subTaskDir if (-e $subTaskDir);
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
	if (-e $self->{completeLog}) {
	    open(FILE,"$self->{completeLog}") || die "can't open log $self->{completeLog}";
	    while (<FILE>) {
		if (/(\d+)/) {
		    $completed{$1} = 1;
		}
	    }
	    close(FILE);
	}
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
