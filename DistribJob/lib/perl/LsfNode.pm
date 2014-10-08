package DJob::DistribJob::LsfNode;

use DJob::DistribJob::Node ":states";
use File::Basename;
use Cwd;
use strict;
use constant DEBUG => 0;

our @ISA = qw(DJob::DistribJob::Node);

sub queueNode {
  my $self = shift;
  if (!$self->getJobid()) {     ##need to run bsub 
    my $localtmpdir;
    ##first create the script...
    my $runFile = $self->{fileName};
    if(!$runFile){
      $localtmpdir = "$ENV{HOME}/djob_tmp_$$";
      mkdir $localtmpdir unless -e $localtmpdir;
      $runFile = "$localtmpdir/djob_nodeScript.$$";
    }elsif($runFile =~ /cancel/){
      $localtmpdir = dirname($runFile);
      $runFile =~ s/cancel/run/;
    }else{
      $runFile = "$runFile.run";
    }
    $self->{localTmpDir} = $localtmpdir;
    if(!-e "$runFile"){
      open(R,">$runFile") || die "unable to create script file '$runFile'\n";
      print R <<"EOF";
#!/bin/sh
function cleanup {
  find $self->{nodeWorkingDirsHome}/\${LSB_JOBID} -user \$LOGNAME -maxdepth 0 -print0 2>/dev/null | xargs -0r rm -rv  >&2
}
trap cleanup SIGINT SIGTERM EXIT

$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort}
EOF
      close R;
      system("chmod +x $runFile");
    }
    $self->{script} = $runFile;

    my $q = $self->getQueue();

    my $bsubcmd = qq^
        bsub \\
        -J DJob_$$ \\
        -o $localtmpdir/djob.%J.out \\
        -e $localtmpdir/djob.%J.err \\
        @{[($q ? " -q $q" : "")]} \\
        @{[($self->{runTime} ? " -W $self->{runTime}" : "")]} \\
        $runFile
    ^;
    DEBUG && warn "DEBUG: bsubcmd: \n$bsubcmd\n\n";
    chomp(my $jid = `$bsubcmd`);
    ($jid) = $jid =~ /^Job <(\d+)>/;
    DEBUG && warn "DEBUG: jobid $jid";    

    if($jid =~ /^\d+$/) { 
      $self->{workingDir} = "$self->{nodeWorkingDirsHome}/$jid";

      $self->setJobid($jid);
      # inject jobids into cancelFile
      if($self->{fileName}){
        open(C,">>$self->{fileName}");
        print C "$self->{jobid} ";
        close C;
      }
    }else{ # else ??
      $self->setState($FAILEDNODE);
      return;
    }	
    
  } 
  $self->setState($QUEUED);
}

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    my $getCmd = "bjobs $self->{jobid} | tail -n1 | awk '{print $6}'";
    my @stat = `$getCmd`;
    return undef if $?;         ##command failed
    $self->{nodeNum} = $stat[0];
  }
  DEBUG && warn "DEBUG: nodeNum $self->{nodeNum}\n";
  return $self->{nodeNum};
}

##over ride this because want to delete those pesky *.OU files
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

  my $task = $self->getTask();
  $task->cleanUpNode($self) if $task;

  
  if($self->{nodeNum} && $self->getState() > $QUEUED && $self->getPort()){
    $self->runCmd("closeAndExit");
    $self->closePort();
  }else{
    system("bkill $self->{jobid} > /dev/null 2>&1");
  }

  if($self->getState() == $FAILEDNODE){ ##don't want to change if is failed node
    $state = $FAILEDNODE;
  }else{
    $self->setState($state ? $state : $COMPLETE); ##complete

  }
}

sub getQueueState {
  my $self = shift;
  return 1 if $self->getState() == $FAILEDNODE || $self->getState() == $COMPLETE;  ##should not be in queue
  my $jobid = $self->getJobid();
  if(!$jobid){
    print STDERR "LsfNode->getQueueState: unable to checkQueueStatus as can't retrieve JobID\n";
    return 0;
  }
  my $checkCmd = "bjobs $jobid 2> /dev/null";
  my $res = `$checkCmd`;
  return $res =~ /RUN/ || $res =~ /PEND/ ? 1 : 0;
  return $? >> 8 ? 0 : 1;  ##returns 0 if error running bjobs with this jobid
}

### delete the log files here ... since it didn't fail won't need them
### also check to see if any log files remain, if don't, then delete the localtmpdir
sub deleteLogFilesAndTmpDir {
  my $self = shift;
  unlink("$self->{localTmpDir}/djob.$self->{jobid}.out") || print STDERR "Unable to unlink '$self->{localTmpDir}/djob.$self->{jobid}.out'\n";
  unlink("$self->{localTmpDir}/djob.$self->{jobid}.err") || print STDERR "Unable to unlink '$self->{localTmpDir}/djob.$self->{jobid}.err'\n";
  my @outfiles = glob("$self->{localTmpDir}/djob.*.out");
  if(scalar(@outfiles) == 0){
    ##remove the script file
    unlink("$self->{script}");
    system("/bin/rm -r $self->{localTmpDir}");
#    print STDERR "Removed $self->{localTmpDir} containing the script and err files\n";
  }
}

# static method
sub getInteractiveShellCommand {
  my ($class, $queue) = @_;
  return "bsub -Is"
}

# static method to extract Job Id from job submitted file text
# used to get job id for distribjob itself
sub getJobIdFromJobSubmittedFile {
  my ($class, $jobInfoString) = @_;

  # Your job 1580354 ("script") has been submitted
  $jobInfoString =~ /Your job (\d+)/;
  return $1;
}

# static method to provide command to run to get status of a job
# used to get status of distribjob itself
sub getCheckStatusCmd {
  my ($class, $jobId) = @_;

  return "bjobs $jobId";
}

# static method to extract status from status file
# used to check status of distribjob itself
# return 1 if still running.
sub checkJobStatus {
  my ($class, $statusFileString, $jobId) = @_;

#JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
#282054  brunkb  EXIT  normal     node062.hpc node057.hpc DJob_18464 Oct  3 14:10

  return $statusFileString =~ /$jobId\s+\S+\s+[RUN|PEND|WAIT]/;
}


1;
