package DJob::DistribJob::SgeNode;

use DJob::DistribJob::Node ":states";
use Cwd;
use strict;
use File::Basename;

our @ISA = qw(DJob::DistribJob::Node);

sub queueNode {
  my $self = shift;
  if (!$self->getJobid()) {     ##need to run qsub 
    ##first create the script...
    my $runFile = "$self->{fileName}";
    if(!$runFile){
      $runFile = "nodeScript.$$";
    }elsif($self->{fileName} =~ /cancel/){
      $runFile =~ s/cancel/run/;
    }else{
      $runFile = "$runFile.run";
    }
    $self->{script} = $runFile;
    if(!-e "$runFile"){
      my $host = `hostname`;
      chomp $host;
      open(R,">$runFile") || die "unable to create script file '$runFile'\n";
      print R <<"EOF";
#!/bin/sh
function cleanup {
  find $self->{nodeWorkingDirsHome}/\${JOB_ID} -user \$LOGNAME -maxdepth 0 -print0 2>/dev/null | xargs -0r rm -rv  >&2
}
trap cleanup SIGINT SIGTERM EXIT

$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort}
EOF
      close R;
      system("chmod +x $runFile");
      if($host =~ /cluster/){
        print STDERR "sleeping to allow nfs to catch up with the script file ... ";
        sleep 60 ;
        print STDERR "done \n";
      }
    }
    my $q = $self->getQueue();

    my $qsubcmd = "qsub -V -cwd -pe DJ $self->{procsPerNode} ". ($q ? "-q $q " : ""). "-l h_vmem=$self->{memPerNode}G $runFile";

    my $tjid = `$qsubcmd`;
    if($tjid =~ /job\s(\d+)/){
      my $jid = $1;
      $self->setJobid($jid);
      print "Node Queued: Jobid = $jid \n";
##NOTE: am changing so that now will use the $TMPDIR for the nodeWorkingDirsHome so that SGE will clean up.
      
      $self->setWorkingDir("$self->{nodeWorkingDirsHome}/$jid");
      if($self->{fileName}){
        open(C,">>$self->{fileName}");
        print C "$self->{jobid} ";
        close C;
      }
    }else{
      print STDERR "\nERROR: unable to determine jobid, scheduler returns '$tjid' ... marking FAILEDNODE\n";
      $self->failNode();
    }
  } 
  $self->setState($QUEUED);
}

sub getQueueState {
  my $self = shift;
  return 1 if $self->getState() == $FAILEDNODE || $self->getState() == $COMPLETE;  ##should not be in queue
  my $jobid = $self->getJobid();
  if(!$jobid){
    print STDERR "SgeNode->getQueueState: unable to checkQueueStatus as can't retrieve JobID\n";
    return 0;
  }
  my $checkCmd = "qstat -j $jobid 2> /dev/null";
  my $res = `$checkCmd`;
  return $? >> 8 ? 0 : 1;  ##returns 0 if error running qstat with this jobid
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



#  my $host = $self->runCmd("hostname");
#  chomp $host;
  
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
  if($state != $FAILEDNODE){  ## if the node has failed don't want to run commands on it ...
  

    
    my $task = $self->getTask();
    $task->cleanUpNode($self) if $task;
    
    
    if($self->{nodeNum} && $self->getState() > $QUEUED && $self->getPort()){
      $self->runCmd("/bin/rm -rf $self->{workingDir}",1);
      $self->runCmd("closeAndExit",1);
      $self->closePort();
    }
  }

  ##now want to get stats and print them:
  if($self->getQueueState()){
    my @stats = `qstat -f -j $self->{jobid}`;
    foreach my $line (@stats){
      if($line =~ /^usage.*?(cpu.*)$/){
        print "  qstat -f -j $self->{jobid}: $1\n";
        last;
      }
    }
    system("qdel $self->{jobid} > /dev/null 2>&1");  
  }else{
    my @stats = `qacct -j $self->{jobid}`;
    foreach my $line (@stats){
      print "qacct -j $self->{jobid}`: $line" if $line =~ /(maxvmem|failed)/i;
    }
  }
  if($self->getState() == $FAILEDNODE){ ##don't want to change if is failed node
    $state = $FAILEDNODE;
  }else{
    $self->setState($state == $FAILEDNODE ? $state : $COMPLETE); ##complete
  }

}


1;
