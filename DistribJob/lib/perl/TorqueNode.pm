package DJob::DistribJob::TorqueNode;

use DJob::DistribJob::Node ":states";
use Cwd;
use strict;

our @ISA = qw(DJob::DistribJob::Node);

sub queueNode {
  my $self = shift;
  if (!$self->getJobid()) {     ##need to run qsub 
    ##first create the script...
    my $runFile = $self->{fileName};
    if(!$runFile){
      $runFile = "nodeScript.$$";
    }elsif($runFile =~ /cancel/){
      $runFile =~ s/cancel/run/;
    }else{
      $runFile = "$runFile.run";
    }
    $self->{script} = $runFile;
    if(!-e "$runFile"){
      open(R,">$runFile") || die "unable to create script file '$runFile'\n";
      print R "#!/bin/sh\n\n$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort}\n";
      close R;
      system("chmod +x $runFile");
    }
    my $qsubcmd = "qsub -N DistribJob -V -j oe -l nodes=1:ppn=$self->{procsPerNode}".($self->{runTime} ? ",walltime=00:$self->{runTime}:00" : "").($self->{queue} ? " -q $self->{queue}" : "")." $runFile";
    print STDERR "\n$qsubcmd\n\n";
    my $jid = `$qsubcmd`;
    chomp $jid;
    $self->{workingDir} = "/$self->{nodeWorkingDirsHome}/$jid";
    $self->setJobid($jid);
    if($self->{fileName}){
      open(C,">>$self->{fileName}");
      print C "$self->{jobid} ";
      close C;
    }
  }
 
  $self->setState($QUEUED);
}

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    my $getCmd = "qstat -n $self->{jobid}";
    my @stat = `$getCmd`;
    return undef if $?;         ##command failed
    if ($stat[-1] =~ /^\s*(\S+)/) {
      $self->{nodeNum} = $1;
    }
  }
  return $self->{nodeNum};
}

# remove this node's job from the queue
sub removeFromQueue {
  my ($self) = @_;
  my $cmd = "qdel $self->{jobid} > /dev/null 2>&1";
  system($cmd);
}

sub getQueueState {
  my $self = shift;
  return 1 if $self->getState() == $FAILEDNODE || $self->getState() == $COMPLETE;  ##should not be in queue
  my $jobid = $self->getJobid();
  if(!$jobid){
    print STDERR "TorqueNode->getQueueState: unable to checkQueueStatus as can't retrieve JobID\n";
    return 0;
  }
  my $checkCmd = "qstat -n $jobid 2> /dev/null";
  my $res = `$checkCmd`;
  return $? >> 8 ? 0 : 1;  ##returns 0 if error running qstat with this jobid
}

sub deleteLogFilesAndTmpDir {
  my $self = shift;
  my $jobidfordel;
  if($self->{jobid} =~ /^(\d+)/){
    $jobidfordel = $1;
  }else{
    print "Unable to remove log files\n";
    return;
  }
  my $delFile = "DistribJob.o$jobidfordel";
#  print "Unlinking $delFile\n";
  unlink("delFile"); ## || print STDERR "Unable to unlink $delFile\n";
  my @outfiles = glob("DistribJob.*");
  if(scalar(@outfiles) == 0){
    ##remove the script file
    unlink("$self->{script}") || print STDERR "Unable to unlink $self->{script}\n";
  }
}

# static method
sub getQueueSubmitCommand {
  my ($class, $queue, $cmdToSubmit) = @_;

  #return "qsub -V -cwd".$queue ? " -q $queue" : "";
  #return "qsub -V -j oe -l nodes=1:ppn=1".($self->{runTime} ? ",walltime=00:$self->{runTime}:00" : "").($self->{queue} ? " -q $self->{queue}" : "");
  return "echo $cmdToSubmit | qsub -V -j oe -l nodes=1:ppn=1,walltime=350:00:00 -q batch";
}

# static method to extract Job Id from the output of the commmand run by getQueueSubmitCommand()
# used to get job id for distribjob itself
sub getJobIdFromJobInfoString {
  my ($class, $jobInfoString) = @_;

  # output message after qsub - 955273.pbs.scm
  #$jobInfoString =~ /Your job (\S+)/;
  $jobInfoString =~ /(\d+).pbs.scm/;
  return $1;
}

# static method to provide command to run to get status of a job, ie, running or done
# used to get status of distribjob itself
sub getCheckStatusCmd {
  my ($class, $jobId) = @_;

  return "qstat -n $jobId";
}

# static method to provide command to run kill jobs
sub getKillJobCmd {
  my ($class, $jobIds) = @_;

  return "qdel $jobIds";
}

# static method to extract status from the output of the command provided by getCheckStatusCmd()
# return 1 if still running.
sub checkJobStatus {
  my ($class, $statusFileString, $jobId) = @_;

#5728531 0.50085 runLiniacJ i_wei        r     10/03/2014

#[hwang@75-108 ~]$ qstat -n 955309
#pbs.scm: 
#                                                                                  Req'd    Req'd       Elap
#Job ID                  Username    Queue    Jobname          SessID  NDS   TSK   Memory   Time    S   Time
#----------------------- ----------- -------- ---------------- ------ ----- ------ ------ --------- - ---------
#955309.pbs.scm          hwang       batch    DistribJob       541093     1      1    --        --  R  00:00:00
#   n74/19

  print STDERR "Status string '$statusFileString' does not contain expected job ID $jobId" unless  $statusFileString =~ /^\s*$jobId/;

  #my $flag = $statusFileString =~ /^\s*$jobId\s+\S+\s+\S+\s+\S+\s+[r|h|w]/;
  my $flag = $statusFileString =~ /^\s*$jobId\s+\S+\s+\S+\s+\S+\s+.*[R|Q]/;
  print STDERR "Found non-running status for job '$jobId' in status string\n $statusFileString\n" if (!$flag);
  return $flag;
}

1;
