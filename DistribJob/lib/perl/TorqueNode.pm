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

    my $num = $self->{memPerNode} ? $self->{memPerNode} : 2; # default is 2GB
    $num = int($num) + ($num > int($num));  # ceil since sapelo cannot take decimal 
 
    my $qsubcmd = "qsub -N DistribJob -V -j oe -l nodes=1:ppn=$self->{procsPerNode}".($self->{queue} ? ":$self->{queue}" : "").($self->{runTime} ? ",walltime=00:$self->{runTime}:00" : "").",mem=$num"."gb"." $runFile";

    # old qsub command "-l nodes=2:ppn=48:batch" won't work on Sapelo2, 
    # Sapelo2 task requires queue name and optionally node type, e.g. -q batch -l nodes=2:ppn=48:AMD
    # need to specifiy node type or leave it empty instead, 
    # e.g -q batch -l nodes=2:ppn=48:AMD or -q batch -l nodes=2:ppn=48
    if ($self->{queue} eq 'batch') {  # if it's Sapelo2, no node type. 
      $qsubcmd = "qsub -N DistribJob -V -j oe -l nodes=1:ppn=$self->{procsPerNode}:AMD".($self->{runTime} ? ",walltime=00:$self->{runTime}:00" : "").",mem=$num"."gb"." $runFile";
    }

#    print STDERR "\n$qsubcmd\n\n";
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

  # old qsub command "-l nodes=2:ppn=48:batch" won't work on Sapelo2, 
  # Sapelo2 task requires queue name and optionally node type, e.g. -q batch -l nodes=2:ppn=48:AMD
  # need to specifiy node type or leave it empty instead, 
  # e.g -q batch -l nodes=2:ppn=48:AMD or -q batch -l nodes=2:ppn=48
  if ($queue eq 'batch') {  # if it's Sapelo2, no node type. 
    
    return "echo $cmdToSubmit | qsub -V -j oe -l nodes=1:ppn=1,walltime=480:00:00";
  }

  return "echo $cmdToSubmit | qsub -V -j oe -l nodes=1:ppn=1".($queue ? ":$queue" : ""). ",walltime=480:00:00";
}

# static method to extract Job Id from the output of the commmand run by getQueueSubmitCommand()
# used to get job id for distribjob itself
sub getJobIdFromJobInfoString {
  my ($class, $jobInfoString) = @_;

  # qsub output - 955273.pbs.scm
  #$jobInfoString =~ /(\d+).pbs.scm/;
  # sapelo2 distribjobJobInfo.txt - 377002.sapelo2
  $jobInfoString =~ /(\d+).sapelo2/;
  return $1;
}

# static method to provide command to run to get status of a job, ie, running or done
# used to get status of distribjob itself
sub getCheckStatusCmd {
  my ($class, $jobId) = @_;

  return "qstat | grep $jobId";
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

#output of getCheckStatusCmd()
#957014.pbs                 STDIN            hwang           00:00:00 R batch

  print STDERR "Status string '$statusFileString' does not contain expected job ID $jobId" unless  $statusFileString =~ /^\s*$jobId/;

  my $flag = $statusFileString =~ /^\s*$jobId.*[R|Q]/;
  my $msg = $flag? "" : "Found non-running status for job '$jobId' in status string\n $statusFileString";
  return ($flag, $msg);
}

1;
