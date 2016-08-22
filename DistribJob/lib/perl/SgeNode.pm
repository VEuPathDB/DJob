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

# remove this node's job from the queue
sub removeFromQueue {
  my ($self) = @_;
  my $cmd = "qdel $self->{jobid} > /dev/null 2>&1";
  system($cmd);
}

# an optional method for subclasses to implement
# called at the end of node->cleanUp
# can query the que to return stats about this run
# print results to stdout
sub reportJobStats {
  my ($self) = @_;
  if($self->getQueueState()){
    my @stats = `qstat -f -j $self->{jobid}`;
    foreach my $line (@stats){
      if($line =~ /^usage.*?(cpu.*)$/){
        print "  qstat -f -j $self->{jobid}: $1\n";
        last;
      }
    }
  }else{
    my @stats = `qacct -j $self->{jobid}`;
    foreach my $line (@stats){
      print "qacct -j $self->{jobid}`: $line" if $line =~ /(maxvmem|failed)/i;
    }
  }

}

sub runJobStatusCheck {
  my ($self, $jobid) = @_;

  my $res = `qstat -j $jobid 2> /dev/null`;
  return $? >> 8 ? 0 : 1;  ##returns 0 if error running qstat with this jobid
}

sub deleteLogFilesAndTmpDir {
  my $self = shift;
  unlink("$self->{script}.e$self->{jobid}") || print STDERR "Unable to unlink $self->{script}.e$self->{jobid}\n";
  unlink("$self->{script}.o$self->{jobid}") || print STDERR "Unable to unlink $self->{script}.o$self->{jobid}\n";
  unlink("$self->{script}.pe$self->{jobid}") || print STDERR "Unable to unlink $self->{script}.pe$self->{jobid}\n";
  unlink("$self->{script}.po$self->{jobid}") || print STDERR "Unable to unlink $self->{script}.po$self->{jobid}\n";
  my @outfiles = glob("$self->{script}.*");
  if(scalar(@outfiles) == 0){
    ##remove the script file
    unlink("$self->{script}") || print STDERR "Unable to unlink $self->{script}\n";
  }
}

# static method
sub getQueueSubmitCommand {
  my ($class, $queue, $cmdToSubmit) = @_;

  return "qsub -V -cwd -q $queue $cmdToSubmit";
}

# static method to extract Job Id from job submitted file text
# used to get job id for distribjob itself
sub getJobIdFromJobInfoString {
  my ($class, $jobInfoString) = @_;

  # Your job 1580354 ("script") has been submitted
  $jobInfoString =~ /Your job (\d+)/;
  return $1;
}

# static method to provide command to run to get status of a job
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

# static method to extract status from status file
# used to check status of distribjob itself
# return 1 if still running.
sub checkJobStatus {
  my ($class, $statusFileString, $jobId) = @_;

#5728531 0.50085 runLiniacJ i_wei        r     10/03/2014

  print STDERR "Status string '$statusFileString' does not contain expected job ID $jobId" unless  $statusFileString =~ /^\s*$jobId/;

  my $flag = $statusFileString =~ /^\s*$jobId\s+\S+\s+\S+\s+\S+\s+[r|h|w]/;
  my $msg = $flag? "" : "Found non-running status for job '$jobId' in status string\n $statusFileString";
  return ($flag, $msg);
}

1;
