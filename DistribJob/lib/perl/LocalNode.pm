package DJob::DistribJob::LocalNode;

use DJob::DistribJob::Node ":states";
use Cwd;
use strict;

our @ISA = qw(DJob::DistribJob::Node);

my $jobID = 1000;
my $localNodeNum = 10;
my $endMatchString = 'FooCmdEnd';
my $endCmdString = "echo \$?.$endMatchString";

sub new {
    my ($class, $nodeNum, $nodeWorkingDirsHome, $slotCount, $runTime, $fileName, $serverHost, $serverPort) = @_;
    my $self = &DJob::DistribJob::Node::new($class, $nodeNum, $nodeWorkingDirsHome, $slotCount, $runTime, $fileName, $serverHost, $serverPort);
    $localNodeNum++;
    $self->{workingDir} = "$nodeWorkingDirsHome/node$localNodeNum";
    print STDERR "$localNodeNum: job dir = ".$self->getWorkingDir()."\n";
    return $self;
}

sub queueNode {
  my $self = shift;
  if (!$self->getJobid()) {    ###...not queued 
    $jobID++;
    $self->setJobid($jobID);
    system("$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort} $jobID &");
  } 
  $self->setState($QUEUED);
}
  

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    $self->{nodeNum} = 'localhost';
  }
  return $self->{nodeNum};
}

# remove this node's job from the queue
sub removeFromQueue {
  my ($self) = @_;
  my $cmd = "kill $self->{jobid} > /dev/null 2>&1";
  system($cmd);
}

sub runJobStatusCheck {
  my ($self, $jobid) = @_;
  return 1;
  my $res = `bjobs $jobid 2> /dev/null`;
  return $res =~ /RUN/ || $res =~ /PEND/ ? 1 : 0;
}

### delete the log files here ... since it didn't fail won't need them
### also check to see if any log files remain, if don't, then delete the localtmpdir
sub deleteLogFilesAndTmpDir {
  my $self = shift;
  unlink("$self->{localTmpDir}/djob.$self->{jobid}.out"); ## || print STDERR "Unable to unlink '$self->{localTmpDir}/djob.$self->{jobid}.out'\n";
  unlink("$self->{localTmpDir}/djob.$self->{jobid}.err"); ## || print STDERR "Unable to unlink '$self->{localTmpDir}/djob.$self->{jobid}.err'\n";
  my @outfiles = glob("$self->{localTmpDir}/djob.*.out");
  if(scalar(@outfiles) == 0){
    ##remove the script file
    unlink("$self->{script}");
    system("/bin/rm -r $self->{localTmpDir}");
#    print STDERR "Removed $self->{localTmpDir} containing the script and err files\n";
  }
}

# static method
sub getQueueSubmitCommand {
  my ($class, $queue, $cmdToSubmit) = @_;
  return "$cmdToSubmit";
}

# static method to extract Job Id from job submitted file text
# used to get job id for distribjob itself
# return job id
sub getJobIdFromJobInfoString {
  my ($class, $jobInfoString) = @_;

  # Job <356327> is submitted to default queue <normal>
  $jobInfoString =~ /Job \<(\d+)\>/;

  return $1;
}

# static method to provide command to run to get status of a job
# used to get status of distribjob itself
sub getCheckStatusCmd {
  my ($class, $jobId) = @_;

  return "$jobId";
}

# static method to provide command to run kill jobs
sub getKillJobCmd {
  my ($class, $jobIds) = @_;

  return "kill $jobIds";
}

# static method to extract status from status file
# used to check status of distribjob itself
# return 1 if still running.
sub checkJobStatus {
  my ($class, $statusFileString, $jobId) = @_;

#JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
#282054  brunkb  EXIT  normal     node062.hpc node057.hpc DJob_18464 Oct  3 14:10

  return 1;
  print STDERR "Status string '$statusFileString' does not contain expected job ID $jobId" unless  $statusFileString =~ /$jobId/;

  my $flag = $statusFileString =~ /$jobId\s+\S+\s+(RUN|PEND|WAIT)/;
  my $msg = $flag? "" : "Found non-running status for job '$jobId' in status string\n $statusFileString";
  return ($flag, $msg);
}

1;
