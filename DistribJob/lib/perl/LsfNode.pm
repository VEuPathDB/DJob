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
	my $mem = int($self->{memPerNode} * 1024);

    my $bsubcmd = qq^
        bsub \\
        -J DJob_$$ \\
        -o $localtmpdir/djob.%J.out \\
        -e $localtmpdir/djob.%J.err \\
        @{[($q ? " -q $q" : "")]} \\
        @{[($mem ? " -M $mem" : "")]} \\
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

# remove this node's job from the queue
sub removeFromQueue {
  my ($self) = @_;
  my $cmd = "bkill $self->{jobid} > /dev/null 2>&1";
  system($cmd);
}

sub runJobStatusCheck {
  my ($self, $jobid) = @_;

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
    my $q = $queue? "-q $queue" : "";
    return "bsub $q $cmdToSubmit";
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

  return "bjobs $jobId";
}

# static method to provide command to run kill jobs
sub getKillJobCmd {
  my ($class, $jobIds) = @_;

  return "bkill $jobIds";
}

# static method to extract status from status file
# used to check status of distribjob itself
# return 1 if still running.
sub checkJobStatus {
  my ($class, $statusFileString, $jobId) = @_;

#JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
#282054  brunkb  EXIT  normal     node062.hpc node057.hpc DJob_18464 Oct  3 14:10

  print STDERR "Status string '$statusFileString' does not contain expected job ID $jobId" unless  $statusFileString =~ /$jobId/;

  my $flag = $statusFileString =~ /$jobId\s+\S+\s+(RUN|PEND|WAIT)/;
  my $msg = $flag? "" : "Found non-running status for job '$jobId' in status string\n $statusFileString";
  return ($flag, $msg);
}


1;
