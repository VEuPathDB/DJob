package DJob::DistribJob::DistribJobUtils;

use strict;

# static methods to support running distribjob

# $mgr needs these methods:
#  - runSshCmdWithRetries
#  - getNodeClass
#  - error
#  - log
#  - logErr

sub runAndMonitorDistribJob {
  my ($mgr, $user, $submitServer, $transferServer, $jobInfoFile, $logFile, $propFile, $numNodes, $time, $queue, $ppn, $maxMemoryGigs) = @_;

  # if not already started, start it up 
  if (!_readInfoFile($mgr, $jobInfoFile, $user, $transferServer) || !distribJobRunning($jobInfoFile, $user, $submitServer, $transferServer, $mgr->getNodeClass())) {

    # first see if by any chance we are already done (would happen if somehow the flow lost track of the job)
    return 1 if _checkClusterTaskLogForDone($logFile, $user, $transferServer);

    # otherwise, start up a new run

  }

  while (1) {
    sleep(10);
    last if !_distribJobRunning($jobInfoFile,$user, $submitServer, $transferServer, $mgr->getNodeClass());
  }

  sleep(1); # wait for log file (though DJob's new wait should be sufficient)

  return _checkClusterTaskLogForDone($mgr, $logFile, $user, $transferServer);
}

sub submitDistribJobToQueue {
  my ($mgr, $user, $submitServer, $transferServer, $jobInfoFile, $logFile, $propFile, $numNodes, $time, $queue, $ppn, $maxMemoryGigs, $rcfile) = @_;

  my $p = $ppn ? "--ppn $ppn " : "";

  # if rcfile, these will be overridden
  my $gushome = '$GUS_HOME';
  my $bashEnv = "-login";  # read bash env from login user's .bashrc

  if ($rcfile) {
    $gushome = '\$GUS_HOME'; # it needs to be escaped for rcfile
    $bashEnv = "-rcfile $rcfile -i";
  }

  my $distribjobCmd = "$gushome/bin/distribjobSubmit $logFile --numNodes $numNodes --runTime $time --propFile $propFile --parallelInit 4 --mpn $maxMemoryGigs --q $queue $p";

  my $submitCmd = $mgr->getNodeClass()->getQueueSubmitCommand($queue, $distribjobCmd);

  my $cmd = "mkdir -p distribjobRuns; cd distribjobRuns; $submitCmd ";

  my $fullCmd = "/bin/bash $bashEnv -c \"$cmd\"";

  # do the submit on submit server, and capture its output
  my $jobInfo = $mgr->runSshCmdWithRetries(0, $fullCmd, "", 1, 0, $user, $submitServer, "");

  $mgr->error("Did not get jobInfo back from command:\n $fullCmd") unless $jobInfo;

  # now write the output of submit into jobInfoFile on transfer server
  my $writeCmd = "cat > $jobInfoFile";

  open(F, "| ssh -2 $user\@$transferServer '/bin/bash -login -c \"$writeCmd\"'") || $mgr->error("Can't open file handle to write job info to transfer server");
  print F $jobInfo;
  close(F);

  # read it back to confirm it got there safely
  my $jobInfoRead = _readInfoFile($mgr, $jobInfoFile, $user, $transferServer);
  chomp $jobInfoRead;
  chomp $jobInfo;
  $mgr->error("Failed writing job info to jobinfo file on cluster.  (Reading it back didn't duplicate what we tried to write)\n expected: '$jobInfo' got: '$jobInfoRead'") unless $jobInfo eq $jobInfoRead;
}


# return 1 if job still running, else 0.  throw errors if we can't figure it out
sub distribJobRunning {
  my ($mgr, $jobInfoFile, $user, $submitServer, $transferServer, $nodeClass) = @_;

  my $jobSubmittedInfo = _readInfoFile($mgr, $jobInfoFile, $user, $transferServer);

  $mgr->error("Job info file on cluster does not exist or is empty: $jobInfoFile") unless $jobSubmittedInfo;

  my $jobId = $nodeClass->getJobIdFromJobInfoString($jobSubmittedInfo);
  $mgr->error("Can't find job id in job submitted file '$jobInfoFile', which contains '$jobSubmittedInfo'") unless $jobId;

  my $checkStatusCmd = $nodeClass->getCheckStatusCmd($jobId);

  my $jobStatusString = $mgr->runSshCmdWithRetries(0, $checkStatusCmd, undef, 1, 1, $user, $submitServer, "2>&1");

  if ($jobStatusString) {
    my ($running, $msg) = $nodeClass->checkJobStatus($jobStatusString, $jobId);
    $mgr->logErr($msg) unless $running;
    return $running;
  } else {
    $mgr->logErr("Empty job status string returned from command '$checkStatusCmd'\n");
    return 0;
  }
}


##############################################################################################
# private methods
##############################################################################################

sub _readInfoFile {
  my ($mgr, $jobInfoFile, $user, $server) = @_;
  my $cmd = "if [ -a $jobInfoFile ];then cat $jobInfoFile; fi";
  my $jobSubmittedInfo = $mgr->runSshCmdWithRetries(0, $cmd, undef, 0, 1, $user, $server, "");
  return $jobSubmittedInfo;
}

sub _checkClusterTaskLogForDone {
  my ($mgr, $logFile, $user, $transferServer) = @_;

  my $cmd = "/bin/bash -login -c \"if [ -a $logFile ]; then tail -1 $logFile; fi\"";

  my $done = $mgr->runSshCmdWithRetries(0, $cmd, undef, 0, 0, $user, $transferServer, "");

  $mgr->logErr("tail of cluster log file is: '$done'");

  return $done && $done =~ /Done/;
}

1;
