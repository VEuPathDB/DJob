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
  find $self->{nodeDir}/\${LSB_JOBID} -user \$LOGNAME -maxdepth 0 -print0 2>/dev/null | xargs -0r rm -rv  >&2
}
trap cleanup SIGINT SIGTERM EXIT

$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort}
EOF
      close R;
      system("chmod +x $runFile");
    }
    $self->{script} = $runFile;
    my $bsubcmd = qq^
        bsub \\
        -J DJob_$$ \\
        -o $localtmpdir/djob.%J.out \\
        -e $localtmpdir/djob.%J.err \\
        @{[($self->{queue} ? " -q $self->{queue}" : "")]} \\
        @{[($self->{runTime} ? " -W $self->{runTime}" : "")]} \\
        $runFile
    ^;
    DEBUG && warn "DEBUG: bsubcmd: \n$bsubcmd\n\n";
    chomp(my $jid = `$bsubcmd`);
    ($jid) = $jid =~ /^Job <(\d+)>/;
    DEBUG && warn "DEBUG: jobid $jid";    

    if($jid =~ /^\d+$/) { 
      $self->{nodeDir} = "$self->{nodeDir}/$jid";

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
#    print STDERR "Cleaning up the node $self->{nodeNum} ... running: /bin/rm -r $self->{nodeDir}\n";
##note that the nodeSocketServer is now deleting the nodeDir
#    $self->runCmd("/bin/rm -r $self->{nodeDir}", 1);  
    ##delete those pesky log files  ... not yet written so comment out ...
#    my $delCmd = "/bin/rm -r $self->{localTmpDir}/djob.".$self->getJobid().".*";
#    print STDERR "deleting log files: '$delCmd'\n";
#    system($delCmd);
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

1;
