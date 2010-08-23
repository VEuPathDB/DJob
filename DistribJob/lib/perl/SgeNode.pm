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
    my $runFile = $self->{fileName};
    if(!$runFile){
      $runFile = "nodeScript.$$";
    }elsif($self->{filename} =~ /cancel/){
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
    my $qsubcmd = "qsub -V -cwd -pe DJ $self->{slotCount} $runFile";
#    print "$qsubcmd\n";
#    my $qsubcmd = "qsub -V -cwd $runFile";
    my $tjid = `$qsubcmd`;
    if($tjid =~ /job\s(\d+)/){
      my $jid = $1;
      $self->setJobid($jid);
##NOTE: am changing so that now will use the $TMPDIR for the nodeDir so that PBS will clean up.
      $self->{nodeDir} = "/tmp/$jid";
      if($self->{fileName}){
        open(C,">>$self->{fileName}");
        print C "$self->{jobid} ";
        close C;
      }
    }else{
      die "unable to determine jobid from $tjid\n";
    }
  } 
  $self->setState($QUEUED);
}

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    my $getCmd = "qstat";
    my @stat = `$getCmd`;
    return undef if $?;         ##command failed
    foreach my $line (@stat){
      if($line =~ /^\s*$self->{jobid}.*(node\d+).q/){
        $self->{nodeNum} = $1;
      }
    }
  }
  return $self->{nodeNum};
}
##over ride this because want to delete those pesky *.OU files
sub cleanUp {
  my ($self,$force, $state) = @_;

  return if $self->getState() >= $COMPLETE; #already cleaned up
    
  if (!$force) {
    foreach my $slot (@{$self->getSlots()}) {
      return unless $slot->isFinished();
    }
  }
  
  ##want to kill any child processes still running to quit cleanly
  if($self->getState() == $INITIALIZINGTASK && $self->{taskPid}){
    kill(1, $self->{taskPid}) unless waitpid($self->{taskPid},1);
  }

  if($state != $FAILEDNODE){  ## if the node has failed don't want to run commands on it ...
  
    print "Cleaning up node $self->{nodeNum}...\n";
    
    my $task = $self->getTask();
    $task->cleanUpNode($self) if $task;
    
    
    if($self->{nodeNum} && $self->getPort()){
      $self->runCmd("/bin/rm -r $self->{nodeDir}", 1);
      $self->runCmd("closeAndExit",1);
      $self->closePort();
      system("qdel $self->{jobid} > /dev/null 2>&1");
    }else{
      system("qdel $self->{jobid} > /dev/null 2>&1");
    }
  }

  $self->setState($state ? $state : $COMPLETE); ##complete
  
  ##delete those pesky files that don't do anything
  my $jid = $self->{jobid};
  if($self->{jobid} =~ /^(\d+)/){
    $jid = $1;
  }
  my $errBase = basename($self->{script});
  my $delCmd = "/bin/rm $errBase.?$jid > /dev/null 2>&1";
  system($delCmd); 
}

1;
