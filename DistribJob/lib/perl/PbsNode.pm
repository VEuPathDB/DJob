package DJob::DistribJob::PbsNode;

use DJob::DistribJob::Node ":states";
use Cwd;
use strict;

our @ISA = qw(DJob::DistribJob::Node);

################################################################################
##NOTE: sites should set the following variable according to how many slots / node
##      their implementation for PBS has
##  NOTE:  THIS HAS BEEN DISCONTINUED...THE NODE CONSTRUCTOR  NOW TAKES AN ARGUMENT (AS DOES DISTRIBJOB)
##     AND SETS THE PROCSPERNODE.  THE DEFAULT IS 2
## my $pbsSlotsPerNode = 2; 
################################################################################

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
    if(!-e "$runFile"){
      open(R,">$runFile") || die "unable to create script file '$runFile'\n";
      print R "#!/bin/sh\n\n$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort}\n";
      close R;
      system("chmod +x $runFile");
    }
    my $qsubcmd = "qsub -N DistribJob -V -j oe -l nodes=1:ppn=$self->{procsPerNode}".($self->{runTime} ? ",walltime=00:$self->{runTime}:00" : "").($self->{queue} ? " -q $self->{queue}" : "")." $runFile";
##    print STDERR "\n$qsubcmd\n\n";
    my $jid = `$qsubcmd`;
    chomp $jid;
    if($jid =~ /^(\d+)/){ 
      $self->{nodeDir} = "$self->{nodeDir}.$1";
    }else{ 
      $self->{nodeDir} = "$self->{nodeDir}.$jid";
    }
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
    if ($stat[-1] =~ /^\s*(\w*\d+)/) {
      $self->{nodeNum} = $1;
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
  
  print "Cleaning up node $self->{nodeNum}...\n";
  
  if($self->{portCon}){
    $self->runCmd("/bin/rm -r $self->{nodeDir}", 1);
    $self->runCmd("closeAndExit");
    $self->closePort();
  }else{
    system("qdel $self->{jobid} > /dev/null 2>&1");
  }
  
  ##delete those pesky OU files that don't do anything
  my $jid = $self->{jobid};
  if($self->{jobid} =~ /^(\d+)/){
    $jid = $1;
  }
  my $delCmd = "/bin/rm $ENV{HOME}/$jid.*OU > /dev/null 2>&1";
  system($delCmd); 
  $self->setState($state ? $state : $COMPLETE); ##complete
}

1;
