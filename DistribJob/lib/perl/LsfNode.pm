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
    if(!-e "$runFile"){
      open(R,">$runFile") || die "unable to create script file '$runFile'\n";
      print R <<"EOF";
#!/bin/sh
function cleanup {
  find $self->{nodeDir}* -user \$LOGNAME -maxdepth 0 -print0 2>/dev/null | xargs -0r rm -rv  >&2
}
trap cleanup SIGINT SIGTERM EXIT
$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort}
EOF
      close R;
      system("chmod +x $runFile");
    }
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
      $self->{nodeDir} = "$self->{nodeDir}.$jid";
    } # else ??

    $self->setJobid($jid);
    # inject jobids into cancelFile
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

  my $task = $self->getTask();
  $task->cleanUpNode($self) if $task;

  
  if($self->{nodeNum} && $self->getPort()){
    $self->runCmd("/bin/rm -r $self->{nodeDir}", 1);
    $self->runCmd("closeAndExit");
    $self->closePort();
  }else{
    system("bkill $self->{jobid} > /dev/null 2>&1");
  }

  $self->setState($state ? $state : $COMPLETE); ##complete
  
}

1;