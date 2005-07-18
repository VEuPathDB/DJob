package DJob::DistribJob::SgeNode;

use DJob::DistribJob::Node ":states";
use CBIL::Util::Utils;
use Cwd;
use IO::Socket;
use strict;

our @ISA = qw(DJob::DistribJob::Node);

################################################################################
##NOTE: sites should set the following variable according to how many slots / node
##      their implementation for PBS has
my $pbsSlotsPerNode = 2;  
################################################################################

my $endMatchString = 'FooCmdEnd';
my $endCmdString = "echo \$?.$endMatchString";

sub queueNode {
  my $self = shift;
  if (!$self->getJobid()) {     ##need to run qsub 
    ##first create the script...
    my $runFile = $self->{fileName};
    $runFile =~ s/cancel/run/;
    if(!-e "$runFile"){
      open(R,">$runFile") || die "unable to create script file '$runFile'\n";
      print R "#!/bin/sh\n\n$ENV{GUS_HOME}/bin/sgeNodeSocketServer.pl $self->{serverHost} $self->{serverPort}\n";
      close R;
      system("chmod +x $runFile");
    }
    my $qsubcmd = "qsub -V -cwd -pe smp 2 $runFile";
    my $tjid = `$qsubcmd`;
    if($tjid =~ /job\s(\d+)/){
      my $jid = $1;
      $self->setJobid($jid);
      open(C,">>$self->{fileName}");
      print C "$self->{jobid} ";
      close C;
    }else{
      die "unable to determine jobid from $tjid\n";
    }
  } 
  $self->setState($QUEUED);
}

sub _init {
  my $self = shift;
  my $ct = 0;
  while (1) {
    last if $self->getNodeAddress();
    sleep $ct < 8 ? 15 : 120;
    $ct++;
  }
  if (!$self->checkNode()) {
    print "Node $self->{nodeNum} is not responding to commands....skipping\n";
    $self->cleanUp(1, $FAILEDNODE);  
    return;
  }
  $self->_initNodeDir(); 
  print "Node $self->{nodeNum} initialized\n";
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

sub _initNodeDir {
  my ($self) = @_;


  if ($self->_fileExists($self->{nodeDir})) {
    $self->runCmd("/bin/rm -r $self->{nodeDir}");
  }

  my $try = 0;
  do {
    die "Can't create $self->{nodeDir} on node $self->{nodeNum}" if ($try++ > 3);
    $self->runCmd("mkdir -p $self->{nodeDir}");
  } until ($self->_fileExists($self->{nodeDir})); 
  return 1;
}

sub runCmd {
  my ($self, $cmd, $ignoreErr) = @_;
  my $sock = $self->getPort();
  print $sock "$cmd\n";
  my $res = "";
  while(<$sock>){
    if(/^(\d+)\.$endMatchString/){
      if($1 && !$ignoreErr){
        print STDERR "Failed with status $1 running $cmd" ;
        print $sock "closeAndExit\n";  ##exits the script on the node..
        close $sock;
        exit(1);
      }
      last;
    }
    $res .= $_;
  }
  return $res;
}

sub getPort {
  my $self = shift;
  if(!$self->{portCon}){
    ##note..need to try a couple of times here because the script may be running but the port is not ready on the nodes to receive connections...seems to take some time
#    print STDERR "Creating new port connection\n";
    my $sock;
    my $ct = 0;
    until($sock){
      $sock = new IO::Socket::INET (
                                    PeerAddr => $self->getNodeAddress(),
                                    PeerPort => '7070',
                                    Proto => 'tcp',
                                   );
      unless($sock){
        die "Could not create socket: $!\n" if $ct++ > 5;
        sleep 3;
        next;
      }
      $self->{portCon} = $sock;
    }
  }
  return $self->{portCon};
}

sub closePort {
  my $self = shift;
  close $self->{portCon};
  undef $self->{portCon};
}

sub execSubTask {
  my ($self, $nodeRunDir, $serverSubtaskDir, $cmd) = @_;
    
  $cmd = "subtaskInvoker $nodeRunDir $serverSubtaskDir/result $self->{jobid} $self->{serverHost} $self->{serverPort} $cmd &";

  return $self->runCmd($cmd);
}

sub _fileExists {
  my($self, $file) = @_;

  for (my $a = 0; $a < 2; $a++) {
    my $test = $self->runCmd("find $file 2> /dev/null", 1);
    return 1 if $test =~ /$file/;
  }
}

sub checkNode {
  my($self) = @_;
  return 1;  ##implement if there are problems..
}

sub tryCommand {
  my($self,$cmd,$retVal,$test) = @_;
  my $pid;
  my $hit = 0;
  eval {
    local $SIG{'ALRM'} = sub { die "alarm\n"; };
    alarm(5);
    $pid = open(BPSH, "$cmd |");
    if (defined($pid)) {
      while (<BPSH>) {
        if (/$retVal/) {
          $hit = 1;
        }
      }
      alarm(0);
    } else {
      die "bpsh error";
    }
  };
  if ($@) {
    if ($@ =~ /alarm/) {
      kill('KILL', $pid);
      #      print STDERR "KILL: $@\n";
    } else {
      print STDERR "Something bad happened: $@\n";
    }
  }
  return $hit ? $test : !$test; 
}

sub cleanUp {
  my ($self,$force, $state) = @_;

  return if $self->getState() >= $COMPLETE; ##already cleaned up
    
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
  system("/bin/rm $ENV{HOME}/$self->{jobid}.*.OU > /dev/null 2>&1"); ##delete that nasty .stdout file
  $self->setState($state ? $state : $COMPLETE); ##complete
}

1;
