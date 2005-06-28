package DJob::DistribJob::PBSNode;

use DJob::DistribJob::Node ":states";
use CBIL::Util::Utils;
use Cwd;
use IO::Socket;
use strict;

our @ISA = qw(DJob::DistribJob::Node);


my $endMatchString = 'FooCmdEnd';
my $endCmdString = "echo \$?.$endMatchString";

sub new {
  my ($class, $nodeNum, $nodeDir, $slotcount, $runTime, $fileName) = @_;
  ##NOTE: need jobid in order to cancel job at end so create an undef one unless have jobid
  my $self = &DJob::DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotcount, $runTime, $fileName);
  $self->{cwd} = getcwd();
  return $self;
}

sub setJobid { my ($self,$jid) = @_; $self->{jobid} = $jid; }
sub getJobid { my $self = shift; return $self->{jobid}; }
sub setNodeDir { my ($self,$nd) = @_; $self->{nodeDir} = $self->{nodeDir} . '/' . $nd; }    

## if don't  have jobid then need to run qsub 
sub initialize {
  my($self) = @_;
  if (!$self->getJobid()) {     ##need to run cmsubmit
    my $qsubcmd = "qsub -V -j oe -l nodes=1:ppn=1".($self->{runTime} ? ",walltime=00:$self->{runTime}:00" : "")." $ENV{GUS_HOME}/bin/nodeSocketServer.pl";
    my $jid = `$qsubcmd`;
    if ($jid =~ s/^(\d+)\..*/$1/) {
      $self->setJobid($1);
      ##need to write this jobid to the cancel file
      open(C,">>$self->{fileName}");
      print C "$self->{jobid} ";
      close C;
    } else {
      die "Unable to  get jobid for this node\n";
    }
  } 
  if ($self->{state} == $QUEUED) {
    if (waitpid($self->{initPid},1)) {
      $self->setState($READYTOINITTASK);
      $self->getNodeAddress() unless $self->{nodeNum};
      $self->setNodeDir($self->{nodeNum});
    } 
    return;
  }
  $self->setState($QUEUED);
  if (!defined $self->{nodeNum}) {  
    my $pid;
  FORK: {
      if ($pid = fork) {
        $self->{initPid} = $pid;
      } elsif (defined $pid) {  
        $self->_init();
        exit;
      } elsif ($! =~ /No more process/) {
        print STDERR "Forking failure: $!\n";
        sleep 1;
        redo FORK;
      } else {
        die "Unable to fork: $!\n";
      }
    } 
  } else {
    $self->setState($READYTOINITTASK);
  }
}

sub _init {
  my $self = shift;
  my $ct = 0;
  while (1) {
    last if $self->getNodeAddress();
    sleep $ct < 6 ? 10 : $ct < 12 ? 60 : 180;
    $ct++;
  }
  if (!$self->checkNode()) {
    print "Node $self->{nodeNum} is not responding to commands....skipping\n";
    $self->cleanUp(1, $FAILEDNODE);  
    return;
  }
  $self->setNodeDir($self->{nodeNum});
  $self->_initNodeDir(); 
  print "Node $self->{nodeNum} initialized\n";
}

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    my $getCmd = "qstat -n $self->{jobid}";
    my @stat = `$getCmd`;
#    print STDERR "getNodeAddress: @stat\n";
    return undef if $?;         ##command failed
    if ($stat[-1] =~ /^\s*(node\d+)/) {
      $self->{nodeNum} = 'm' . $1;
#      print STDERR "getNodeAddress: '$1'\n";
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
  until ($self->_fileExists($self->{nodeDir})) {
    die "Can't create $self->{nodeDir} on node $self->{nodeNum}" if ($try++ > 3);
    $self->runCmd("mkdir -p $self->{nodeDir}");
  }
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
#    print STDERR "Creating new port connection  to $self->{nodeNum}\n";
    my $sock;
    my $ct = 0;
    until($sock){
      $sock = new IO::Socket::INET (
                                    PeerAddr => $self->getNodeAddress(),
                                    PeerPort => '7070',
                                    Proto => 'tcp',
                                   );
      unless($sock){
        die "Could not create socket: $!\n" if $ct > 10;
        sleep 5;
        $ct++;
        next;
      }
      $self->{portCon} = $sock;
    }
  }
  return $self->{portCon};
}

sub closePort {
  my $self = shift;
  close $self->{portCon} if $self->{portCon};
  undef $self->{portCon};
}

sub execSubTask {
  my ($self, $nodeRunDir, $serverSubtaskDir, $cmd) = @_;
    
  $cmd = "subtaskInvoker $nodeRunDir $serverSubtaskDir/result $cmd &";

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
  if ($self->getState() == $QUEUED && $self->{initPid}) {
    kill(1, $self->{initPid}) unless waitpid($self->{initPid},1);
  }

  print "Cleaning up node $self->{nodeNum}...\n";
   
  if($self->{portCon}){
    $self->runCmd("/bin/rm -r $self->{nodeDir}", 1);
    $self->runCmd("closeAndExit");
  }else{
    system("qdel $self->{jobid} > /dev/null 2>&1");
  }
  system("/bin/rm $ENV{HOME}/$self->{jobid}.*.OU > /dev/null 2>&1"); ##delete that nasty .stdout file
#  system("/bin/rm $ENV{HOME}/$self->{jobid}.*.ER > /dev/null 2>&1"); ##delete that nasty .stdout file
  $self->setState($state ? $state : $COMPLETE); ##complete
}

1;
