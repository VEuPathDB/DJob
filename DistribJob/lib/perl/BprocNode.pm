package DJob::DistribJob::BprocNode;

use DJob::DistribJob::Node ":states";
use CBIL::Util::Utils;
use Cwd;

@ISA = qw(DJob::DistribJob::Node);

use strict;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotcount, $runTime, $fileName) = @_;
    ##NOTE: need jobid in order to cancel job at end so create an undef one unless have jobid
    die "ERROR: BprocNode must be created with a jobid or undef nodenum\n" if (defined $nodeNum && length($nodeNum) <= 3); 
    my $self = &DJob::DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotcount, $runTime, $fileName);
    $self->setJobid($nodeNum);   ##job ids are always long...
    $self->{cwd} = getcwd();
    return $self;
}

sub setJobid { my ($self,$jid) = @_; $self->{jobid} = $jid; }
sub getJobid { my $self = shift; return $self->{jobid}; }
    

## if don't  have jobid then need to run cmsubmit
## then cmgetnodes to get nodes...
sub initialize {
  my($self) = @_;
#  print STDERR "Initializing node: '$self->{nodeNum}'\n";
  if(!$self->getJobid()){  ##need to run cmsubmit
    my $cmd = "cmsubmit -p 2 -t $self->{runTime} -d $self->{cwd} -y $ENV{GUS_HOME}/bin/liniachold";
    my $cm = `$cmd`;
    if($cm =~ /Job Identification Number:\s+(\S+)/m){
      $self->setJobid($1);
      ##need to write this jobid to the cancel file
      open(C,">>$self->{fileName}");
      print C "$self->{jobid} ";
      close C;
    }else{
      die "ERROR: unable to get jobid for this node\n";
    }
  } 
  if($self->{state} == $QUEUED){
    if(waitpid($self->{initPid},1)){
      $self->setState($READYTOINITTASK);
      $self->runCmgetnodes() unless $self->{nodeNum};
    } 
    return;
  }
  $self->setState($QUEUED);
  if(!defined $self->{nodeNum}){  
    my $pid;
  FORK: {
      if($pid = fork){
        $self->{initPid} = $pid;
      }elsif(defined $pid) {  
        $self->_init();
        exit;
      }elsif($! =~ /No more process/){
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
  my $sleep = int(rand(15));
  $sleep += 15;
  while(1){
    sleep $ct < 4 ? $sleep : $ct < 12 ? 120 : 300;
    last if $self->runCmgetnodes();
    $ct++;
  }
  if(!$self->checkNode()){
    print "Node $self->{nodeNum} is not responding to commands....skipping\n";
    $self->cleanUp(1, $FAILEDNODE);  
    ##send email to admins..
    open(MAIL, "|mail brunkb\@cpbi.upenn.edu -s 'Genomics Liniac Node $self->{nodeNum} not responding'");
    print MAIL "Genomics node $self->{nodeNum} is not responding to bpsh commands\n";
    close MAIL;
    return;
  }
  $self->_initNodeDir(); 
  print "Node $self->{nodeNum} initialized\n";
}

sub runCmgetnodes {
  my $self = shift;
  my $getCmd = "cmgetnodes $self->{jobid}";
#  print STDERR "Running cmgetnodes...$getCmd\n";
  my $n = `$getCmd`;
  return undef if $?; ##command failed
  my @node = split(/\s+/,$n);
  if(scalar(@node) != 2){
    system("canceljob $self->{jobid}");  ##cancel the job..
    print STDERR "ERROR: Bprocnodes must only have a single node num: jobid '$self->{jobid}', nodes '$n'...reinitializiing\n" ;
    system("canceljob $self->{jobid} > /dev/null");  ##cancel the existing job...
    $self->{jobid} = undef;
    $self->{nodeNum} = undef;
    $self->setState(0);
    $self->initialize();
  }
  $self->{nodeNum} = $node[0];
  return 1;
}

sub _initNodeDir {
    my ($self) = @_;


    if ($self->_fileExists($self->{nodeDir})) {
	$self->runCmd("/bin/rm -r $self->{nodeDir}");
    }

    my $try = 0;
    until($self->_fileExists($self->{nodeDir})){
	die "Can't create $self->{nodeDir} on node $self->{nodeNum}" 
	    if ($try++ > 3);
	$self->runCmd("mkdir -p $self->{nodeDir}");
    }
    return 1;
}

sub runCmd {
    my ($self, $cmd) = @_;
    
    return &CBIL::Util::Utils::runCmd("bpsh -n $self->{nodeNum} $cmd");
}

sub execSubTask {
    my ($self, $nodeRunDir, $serverSubtaskDir, $cmd) = @_;
    
    $cmd = "subtaskInvoker $nodeRunDir $serverSubtaskDir/result " . $cmd;

    my $rc = system("bpsh -n $self->{nodeNum} $cmd &");
    my $status = $rc >> 8;
    die "Failed with status $status running '$cmd'\n" if $status;
}

sub _fileExists {
  my($self, $file) = @_;

  for(my $a = 0; $a < 2; $a++){
    my $test = `bpsh -n $self->{nodeNum} find $file 2> /dev/null`;
    return 1 if $test =~ /$file/;
  }
}

sub checkNode {
  my($self) = @_;
  my $node = $self->getNum();
  ##First run diagnose to see if load is too high
#  my @diag = `diagnose -n node$node`;
#  foreach my $line (@diag){
#    if($line =~ /node$node.*DEF\s*(\d+\.\d+)/){
#      if($1 > 3.2){
#       print STDERR $line;
#       print STDERR "Node $self->{nodeNum}  Load: $1\n";
#       return 0;
#     }
#    }
#  }
  ##need to try these things and time out!!
  return 0 unless $self->tryCommand("bpsh -n $node hostname","node$node",1);
  return $self->tryCommand("bpsh -n $node ls $ENV{HOME}","No such",0);
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
        if(/$retVal/){
          $hit = 1;
        }
      }
      alarm(0);
    } else {
      die "bpsh error";
    }
  };
  if($@){
    if ($@ =~ /alarm/) {
      kill('KILL', $pid);
#      print STDERR "KILL: $@\n";
    }else{
      print STDERR "Something bad happened: $@\n";
    }
  }
  return $hit ? $test : !$test; 
}

sub cleanUp {
    my ($self,$force, $state) = @_;

    return if $self->getState() >= $COMPLETE;  ##already cleaned up
    
    if(!$force){
      foreach my $slot (@{$self->getSlots()}){
        return unless $slot->isFinished();
      }
    }

    ##want to kill any child processes still running to quit cleanly
    if($self->getState() == $QUEUED && $self->{initPid}){
      kill(1, $self->{initPid}) unless waitpid($self->{initPid},1);
    }

    print "Cleaning up node $self->{nodeNum}...\n";
   
    $self->runCmd("/bin/rm -r $self->{nodeDir}") if $self->getState() > $QUEUED;
    unlink("$self->{jobid}.stdout");  ##delete that nasty .stdout file
    system("canceljob $self->{jobid} > /dev/null");  ##release the node back into scheduler
    $self->setState($state ? $state : $COMPLETE);  ##complete
}

1;
