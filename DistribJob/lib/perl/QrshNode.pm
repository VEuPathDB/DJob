package DJob::DistribJob::QrshNode;

use DJob::DistribJob::Node;
use CBIL::Util::Utils;
use IPC::Open2;
use Symbol;

@ISA = qw(DJob::DistribJob::Node);

use strict;

my $endMatchString = "FooCmdEnd";
my $endCmdString = "echo $?.$endMatchString";

sub new {
    my ($class, $nodeNum, $nodeDir, $slotCount,$initDir) = @_;
    my $self = {};  ##&DJob::DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotCount);
    bless($self, $class);
    $self->{nodeDir} = $nodeDir;
    $self->{slotCount} = $slotCount;
    return $self;
}

sub initialize {
  my ($self,$initDir) = shift;
  return undef if $self->{cmdPid} > 1;  ##already have a valid open connection

  ##test to  see if can get connection...
  my $hn = `qrsh -cwd -pe smp 2 -nostdin -noshell hostname 2> /dev/null`;
  return unless $hn =~ /^node/;
  $self->{cmdPid} = open2($self->{readH}, $self->{writeH}, "qrsh -V -cwd -now no -pe smp 2 bash -s");
  if($self->{cmdPid} =~ /^open2/){
    print "ERROR: $@\n";
    print STDERR "\nOPEN FAILED:\n$self->{cmdPid}\n"; ##open failed...
    close $self->{readH};
    close $self->{writeH};
    return undef;
  }
  ##set file handles to flush immediately
  my $oldfh = select($self->{readH});
  $| = 1;
  select($self->{writeH});
  $| = 1;
  select $oldfh;
  
  if(!$initDir){
    $self->_initNodeDir();  
  }
  $self->setState(2);
  print "Node ".$self->getNum()." is initialized and ready to initialize Task\n";
}

sub _initNodeDir {
    my ($self) = @_;

    return 0 unless ($self->checkNode());

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

##over ride so make specific to QrshNode
sub getNum {
  my $self = shift;
  if(!defined $self->{nodeNum} && $self->getState() >= 2){
    my $hostname = $self->runCmd('hostname');
    if($hostname =~ /^(node\d+)/){
      $self->{nodeNum} = "$1.q";
    }else{
      print STDERR "Unable to determine the Node number: $hostname\n";
      $self->setState(6);
    }
  }
  return $self->{nodeNum};
}

sub getSshNode {
    my $self = shift;
    unless($self->{sshnode}){
	my $tmp = $self->getNum();
	$tmp =~ s/\.q//;
	$self->{sshnode} = $tmp;
    }
    return $self->{sshnode};
}

sub runCmd {
    my ($self, $cmd) = @_;
    my $cmdHandle = $self->{writeH}; 
#    print $cmdHandle "$cmd$endCmdString";
    print $cmdHandle "$cmd\n";
    print $cmdHandle "$endCmdString\n";
    my $resHandle = $self->{readH};
    my $res = "";
    while(<$resHandle>){
	if(/^(\d+)\.$endMatchString/){
	    die "Failed with status $1 running '$cmd'\n" if $1;
	    last;
	}
	$res .= $_;
    }
    return $res;
}

sub execSubTask {
    my ($self, $nodeRunDir, $serverSubtaskDir, $cmd) = @_;
    
    $cmd = "$ENV{GUS_HOME}/bin/subtaskInvoker $nodeRunDir $serverSubtaskDir/result $cmd &";

    my $output = $self->runCmd($cmd);

    return $output;
}

sub _fileExists {
  my($self, $file) = @_;

  for(my $a = 0; $a < 2; $a++){
    my $test = $self->runCmd("find $file 2> /dev/null");
    return 1 if $test =~ /$file/;
  }
}

sub checkNode {
  my($self) = @_;
  return 1; ##don't check unless there are problems
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
      die "qrsh error";
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

# over ride the cleanUp method so do the waitpid thing...
sub cleanUp {
    my ($self,$force) = @_;
    return if $self->getState() == 5;
    ##want to only clean up if both slots are finished
    if(!$force){
      foreach my $slot (@{$self->getSlots()}){
	return unless $slot->isFinished();
      }
    }
    print "Cleaning up node ".$self->getNum()."\n";
    if($self->getState() >= 2){
      $self->runCmd("/bin/rm -r $self->{nodeDir}");
      close $self->{readH};
      close $self->{writeH};
      $self->{clean} = 1;
      waitpid($self->{cmdPid}, 0);
    }
    $self->setState(5);
}


1;
