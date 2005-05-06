package DJob::DistribJob::PBSNode;

use DJob::DistribJob::Node;
use CBIL::Util::Utils;
use IPC::Open2;
use Symbol;
use lib '/scratch/mheiges/testbed';

@ISA = qw(DJob::DistribJob::Node);

use Expect;
use strict;

my $endMatchString = "FooCmdEnd";
my $endCmdString = "echo $?.$endMatchString";

my $log = $ENV{HOME} . "/distribjob.log";

sub new {
    my ($class, $nodeNum, $nodeDir, $slotCount, $initDir) = @_;
    my $self = {};
    bless($self, $class);
    $self->{expect} = Expect->new();

    $self->{nodeDir} = "$nodeDir/$nodeNum";
    $self->{nodeNum} = $nodeNum;
    $self->{slotCount} = $slotCount;
    ##MARK, change the next line of code to work with PBS.  You should reserve all processors on a single machine
    ##opens channel to node
    
    $self->{expect}->spawn("qsub -I -l nodes=1") or die "Could not start qsub: $!\n";
    unless ($self->{expect}->expect(undef, "mnode")) {
    }
    unless ($e->expect(10, '-re', '\$ $')) {
        print "ERROR";
        die "qsub did not return a prompt in a timely manner. I'm not waiting any longer.\n";
    }
    if(!$initDir){
       return undef unless $self->_initNodeDir();  
    }
  
    #$self->_log(__FILE__ . ":" . __LINE__ . " PBSNode created."); # DEBUG
    $self->runCmd("hostname > /scratch/mheiges/distribjob.test.$nodeNum");# DEBUG
  
  return $self;
}

sub runCmd {
    my ($self, $cmd) = @_;
    my $cmdHandle = $self->{expect}; 
    $cmdHandle->clear_accum();
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


sub _log {
    my $self = shift;
    my $msg = shift;
    open(LOG, ">>$log");
    print LOG $msg . "\n";
    close LOG;
}

##is essentially like local node for you
sub _initNodeDir {
    my($self) = @_;
    if (-e $self->{nodeDir}) {
	$self->runCmd("/bin/rm -r $self->{nodeDir}");
    $self->_log(__FILE__ . ":" . __LINE__ . " removing : " . $self->{nodeDir}); # DEBUG
    }
    $self->_log(__FILE__ . ":" . __LINE__ . " creating: " . $self->{nodeDir}); # DEBUG
    $self->runCmd("mkdir -p $self->{nodeDir}");
}

##over ride so make specific to QrshNode
##MARK, you will need to alter this to return the node number when using PBS
sub getNum {
  my $self = shift;
  if(!defined $self->{nodeNum}){
    my $hostname = $self->runCmd('hostname');
    if($hostname =~ /^(mnode\d+)/){
      $self->{nodeNum} = $1;
    }else{
      print STDERR "Unable to determine the Node number\n";
    }
  }
  return $self->{nodeNum};
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
    my ($self) = @_;
    unless ($self->{clean}) {
      print "Cleaning up node ".$self->getNum()."\n";
      $self->runCmd("/bin/rm -r $self->{nodeDir}");
      close $self->{readH};
      close $self->{writeH};
      $self->{clean} = 1;
      #qx(qdel $self->{qJobID});
      $self->runCmd("killall sleep");
      waitpid($self->{cmdPid}, 0);
    }
}


1;
