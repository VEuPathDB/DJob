package DJob::DistribJob::PBSNode;

use DJob::DistribJob::Node;
use CBIL::Util::Utils;
use IPC::Open2;
use Symbol;

@ISA = qw(DJob::DistribJob::Node);

use strict;

my $endMatchString = "FooCmdEnd";
my $endCmdString = "echo $?.$endMatchString";
my $log = $ENV{HOME} . "/distribjob.log";

sub new {
    my ($class, $nodeNum, $nodeDir, $slotCount, $initDir) = @_;
    my $self = {};  ##&DJob::DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotCount);
    bless($self, $class);
    $self->{nodeDir} = "$nodeDir/$nodeNum";
    $self->{nodeNum} = $nodeNum;
    $self->{slotCount} = $slotCount;
    my $callstr = "ssh -qx $self->{nodeNum} bash --login -s";

    $self->{cmdPid} = open2($self->{readH}, $self->{writeH}, $callstr);  

    return undef if $self->{cmdPid} =~ /^open2/;
    ##set file handles to flush immediately
    my $oldfh = select($self->{readH});
    $| = 1;
    select($self->{writeH});
    $| = 1;
    select $oldfh;

    return undef unless $self->checkNode();

    $self->_initNodeDir(); 
    
    return $self;
}

sub _log {
    my $self = shift;
    my $msg = shift;
    open(LOG, ">>$log");
    print LOG $msg . "\n";
    close LOG;
}

## following the setup of a local node 
sub _initNodeDir {
    my($self) = @_;

    if (-e $self->{nodeDir}) {
	$self->runCmd("/bin/rm -rf $self->{nodeDir}");
    }
    $self->runCmd("mkdir -p $self->{nodeDir}");
}

## return the node number when using PBS
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
  # a sleep process should be running through the PBS queue
  return $self->runCmd( "ps  -C sleep --no-heading");
}

# over ride the cleanUp method so do the waitpid thing...
sub cleanUp {
    my ($self) = @_;
    unless ($self->{clean}) {
      if (-e $self->{nodeDir}) {
          print "removing " . $self->{nodeDir} . "\n";
          $self->runCmd("/bin/rm -rf $self->{nodeDir}");
      }
      close $self->{readH};
      close $self->{writeH};
      $self->{clean} = 1;
      #qx(qdel $self->{qJobID});
      $self->runCmd("killall sleep");
      waitpid($self->{cmdPid}, 0);
    }
}

sub DESTROY {
    &cleanup();
}
1;
