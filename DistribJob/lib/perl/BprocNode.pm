package DJob::DistribJob::BprocNode;

use DJob::DistribJob::Node;
use CBIL::Util::Utils;

@ISA = qw(DJob::DistribJob::Node);

use strict;

sub new {
    my ($class, $nodeNum, $nodeDir, $slotcount) = @_;
    my $self = &DJob::DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotcount);
    return undef unless $self->checkNode();
    $self->_initNodeDir();  
    $self->{nodeDir} = $nodeDir;
    return $self;
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
  my @diag = `diagnose -n node$node`;
  foreach my $line (@diag){
    if($line =~ /node$node.*DEF\s*(\d+\.\d+)/){
      if($1 > 3.2){
#       print STDERR $line;
       print STDERR "  Load: $1\n";
       return 0;
     }
    }
  }
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

1;
