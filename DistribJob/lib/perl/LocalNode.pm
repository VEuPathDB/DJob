package DJob::DistribJob::LocalNode;

use DJob::DistribJob::Node ":states";
use Cwd;
use strict;

our @ISA = qw(DJob::DistribJob::Node);

my $jobID = 1000;
my $endMatchString = 'FooCmdEnd';
my $endCmdString = "echo \$?.$endMatchString";

sub new {
    my ($class, $nodeNum, $nodeDir, $slotCount, $runTime, $fileName, $serverHost, $serverPort) = @_;
    my $self = &DJob::DistribJob::Node::new($class, $nodeNum, $nodeDir, $slotCount, $runTime, $fileName, $serverHost, $serverPort);
    $self->{nodeDir} = "$nodeDir/$nodeNum";
    return $self;
}

sub queueNode {
  my $self = shift;
  if (!$self->getJobid()) {    ###...not queued 
    $jobID++;
    $self->setJobid($jobID);
    system("$ENV{GUS_HOME}/bin/nodeSocketServer.pl $self->{serverHost} $self->{serverPort} $jobID &");
  } 
  $self->setState($QUEUED);
}
  

sub getNodeAddress {
  my $self = shift;
  if (!defined $self->{nodeNum}) {
    $self->{nodeNum} = 'localhost';
  }
  return $self->{nodeNum};
}

1;
