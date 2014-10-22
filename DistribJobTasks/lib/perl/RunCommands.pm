package DJob::DistribJobTasks::RunCommands;

## simple task to take in a file where each line is a command that is to be run
## Can specify list of files that need to be copied to the node initially
## and the name of an output file that the command produces.  If no output file
## is specified then the assumption is that output will go to stdout.

## Brian Brunk (Penn Bioinformatics Core)

use DJob::DistribJob::Task;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["inputFile",    "",     "full path to the input file with commands one / line"],
 ["outputFile",   "none",   "eg, myRun.out"],
 ["debug",  "false", ""],
 );
my @commands;

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    my $debug = $self->getProperty("debug");
    $self->{'debug'} = 1 if $debug =~ /true/i;
    return $self;
}

# called once
sub initServer {
    my ($self, $inputDir) = @_;
    my $ct = 0;
    open(F, $self->getProperty("inputFile"));
    while(<F>){
      next if /^\s*\#/;
      $ct++;
      push(@commands,$_);
    }
    $self->{commands} =  \@commands;
    $self->{size} = scalar(@commands);
    print "Running ".$self->getInputSetSize()." commands from input file\n\n";
    # do nothing
}

sub initNode {
    my ($self, $node, $inputDir) = @_;
    return 1;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;
    return $self->{size};
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir,$subTask) = @_;
    open(O,">$subTaskDir/run.sh");
    my $num = $start + 1;
    my @cmds = @commands[$start..$end - 1];
    print O "@cmds\n";
    close O; 
    $node->runCmd("cp $subTaskDir/run.sh $nodeSlotDir/");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $cmd = "/bin/sh run.sh";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    my $out = $self->getProperty('outputFile') eq 'none' ? "subtask.output" : $self->getProperty('outputFile'); 
    my $cmd = "cat $nodeExecDir/$out >> $mainResultDir/$out";
    $node->runCmd($cmd); 
}


# a node is now passed in that can be used to run commands on a node using $node->runCmd("cmd")
# NOTE that in order to use this must add keepNodeForPostProcessing=yes to controller.prop file
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;
  return 1;
}

1;
