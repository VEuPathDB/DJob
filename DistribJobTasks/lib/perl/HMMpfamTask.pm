package DJob::DistribJobTasks::HMMpfamTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["hmmpfamDir",   "",   "path to directory that contains the hmmpfam script"],
 ["dbFilePath",      "",     "full path to database file"],
 ["inputFilePath",   "",     "full path to input file"],
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;
    return 1;
}

##copy the pfam database file(s) to the nodedir
sub initNode {
    my ($self, $node, $inputDir) = @_;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->{props}->getProp("inputFilePath");

    if (-e "$fastaFileName.gz") {
	&runCmd("gunzip $fastaFileName.gz");
    }

    print "Counting sequences in $fastaFileName\n";
    $self->{fastaFile} = CBIL::Bio::FastaFileSequential->new($fastaFileName);
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir,$subTask) = @_;

    if(!$subTask->getRedoSubtask()){
      $self->{fastaFile}->writeSeqsToFile($start, $end, "$subTaskDir/seqsubset.fsa");
    }

    $self->runCmdOnNode($node, "cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $dbFilePath = $self->{props}->getProp("dbFilePath");
    my $hmmpfamDir = $self->{props}->getProp("hmmpfamDir");

    my $cmd =  "$hmmpfamDir/hmmpfam -A 0 --acc --cut_ga $dbFilePath $nodeExecDir/seqsubset.fsa";

    return $cmd;
}


sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $self->runCmdOnNode($node, "cat $nodeExecDir/subtask.output >> $mainResultDir/hmmpfam.out");
}
1;
