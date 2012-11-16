   package DJob::DistribJobTasks::tRNAscanTask;

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
 ["tRNAscanDir",   "default",   "path to directory that contains the tRNAscan script"],
 ["inputFilePath",   "",     "full path to input file"],
 ["trainingOption", "", "training set used, eg. -G (general) or -C (covariance only)"]
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
    return 1;
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
    $node->runCmd("touch $subTaskDir/seqsubset.fsa");

    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $tRNAscanDir = $self->{props}->getProp("tRNAscanDir");
    $tRNAscanDir = $tRNAscanDir eq 'default' ? getProgramDir("tRNAscan-SE") : $tRNAscanDir;

    my $trainingOption = $self->{props}->getProp("trainingOption"); 

    # no longer support training option.  out of date and too slow
#    my $cmd =  "$tRNAscanDir/tRNAscan-SE -$trainingOption $nodeExecDir/seqsubset.fsa";
    my $cmd =  "$tRNAscanDir/tRNAscan-SE  $nodeExecDir/seqsubset.fsa";

    return $cmd;
}


sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $node->runCmd("cat $nodeExecDir/subtask.output >> $mainResultDir/trnascan.out");
}
1;
