package DJob::DistribJobTasks::RepeatMaskerTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["rmPath",        "",  "eg: /export/Bioinformatics/usr/local/src/bio/RepeatMasker/latest"],
 ["inputFilePath", "",  "full path to input file"],
 ["trimDangling",  "",  "y or n"],
 ["rmParamsFile",     "",  "full path to file containing RepeatMasker command line options. Can be all on one line or one line / option."],
 ["dangleMax", "", "trim this or fewer bases"]

 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

sub initServer {
    my ($self, $inputDir) = @_;

}

sub initNode {
    my ($self, $node, $inputDir) = @_;

}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->getProperty("inputFilePath");

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
      $self->{fastaFile}->writeSeqsToFile($start, $end, 
					"$subTaskDir/seqsubset.fsa");

      my $rmParamsFile = $self->getProperty("rmParamsFile");
      &runCmd("cp $inputDir/$rmParamsFile $subTaskDir");
    }
    $node->runCmd("touch $subTaskDir/seqsubset.fsa");
    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $rmPath = $self->getProperty("rmPath");
    my $rmParamsFile = $self->getProperty("rmParamsFile");
    my $trimDangling = $self->getProperty("trimDangling") eq "y"? "--trimDangling" : "";
    my $dangleMax = $self->getProperty("dangleMax");

    my $cmd = "repeatMasker --rmPath $rmPath --seqFile seqsubset.fsa --outFile blocked.seq --errorFile blocked.err $trimDangling --dangleMax $dangleMax --rmParamsFile $rmParamsFile";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $node->runCmd("cat $nodeExecDir/blocked.seq >> $mainResultDir/blocked.seq");
    $node->runCmd("cat $nodeExecDir/blocked.err >> $mainResultDir/blocked.err");
}
1;
