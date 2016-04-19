package DJob::DistribJobTasks::XSignalsTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["prefixFilesPath", "",  "prefix to RNA-seq junction and coverage files"],
 ["inputFilePath", "",  "full path to input FASTA file"],
 ["blockLen", "20",  "block length used in discarding permutations"],
 ["stranded", "no", "whether input RNA-Seq data is stranded or not\n"],
 ["numPermutations", "1000", "number of permutations for computing change points"]
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

sub initServer {
    my ($self, $inputDir) = @_;
    my $prefixOutputFiles = $self->getProperty("prefixFilesPath");
    $self->{nodeForInit}->runCmd("rm -f $prefixOutputFiles.xsignal");
    $self->{nodeForInit}->runCmd("rm -f $prefixOutputFiles.xsigscores");
    $self->{nodeForInit}->runCmd("rm -f ;$prefixOutputFiles.stderr");
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
					"$subTaskDir/seqsubset.fa");
    }
    $self->runCmdOnNode("touch $subTaskDir/seqsubset.fa.touch",1);
    $self->runCmdOnNode("/bin/rm $subTaskDir/seqsubset.fa.touch",1);
    $self->runCmdOnNode("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $blockLen = $self->getProperty("blockLen");
    my $prefixFilesPath = $self->getProperty("prefixFilesPath");
    my $stranded = ($self->getProperty("stranded") eq "no" ? "" : "--stranded");
    my $num_perms = $self->getProperty("numPermutations");
    my $cmd = "get_covxsignals $stranded --block-len=$blockLen --prefix-files=$prefixFilesPath --prefix-output=seqsubset --num-permutations=$num_perms seqsubset.fa";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    my $prefixOutputFiles = $self->getProperty("prefixFilesPath");
    $self->runCmdOnNode("cat $nodeExecDir/seqsubset.xsignal >> $prefixOutputFiles.xsignal");
    $self->runCmdOnNode("cat $nodeExecDir/seqsubset.xsigscores >> $prefixOutputFiles.xsigscores");
    $self->runCmdOnNode("cat $nodeExecDir/subtask.stderr >> $prefixOutputFiles.stderr");
}
1;
