package DJob::DistribJobTasks::CRAIGPredTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["paramsFile", "", "Name of the file containing the gene model parameters"],
 ["prefixFilesPath", "",  "prefix to resource files, such as RNA-seq junction, coverage files and/or evidence signal and state files for ensemble models"],
 ["inputFilePath", "",  "full path to input FASTA file"],
 ["predictUTRs", "yes",  "Predictions will contain UTR regions. If not specified and paramters model UTRs, expect worse performance, so choose carefully"],
 ["format", "locs", "output format"],
 ["START-resource", " ", "resource name for Start signals"],
 ["STOP-resource", " ", "resource name for Stop signals"],
 ["DONOR-resource", " ", "resource name for Donor signals"],
 ["ACCEPTOR-resource", " ", "resource name for Acceptor signals"]
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
					"$subTaskDir/seqsubset.fa");
    }
    $self->runCmdOnNode($node, "touch $subTaskDir/seqsubset.fa.touch",1);
    $self->runCmdOnNode($node, "/bin/rm $subTaskDir/seqsubset.fa.touch",1);
    $self->runCmdOnNode($node, "cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $format = $self->getProperty("format");
    my $params = $self->getProperty("paramsFile");
    my $predictUTRs = ($self->getProperty("predictUTRs") eq "yes" ? "--predict-utrs" : "");
    my $prefixFilesPath = $self->getProperty("prefixFilesPath");
    my $startResource = ($self->getProperty("START-resource") eq " " ? "" : "--START-resource=".$self->getProperty("START-resource"));
    my $stopResource = ($self->getProperty("STOP-resource") eq " " ? "" : "--STOP-resource=".$self->getProperty("STOP-resource"));
    my $donorResource = ($self->getProperty("DONOR-resource") eq " " ? "" : "--DONOR-resource=".$self->getProperty("DONOR-resource"));
    my $acceptorResource = ($self->getProperty("ACCEPTOR-resource") eq " " ? "" : "--ACCEPTOR-resource=".$self->getProperty("ACCEPTOR-resource"));
    my $cmd = "craigPredict $predictUTRs  $startResource $stopResource $donorResource $acceptorResource --output-file=seqsubset.locs --format=$format --best=1 --strand=both --prefix-files=$prefixFilesPath $params seqsubset.fa";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    $self->runCmdOnNode($node, "cat $nodeExecDir/seqsubset.locs >> $mainResultDir/output0.locs");
    $self->runCmdOnNode($node, "cat $nodeExecDir/subtask.stderr >> $mainResultDir/stderr0", 1);
}
1;
