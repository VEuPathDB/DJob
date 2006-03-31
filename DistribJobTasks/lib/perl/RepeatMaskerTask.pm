package DJob::DistribJobTasks::RepeatMaskerTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFile;
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
 ["rmOptions",     "NONE",  ""],
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

    print "Creating index for $fastaFileName (may take a while)\n";
    $self->{fastaFile} = CBIL::Bio::FastaFile->new($fastaFileName);
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir) = @_;

    $self->{fastaFile}->writeSeqsToFile($start, $end, 
					"$subTaskDir/seqsubset.fsa");

    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $rmPath = $self->getProperty("rmPath");
    my $rmOptions = $self->getProperty("rmOptions");
    my $trimDangling = $self->getProperty("trimDangling") eq "y"? "--trimDangling" : "";
    my $dangleMax = $self->getProperty("dangleMax");

    my $options = $rmOptions eq "NONE"? "" : "--rmOptions '$rmOptions'";
#    my $cmd = "repeatMasker --rmPath $rmPath $options --seqFile $nodeExecDir/seqsubset.fsa --outFile blocked.seq --errorFile blocked.err $trimDangling --dangleMax $dangleMax";
    my $cmd = "$rmPath/RepeatMasker $rmOptions seqsubset.fsa";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $node->runCmd("cat $nodeExecDir/seqsubset.fsa.masked >> $mainResultDir/blocked.seq");
    $node->runCmd("cat $nodeExecDir/subtask.stderr >> $mainResultDir/blocked.err");
}
1;
