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

    my $fastaFileName = $self->{props}->getProp("inputFilePath");

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

sub runSubTask { 
    my ($self, $node, $inputDir, $subTaskDir, $nodeSlotDir) = @_;

    my $rmPath = $self->{props}->getProp("rmPath");
    my $rmOptions = $self->{props}->getProp("rmOptions");
    my $trimDangling = $self->{props}->getProp("trimDangling") eq "y"? "--trimDangling" : "";

    my $options = $rmOptions eq "NONE"? "" : "--rmOptions $rmOptions";
    my $cmd = "repeatMasker --rmPath $rmPath $options --seqFile $nodeSlotDir/seqsubset.fsa --outFile blocked.seq --errorFile blocked.err $trimDangling";

    $node->execSubTask("$nodeSlotDir/result", "$subTaskDir/result", $cmd);
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $subTaskResultDir, $mainResultDir) = @_;

    &runCmd("cat $subTaskResultDir/blocked.seq >> $mainResultDir/blocked.seq");
    &runCmd("cat $subTaskResultDir/blocked.err >> $mainResultDir/blocked.err");
}
1;
