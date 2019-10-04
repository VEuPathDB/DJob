package DJob::DistribJobTasks::CNVTasks;

use lib "$ENV{GUS_HOME}/lib/perl";
use DJob::DistribJob::Task;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
    ["genomicSeqsFile", "none", "full path to reference genome fasta file"],
    ["bamFile", "none", "full path to alignment file"],
    ["gtfFile", "none", "full path to reference genome annotation in gtf format"],
    ["geneFootprintFile", "none", "full path to gene footprint file"],
    ["samtoolsIndex", "none", "full path to samtools index of reference fasta file"],
    ["sampleName", "", "sample or strain to be analysed"],
    ["window", 1000, "window size for binned coverage"],
    ["snpsClusterDir", "", "Dir containing snps data on cluster"],
);


sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
  my ($self, $inputDir) = @_;
}

sub initNode {
    my ($self, $node, $inputDir) = @_;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;
    return 1;
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $genomicSeqsFile = $self->getProperty("genomicSeqsFile");
    my $bamFile = $self->getProperty("bamFile");
    my $gtfFile = $self->getProperty("gtfFile");
    my $geneFootprintFile = $self->getProperty("geneFootprintFile");
    my $samtoolsIndex = $self->getProperty("samtoolsIndex");
    my $sampleName = $self->getProperty("sampleName");
    my $window = $self->getProperty("window");
    my $snpsClusterDir = $self->getProperty("snpsClusterDir");
    my $workingDir = "$node->{masterDir}/mainresult";
       
    
    my $cmd = "runCNVTasks.pl --genomicSeqsFile $genomicSeqsFile --bamFile $bamFile --gtfFile $gtfFile --geneFootprintFile $geneFootprintFile --samtoolsIndex $samtoolsIndex --workingDir $workingDir --sampleName $sampleName --window $window --snpsClusterDir $snpsClusterDir";
      
    return $cmd;
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
}

sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;
  return 1;  
}



1;
