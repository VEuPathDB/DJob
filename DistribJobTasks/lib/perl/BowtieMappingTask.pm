package DJob::DistribJobTasks::BowtieMappingTask;

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
	["mateA", "none", "full path to file of reads"],
	["mateB", "none", "full path to file of paired ends reads"],
	["bowtieIndex", "none", "full path of the bowtie indices .. likely same as fasta file"],
	["isColorspace", "false", "input sequence reads are in SOLiD colorspace.  Quality files must be exactly matename.qual"],
    ["removePCRDuplicates", "true", "remove PCR duplicates for any analysis involving read depth, e.g., ploidy, CNV, mapping replication origins"],
    ["extraBowtieParams", "none", "Bowtie parameters other than default"],
	["bowtie2", "default", "full path to the bowtie2 bin dir"],
	["sampleName", "", "strain to be put into output"],
	["deleteIntermediateFiles", "true", "[true]|false: if true then deletes intermediate files to save space"],

 ["writeBedFile", "true", "[true]|false: if true then runs bamToBed on unique and non unique mappers"],
 ["topLevelSeqSizeFile", "none", "required if writeBedFile turned on"],

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

    my $mateA = $self->getProperty ("mateA");
    my $mateB = $self->getProperty ("mateB");
    my $bowtieIndex = $self->getProperty ("bowtieIndex");
    my $isColorspace = $self->getProperty ("isColorspace");
    my $removePCRDuplicates = $self->getProperty ("removePCRDuplicates");
    my $extraBowtieParams = $self->getProperty ("extraBowtieParams");
    my $sampleName = $self->getProperty ("sampleName");
    my $wDir = "$node->{masterDir}/mainresult";
    my $bowtie2 = $self->getProperty ("bowtie2");


       
    
    my $cmd = "runBowtieMapping.pl --mateA $mateA".(-e "$mateB" ? " --mateB $mateB" : "");
    $cmd .= " --bowtieIndex $bowtieIndex";
    $cmd .= " --bowtie2 $bowtie2";
    if($self->getProperty('extraBowtieParams') ne 'none'){
      # The swaps are so getOpt::Long can take bowtieParams as a string
      $extraBowtieParams =~ s/-/_/g;
      $extraBowtieParams =~ s/ /#/g;
      $cmd .= " --extraBowtieParams $extraBowtieParams";
    }
    if($self->getProperty('isColorspace') eq 'true'){
      $cmd .= " --isColorspace";
    }
    if ($self->getProperty('removePCRDuplicates') eq 'true'){
      $cmd.= " --removePCRDuplicates";
    }
    $cmd .= " --sampleName $sampleName";
    $cmd .= " --workingDir $wDir" . ($self->getProperty('deleteIntermediateFiles') eq 'true' ? " --deleteIntermediateFiles" : "");
      
    return $cmd;
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
}

sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;

  my $writeBedFile = $self->getProperty("writeBedFile");

  # BED 
  if($writeBedFile && lc($writeBedFile) eq 'true') {

      my $topLevelSeqSizeFile = $self->getProperty("topLevelSeqSizeFile");
      my $sampleName = $self->getProperty ("sampleName");

      unless(-e $topLevelSeqSizeFile) {
	  die "Top Level Seq Size FIle $topLevelSeqSizeFile does not exist";
      }

      $self->runCmnOnNode($node, "samtools index $mainResultDir/${sampleName}.bam");
      $self->runCmdOnNode($node, "bamutils tobedgraph $mainResultDir/${sampleName}.bam >$mainResultDir/${sampleName}.bed");
  }


  }

1;
