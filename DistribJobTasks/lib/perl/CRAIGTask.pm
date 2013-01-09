package DJob::DistribJobTasks::CRAIGTask;
#use lib "/home/abernal/GUS/gus_home/lib/perl";

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
 ["RSqType", "rum", "RNASeq type rum|gsnap|epdb. Only rum currently supported"],
 ["RSqDataSetID", "", "RNASeq dataset id"],
 ["RSqSampleID", "", "RNASeq sample id"],
 ["RSqInputDir", "", "Full path to RNASeq's alignments directory"],
 ["RSqStranded", "", "Whether RNASeq is stranded or not"],
 ["RSqWrapperIn", "cat", "Script used to read RNASeq coverage information"],
 ["dataDirOutput", "", "Full path to Preprocessing output directory"],
 ["InputAnnotFile", "", "Full path to genome isoform annotation"],
 ["InputAnnotFormat","gff3","Format for genome isoform annotation locs|gff3|gtf|genbank"],
 ["InputContigFile","","Full path to genome contig sequences"],
 ["InputContigFormat","fasta","Format for genome isoform annotation fasta|genbank"],
 ["transcriptTag", "exon", "Tag for recognizing transcribed regions in InputAnnotFile"],
 ["CDSTag", "CDS", "Tag for recognizing coding regions in InputAnnotFile"],
 ["gcClasses", "100", "gc isochores present in the genome. f.e. human would be 43,51,57"],
 ["species", "", "name of species tgondii|plasmo"],
 ["coverageDensity", "5", "Sample's coverage density. A number between 1(very sparse, like gregory datasets) to 9(very dense, like reid day4 dataset)."],
 ["geneModelType", "ngscraig", "model type ngscraig/ecraig/craig"],
 ["blockLen", "20",  "block length used in discarding permutations"],
 ["numPermutations", "1000", "number of permutations for computing change points"],
 ["modelUTRs", "yes", "models and predicts UTRs as part of genes"],
 ["modelDirOutput", "", "Full path to Learning and prediction directory"],
 ["genUTROnlyModel", "yes", "generates utr_only predictions"],
 ["numDJobNodes", "20", "Number of nodes used by internal distribjob calls"],
 ["DJobNodeClass", "LocalNode", "Node Class for the internal DJob calls, either LocalNode or SgeNode"],
 ["DJobInputBaseDir", "", "Base Input dir for internal distribjob calls"],
 ["minPercJunctionAlignments", "", "Minimum percentage of genes with introns that have good junction alignments that are needed in order to use the junction information as the sole source for splicing"]
);

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
  my ($self, $inputDir) = @_;
  # creating preconfig file
  open(PRECONFIG,">$inputDir/preconfig");
  my $rsqType = $self->getProperty('RSqType');
  my $rsqDataSetId = $self->getProperty('RSqDataSetID');
  my $rsqSampleId = $self->getProperty('RSqSampleID');
  my $rsqInputDir = $self->getProperty('RSqInputDir');
  my $rsqStranded = $self->getProperty('RSqStranded');
  my $rsqWrapperIn = $self->getProperty('RSqWrapperIn');
  print PRECONFIG "rnaseq\t$rsqType\t$rsqDataSetId\t$rsqSampleId\t$rsqInputDir\t$rsqStranded\t$rsqWrapperIn\n";
  close(PRECONFIG);

  my $outputDir = $self->getProperty('dataDirOutput');
  my $inputAnnotFmt = $self->getProperty('InputAnnotFormat');
  my $inputContig = $self->getProperty('InputContigFile');
  my $inputContigFmt = $self->getProperty('InputContigFormat');
  my $transcriptTag = $self->getProperty('transcriptTag');
  my $cdsTag = $self->getProperty('CDSTag');
  my $coverageDensity = $self->getProperty('coverageDensity');
  my $gcClasses = $self->getProperty('gcClasses');
  my $modelType = $self->getProperty('geneModelType');
  my $numPerms = $self->getProperty('numPermutations');
  my $blockLen = $self->getProperty('blockLen');
  my $numDJobNodes = $self->getProperty('numDJobNodes');
  my $djobNodeClass = $self->getProperty('DJobNodeClass'); 
  my $DJobInputBaseDir = $self->getProperty('DJobInputBaseDir');
  my $species = $self->getProperty('species');
  my $inputAnnotFile = $self->getProperty('InputAnnotFile');
  my $inputContigFile = $self->getProperty('InputContigFile');

  if (!-e $DJobInputBaseDir) {
      `mkdir $DJobInputBaseDir`;
  }

  if (!-e "$DJobInputBaseDir/preprocess") {
      `mkdir $DJobInputBaseDir/preprocess`;
  }
  else {
      `rm -rf $DJobInputBaseDir/preprocess/*`; 
  }    

  my $cmd = "craigPreprocess.py --pre-config $inputDir/preconfig --out-dir $outputDir --annot-fmt $inputAnnotFmt --contig-fmt $inputContigFmt --transcript-tag $transcriptTag --cds-tag $cdsTag --gc-classes $gcClasses --model $modelType --num-permutations $numPerms --block-length $blockLen --coverage-density $coverageDensity --djob-num-nodes $numDJobNodes --djob-node-class $djobNodeClass --djob-input $DJobInputBaseDir/preprocess $species $inputAnnotFile $inputContigFile config";
  
  $self->{nodeForInit}->runCmd($cmd);
  
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
    my $modelDirOutput = $self->getProperty('modelDirOutput');
    my $geneModelType = $self->getProperty('geneModelType');
    my $minPercJunctionAlignments = $self->getProperty('minPercJunctionAlignments');
    my $numDJobNodes = $self->getProperty('numDJobNodes');
    my $djobNodeClass = $self->getProperty('DJobNodeClass'); 
    my $DJobInputBaseDir = $self->getProperty('DJobInputBaseDir');
    my $genUTROnlyModel = $self->getProperty('genUTROnlyModel');
    my $modelUTRs = $self->getProperty('modelUTRs');
    my $species = $self->getProperty('species');
    my $dataDirOutput = $self->getProperty('dataDirOutput');

    if (!-e "$DJobInputBaseDir/model") {
	`mkdir $DJobInputBaseDir/model`;
    }
    else {
	`rm -rf $DJobInputBaseDir/model/*`; 
    }    

    my $cmd = "craig4eupath.py --force-train --out-dir $modelDirOutput --model $geneModelType --min-perc-junc-aligns $minPercJunctionAlignments --djob-num-nodes  $numDJobNodes --djob-node-class $djobNodeClass --djob-input $DJobInputBaseDir/model";
    
    if($genUTROnlyModel == "yes") {
	$cmd = $cmd." --utr-only-model";
    }
    
    if($modelUTRs == "yes") {
	$cmd = $cmd." --model-utrs";
    }
    
    $cmd = $cmd." $species $dataDirOutput/config";
#    print "Returning command: $cmd\n";
#    exit(0);  ##for testing
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
