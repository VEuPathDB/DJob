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
 ["RSqOrientation", "", "Orientation of RNASeq reads: N|F|R. Use 'N' for non-stranded RNA reads. For paired reads you need to provide with orientation for each individual read, i.e. RR, FF, FR (illumina)"],
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
# ["coverageDensity", "-1", "Sample's coverage density. A number between 1(very sparse, like gregory datasets) to 9(very dense, like reid day4 dataset). This is roughyl a linear factor, but most datasets with good mapping should always be above 4 or 5. For example, reid3 should around 7, sibley around 5, oocyst(the bottom line for good datasets) around 4. Gregory data sets should all be 1 or 2. This can be extracted from the mapping_stats.txt file in rum. Just check the number of total mapped reads, consider paired reads to be better than single reads and land in a number. This estimated number is not so much important as long as it's a relatively good guess. If absent the program will try to estimate one by default, using the total number of junctions that have been covered"],
 ["geneModelType", "ngscraig", "model type ngscraig/ecraig/craig"],
 ["blockLen", "20",  "block length used in discarding permutations"],
 ["numPermutations", "1000", "number of permutations for computing change points"],
 ["trainingIterations", "10", "number of training iterations for learning convergence"],
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
    my $rsqOrientation = $self->getProperty('RSqOrientation');
    my $rsqWrapperIn = $self->getProperty('RSqWrapperIn');
    print PRECONFIG "rnaseq\t$rsqType\t$rsqDataSetId\t$rsqSampleId\t$rsqInputDir\t$rsqOrientation\t$rsqWrapperIn\n";
    close(PRECONFIG);
    
    my $dataDirOutput = $self->getProperty('dataDirOutput');
    my $inputAnnotFmt = $self->getProperty('InputAnnotFormat');
    my $inputContig = $self->getProperty('InputContigFile');
    my $inputContigFmt = $self->getProperty('InputContigFormat');
    my $transcriptTag = $self->getProperty('transcriptTag');
    my $cdsTag = $self->getProperty('CDSTag');
#    my $coverageDensity = $self->getProperty('coverageDensity');
    my $gcClasses = $self->getProperty('gcClasses');
    my $geneModelType = $self->getProperty('geneModelType');
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
    
    my $cmd1 = "craigPreprocess.py --pre-config $inputDir/preconfig --out-dir $dataDirOutput --annot-fmt $inputAnnotFmt --contig-fmt $inputContigFmt --transcript-tag $transcriptTag --cds-tag $cdsTag --gc-classes $gcClasses --model $geneModelType --num-permutations $numPerms --block-length $blockLen --djob-num-nodes $numDJobNodes --djob-node-class $djobNodeClass --djob-input $DJobInputBaseDir/preprocess $species $inputAnnotFile $inputContigFile config";

    my $modelDirOutput = $self->getProperty('modelDirOutput');
    my $minPercJunctionAlignments = $self->getProperty('minPercJunctionAlignments');

    my $genUTROnlyModel = $self->getProperty('genUTROnlyModel');
    my $modelUTRs = $self->getProperty('modelUTRs');
    my $dataDirOutput = $self->getProperty('dataDirOutput');
    my $trainingIterations = $self->getProperty('trainingIterations');
    
    if (!-e "$DJobInputBaseDir/model") {
	`mkdir $DJobInputBaseDir/model`;
    }
    else {
	`rm -rf $DJobInputBaseDir/model/*`; 
    }    
    
    my $cmd2 = "craig4eupath.py --training-iterations $trainingIterations --force-train --out-dir $modelDirOutput --model $geneModelType --min-perc-junc-aligns $minPercJunctionAlignments --djob-num-nodes  $numDJobNodes --djob-node-class $djobNodeClass --djob-input $DJobInputBaseDir/model";
    
    if($genUTROnlyModel == "yes") {
	$cmd2 = $cmd2." --utr-only-model";
    }
    
    if($modelUTRs == "yes") {
	$cmd2 = $cmd2." --model-utrs";
    }
    
    $cmd2 = $cmd2." $species $dataDirOutput/config";
#    exit(0);  ##for testin
    print "CRAIG Run Command: $cmd1 && $cmd2\n";
    $self->{nodeForInit}->runCmd($cmd1." && ".$cmd2);
#    $self->{nodeForInit}->runCmd($cmd1);
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
    my $cmd = "echo hola";
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
