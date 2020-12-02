package DJob::DistribJobTasks::Hisat2Task;

use lib "$ENV{GUS_HOME}/lib/perl";
use DJob::DistribJob::Task;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;
use Data::Dumper;

@ISA = (DJob::DistribJob::Task);
use strict;
# [name, default (or null if reqd), comment]
my @properties = 
    (
     ["mateA",   "",     "full path to reads file"],
     ["mateB",   "none",     "full path to paired reads file (optional)"],
     ["genomeDatabase",   "",     "full path to the genome database"],# hisat2 database
     ["maskFile",   "none",     "full path to the gtf masked file (rRNAs removed); required for HTseq"],
     ["sraSampleIdQueryList", "none", "Comma delimited list of identifiers that can be used to retrieve SRS samples"], # these should go to EBI, but I'll leave here just in case
     ["hisat2", "default", "full path to the bowtie2 bin dir"],
     ["extraHisatParams", "none", "Hisat2 parameters other than default"], #hisat2 params
     ["deleteIntermediateFiles", "true", "[true]|false: if true then deletes intermediate files to save space"],
     ["quantify", "true", "[true]|false: if true then runs HTSeq"],
     ["writeCovFiles", "true", "[true]|false: if true then runs bamutils"],
     ["isStrandSpecific", "false", "[true]|false"],
     ["quantifyJunctions", "true", "[true]|false: if true then runs gsnapSam2Junctions"],
     ["topLevelGeneFootprintFile", "none", "required if quantify is true"],
     ["hasPairedEnds", "false", "true|[false]"],
     ["ppn", 4, "number of processors to use"],
     ["maxIntronLen", 500000, "maximum intron length"]
    );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;
    my $baseName;
    ##need to download fastq from sra if sample ids passed in.
    my $sidlist = $self->getProperty('sraSampleIdQueryList');
    my $isPairedEnd = $self->getProperty('hasPairedEnds');

    my $mateA = $self->getProperty('mateA');
    my $mateB = $self->getProperty('mateB');

    # data from SRA - these should go to EBI but leaving here just in case
    if($sidlist && $sidlist ne 'none'){ 
        
        if(!$mateA || $mateA eq 'none'){
            $mateA = "$inputDir/reads_1.fastq";
            $self->setProperty('mateA',"$mateA");
            
            $mateB = "$inputDir/reads_2.fastq";
            $self->setProperty('mateB',"$mateB");
            $baseName = "reads";
        }
        if(-e "$mateA"){
            print "reads file $mateA already present so not retrieving from SRA\n";
        }
        else{  ##need to retrieve here
            print "retrieving reads from SRA for '$sidlist'\n";
            &runCmd("getFastqFromSra.pl --workingDir $inputDir --readsOne $mateA --readsTwo $mateB --sampleIdList '$sidlist' --pairs $isPairedEnd");
        }
    }
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
    my ($self, $node, $inputDir, $nodeExecDir,$subtaskNumber,$mainResultDir) = @_;

    my $mateA = $self->getProperty ("mateA");
    my $mateB = $self->getProperty ("mateB");
    $mateB = undef if(lc($mateB) eq 'none');
    my $genomeDatabase = $self->getProperty("genomeDatabase");
    my $ppn = $self->getProperty("ppn");
    my $maxIntronLen = $self->getProperty("maxIntronLen");
    my $maskedFile = $self->getProperty("maskFile");
    my $topLevelGeneFootprintFile = $self->getProperty("topLevelGeneFootprintFile");
    my $hisat2 = $self->getProperty("hisat2");

    my $cmd = "runHisat2Pipeline.pl --mateA $mateA --hisat2Index $genomeDatabase --sampleName results --workingDir $mainResultDir --ppn $ppn --maxIntronLen $maxIntronLen --maskedFile $maskedFile --topLevelGeneFootprintFile $topLevelGeneFootprintFile --hisat2 $hisat2";


    if ($mateB && -e $mateB) {
        $cmd = $cmd . " --mateB $mateB";
    }

    my $extraHisatParams = $self->getProperty("extraHisatParams");
    if($extraHisatParams ne 'none'){
        # The swaps are so getOpt::Long can take bowtieParams as a string
        $extraHisatParams =~ s/-/_/g;
        $extraHisatParams =~ s/ /#/g;
    } else {
        $extraHisatParams = '';
    }
    $cmd = $cmd . " --extraHisatParams $extraHisatParams";
      

    #add boolean properties
    my $deleteIntermediateFiles = $self->getProperty("deleteIntermediateFiles");
    my $quantify = $self->getProperty("quantify");
    my $quantifyJunctions = $self->getProperty("quantifyJunctions");
    my $writeCovFiles = $self->getProperty("writeCovFiles");
    my $isStrandSpecific = $self->getProperty ("isStrandSpecific");

    if ($deleteIntermediateFiles && lc($deleteIntermediateFiles) eq 'true') {
        $cmd = $cmd . " --deleteIntermediateFiles";
    }
    if ($quantify && lc($quantify) eq 'true') {
        $cmd = $cmd . " --quantify";
    }
    if ($quantifyJunctions && lc($quantifyJunctions) eq 'true') {
        $cmd = $cmd . " --junctions";
    }
    if ($writeCovFiles && lc($writeCovFiles) eq 'true') {
        $cmd = $cmd . " --writeCov";
    }
    if ($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
        $cmd = $cmd . " --isStranded";
    }

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;
  return 1;
}

1;
