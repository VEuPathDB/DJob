package DJob::DistribJobTasks::GSNAPTask;

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
 ["mateA",   "",     "full path to reads file"],
 ["mateB",   "none",     "full path to paired reads file (optional)"],
 ["genomeDatabase",   "",     "full path to the genome database"],
 ["iitFile",   "none",     "full path to the iit file for splice sites"],
 ["maskedFile",   "none",     "full path to the gtf masked file (rRNAs removed); required for Cufflinks and HTseq"],
 ["sraSampleIdQueryList", "none", "Comma delimited list of identifiers that can be used to retrieve SRS samples"],
 ["extraGsnapParams", "none", "GSNAP parameters other than default"],
 ["outputFileBasename", "results", "Base name for the results file"],
 ["nPaths",   "30",     "Limits the number of nonunique mappers printed to a max of [30]"],
 ["deleteIntermediateFiles", "true", "[true]|false: if true then deletes intermediate files to save space"],
 ["quantify", "true", "[true]|false: if true then runs Cufflinks and HTSeq"],
 ["writeCovFiles", "true", "[true]|false: if true then runs sam2cov"],
 ["isStrandSpecific", "false", "[true]|false"],
 ["quantifyJunctions", "true", "[true]|false: if true then runs gsnapSam2Junctions"],
 ["topLevelFastaFaiFile", "none", "required if writeCovFiles is turned on"],
 ["topLevelGeneFootprintFile", "none", "required if quantify is true"]
);

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
  my ($self, $inputDir) = @_;
  ##need to download fastq from sra if sample ids passed in.
  my $sidlist = $self->getProperty('sraSampleIdQueryList');

  if($sidlist && $sidlist ne 'none'){ ##have a value and other than default
    my $mateA = $self->getProperty('mateA');
    my $mateB = $self->getProperty('mateB');

    if(!$mateA || $mateA eq 'none'){
      $mateA = "$inputDir/reads_1.fastq";
      $self->setProperty('mateA',"$mateA");

      $mateB = "$inputDir/reads_2.fastq";
      $self->setProperty('mateB',"$mateB");
    }

    if(-e "$mateA"){
      print "reads file $mateA already present so not retrieving from SRA\n";
    }else{  ##need to retrieve here
      print "retrieving reads from SRA for '$sidlist'\n";
      $self->{nodeForInit}->runCmd("getFastqFromSra.pl --workingDir $inputDir --readsOne $mateA --readsTwo $mateB --sampleIdList '$sidlist'");
    }
  } 

}

sub initNode {
    my ($self, $node, $inputDir) = @_;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $reads = $self->getProperty('mateA');
    my $readLineCount = `wc -l $reads`;

    # Try to process ~ 1 Million reads or less per node.  div by 4 because of fastq.  Fasta files will have up to 2 million reads per process
    my $guessInputSize = int($readLineCount / 4);
    if($guessInputSize < 1) {
      $guessInputSize = 1;
    }

    $self->{inputSetSize} = $guessInputSize;

    return $self->{inputSetSize};
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
    my $iitFile = $self->getProperty("iitFile");

    my $databaseDirectory = dirname($genomeDatabase);
    my $databaseName = basename($genomeDatabase);

    my $nPaths = $self->getProperty("nPaths");

    my $wDir = "$node->{masterDir}/mainresult";

    my $totalSubtasks = int($self->{size} / $self->{subTaskSize});
    $totalSubtasks += 1 if $self->{size} % $self->{subTaskSize};

    my $q = $subtaskNumber - 1 . "/" . $totalSubtasks;

    my $extraGsnapParams = $self->getProperty("extraGsnapParams") eq "none" ? undef : $self->getProperty("extraGsnapParams");

    my $cmd = "gsnap $extraGsnapParams --force-xs-dir -q $q  --quiet-if-excessive -N 1 -s $iitFile -A sam -n $nPaths -D $databaseDirectory -d $databaseName  $mateA $mateB";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $node->runCmd("samtools view -Sb $nodeExecDir/subtask.output > $mainResultDir/${subTaskNum}_node.bam 2>>$nodeExecDir/subtask.stderr");

    return $node->getErr();
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;

  my $outputFileBasename = $self->getProperty("outputFileBasename");

  my @bams = glob "$mainResultDir/*_node.bam";

  die "Did not find  bam files in $mainResultDir/*_node.bam" unless(scalar @bams > 0);

  if(scalar @bams > 1) {
    $node->runCmd("samtools merge $mainResultDir/${outputFileBasename}.bam $mainResultDir/*_node.bam");
  }
  else {
    $node->runCmd("cp $bams[0] $mainResultDir/${outputFileBasename}.bam");
  }

  # sort bams by location
  $node->runCmd("samtools sort $mainResultDir/${outputFileBasename}.bam $mainResultDir/${outputFileBasename}_sorted");

  # clean up some extra files
  unlink glob "$mainResultDir/*_node.bam";

  my $sidlist = $self->getProperty('sraSampleIdQueryList');

  if($sidlist && $sidlist ne 'none' && lc($self->getProperty('deleteIntermediateFiles')) eq 'true'){ ##have a value and other than default so reads were retrieved from sra
    my $mateA = $self->getProperty('mateA');
    my $mateB = $self->getProperty('mateB');
    unlink($mateA) if -e "$mateA";
    unlink($mateB) if -e "$mateB";
  }

  my $runQuant = $self->getProperty("quantify");
  my $writeCovFiles = $self->getProperty("writeCovFiles");
  my $quantifyJunctions = $self->getProperty("quantifyJunctions");
  my $isStrandSpecific = $self->getProperty("isStrandSpecific");

  # Quantification
  if($runQuant && lc($runQuant) eq 'true') {
    my $maskedFile = $self->getProperty("maskedFile");
    my $topLevelGeneFootprintFile = $self->getProperty("topLevelGeneFootprintFile");

   # Cufflinks

    if($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
	$node->runCmd("cufflinks --no-effective-length-correction --compatible-hits-norm --library-type fr-firststrand -o $mainResultDir -G $maskedFile $mainResultDir/${outputFileBasename}_sorted.bam");
	rename "$mainResultDir/genes.fpkm_tracking", "$mainResultDir/genes.cuff.firststrand.fpkm_tracking";
	rename "$mainResultDir/isoforms.fpkm_tracking", "$mainResultDir/isoforms.cuff.firststrand.fpkm_tracking";

	$node->runCmd("cufflinks --no-effective-length-correction --compatible-hits-norm --library-type fr-secondstrand -o $mainResultDir -G $maskedFile $mainResultDir/${outputFileBasename}_sorted.bam");
	rename "$mainResultDir/genes.fpkm_tracking", "$mainResultDir/genes.cuff.secondstrand.fpkm_tracking";
	rename "$mainResultDir/isoforms.fpkm_tracking", "$mainResultDir/isoforms.cuff.secondstrand.fpkm_tracking";
    }

    else {
	$node->runCmd("cufflinks --no-effective-length-correction --compatible-hits-norm --library-type fr-unstranded -o $mainResultDir -G $maskedFile $mainResultDir/${outputFileBasename}_sorted.bam");
	rename "$mainResultDir/genes.fpkm_tracking", "$mainResultDir/genes.cuff.unstranded.fpkm_tracking";
	rename "$mainResultDir/isoforms.fpkm_tracking", "$mainResultDir/isoforms.cuff.unstranded.fpkm_tracking";
    }
    
    # HTSeq
    my @modes = ('union', 'intersection-nonempty', 'intersection-strict');
    if ($isStrandSpecific && lc($isStrandSpecific) eq 'true') {

	for (my $i=0; $i<@modes; $i++) {
	    my $mode = $modes[$i];
	    $node->runCmd("python -m HTSeq.scripts.count --format=bam --order=pos --stranded=reverse --type=exon --idattr=gene_id --mode=$mode $mainResultDir/${outputFileBasename}_sorted.bam $maskedFile > $mainResultDir/genes.htseq-$mode.firststrand.counts");
	    $node->runCmd("python -m HTSeq.scripts.count --format=bam --order=pos --stranded=yes --type=exon --idattr=gene_id --mode=$mode $mainResultDir/${outputFileBasename}_sorted.bam $maskedFile > $mainResultDir/genes.htseq-$mode.secondstrand.counts");

	    $node->runCmd("makeFpkmFromHtseqCounts.pl --geneFootprintFile $topLevelGeneFootprintFile --countFile $mainResultDir/genes.htseq-$mode.firststrand.counts --fpkmFile $mainResultDir/genes.htseq-$mode.firststrand.fpkm --antisenseCountFile $mainResultDir/genes.htseq-$mode.secondstrand.counts --antisenseFpkmFile $mainResultDir/genes.htseq-$mode.secondstrand.fpkm");
	}
    }
    else {
      for (my $i=0; $i<@modes; $i++) {
        my $mode = $modes[$i];
        $node->runCmd("python -m HTSeq.scripts.count --format=bam --order=pos --stranded=no --type=exon --idattr=gene_id --mode=$mode $mainResultDir/${outputFileBasename}_sorted.bam $maskedFile > $mainResultDir/genes.htseq-$mode.unstranded.counts");
	$node->runCmd("makeFpkmFromHtseqCounts.pl --geneFootprintFile $topLevelGeneFootprintFile --countFile $mainResultDir/genes.htseq-$mode.unstranded.counts --fpkmFile $mainResultDir/genes.htseq-$mode.unstranded.fpkm");
      }
    }
  }

  # Junctions
  if($quantifyJunctions && lc($quantifyJunctions eq 'true')) {
    $node->runCmd("gsnapSam2Junctions.pl  --is_bam  --input_file $mainResultDir/${outputFileBasename}_sorted.bam --output_file $mainResultDir/junctions.tab");
  }

  # COVERAGE PLOTS
  if($writeCovFiles && lc($writeCovFiles) eq 'true') {
    my $topLevelFastaFaiFile = $self->getProperty("topLevelFastaFaiFile");
    unless(-e $topLevelFastaFaiFile) {
      die "Top Level Genome fa.fai File $topLevelFastaFaiFile does not exist";
    }
    $node->runCmd("samtools view $mainResultDir/${outputFileBasename}.bam > $mainResultDir/${outputFileBasename}.sam");
    my $mateB = $self->getProperty('mateB');
    my $isPairedEnd = 1;
    $isPairedEnd = 0 if(lc($mateB) eq 'none');
 
    if($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
      $node->runCmd("sam2cov -s 1 -e $isPairedEnd -p $mainResultDir/${outputFileBasename}.firststrand. $topLevelFastaFaiFile $mainResultDir/${outputFileBasename}.sam");
      $node->runCmd("sam2cov -s 2 -e $isPairedEnd -p $mainResultDir/${outputFileBasename}.secondstrand. $topLevelFastaFaiFile $mainResultDir/${outputFileBasename}.sam");
    }
    else {
      $node->runCmd("sam2cov -s 0 -e $isPairedEnd -p $mainResultDir/${outputFileBasename}.unstranded. $topLevelFastaFaiFile $mainResultDir/${outputFileBasename}.sam");
    }
    unlink("$mainResultDir/${outputFileBasename}.sam");
  }

  return 1;
}

1;
