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
 ["gtfFile",   "none",     "full path to the gtf file; Only required for cufflinks and junctions quantification"],
 ["maskFile",   "none",     "full path to the gtf mask file; Only required for cufflinks and junctions quantification"],
 ["sraSampleIdQueryList", "none", "Comma delimited list of identifiers that can be used to retrieve SRS samples"],
 ["extraGsnapParams", "none", "GSNAP parameters other than default"],
 ["outputFileBasename", "results", "Base name for the results file"],
 ["nPaths",   "30",     "Limits the number of nonunique mappers printed to a max of [30]"],
 ["deleteIntermediateFiles", "true", "[true]|false: if true then deletes intermediate files to save space"],
 ["quantifyWithCufflinks", "true", "[true]|false: if true then runs cufflinks on Unique and Multi Mappers"],
 ["writeBedFile", "true", "[true]|false: if true then runs bamToBed on unique and non unique mappers"],
 ["isStrandSpecific", "false", "[true]|false: if true then runs bamToBed on unique and non unique mappers"],
 ["quantifyJunctions", "true", "[true]|false: if true then runs cufflinks on Unique and Multi Mappers"],
 ["topLevelSeqSizeFile", "none", "required if writeBedFile turned on"],
 ["topLevelFastaFile", "none", "required if quantifyWithCufflinks is true"],
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

    my $cmd = "gsnap $extraGsnapParams --force-xs-dir -q $q --split-output 'split_output' --quiet-if-excessive -N 1 -s $iitFile -A sam -n $nPaths -D $databaseDirectory -d $databaseName  $mateA $mateB";
    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    my $outputBase = "$nodeExecDir/split_output";

    my $nes = `grep ^\@ -vl $outputBase*`;

    my @nonEmptySam = split(/\n/, $nes);

    # transform split output into bam files
    foreach my $sam (@nonEmptySam) {
      $node->runCmd("samtools view -Sb $sam > ${sam}.bam 2>>$nodeExecDir/subtask.stderr");
    }


    my @unique = map { $outputBase . "." . $_ } ("concordant_transloc",
                                                 "concordant_uniq",
                                                "halfmapping_transloc",
                                                 "halfmapping_uniq",
                                                 "paired_uniq_circular",
                                                 "paired_uniq_inv",
                                                 "paired_uniq_long",
                                                 "paired_uniq_scr",
                                                 "unpaired_transloc",
                                                 "unpaired_uniq",
    );

    my @all = map { $outputBase . "." . $_ } ("concordant_mult",
                                             "halfmapping_mult",
                                             "paired_mult",
                                             "unpaired_mult",
					     "concordant_transloc",
					     "concordant_uniq",
					     "halfmapping_transloc",
					     "halfmapping_uniq",
					     "paired_uniq_circular",
					     "paired_uniq_inv",
					     "paired_uniq_long",
					     "paired_uniq_scr",
					     "unpaired_transloc",
					     "unpaired_uniq",
    );


    my @uniqueBams;
    foreach my $bam (@unique) {
      push @uniqueBams, "${bam}.bam" if(-e "${bam}.bam");
    }
    my $uniqueBams = join(" ", @uniqueBams);

    my @allBams;
    foreach my $bam (@all) {
      push @allBams, "${bam}.bam" if(-e "${bam}.bam");
    }
    my $allBams = join(" ", @allBams);

    # merge into Unique and non unique files
    if(scalar @uniqueBams > 1) {
	$node->runCmd("samtools merge $nodeExecDir/unique.bam $uniqueBams");
    } 
    else {
	$node->runCmd("cp $uniqueBams $nodeExecDir/unique.bam" );
    }

    if(scalar @allBams > 1) {
	$node->runCmd("samtools merge $nodeExecDir/all.bam $allBams");
    }
    else {
	$node->runCmd("cp $allBams $nodeExecDir/all.bam");
    }

    # copy unique and non unique to mainresult dir;  Cannot merge into one file yet
    $node->runCmd("cp $nodeExecDir/unique.bam  $mainResultDir/$subTaskNum.unique.bam");
    $node->runCmd("cp $nodeExecDir/all.bam  $mainResultDir/$subTaskNum.all.bam");

    unless(-e "$mainResultDir/$subTaskNum.unique.bam" && -e "$mainResultDir/$subTaskNum.all.bam") {
      print "unique or all file does not exist in the mainResultDir\n";
      return 1;
    }

    return $node->getErr();
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;

  my $outputFileBasename = $self->getProperty("outputFileBasename");


  my @masterUnique = glob "$mainResultDir/*.unique.bam";
  my @masterAll = glob "$mainResultDir/*.all.bam";


  die "Did not find *unique* bam files in $mainResultDir/*.unique.bam" unless(scalar @masterUnique > 0);
  die "Did not find *all* bam files in $mainResultDir/*.all.bam" unless(scalar @masterAll > 0);

  # merge subtasks into unique and nonunique bams 
   if(scalar @masterUnique > 1) {
    $node->runCmd("samtools merge $mainResultDir/${outputFileBasename}_unique.bam $mainResultDir/*.unique.bam");
  }
  else {
    $node->runCmd("cp $masterUnique[0] $mainResultDir/${outputFileBasename}_unique.bam");
  }

  if(scalar @masterAll > 1) {
    $node->runCmd("samtools merge $mainResultDir/${outputFileBasename}_all.bam $mainResultDir/*.all.bam");
  }
  else {
    $node->runCmd("cp $masterAll[0] $mainResultDir/${outputFileBasename}_all.bam");
  }

  # sort bams
  $node->runCmd("samtools sort $mainResultDir/${outputFileBasename}_unique.bam $mainResultDir/${outputFileBasename}_unique_sorted");
  $node->runCmd("samtools sort $mainResultDir/${outputFileBasename}_all.bam $mainResultDir/${outputFileBasename}_all_sorted");

  # clean up some extra files
  unlink glob "$mainResultDir/*all.bam";
  unlink glob "$mainResultDir/*unique.bam";

  my $sidlist = $self->getProperty('sraSampleIdQueryList');

  if($sidlist && $sidlist ne 'none' && lc($self->getProperty('deleteIntermediateFiles')) eq 'true'){ ##have a value and other than default so reads were retrieved from sra
    my $mateA = $self->getProperty('mateA');
    my $mateB = $self->getProperty('mateB');
    unlink($mateA) if -e "$mateA";
    unlink($mateB) if -e "$mateB";
  }

  die "UNIQUE SORTED File [$mainResultDir/${outputFileBasename}_unique_sorted.bam] does not exist" unless(-e "$mainResultDir/${outputFileBasename}_unique_sorted.bam");
  die "ALL SORTED File [$mainResultDir/${outputFileBasename}_all_sorted.bam] does not exist" unless(-e "$mainResultDir/${outputFileBasename}_all_sorted.bam");


  # make bam index files
  $node->runCmd("samtools index $mainResultDir/${outputFileBasename}_unique_sorted.bam");
  $node->runCmd("samtools index $mainResultDir/${outputFileBasename}_all_sorted.bam");


  my $runCufflinks = $self->getProperty("quantifyWithCufflinks");
  my $writeBedFile = $self->getProperty("writeBedFile");
  my $quantifyJunctions = $self->getProperty("quantifyJunctions");


  my $isStrandSpecific = $self->getProperty("isStrandSpecific");

  # FPKM From Cufflinks
  if($runCufflinks && lc($runCufflinks) eq 'true') {
    my $gtfFile = $self->getProperty("gtfFile");
    my $maskFile = $self->getProperty("maskFile");
    my $topLevelFastaFile = $self->getProperty("topLevelFastaFile");

    my @libraryTypes = ("fr-unstranded");

    if($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
      @libraryTypes = ("fr-firststrand", "fr-secondstrand");
    }

    foreach my $lt (@libraryTypes) {
      $node->runCmd("cufflinks --no-effective-length-correction --compatible-hits-norm -M $maskFile --library-type '$lt' -o $mainResultDir -G $gtfFile $mainResultDir/${outputFileBasename}_unique_sorted.bam");
      rename "$mainResultDir/genes.fpkm_tracking", "$mainResultDir/genes.unique.fpkm_tracking.$lt";
      rename "$mainResultDir/isoforms.fpkm_tracking", "$mainResultDir/isoforms.unique.fpkm_tracking.$lt";

      $node->runCmd("cufflinks --no-effective-length-correction --compatible-hits-norm -M $maskFile --library-type '$lt' -o $mainResultDir -G $gtfFile $mainResultDir/${outputFileBasename}_all_sorted.bam");
      rename "$mainResultDir/genes.fpkm_tracking", "$mainResultDir/genes.all.fpkm_tracking.$lt";
      rename "$mainResultDir/isoforms.fpkm_tracking", "$mainResultDir/isoforms.all.fpkm_tracking.$lt";
    }
  }

  # Junctions
  if($quantifyJunctions && lc($quantifyJunctions eq 'true')) {
    $node->runCmd("gsnapSam2Junctions.pl  --is_bam  --input_file $mainResultDir/${outputFileBasename}_all_sorted.bam --output_file $mainResultDir/junctions.tab");
  }

  # BED 
  if($writeBedFile && lc($writeBedFile) eq 'true') {

    my $topLevelSeqSizeFile = $self->getProperty("topLevelSeqSizeFile");
    unless(-e $topLevelSeqSizeFile) {
      die "Top Level Seq Size FIle $topLevelSeqSizeFile does not exist";
    }

    # For strand specific datasets ... write out for and rev bed files
    if($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
      $node->runCmd("bamutils tobedgraph -plus $mainResultDir/${outputFileBasename}_unique_sorted.bam >$mainResultDir/${outputFileBasename}_unique_sorted_forward.bed");
      $node->runCmd("bamutils tobedgraph -plus $mainResultDir/${outputFileBasename}_all_sorted.bam >$mainResultDir/${outputFileBasename}_all_sorted_forward.bed");

      $node->runCmd("bamutils tobedgraph -minus $mainResultDir/${outputFileBasename}_unique_sorted.bam >$mainResultDir/${outputFileBasename}_unique_sorted_reverse.bed");
      $node->runCmd("bamutils tobedgraph -minus $mainResultDir/${outputFileBasename}_all_sorted.bam >$mainResultDir/${outputFileBasename}_all_sorted_reverse.bed");
    }
    else {
      $node->runCmd("bamutils tobedgraph $mainResultDir/${outputFileBasename}_unique_sorted.bam >$mainResultDir/${outputFileBasename}_unique_sorted.bed");
      $node->runCmd("bamutils tobedgraph $mainResultDir/${outputFileBasename}_all_sorted.bam >$mainResultDir/${outputFileBasename}_all_sorted.bed");
    }
  }

  return 1;
}

1;
