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
 ["sraSampleIdQueryList", "none", "Comma delimited list of identifiers that can be used to retrieve SRS samples"],
 ["extraGsnapParams", "none", "GSNAP parameters other than default"],
 ["outputFileBasename", "results", "Base name for the results file"],
 ["nPaths",   "30",     "Limits the number of nonunique mappers printed to a max of [30]"],
 ["deleteIntermediateFiles", "true", "[true]|false: if true then deletes intermediate files to save space"],
 ["quantifyWithCufflinks", "true", "[true]|false: if true then runs cufflinks on Unique and Multi Mappers"],
 ["writeBedFile", "true", "[true]|false: if true then runs bamToBed on unique and non unique mappers"],
 ["isStrandSpecific", "false", "[true]|false: if true then runs bamToBed on unique and non unique mappers"],

 ["quantifyJunctions", "true", "[true]|false: if true then runs cufflinks on Unique and Multi Mappers"],
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
    my $genomeDatabase = $self->getProperty("genomeDatabase");
    my $iitFile = $self->getProperty("iitFile");

    my $databaseDirectory = dirname($genomeDatabase);
    my $databaseName = basename($genomeDatabase);

    my $nPaths = $self->getProperty("nPaths");

    my $wDir = "$node->{masterDir}/mainresult";

    my $totalSubtasks = int($self->{size} / $self->{subTaskSize});
    $totalSubtasks += 1 if $self->{size} % $self->{subTaskSize};

    my $q = $subtaskNumber - 1 . "/" . $totalSubtasks;

    my $extraGsanpParams = $self->getProperty("extraGsnapParams") eq "none" ? undef : $self->getProperty("extraGsnapParams");

    my $cmd = "gsnap $extraGsnapParams --force-xs-dir -q $q --split-output 'split_output' --quiet-if-excessive --nofails -N 1 -s $iitFile -A sam -n $nPaths -D $databaseDirectory -d $databaseName  $mateA $mateB";
    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    my $outputBase = "$nodeExecDir/split_output";

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

    my @nu = map { $outputBase . "." . $_ } ("concordant_mult",
                                             "halfmapping_mult",
                                             "paired_mult",
                                             "unpaired_mult",
    );

    # transform split output into bam files
    foreach my $bam (@unique, @nu) {
      $node->runCmd("samtools view -Sb $bam > ${bam}.bam 2>>$nodeExecDir/subtask.stderr; echo done");
    }

    my @uniqueBams = map { "${_}.bam" } @unique ;
    my $uniqueBams = join(" ", @uniqueBams);

    my @nuBams = map { "${_}.bam" } @nu;
    my $nuBams = join(" ", @nuBams);

    # merge into Unique and non unique files
    $node->runCmd("samtools merge $nodeExecDir/unique.bam $uniqueBams" );
    $node->runCmd("samtools merge $nodeExecDir/nu.bam $nuBams");

    # copy unique and non unique to mainresult dir;  Cannot merge into one file yet
    $node->runCmd("cp $nodeExecDir/unique.bam  $mainResultDir/$subTaskNum.unique.bam");
    $node->runCmd("cp $nodeExecDir/nu.bam  $mainResultDir/$subTaskNum.nu.bam");
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;

  my $outputFileBasename = $self->getProperty("outputFileBasename");

  mkdir "$mainResultDir/unique";
  mkdir "$mainResultDir/nu";

  # merge subtasks into unique and nonunique bams 
  $node->runCmd("samtools merge $mainResultDir/unique/${outputFileBasename}_unique.bam $mainResultDir/*.unique.bam");
  $node->runCmd("samtools sort $mainResultDir/unique/${outputFileBasename}_unique.bam $mainResultDir/unique/${outputFileBasename}_unique_sorted");

  # sort bams
  $node->runCmd("samtools merge $mainResultDir/nu/${outputFileBasename}_nu.bam $mainResultDir/*.nu.bam");
  $node->runCmd("samtools sort $mainResultDir/nu/${outputFileBasename}_nu.bam $mainResultDir/nu/${outputFileBasename}_nu_sorted");

  # clean up 
  unlink glob "$mainResultDir/*.bam";

  my $sidlist = $self->getProperty('sraSampleIdQueryList');

  if($sidlist && $sidlist ne 'none' && lc($self->getProperty('deleteIntermediateFiles')) eq 'true'){ ##have a value and other than default so reads were retrieved from sra
    my $mateA = $self->getProperty('mateA');
    my $mateB = $self->getProperty('mateB');
    unlink($mateA) if -e "$mateA";
    unlink($mateB) if -e "$mateB";
  }

  my $runCufflinks = $self->getProperty("quantifyWithCufflinks");
  my $writeBedFile = $self->getProperty("writeBedFile");
  my $isStrandSpecific = $self->getProperty("isStrandSpecific");

  # FPKM From Cufflinks
  if($runCufflinks && lc($runCufflinks) eq 'true') {
    my $gtfFile = $self->getProperty("gtfFile");

    my @libraryTypes = ("fr-unstranded");

    if($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
      @libraryTypes = ("fr-firststrand", "fr-secondstrand");
    }

    foreach my $lt (@libraryTypes) {
      $node->runCmd("cufflinks --library-type '$lt' -o $mainResultDir/unique -G $gtfFile $mainResultDir/unique/${outputFileBasename}_unique_sorted.bam");
      rename "$mainResultDir/unique/genes.fpkm_tracking" "$mainResultDir/unique/genes.fpkm_tracking.$lt";

      $node->runCmd("cufflinks --library-type '$lt' -o $mainResultDir/nu -G $gtfFile $mainResultDir/nu/${outputFileBasename}_nu_sorted.bam");
      rename "$mainResultDir/nu/genes.fpkm_tracking" "$mainResultDir/nu/genes.fpkm_tracking.$lt";
    }
  }

  # BED 
  if($writeBedFile && lc($writeBedFile) eq 'true') {

    # For strand specific datasets ... write out for and rev bed files
    if($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
      $node->runCmd("samtools view -b -F 16 $mainResultDir/unique/${outputFileBasename}_unique_sorted.bam >$mainResultDir/unique/${outputFileBasename}_unique_sorted_forward.bam"  );
      $node->runCmd("samtools view -b -f 16 $mainResultDir/unique/${outputFileBasename}_unique_sorted.bam >$mainResultDir/unique/${outputFileBasename}_unique_sorted_reverse.bam"  );

      $node->runCmd("bamToBed -i $mainResultDir/unique/${outputFileBasename}_unique_sorted_forward.bam >$mainResultDir/unique/${outputFileBasename}_unique_sorted_forward.bed");
      $node->runCmd("bamToBed -i $mainResultDir/nu/${outputFileBasename}_nu_sorted_forward.bam >$mainResultDir/nu/${outputFileBasename}_nu_sorted_forward.bed");

      $node->runCmd("bamToBed -i $mainResultDir/unique/${outputFileBasename}_unique_sorted_reverse.bam >$mainResultDir/unique/${outputFileBasename}_unique_sorted_reverse.bed");
      $node->runCmd("bamToBed -i $mainResultDir/nu/${outputFileBasename}_nu_sorted_reverse.bam >$mainResultDir/nu/${outputFileBasename}_nu_sorted_reverse.bed");
    }
    else {
      $node->runCmd("bamToBed -i $mainResultDir/unique/${outputFileBasename}_unique_sorted.bam >$mainResultDir/unique/${outputFileBasename}_unique_sorted.bed");
      $node->runCmd("bamToBed -i $mainResultDir/nu/${outputFileBasename}_nu_sorted.bam >$mainResultDir/nu/${outputFileBasename}_nu_sorted.bed");
    }
  }


}

1;
