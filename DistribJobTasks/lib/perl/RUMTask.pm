package DJob::DistribJobTasks::RUMTask;

use DJob::DistribJob::Task;
#use CBIL::Bio::FastaFile;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["readFilePath",   "",     "full path to read file"],
 ["pairedReadFilePath",   "none",     "full path to paired read file (optional)"],
 ["genomeFastaFile",   "",     "genome file in fasta format"],
 ["genomeBowtieIndex",   "none",     "genome bowtie index file"],
 ["transcriptFastaFile",   "none",     "transcript file in fasta format"],
 ["transcriptBowtieIndex",   "none",     "transcript bowtie index files"],
 ["geneAnnotationFile",   "",     "geneAnnotationFile in Gregs format (ucsc format)"],
 ["bowtieBinDir",   "",     "bowtie bin directory"],
 ["blatExec",   "",     "blat executable, must be full path unless defined in your path"],
 ["mdustExec",   "",     "mdust executable, must be full path unless defined in your path"],
 ["perlScriptsDir",   "",     "directory for RUM scripts"],
 ["limitNU",   "30",     "Limits the number of ambiguous mappers to a max of [30]"],
 ["numInsertions",   "1",     "number of instertions to allow when parsing BLAT [1] (note, if paired end data the number of insertions is constrained to 0 or 1"],
 ["minBlatIdentity",   "93",     "run BLAT with minimum identity of [93]"],
 ["createSAMFile",   "false",     "create SAM file if 1 ([0] | 1)"],
 ["countMismatches",   "false",     "report in the last column the number of mismatches, ignoring insertions ([0] | 1)"]
 );

## note passing in additional param that skips computing the size in constructor
sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties, 1); 
    return $self;
}

# called once .. do everything here upstream of sending chunks to nodes
# if have pairedReadFile check to see that size of chunk is even number
# combine files into single fasta file
# get qual values?
# make bowtie indices (genome and transcript)
# either create fasta index here or in getCount method
sub initServer {
  my ($self, $inputDir) = @_;
  
  ## NOTE: need to make certain that subtask size is an even number.
  die "subTaskSize must be an even number!\n" if ($self->{subTaskSize} % 2);

  my $readFilePath = $self->getProperty("readFilePath");
  my $pairedReadFilePath = $self->getProperty("pairedReadFilePath");
  my $genomeFastaFile = $self->getProperty("genomeFastaFile");
  my $transcriptFastaFile = $self->getProperty("transcriptFastaFile");
  my $bowtieBinDir = $self->getProperty("bowtieBinDir");
  my $perlScriptsDir = $self->getProperty("perlScriptsDir");
  my $genomeBowtieIndex = $self->getProperty("genomeBowtieIndex");
  my $transcriptBowtieIndex = $self->getProperty("transcriptBowtieIndex");
  my $createSAMFile = $self->getProperty("createSAMFile");
  
  
  die "readFilePath $readFilePath does not exist" unless -e "$readFilePath";
  die "pairedReadFilePath $pairedReadFilePath does not exist" if $pairedReadFilePath != 'none' && !(-e "$pairedReadFilePath");
  die "readFilePath equals pairedReadFilePath" if $pairedReadFilePath != 'none' && $pairedReadFilePath eq $readFilePath;
  die "genomeFastaFile $genomeFastaFile does not exist" unless -e "$genomeFastaFile";
#  die "--transcriptFastaFile or --transcriptBowtieIndex must be provided" unless -e "$transcriptFastaFile" || -e "$transcriptBowtieIndex.1.ebwt";
  die "bowtieBinDir $bowtieBinDir does not exist" unless -e "$bowtieBinDir";
  die "perlScriptsDir $perlScriptsDir does not exist" unless -e "$perlScriptsDir";

  ##if aren't creating SAM file then don't need to create qual file or indices
  if($createSAMFile =~ /true/i){
    $self->{quals} = 1;
  }else{
    $self->{quals} = 0;
  }
  
  ##make fasta file and qual files .. note that will chunk them up here as well.
  if(!(-d "$inputDir/subtasks")){
    mkdir("$inputDir/subtasks");
  }
  if(!(-e "$inputDir/subtasks/sequenceCount")){ ##file will exist if have already done this
    print "parsing reads file(s) to fasta and qual files\n";
    &runCmd("perl $perlScriptsDir/parse2fasta.pl $readFilePath".(-e $pairedReadFilePath ? " $pairedReadFilePath" : "")." | makeSubtasksFromFAFile.pl --stdin --outStem reads --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize}");
    &runCmd("perl $perlScriptsDir/fastq2qualities.pl $readFilePath".(-e $pairedReadFilePath ? " $pairedReadFilePath" : "")." | makeSubtasksFromFAFile.pl --stdin --outStem quals --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize}") if $self->{quals};
  }

  $self->{pairedEnd} = (-e "$pairedReadFilePath") ? "paired" : "single";

  ##set the count .. inputSetSize
  $self->{inputSetSize} = `cat $inputDir/subtasks/sequenceCount`;
  chomp $self->{inputSetSize};

  ##now check to see if have valid quals file
  if($self->{quals}){
    my $X = `head -2 $inputDir/subtasks/quals.1 | tail -1`;
    if($X =~ /Sorry, can't figure/) {
      $self->{quals} = 0;
    }else{
      $self->{quals} = 1;
    }
  }

  ##get read length
  my $read = `head -2 $inputDir/subtasks/reads.1  | tail -1`;
  chomp $read;
  $self->{readLength} = length($read);

  ## make bowtie genome index .. could check to see if it exists first
  ## obviously don't need to do this if pass in bowtie index.
  if(-e "$genomeBowtieIndex.1.ebwt"){
    $self->{bowtie_genome} = $genomeBowtieIndex;
    $self->{bowtie_genome_name} = basename($genomeBowtieIndex);
  }else{
    my $gdir = dirname($genomeFastaFile);
    $self->{bowtie_genome_name} = basename($genomeFastaFile);
    $self->{bowtie_genome} = "$gdir/$self->{bowtie_genome_name}";
    my @lsg = `ls -rt $genomeFastaFile $self->{bowtie_genome}.1.ebwt`;
    map { chomp } @lsg;
    if (scalar(@lsg) != 2 || $lsg[0] ne $genomeFastaFile) { 
      print "Building bowtie index for genome\n";
      &runCmd("$bowtieBinDir/bowtie-build $genomeFastaFile $self->{bowtie_genome}");
    }
  }

  ## make bowtie transcript index if need be .. could check to see if it exists first
  if(-e "$transcriptFastaFile" || -e "$transcriptBowtieIndex.1.ebwt"){
    if(-e "$transcriptBowtieIndex.1.ebwt"){
      $self->{bowtie_transcript} = $transcriptBowtieIndex;
      $self->{bowtie_transcript_name} = basename($transcriptBowtieIndex);
    }else{
      my $tdir = dirname($transcriptFastaFile);
      $self->{bowtie_transcript_name} = basename($transcriptFastaFile);
      $self->{bowtie_transcript} = "$tdir/$self->{bowtie_transcript_name}";
      my @lst = `ls -rt $transcriptFastaFile $self->{bowtie_transcript}.1.ebwt`;
      map { chomp } @lst;
      if (scalar(@lst) != 2 || $lst[0] ne $transcriptFastaFile) { 
        print "Building bowtie index for transcriptome\n";
        &runCmd("$bowtieBinDir/bowtie-build $transcriptFastaFile $self->{bowtie_transcript}");
      }
    }
  }

  ##get and report size ....
   print "Finding input set size\n";
  $self->{size} = $self->getInputSetSize($inputDir);
  if ($self->{size} == 0) {
    print "Error: Input set size is 0\n";
    exit 1;
  }
  print "Input set size is $self->{size}\n";
  my $c = int($self->{size} / $self->{subTaskSize});
  $c += 1 if $self->{size} % $self->{subTaskSize};
  print "Subtask set size is $self->{subTaskSize} ($c subtasks)\n";
}

sub initNode {
  my ($self, $node, $inputDir) = @_;
  
  my $genomeFastaFile = $self->getProperty("genomeFastaFile");
  my $nodeDir = $node->getDir();
  $node->runCmd("cp $genomeFastaFile $nodeDir");
  $node->runCmd("cp $self->{bowtie_genome}.* $nodeDir");
  $node->runCmd("cp $self->{bowtie_transcript}.* $nodeDir");
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    return $self->{inputSetSize};
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir) = @_;

    my $subtaskNum = int($start / $self->{subTaskSize}) + 1;
    $node->runCmd("cp $inputDir/subtasks/reads.$subtaskNum $nodeExecDir/seqSubset.fa");
    $node->runCmd("cp $inputDir/subtasks/quals.$subtaskNum $nodeExecDir/qualsSubset.fa") if $self->{quals};

}

sub makeSubTaskCommand { 
  my ($self, $node, $inputDir, $nodeExecDir) = @_;
    
  if($self->{subtaskCmd}){
    return $self->{subtaskCmd};
  }else{
    my $genomeFastaFile = $node->getDir() . "/" . basename($self->getProperty("genomeFastaFile"));
    my $bowtieBinDir = $self->getProperty("bowtieBinDir");
    my $perlScriptsDir = $self->getProperty("perlScriptsDir");
    my $genomeBowtieIndex =  $node->getDir() . "/" . basename($self->{bowtie_genome_name});
    my $transcriptBowtieIndex = $node->getDir() . "/" . basename($self->{bowtie_transcript_name});
    my $geneAnnotationFile = $self->getProperty("geneAnnotationFile");
    my $blatExec = $self->getProperty("blatExec");
    my $mdustExec = $self->getProperty("mdustExec");
    my $limitNU = $self->getProperty("limitNU");
    my $minBlatIdentity = $self->getProperty("minBlatIdentity");
    my $numInsertions = $self->getProperty("numInsertions");
    my $createSAMFile = $self->getProperty("createSAMFile");
    my $countMismatches = $self->getProperty("countMismatches");
    
    my $cmd =  "runRUMOnNode.pl --readsFile seqSubset.fa --qualFile ".($self->{quals} ? "qualsSubset.fa" : "none")." --genomeFastaFile $genomeFastaFile --genomeBowtieIndex $genomeBowtieIndex".($self->{bowtie_transcript} ? " --transcriptBowtieIndex $transcriptBowtieIndex" : "")." --geneAnnotationFile $geneAnnotationFile --bowtieExec $bowtieBinDir/bowtie --blatExec $blatExec --mdustExec $mdustExec --perlScriptsDir $perlScriptsDir --limitNU $limitNU --pairedEnd $self->{pairedEnd} --minBlatIdentity $minBlatIdentity --numInsertions $numInsertions --createSAMFile ".($createSAMFile =~ /true/i ? "1" : "0")." --countMismatches ".($countMismatches =~ /true/i ? "1" : "0");
    $self->{subtaskCmd} = $cmd;
    return $cmd;
  }
}

## bring back alignment and sam files and append subtasknum.  then at end concatenate.
sub integrateSubTaskResults {
  my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
  $node->runCmd("cp $nodeExecDir/RUM_Unique $mainResultDir/RUM_Unique.$subTaskNum");
  $node->runCmd("cp $nodeExecDir/RUM_NU $mainResultDir/RUM_NU.$subTaskNum");
  $node->runCmd("cp $nodeExecDir/RUM_sam $mainResultDir/RUM_sam.$subTaskNum") if $self->getProperty("createSAMFile") =~ /true/i;
}

## concatenate files here so are in order
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir) = @_;
  my $currDir = `pwd`;
  chomp $currDir;
  chdir("$mainResultDir") || die "$!";
  print STDERR "Concatenating RUM_Unique files\n";
  foreach my $f ($self->sortResultFiles('RUM_Unique.*')){
#    print STDERR "  Addng $f\n";
    &runCmd("cat $f >> RUM_Unique.all");
    unlink($f);
  }

  print STDERR "Concatenating RUM_NU files\n";
  foreach my $f ($self->sortResultFiles('RUM_NU.*')){
    &runCmd("cat $f >> RUM_NU.all");
    unlink($f);
  }

  if($self->getProperty("createSAMFile") =~ /true/i){
    print STDERR "Concatenating RUM_sam files\n";
    foreach my $f ($self->sortResultFiles('RUM_sam.*')){
      &runCmd("cat $f >> RUM_sam.all");
      unlink($f);
    }
  }
  chdir("$currDir") || die "$!";
  

  return 1;
}

sub sortResultFiles {
  my ($self,$fn) = @_;  ##takes in string to do ls with
  my @files = `ls $fn`;
  my @tmp;
  foreach my $f (@files){
    next if $f  =~ /all$/;
    chomp $f;
    next if $f  =~ /all$/;
    if($f =~ /\.(\d+)$/){
      push(@tmp,[$f,$1]);
    }else{
      die "ERROR sorting result files";
    }
  }
  my @sort;
  foreach my $a (sort{$a->[1] <=> $b->[1]}@tmp){
    push(@sort,$a->[0]);
  }
  return @sort;
}

1;
