package DJob::DistribJobTasks::RUMTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFile;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["readFilePath",   "",     "full path to read file"],
 ["pairedReadFilePath",   "",     "full path to paired read file (optional)"],
 ["genomeFastaFile",   "",     "genome file in fasta format"],
 ["genomeBowtieIndex",   "",     "genome bowtie index file"],
 ["transcriptFastaFile",   "",     "transcript file in fasta format"],
 ["transcriptBowtieIndex",   "",     "transcript bowtie index files"],
 ["geneAnnotationFile",   "",     "geneAnnotationFile in Gregs format (ucsc format)"],
 ["bowtieBinDir",   "",     "bowtie bin directory"],
 ["blatExec",   "",     "blat executable, must be full path unless defined in your path"],
 ["mdustExec",   "",     "mdust executable, must be full path unless defined in your path"],
 ["perlScriptsDir",   "",     "directory for RUM scripts"],
 ["limitNU",   "30",     "Limits the number of ambiguous mappers to a max of [30]"],
 ["numInsertions",   "1",     "number of instertions to allow when parsing BLAT [1] (note, if paired end data the number of insertions is constrained to 0 or 1"],
 ["minBlatIdentity",   "93",     "run BLAT with minimum identity of [93]"],
 ["createSAMfile",   "0",     "create SAM file if 1 ([0] | 1)"],
 ["countMismatches",   "0",     "report in the last column the number of mismatches, ignoring insertions ([0] | 1)"]
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
  
  
  die "readFilePath $readFilePath does not exist" unless -e "$readFilePath";
  die "pairedReadFilePath $pairedReadFilePath does not exist" if $pairedReadFilePath && !(-e "$pairedReadFilePath");
  die "readFilePath equals pairedReadFilePath" if $pairedReadFilePath && $pairedReadFilePath eq $readFilePath;
  die "genomeFastaFile $genomeFastaFile does not exist" unless -e "$genomeFastaFile";
  die "--transcriptFastaFile or --transcriptBowtieIndex must be provided" unless -e "$transcriptFastaFile" || -e "$transcriptBowtieIndex.1.ebwt";
  die "bowtieBinDir $bowtieBinDir does not exist" unless -e "$bowtieBinDir";
  die "perlScriptsDir $perlScriptsDir does not exist" unless -e "$perlScriptsDir";
  
  ##make fasta file and qual file
  $self->{"reads_fa"} = "$inputDir/reads.fa";
  $self->{"quals_fa"} = "$inputDir/quals.fa";
  ## NOTE: may have already done this if restarting so check to see if have files and they are newer than input
  my @ls = `ls -rt $readFilePath $self->{reads_fa}`;
  map { chomp } @ls;
  if (scalar(@ls) != 2 || $ls[0] ne $readFilePath) { 
    if($pairedReadFilePath){  ## doing paired reads
      &runCmd("perl $perlScriptsDir/parse2fasta.pl $readFilePath $pairedReadFilePath > $self->{reads_fa}");
      &runCmd("perl $perlScriptsDir/fastq2qualities.pl $readFilePath $pairedReadFilePath > $self->{quals_fa}");
    }else{
      &runCmd("perl $perlScriptsDir/parse2fasta.pl $readFilePath > $self->{reads_fa}");
      &runCmd("perl $perlScriptsDir/fastq2qualities.pl $readFilePath > $self->{quals_fa}");
    }
  }

  $self->{pairedEnd} = $pairedReadFilePath ? "paired" : "single";

  ##now check to see if have valid quals file
  my $X = `head -2 $self->{quals_fa} | tail -1`;
  if($X =~ /Sorry, can't figure/) {
    $self->{quals_fa} = 'none';
  }    

  ##get read length
  my $read = `head -2 $self->{reads_fa} | tail -1`;
  chomp $read;
  $self->{readLength} = length($read);

  ## make bowtie genome index .. could check to see if it exists first
  ## obviously don't need to do this if pass in bowtie index.
  if(-e "$genomeBowtieIndex.1.ebwt"){
    $self->{bowtie_genome} = $genomeBowtieIndex;
  }else{
    $self->{bowtie_genome} = "$inputDir/bowtie_genome";
    my @lsg = `ls -rt $genomeFastaFile $self->{bowtie_genome}.1.ebwt`;
    map { chomp } @lsg;
    if (scalar(@lsg) != 2 || $lsg[0] ne $genomeFastaFile) { 
      &runCmd("$bowtieBinDir/bowtie_build $genomeFastaFile $self->{bowtie_genome}");
    }
  }

  ## make bowtie transcript index if need be .. could check to see if it exists first
  if(-e "$transcriptFastaFile" || -e "$transcriptBowtieIndex.1.ebwt"){
    if(-e "$transcriptGenomeIndex.1.ebwt"){
      $self->{bowtie_transcript} = $transcriptBowtieIndex;
    }else{
      $self->{bowtie_transcript} = "$inputDir/bowtie_transcript";
      my @lst = `ls -rt $transcriptFastaFile $self->{bowtie_transcript}.1.ebwt`;
      map { chomp } @lst;
      if (scalar(@lst) != 2 || $lst[0] ne $transcriptFastaFile) { 
        &runCmd("$bowtieBinDir/bowtie_build $transcriptFastaFile $self->{bowtie_transcript}");
      }
    }
  }

  ##create indices needed to write subsets to work on
  $self->{readsIndex} = CBIL::Bio::FastaFile->new($self->{reads_fa});
  $self->{qualsIndex} = CBIL::Bio::FastaFile->new($self->{quals_fa}) if $self->{quals};
  if($self->{quals} && $self->{readsIndex}->getCount() != $self->{qualsIndex}->getCount()){
    die "ERROR: different number of quals from reads\n";
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

    ## could leave indices on server .. if so do nothing
    ## otherwise should copy all indices to nodeDir

}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    return $self->{readsIndex}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir) = @_;

    $self->{readsIndex}->writeSeqsToFile($start, $end, "$serverSubTaskDir/seqSubset.fa");
    ##also need to write the quals if using  these
    $self->{qualsIndex}->writeSeqsToFile($start, $end, "$serverSubTaskDir/qualsSubset.fa") if $self->{quals};

    $node->runCmd("cp -r $serverSubTaskDir/* $nodeExecDir");
}

sub makeSubTaskCommand { 
  my ($self, $node, $inputDir, $nodeExecDir) = @_;
    
  my $genomeFastaFile = $self->getProperty("genomeFastaFile");
  my $transcriptFastaFile = $self->getProperty("transcriptFastaFile");
  my $bowtieBinDir = $self->getProperty("bowtieBinDir");
  my $perlScriptsDir = $self->getProperty("perlScriptsDir");
  my $genomeBowtieIndex = $self->getProperty("genomeBowtieIndex");
  my $transcriptBowtieIndex = $self->getProperty("transcriptBowtieIndex");
  my $geneAnnotationFile = $self->getProperty("geneAnnotationFile");
  my $blatExec = $self->getProperty("blatExec");
  my $mdustExec = $self->getProperty("mdustExec");
  my $limitNU = $self->getProperty("limitNU");
  my $minBlatIdentity = $self->getProperty("minBlatIdentity");
  my $numInsertions = $self->getProperty("numInsertions");
  my $createSAMFile = $self->getProperty("createSAMFile");
  my $countMismatches = $self->getProperty("countMismatches");

    my $cmd =  "runRUMOnNode.pl --readsFile seqSubset.fa --qualFile qualsSubset.fa --genomeFastaFile $genomeFastaFile --genomeBowtieIndex $self->{bowtie_genome}".($self->{bowtie_transcript} ? " --transcriptBowtieIndex $self->{bowtie_transcript}" : "")." --geneAnnotationFile $geneAnnotationFile --bowtieExec $bowtieBinDir/bowtie --blatExec $blatExec --mdustExec $mdustExec --perlScriptsDir $perlScriptsDir --limitNU $limitNU --pairedEnd $self->{pairedEnd} --minBlatIdentity $minBlatIdentity --numInsertions $numInsertions --createSAMFile $createSAMFile --countMismatches $countMismatches";

    return $cmd;
}

## bring back alignment and sam files and append subtasknum.  then at end concatenate.
sub integrateSubTaskResults {
  my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
  $node->runCmd("mv $nodeExecDir/RUM_Unique $mainResultDir/RUM_Unique.$subTaskNum");
  $node->runCmd("mv $nodeExecDir/RUM_NU $mainResultDir/RUM_NU.$subTaskNum");
  $node->runCmd("mv $nodeExecDir/RUM_sam $mainResultDir/RUM_sam.$subTaskNum") if $self->getProperty("createSAMFile");
}

## concatenate files here so are in order
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir) = @_;
  print "Concatenating RUM_Unique files\n";
  foreach my $f ($self->sortResultFiles('RUM_Unique.*')){
    &runCmd("cat $f >> RUM_Unique.all");
    &runCmd("rm $f");
  }

  print "Concatenating RUM_NU files\n";
  foreach my $f ($self->sortResultFiles('RUM_NU.*')){
    &runCmd("cat $f >> RUM_NU.all");
    &runCmd("rm $f");
  }

  if($self->getProperty("createSAMFile"){
    print "Concatenating RUM_sam files\n";
    foreach my $f ($self->sortResultFiles('RUM_sam.*')){
      &runCmd("cat $f >> RUM_sam.all");
      &runCmd("rm $f");
    }
  }

  return 1;
}

sub sortResultFiles {
  my ($self,$fn) = @_;  ##takes in string to do ls with
  my @files = `ls $fn`;
  my @tmp;
  foreach my $f (@files){
    chomp;
    if($f =~ /\.(\d+)$/){
      push(@tmp,[$f,$1]);
    }else{
      die "ERROR sorting result files";
    }
  }
  my @sort;
  foreach my $a (sort{$a->[1] <=> $b->[1]}@tmp){
    push @sort($a->[0]);
  }
  return @sort;
}

1;
