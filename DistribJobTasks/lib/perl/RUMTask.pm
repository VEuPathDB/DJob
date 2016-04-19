package DJob::DistribJobTasks::RUMTask;
use DJob::DistribJob::Task;
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
 ["geneAnnotationFile",   "none",     "geneAnnotationFile in Greg's format (ucsc format)"],
 ["bowtieBinDir",   "default",     "bowtie bin directory"],
 ["blatExec",   "",     "blat executable, must be full path unless defined in your path"],
 ["mdustExec",   "",     "mdust executable, must be full path unless defined in your path"],
 ["perlScriptsDir",   "$ENV{GUS_HOME}/bin",     "directory for RUM scripts"],
 ["limitNU",   "30",     "Limits the number of ambiguous mappers to a max of [30]"],
 ["numInsertions",   "1",     "number of instertions to allow when parsing BLAT [1] (note, if paired end data the number of insertions is constrained to 0 or 1"],
 ["minBlatIdentity",   "93",     "run BLAT with minimum identity of [93]"],
 ["minlength",   "false",     "report only alignments this long or longer [false]"],
 ["createSAMFile",   "false",     "create SAM file if 1 ([false] | true)"],
 ["strandSpecific",   "false",     "data is strand specific if 1 ([false] | true)"],
 ["createJunctionsFile", "false", "create Juncctions file if 1 ([false] | true)"],
 ["countMismatches",   "false",     "report in the last column the number of mismatches, ignoring insertions ([false] | true)"],
 ["SNPs",   "false",     "run snp finder ([false] | true)"], ##not currently supported
 ["variableLengthReads",   "false",     "reads have variable lengths [false] | true)"],
 ["postProcess",   "false",     "run Post Processing steps ([false] | true)"],
 ["saveIntermediateFiles",   "false",     "copy back all intermediate files to mainResultDir ([false] | true)"]
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
    my $variableLengthReads = $self->getProperty("variableLengthReads");
    
    die "readFilePath $readFilePath does not exist" unless -e "$readFilePath";
    die "pairedReadFilePath $pairedReadFilePath does not exist" if $pairedReadFilePath ne 'none' && !(-e "$pairedReadFilePath");
    die "transcriptFastaFile $transcriptFastaFile does not exist" if $transcriptFastaFile != 'none' && !(-e "$transcriptFastaFile");
    die "readFilePath equals pairedReadFilePath" if $pairedReadFilePath != 'none' && $pairedReadFilePath eq $readFilePath;
    die "genomeFastaFile $genomeFastaFile does not exist" unless -e "$genomeFastaFile";
#  die "--transcriptFastaFile or --transcriptBowtieIndex must be provided" unless -e "$transcriptFastaFile" || -e "$transcriptBowtieIndex.1.ebwt";

    if($bowtieBinDir ne 'default'){
      die "bowtieBinDir $bowtieBinDir does not exist" unless -e "$bowtieBinDir";
      $ENV{PATH} .= ":$bowtieBinDir";
    }
    die "perlScriptsDir $perlScriptsDir does not exist" unless -e "$perlScriptsDir";
    
    ##if aren't creating SAM file then don't need to create qual file
    if($createSAMFile =~ /true/i){
	$self->{quals} = 1;
    }else{
	$self->{quals} = 0;
    }
    
    ##make fasta file and qual files .. note that will chunk them up here as well.
    my $date;
    if(!(-d "$inputDir/subtasks")){
	mkdir("$inputDir/subtasks");
    }
    my $madeReadFile = 0;
  REDOSUBTASKS:
    if(!(-e "$inputDir/subtasks/sequenceCount")){ ##file will exist if have already done this
      $date = `date`; chomp $date;
      print "[$date] parsing reads file(s) to fasta and qual files\n";
      $self->{nodeForInit}->runCmd("perl $perlScriptsDir/fastq2qualities.pl $readFilePath".(-e $pairedReadFilePath ? " $pairedReadFilePath" : "")." | makeSubtasksFromFAFileForRUM.pl --stdin --outStem quals --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize}") if $self->{quals};
      $self->{nodeForInit}->runCmd("perl $perlScriptsDir/parse2fasta.pl $readFilePath".(-e $pairedReadFilePath ? " $pairedReadFilePath" : "")." | makeSubtasksFromFAFileForRUM.pl --stdin --outStem reads --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize}");
      $madeReadFile = 1;
    }

    ##set the count .. inputSetSize
    $self->{inputSetSize} = `cat $inputDir/subtasks/sequenceCount`;
    chomp $self->{inputSetSize};

    ##get and report size ....
    $self->{size} = $self->getInputSetSize($inputDir);
    if ($self->{size} == 0) {
	print "Error: Input set size is 0\n";
	exit 1;
    }
    my $c = int($self->{size} / $self->{subTaskSize});
    $c += 1 if $self->{size} % $self->{subTaskSize};
    $self->{numSubtasks} = $c;

    ##check to see if the the number of files created is consistent with number of subtasks
    ## if not, make the files above because subtask size has changed;
    ## not going to test for change in inputfile as if user is silly enough to change insitu they'll have to live with it
    ## intended to catch case where subtask size has been changed.
    my @sFiles = glob("$inputDir/subtasks/reads.*");
    if(scalar(@sFiles) != $self->{numSubtasks}){
      die "Unable to create subtask reads files properly\n" if $madeReadFile;
      system("/bin/rm $inputDir/subtasks/*");
      goto REDOSUBTASKS;
    }
    print "Input set size is $self->{size}\n";
    print "Subtask set size is $self->{subTaskSize} ($self->{numSubtasks} subtasks)\n";
    
    $self->{pairedEnd} = (-e "$pairedReadFilePath") ? "paired" : "single";
    
    
    ##now check to see if have valid quals file
    if($self->{quals}){
	my $X = `head -2 $inputDir/subtasks/quals.1 | tail -1`;
	if($X =~ /Sorry, can.t figure/) {
	    $self->{quals} = 0;
	} else {
	    $self->{quals} = 1;
	}
    }

    ##get read length
    my $variableLengthReads = $self->getProperty("variableLengthReads");
    if($variableLengthReads eq 'true') {
	$self->{readLength} = "v";
    } else {
	my $read = `head -2 $inputDir/subtasks/reads.1  | tail -1`;
	chomp $read;
	$self->{readLength} = length($read);
    }
    
    if($self->getProperty("minlength") == 'false') {
	if($self->{readLength} < 80) {
	    if($self->{match_length_cutoff} == 0) {
		$self->{match_length_cutoff} = 35;
	    }
	} else {
	    if($self->{match_length_cutoff} == 0) {
		$self->{match_length_cutoff} = 50;
	    }
	}
	if($self->{match_length_cutoff} >= .8 * $self->{readLength}) {
	    if($self->{match_length_cutoff} == 0) {
		$self->{match_length_cutoff} = int(.6 * $self->{readLength});
	    }
	}
    } else {
        $self->{match_length_cutoff} = $self->{minlength};
    }
    my $match_length_cutoff = $self->{match_length_cutoff};
    
    
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
          $date = `date`; chomp $date;
	    print "$date: Building bowtie index for genome\n";
	    $self->{nodeForInit}->runCmd("bowtie-build $genomeFastaFile $self->{bowtie_genome}");
	}
    }
    
    ## make bowtie transcript index if need be .. could check to see if it exists first
    ## obviously don't need to do this if pass in bowtie index.
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
              $date = `date`; chomp $date;
		print "$date: Building bowtie index for transcriptome\n";
		$self->{nodeForInit}->runCmd("bowtie-build $transcriptFastaFile $self->{bowtie_transcript}");
	    }
	}
    }
    

    
    $date = `date`;
    open(LOGFILE, ">rum.log_master");
    print LOGFILE "\nstart: $date\n";
    if($variableLengthReads eq "false") {
        my $readlength = $self->{readLength};
	print LOGFILE "readlength: $readlength\n";
    } else {
	print LOGFILE "readlength: variable\n";
    }
    
    my $minlength = $self->{minlength};
    if($variableLengthReads eq "false" || $minlength > 0) {
	if($minlength == 0) {
	    print LOGFILE "minimum length alignment to report: $match_length_cutoff\n  *** NOTE: If you want shorter alignments reported, use the -minlength option.\n";
	} else {
	    print LOGFILE "minimum length alignment to report: $match_length_cutoff.\n";
	}
	print "\n *** NOTE: I am going to report alginments of length $match_length_cutoff.\n";
	print "If you want shorter alignments to be reported, use the -minlength option.\n\n";
    } else {
	print LOGFILE "minimum length alignment to report: NA since read length is variable\n";
    }

    print LOGFILE "paired- or single-end: $self->{pairedEnd}\n";
    my $x = $self->getProperty("limitNU");
    print LOGFILE "limitNU: $x\n";
#    print LOGFILE "dna: $dna\n";

    
}

sub initNode {
  my ($self, $node, $inputDir) = @_;
  return 1;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    return $self->{inputSetSize};
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;

    my $subtaskNum = int($start / $self->{subTaskSize}) + 1;
    $self->runCmdOnNode("touch $inputDir/subtasks/reads.$subtaskNum.touch",1);
    $self->runCmdOnNode("/bin/rm $inputDir/subtasks/reads.$subtaskNum.touch",1);
    $self->runCmdOnNode("cp $inputDir/subtasks/reads.$subtaskNum $nodeExecDir/seqSubset.fa");
    $self->runCmdOnNode("cp $inputDir/subtasks/quals.$subtaskNum $nodeExecDir/qualsSubset.fa") if $self->{quals};

}

sub makeSubTaskCommand { 
  my ($self, $node, $inputDir, $nodeExecDir,$subtaskNumber,$mainResultDir) = @_;
    
  if(!$self->{subtaskCmd}){
    my $genomeFastaFile = $self->getProperty("genomeFastaFile");
    my $bowtieBinDir = $self->getProperty("bowtieBinDir");
    $bowtieBinDir = $bowtieBinDir eq 'default' ? "" : "$bowtieBinDir/";
    my $perlScriptsDir = $self->getProperty("perlScriptsDir");
    my $genomeBowtieIndex =  $self->{bowtie_genome};
    my $transcriptBowtieIndex = $self->{bowtie_transcript};
    my $geneAnnotationFile = $self->getProperty("geneAnnotationFile");
    my $blatExec = $self->getProperty("blatExec");
    my $mdustExec = $self->getProperty("mdustExec");
    my $limitNU = $self->getProperty("limitNU");
    my $minBlatIdentity = $self->getProperty("minBlatIdentity");
    my $numInsertions = $self->getProperty("numInsertions");
    my $createSAMFile = $self->getProperty("createSAMFile");
    my $countMismatches = $self->getProperty("countMismatches");
    my $minlength = $self->getProperty("minlength");
    my $matchLengthCutoff = $self->{"match_length_cutoff"};
    my $readlength = $self->{"readLength"};
    my $strandspecific = $self->getProperty("strandSpecific");

    my $cmd =  "runRUMOnNode.pl --readlength $readlength --readsFile seqSubset.fa --qualFile ".($self->{quals} ? "qualsSubset.fa" : "none")." --genomeFastaFile $genomeFastaFile --genomeBowtieIndex $genomeBowtieIndex".($self->{bowtie_transcript} ? " --transcriptBowtieIndex $transcriptBowtieIndex" : "").($geneAnnotationFile =~ /none$/ ? "" : " --geneAnnotationFile $geneAnnotationFile")." --bowtieExec $bowtieBinDir"."bowtie --blatExec $blatExec --mdustExec $mdustExec --perlScriptsDir $perlScriptsDir --limitNU $limitNU --pairedEnd $self->{pairedEnd} --minlength $minlength --minBlatIdentity $minBlatIdentity --numInsertions $numInsertions --createSAMFile ".($createSAMFile =~ /true/i ? "1" : "0")." --countMismatches ".($countMismatches =~ /true/i ? "1" : "0")." --mainResultDir $mainResultDir --matchLengthCutoff $matchLengthCutoff --strandSpecific $strandspecific";

    $self->{subtaskCmd} = $cmd;
  }
  return $self->{subtaskCmd} . " --subtaskNumber $subtaskNumber";
}

sub integrateSubTaskResults {
  my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
  if($self->getProperty("saveIntermediateFiles") =~ /true/i){
    mkdir("$mainResultDir/subtask.$subTaskNum");
    $self->runCmdOnNode("cp $nodeExecDir/* $mainResultDir/subtask.$subTaskNum");
  }
}

# do all post-processing here:
# a node is now passed in that can be used to run commands on a node using $self->runCmdOnNode("cmd")
# NOTE that in order to use this must add keepNodeForPostProcessing=yes to controller.prop file

sub cleanUpServer {
    my($self, $inputDir, $mainResultDir, $node) = @_;
    return 1 if $self->getProperty('postProcess') ne 'true';


    my $genomeFastaFile = $self->getProperty("genomeFastaFile");
    my $geneAnnotationFile = $self->getProperty("geneAnnotationFile");
    my $perlScriptsDir = $self->getProperty("perlScriptsDir");

    $self->runCmdOnNode("postProcessRUMTask --genomeFastaFile $genomeFastaFile --geneAnnotationFile $geneAnnotationFile --mainResultDir $mainResultDir --perlScriptsDir $perlScriptsDir --createJunctions".($self->{bowtie_transcript} ? " --haveTranscripts" : ""));
    
}

1;
