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
 ["bowtieBinDir",   "",     "bowtie bin directory"],
 ["blatExec",   "",     "blat executable, must be full path unless defined in your path"],
 ["mdustExec",   "",     "mdust executable, must be full path unless defined in your path"],
 ["perlScriptsDir",   "",     "directory for RUM scripts"],
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
 ["postProcess",   "true",     "run Post Processing steps false | [true])"],
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
    die "pairedReadFilePath $pairedReadFilePath does not exist" if $pairedReadFilePath != 'none' && !(-e "$pairedReadFilePath");
    die "readFilePath equals pairedReadFilePath" if $pairedReadFilePath != 'none' && $pairedReadFilePath eq $readFilePath;
    die "genomeFastaFile $genomeFastaFile does not exist" unless -e "$genomeFastaFile";
#  die "--transcriptFastaFile or --transcriptBowtieIndex must be provided" unless -e "$transcriptFastaFile" || -e "$transcriptBowtieIndex.1.ebwt";
    die "bowtieBinDir $bowtieBinDir does not exist" unless -e "$bowtieBinDir";
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
      &runCmd("perl $perlScriptsDir/fastq2qualities.pl $readFilePath".(-e $pairedReadFilePath ? " $pairedReadFilePath" : "")." | makeSubtasksFromFAFileForRUM.pl --stdin --outStem quals --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize}") if $self->{quals};
      &runCmd("perl $perlScriptsDir/parse2fasta.pl $readFilePath".(-e $pairedReadFilePath ? " $pairedReadFilePath" : "")." | makeSubtasksFromFAFileForRUM.pl --stdin --outStem reads --directory $inputDir/subtasks --subtaskSize $self->{subTaskSize}");
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
	    &runCmd("$bowtieBinDir/bowtie-build $genomeFastaFile $self->{bowtie_genome}");
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
		&runCmd("$bowtieBinDir/bowtie-build $transcriptFastaFile $self->{bowtie_transcript}");
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
  my $genomeFastaFile = $self->getProperty("genomeFastaFile");
  my $geneAnnotationFile = $self->getProperty("geneAnnotationFile");
  my $nodeDir = $node->getDir();
  $node->runCmd("cp $genomeFastaFile $nodeDir");
  $node->runCmd("cp $self->{bowtie_genome}.* $nodeDir");
  $node->runCmd("cp $self->{bowtie_transcript}.* $nodeDir") if $self->{bowtie_transcript};
  $node->runCmd("cp $geneAnnotationFile $nodeDir") if (-e "$geneAnnotationFile");
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    return $self->{inputSetSize};
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;

    my $subtaskNum = int($start / $self->{subTaskSize}) + 1;
    $node->runCmd("cp $inputDir/subtasks/reads.$subtaskNum $nodeExecDir/seqSubset.fa");
    $node->runCmd("cp $inputDir/subtasks/quals.$subtaskNum $nodeExecDir/qualsSubset.fa") if $self->{quals};

}

sub makeSubTaskCommand { 
  my ($self, $node, $inputDir, $nodeExecDir,$subtaskNumber,$mainResultDir) = @_;
    
  if(!$self->{subtaskCmd}){
    my $genomeFastaFile = $self->getProperty("genomeFastaFile");
    my $bowtieBinDir = $self->getProperty("bowtieBinDir");
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

    my $cmd =  "runRUMOnNode.pl --readlength $readlength --readsFile seqSubset.fa --qualFile ".($self->{quals} ? "qualsSubset.fa" : "none")." --genomeFastaFile $genomeFastaFile --genomeBowtieIndex $genomeBowtieIndex".($self->{bowtie_transcript} ? " --transcriptBowtieIndex $transcriptBowtieIndex" : "").($geneAnnotationFile =~ /none$/ ? "" : " --geneAnnotationFile $geneAnnotationFile")." --bowtieExec $bowtieBinDir/bowtie --blatExec $blatExec --mdustExec $mdustExec --perlScriptsDir $perlScriptsDir --limitNU $limitNU --pairedEnd $self->{pairedEnd} --minlength $minlength --minBlatIdentity $minBlatIdentity --numInsertions $numInsertions --createSAMFile ".($createSAMFile =~ /true/i ? "1" : "0")." --countMismatches ".($countMismatches =~ /true/i ? "1" : "0")." --mainResultDir $mainResultDir --matchLengthCutoff $matchLengthCutoff --strandSpecific $strandspecific";

    $self->{subtaskCmd} = $cmd;
  }
  return $self->{subtaskCmd} . " --subtaskNumber $subtaskNumber";
}

sub integrateSubTaskResults {
  my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
  if($self->getProperty("saveIntermediateFiles") =~ /true/i){
    mkdir("$mainResultDir/subtask.$subTaskNum");
    $node->runCmd("cp $nodeExecDir/* $mainResultDir/subtask.$subTaskNum");
  }
}

# do all post-processing here:
# a node is now passed in that can be used to run commands on a node using $node->runCmd("cmd")
# NOTE that in order to use this must add keepNodeForPostProcessing=yes to controller.prop file

sub cleanUpServer {
    my($self, $inputDir, $mainResultDir, $node) = @_;
    return 1 if $self->getProperty('postProcess') ne 'true';

    my $genomeFastaFile = $self->getProperty("genomeFastaFile");
    my $geneAnnotationFile = $self->getProperty("geneAnnotationFile");
    
    my $currDir = `pwd`;
    chomp $currDir;
    chdir("$mainResultDir") || die "$!";
    # concatenate files so are in order
    my $date;
    $date = `date`; chomp $date;
    print STDERR "[$date] Concatenating RUM_Unique files\n";
    my $numchunks = $self->{numSubtasks};
    if(!-e "$mainResultDir/RUM_Unique"){
      unlink("$mainResultDir/RUM_Unique.tmp");
      my @unique = $self->sortResultFiles('RUM_Unique.* | grep -v sorted');
      my $cnt_x = 0;
      if(scalar(@unique) != $self->{numSubtasks} && $cnt_x < 20){
	sleep(5);
	@unique = $self->sortResultFiles('RUM_Unique.* | grep -v sorted');
	$cnt_x++;
      }    
      if(scalar(@unique) != $self->{numSubtasks}){
	print STDERR "Not concatenating files: number of subtasks ($self->{numSubtasks}) differs from number of files (".scalar(@unique).") ($cnt_x)\n";
	return 1;
      }
      foreach my $f (@unique){
	print STDERR "  Adding $f\n";
	$node->runCmd("cat $mainResultDir/$f >> $mainResultDir/RUM_Unique.tmp");
      }
      sleep 10;
      &runCmd("/bin/mv $mainResultDir/RUM_Unique.tmp $mainResultDir/RUM_Unique");
      foreach my $f (@unique){
        unlink($f);
      }
    }else{
      print STDERR "  DONE\n";
    }

    $date = `date`; chomp $date;
    print STDERR "[$date] Concatenating RUM_NU files\n";
    if(!-e "$mainResultDir/RUM_NU"){
      unlink("$mainResultDir/RUM_NU.tmp");
      my @nu = $self->sortResultFiles('RUM_NU.* | grep -v sorted');
      foreach my $f (@nu){
	print STDERR "  Adding $f\n";
	$node->runCmd("cat $mainResultDir/$f >> $mainResultDir/RUM_NU.tmp");
	unlink($f);
      }
      sleep 10;
      &runCmd("/bin/mv $mainResultDir/RUM_NU.tmp $mainResultDir/RUM_NU");
      foreach my $f (@nu){
        unlink($f);
      }
    }else{
      print STDERR "  DONE\n";
    }

    $date = `date`; chomp $date;
    print STDERR "[$date] Concatenating RUM.sam files\n";
    if(-e "$mainResultDir/RUM_sam.1"){  #input file still exists so need to run.
      unlink("$mainResultDir/RUM.sam_tmp");
      my %samheader;
      my @samheaders = $self->sortResultFiles('sam_header.*');
      foreach my $f (@samheaders){
        open(SAMHEADER, $f);
        while(my $line = <SAMHEADER>) {
          chomp($line);
          $line =~ /SN:([^\s]+)\s/;
          $samheader{$1}=$line;
        }
        close(SAMHEADER);
      }
      open(SAMOUT, ">RUM.sam_tmp");
      foreach my $key (sort {cmpChrs($a,$b)} keys %samheader) {
        my $shout = $samheader{$key};
        print SAMOUT "$shout\n";
      }
      close(SAMOUT);
      foreach my $f ($self->sortResultFiles('RUM_sam.*')){
        print STDERR "  Adding $f\n";
        $node->runCmd("cat $mainResultDir/$f >> $mainResultDir/RUM.sam_tmp");
      }
      sleep 10;
      &runCmd("/bin/mv $mainResultDir/RUM.sam_tmp $mainResultDir/RUM.sam");
      foreach my $f ($self->sortResultFiles('RUM_sam.*')){
        unlink($f);
      }
    }
    # I did the following because sometimes the file RUM.sam isn't written
# by the time it reaches this, because of occasional turgid file system on the cluster.
    my $flag_x = 0;
    my $cnt_x = 0;
    my $NumSeqs;
    while($flag_x == 0) {
      $cnt_x++;
      if($cnt_x > 100) {
        $flag_x = 1;
      }
      my $X = `tail -1 $mainResultDir/RUM.sam`;
      $X =~ /^seq.(\d+)/;
      $NumSeqs = $1;
      if($NumSeqs =~ /^\d+$/) {
        $flag_x = 1;
      } else {
        sleep(5);
      }
    }
    $date = `date`; chomp $date;
    print STDERR "[$date] running samtools on file to create sorted bam file\n";
    if(!-e "$mainResultDir/RUM.bam"){
      unlink("$mainResultDir/RUM.bam_tmp");
      $node->runCmd("samtools view -b -S $mainResultDir/RUM.sam > $mainResultDir/RUM.bam_tmp");
      &runCmd("/bin/mv $mainResultDir/RUM.bam_tmp $mainResultDir/RUM.bam");
    }
    if(!-e "$mainResultDir/RUM.sorted.bam"){
      unlink("$mainResultDir/RUM.sorted_tmp.bam");
      $node->runCmd("samtools sort $mainResultDir/RUM.bam $mainResultDir/RUM.sorted_tmp");
      &runCmd("/bin/mv $mainResultDir/RUM.sorted_tmp.bam $mainResultDir/RUM.sorted.bam");
    }
    if(!-e "$mainResultDir/RUM.sorted.bam.bai"){
      unlink("$mainResultDir/RUM.sorted.bam_tmp.bai");
      $node->runCmd("samtools index $mainResultDir/RUM.sorted.bam $mainResultDir/RUM.sorted.bam_tmp.bai");
      &runCmd("/bin/mv $mainResultDir/RUM.sorted.bam_tmp.bai $mainResultDir/RUM.sorted.bam.bai");
    }

    my $perlScriptsDir = $self->getProperty("perlScriptsDir");
    $date = `date`; chomp $date;
    print STDERR "[$date] Counting mapped reads\n";
    if(!-e "$mainResultDir/mapping_stats.txt"){
      unlink("$mainResultDir/mapping_stats.txt.tmp");
      $node->runCmd("perl $perlScriptsDir/count_reads_mapped.pl $mainResultDir/RUM_Unique $mainResultDir/RUM_NU -minseq 1 -maxseq $NumSeqs > $mainResultDir/mapping_stats.txt.tmp 2>> $mainResultDir/PostProcessing-errorlog");
      &runCmd("/bin/mv $mainResultDir/mapping_stats.txt.tmp $mainResultDir/mapping_stats.txt");
    }

## --- good to here --- ###
    my $transcriptBowtieIndex = $self->{bowtie_transcript};
    if($transcriptBowtieIndex ne 'none') {
      $date = `date`; chomp $date;
      print STDERR "[$date] Merging blat and bowtie quantifications\n";
      if($self->getProperty("strandSpecific") eq 'true') {
        if(!-e "$mainResultDir/feature_quantifications.ps"){
          unlink("$mainResultDir/feature_quantifications.ps.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications.ps.tmp -strand ps -countsonly 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications.ps.tmp $mainResultDir/feature_quantifications.ps");
        }
        if(!-e "$mainResultDir/feature_quantifications_normalized.ps"){
          unlink("$mainResultDir/feature_quantifications_normalized.ps.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications_normalized.ps.tmp -strand ps 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications_normalized.ps.tmp $mainResultDir/feature_quantifications_normalized.ps");
        }
        if(!-e "$mainResultDir/feature_quantifications.ms"){
          unlink("$mainResultDir/feature_quantifications.ms.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications.ms.tmp -strand ms -countsonly 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications.ms.tmp $mainResultDir/feature_quantifications.ms");
        }
        if(!-e "$mainResultDir/feature_quantifications_normalized.ms"){
          unlink("$mainResultDir/feature_quantifications_normalized.ms.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications_normalized.ms.tmp -strand ms 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications_normalized.ms.tmp $mainResultDir/feature_quantifications_normalized.ms");
        }
          
        if(!-e "$mainResultDir/feature_quantifications.pa"){
          unlink("$mainResultDir/feature_quantifications.pa.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications.pa.tmp -strand pa -countsonly 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications.pa.tmp $mainResultDir/feature_quantifications.pa");
        }
        if(!-e "$mainResultDir/feature_quantifications_normalized.pa"){
          unlink("$mainResultDir/feature_quantifications_normalized.pa.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications_normalized.pa.tmp -strand pa 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications_normalized.pa.tmp $mainResultDir/feature_quantifications_normalized.pa");
        }
          
        if(!-e "$mainResultDir/feature_quantifications.ma"){
          unlink("$mainResultDir/feature_quantifications.ma.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications.ma.tmp -strand ma -countsonly 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications.ma.tmp $mainResultDir/feature_quantifications.ma");
        }
        if(!-e "$mainResultDir/feature_quantifications_normalized.ma"){
          unlink("$mainResultDir/feature_quantifications_normalized.ma.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications_normalized.ma.tmp -strand ma 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications_normalized.ma.tmp $mainResultDir/feature_quantifications_normalized.ma");
        }
          
        if(!-e "$mainResultDir/feature_counts"){
          unlink("$mainResultDir/feature_counts.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants_strandspecific.pl $mainResultDir/feature_quantifications.ps $mainResultDir/feature_quantifications.ms $mainResultDir/feature_quantifications.pa $mainResultDir/feature_quantifications.ma $geneAnnotationFile $mainResultDir/feature_counts.tmp 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_counts.tmp $mainResultDir/feature_counts");
        }
        if(!-e "$mainResultDir/feature_quantifications"){
          unlink("$mainResultDir/feature_quantifications.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants_strandspecific.pl $mainResultDir/feature_quantifications_normalized.ps $mainResultDir/feature_quantifications_normalized.ms $mainResultDir/feature_quantifications_normalized.pa $mainResultDir/feature_quantifications_normalized.ma $geneAnnotationFile $mainResultDir/feature_quantifications.tmp 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications.tmp $mainResultDir/feature_quantifications");
        }
      } else {
        if(!-e "$mainResultDir/feature_counts"){
          unlink("$mainResultDir/feature_counts.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_counts.tmp -countsonly 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_counts.tmp $mainResultDir/feature_counts");
        }
        if(!-e "$mainResultDir/feature_quantifications"){
          unlink("$mainResultDir/feature_quantifications.tmp");
          $node->runCmd("perl $perlScriptsDir/merge_quants.pl $mainResultDir $numchunks $mainResultDir/feature_quantifications.tmp 2>> $mainResultDir/PostProcessing-errorlog");
          &runCmd("/bin/mv $mainResultDir/feature_quantifications.tmp $mainResultDir/feature_quantifications");
        }
      }
    }

    my $string;
    $date = `date`; chomp $date;
    print STDERR "[$date] Merging sorted RUM files\n";
    if(!-e "$mainResultDir/RUM_Unique.sorted"){
      unlink("$mainResultDir/RUM_Unique.sorted_tmp");
      $string = "$mainResultDir/RUM_Unique.sorted_tmp";
      for(my $j=1; $j<$numchunks+1; $j++) {
	$string = $string . " $mainResultDir/RUM_Unique.sorted.$j";
      }
      $node->runCmd("perl $perlScriptsDir/merge_sorted_RUM_files.pl $string 2>> $mainResultDir/PostProcessing-errorlog");
      &runCmd("/bin/mv $mainResultDir/RUM_Unique.sorted_tmp $mainResultDir/RUM_Unique.sorted");
    }
    if(!-e "$mainResultDir/RUM_NU.sorted"){
      unlink("$mainResultDir/RUM_NU.sorted_tmp");
      $string = "$mainResultDir/RUM_NU.sorted_tmp";
      for(my $j=1; $j<$numchunks+1; $j++) {
	$string = $string . " $mainResultDir/RUM_NU.sorted.$j";
      }
      $node->runCmd("perl $perlScriptsDir/merge_sorted_RUM_files.pl $string 2>> $mainResultDir/PostProcessing-errorlog");
      &runCmd("/bin/mv $mainResultDir/RUM_NU.sorted_tmp $mainResultDir/RUM_NU.sorted");
    }

    # add counts of reads per chromosome to the mapping_stats.txt file:
    $date = `date`; chomp $date;
    print STDERR "[$date] Colating mapping stats\n";
    if(!-e "$mainResultDir/mapping_stats.txt"){
      unlink("$mainResultDir/mapping_stats.txt.tmp");
      $string = "$mainResultDir/mapping_stats.txt.tmp";
      for(my $j=1; $j<$numchunks+1; $j++) {
	$string = $string . " $mainResultDir/chr_counts_u.$j";
      }
      $node->runCmd("perl $perlScriptsDir/merge_chr_counts.pl $string 2>> $mainResultDir/PostProcessing-errorlog");
      $string = "$mainResultDir/mapping_stats.txt.tmp";
      for(my $j=1; $j<$numchunks+1; $j++) {
	$string = $string . " $mainResultDir/chr_counts_nu.$j";
      }
      $node->runCmd("perl $perlScriptsDir/merge_chr_counts.pl $string 2>> $mainResultDir/PostProcessing-errorlog");
      &runCmd("/bin/mv $mainResultDir/mapping_stats.txt.tmp $mainResultDir/mapping_stats.txt");
    }

    if($self->getProperty("createJunctionsFile") eq 'true') {
      $date = `date`; chomp $date;
      print STDERR "[$date] Creating junctions file\n";
      if(!-e "$mainResultDir/junctions_all.rum"){
        unlink("$mainResultDir/junctions_all.rum.tmp");
        $node->runCmd("perl $perlScriptsDir/make_RUM_junctions_file.pl $mainResultDir/RUM_Unique $mainResultDir/RUM_NU $genomeFastaFile $geneAnnotationFile $mainResultDir/junctions_all.rum.tmp $mainResultDir/junctions_all.bed $mainResultDir/junctions_high-quality.bed -faok 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/junctions_all.rum.tmp $mainResultDir/junctions_all.rum");
      }
    }

    $date = `date`; chomp $date;
    print STDERR "[$date] Creating coverage files\n";
    if(!-e "$mainResultDir/RUM_Unique.cov"){
      unlink("$mainResultDir/RUM_Unique.cov.tmp");
      $node->runCmd("perl $perlScriptsDir/rum2cov.pl $mainResultDir/RUM_Unique.sorted $mainResultDir/RUM_Unique.cov.tmp -name \"Unique Mappers\" 2>> $mainResultDir/PostProcessing-errorlog");
      &runCmd("/bin/mv $mainResultDir/RUM_Unique.cov.tmp $mainResultDir/RUM_Unique.cov");
    }

    if(!-e "$mainResultDir/RUM_Unique.normalized.cov"){
      unlink("$mainResultDir/RUM_Unique.normalized.cov.tmp");
      $node->runCmd("perl $perlScriptsDir/normalizeCov.pl $mainResultDir/RUM_Unique.cov -open > $mainResultDir/RUM_Unique.normalized.cov.tmp 2>> $mainResultDir/PostProcessing-errorlog");
      &runCmd("/bin/mv $mainResultDir/RUM_Unique.normalized.cov.tmp $mainResultDir/RUM_Unique.normalized.cov");
    }

    if(!-e "$mainResultDir/RUM_NU.cov"){
      unlink("$mainResultDir/RUM_NU.cov.tmp");
      $node->runCmd("perl $perlScriptsDir/rum2cov.pl $mainResultDir/RUM_NU.sorted $mainResultDir/RUM_NU.cov.tmp -name \"Non-Unique Mappers\" 2>> $mainResultDir/PostProcessing-errorlog");
      &runCmd("/bin/mv $mainResultDir/RUM_NU.cov.tmp $mainResultDir/RUM_NU.cov");
    }
    if(!-e "$mainResultDir/RUM_NU.normalized.cov"){
      unlink("$mainResultDir/RUM_NU.normalized.cov.tmp");
      $node->runCmd("perl $perlScriptsDir/normalizeCov.pl $mainResultDir/RUM_NU.cov -open > $mainResultDir/RUM_NU.normalized.cov.tmp 2>> $mainResultDir/PostProcessing-errorlog");
      &runCmd("/bin/mv $mainResultDir/RUM_NU.normalized.cov.tmp $mainResultDir/RUM_NU.normalized.cov");
    }

    if($self->getProperty("strandSpecific") eq 'true') {
      $date = `date`; chomp $date;
      print STDERR "[$date] Making strand specific files\n";
	# breakup RUM_Unique and RUM_NU files into plus and minus
      if(!-e "$mainResultDir/RUM_Unique.sorted.plus"){
        unlink("$mainResultDir/RUM_Unique.sorted.plus.tmp");
	$node->runCmd("perl $perlScriptsDir/breakup_RUM_files_by_strand.pl $mainResultDir/RUM_Unique.sorted $mainResultDir/RUM_Unique.sorted.plus.tmp $mainResultDir/RUM_Unique.sorted.minus 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_Unique.sorted.plus.tmp $mainResultDir/RUM_Unique.sorted.plus");
      }
      if(!-e "$mainResultDir/RUM_NU.sorted.plus"){
        unlink("$mainResultDir/RUM_NU.sorted.plus.tmp");
	$node->runCmd("perl $perlScriptsDir/breakup_RUM_files_by_strand.pl $mainResultDir/RUM_NU.sorted $mainResultDir/RUM_NU.sorted.plus.tmp $mainResultDir/RUM_NU.sorted.minus 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_NU.sorted.plus.tmp $mainResultDir/RUM_NU.sorted.plus");
      }
	# run rum2cov on all four files
      if(!-e "$mainResultDir/RUM_Unique.plus.cov"){
        unlink("$mainResultDir/RUM_Unique.plus.cov.tmp");
	$node->runCmd("perl $perlScriptsDir/rum2cov.pl $mainResultDir/RUM_Unique.sorted.plus $mainResultDir/RUM_Unique.plus.cov.tmp -name \"Unique Mappers Plus Strand\" 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_Unique.plus.cov.tmp $mainResultDir/RUM_Unique.plus.cov");
      }
      if(!-e "$mainResultDir/RUM_Unique.plus.normalized.cov"){
        unlink("$mainResultDir/RUM_Unique.plus.normalized.cov.tmp");
	$node->runCmd("perl $perlScriptsDir/normalizeCov.pl $mainResultDir/RUM_Unique.plus.cov -open > $mainResultDir/RUM_Unique.plus.normalized.cov.tmp 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_Unique.plus.normalized.cov.tmp $mainResultDir/RUM_Unique.plus.normalized.cov");
      }
      if(!-e "$mainResultDir/RUM_Unique.minus.cov"){
        unlink("$mainResultDir/RUM_Unique.minus.cov.tmp");
 	$node->runCmd("perl $perlScriptsDir/rum2cov.pl $mainResultDir/RUM_Unique.sorted.minus $mainResultDir/RUM_Unique.minus.cov.tmp -name \"Unique Mappers Minus Strand\" 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_Unique.minus.cov.tmp $mainResultDir/RUM_Unique.minus.cov");
      }
      if(!-e "$mainResultDir/RUM_Unique.minus.normalized.cov"){
        unlink("$mainResultDir/RUM_Unique.minus.normalized.cov.tmp");
	$node->runCmd("perl $perlScriptsDir/normalizeCov.pl $mainResultDir/RUM_Unique.minus.cov -open > $mainResultDir/RUM_Unique.minus.normalized.cov.tmp 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_Unique.minus.normalized.cov.tmp $mainResultDir/RUM_Unique.minus.normalized.cov");
      }
      if(!-e "$mainResultDir/RUM_NU.plus.cov"){
        unlink("$mainResultDir/RUM_NU.plus.cov.tmp");
	$node->runCmd("perl $perlScriptsDir/rum2cov.pl $mainResultDir/RUM_NU.sorted.plus $mainResultDir/RUM_NU.plus.cov.tmp -name \"Non-Unique Mappers Plus Strand\" 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_NU.plus.cov.tmp $mainResultDir/RUM_NU.plus.cov");
      }
      if(!-e "$mainResultDir/RUM_NU.plus.normalized.cov"){
        unlink("$mainResultDir/RUM_NU.plus.normalized.cov.tmp");
	$node->runCmd("perl $perlScriptsDir/normalizeCov.pl $mainResultDir/RUM_NU.plus.cov -open > $mainResultDir/RUM_NU.plus.normalized.cov.tmp 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_NU.plus.normalized.cov.tmp $mainResultDir/RUM_NU.plus.normalized.cov");
      }
      if(!-e "$mainResultDir/RUM_NU.minus.cov"){
        unlink("$mainResultDir/RUM_NU.minus.cov.tmp");
	$node->runCmd("perl $perlScriptsDir/rum2cov.pl $mainResultDir/RUM_NU.sorted.minus $mainResultDir/RUM_NU.minus.cov.tmp -name \"Non-Unique Mappers Minus Strand\" 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_NU.minus.cov.tmp $mainResultDir/RUM_NU.minus.cov");
      }
      if(!-e "$mainResultDir/RUM_NU.minus.normalized.cov"){
        unlink("$mainResultDir/RUM_NU.minus.normalized.cov.tmp");
	$node->runCmd("perl $perlScriptsDir/normalizeCov.pl $mainResultDir/RUM_NU.minus.cov -open > $mainResultDir/RUM_NU.minus.normalized.cov.tmp 2>> $mainResultDir/PostProcessing-errorlog");
        &runCmd("/bin/mv $mainResultDir/RUM_NU.minus.normalized.cov.tmp $mainResultDir/RUM_NU.minus.normalized.cov");
      }
    }

# cleanup temp files
    $date = `date`; chomp $date;
    print STDERR "[$date] Cleaning up temp files\n";

    #`rm RUM_Unique`;
    #`rm RUM_NU`;
    #`rm RUM_Unique.sorted`;
    #`rm RUM_NU.sorted`;
    if($self->getProperty("strandSpecific") eq 'true') {
	`rm RUM_Unique.sorted.plus`;
	`rm RUM_Unique.sorted.minus`;
	`rm RUM_NU.sorted.plus`;
	`rm RUM_NU.sorted.minus`;
    }
    for(my $j=1; $j<$numchunks+1; $j++) {
	`rm RUM_Unique.sorted.$j`;
	`rm RUM_NU.sorted.$j`;
    }    
    `rm RUM.sam*`;
    `rm chr_counts*`;
    `rm sam_header*`;
    if($transcriptBowtieIndex ne 'none') {
	`rm quant*`;
    }
    if($transcriptBowtieIndex ne 'none' && $self->getProperty("strandSpecific") eq 'true') {
	`rm feature_quantifications.ps`;
	`rm feature_quantifications.ms`;
	`rm feature_quantifications.pa`;
	`rm feature_quantifications.ma`;
	`rm feature_quantifications_normalized.ps`;
	`rm feature_quantifications_normalized.ms`;
	`rm feature_quantifications_normalized.pa`;
	`rm feature_quantifications_normalized.ma`;
    }

# do SNP calling, if requested
# not currently supported!!

#    my $SNPs = $self->getProperty("SNPs");
#    if($SNPs eq 'true') {
#      $date = `date`; chomp $date;
#      print STDERR "[$date] Calling SNPs\n";
#	$node->runCmd("identifySNPsFromBamFile.pl --genomeFastaFile $genomeFastaFile --varScanJarFile $ENV{GUS_HOME}/lib/java/VarScan.v2.2.5.jar --samtoolsPath /gpfs/fs0/share/apps/bs/bioscope/bin/samtools --bamFile $mainResultDir/RUM.bam");
#    }
    
    $date = `date`; chomp $date;
    print STDERR "[$date] Post processing steps complete\n";
    chdir("$currDir") || die "$!";
    return 1;
}

sub sortResultFiles {
  my ($self,$fn) = @_;  ##takes in string to do ls with
  my @files = `ls $fn`;
  my @tmp;
  foreach my $f (@files){
    next if !($f  =~ /\d$/);
    chomp $f;
    next if !($f  =~ /\d$/);
    if($f =~ /[a-zA-Z]\.(\d+)/){
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

sub isroman($) {
    my $arg = shift;
    $arg ne '' and
      $arg =~ /^(?: M{0,3})
                (?: D?C{0,3} | C[DM])
                (?: L?X{0,3} | X[LC])
                (?: V?I{0,3} | I[VX])$/ix;
}

sub arabic($) {
    my $arg = shift;
    my %roman2arabic = qw(I 1 V 5 X 10 L 50 C 100 D 500 M 1000);
    my %roman_digit = qw(1 IV 10 XL 100 CD 1000 MMMMMM);
    my @figure = reverse sort keys %roman_digit;
    $roman_digit{$_} = [split(//, $roman_digit{$_}, 2)] foreach @figure;
    isroman $arg or return undef;
    my $last_digit = 1000;
    my $arabic=0;
    foreach (split(//, uc $arg)) {
        my ($digit) = $roman2arabic{$_};
        $arabic -= 2 * $last_digit if $last_digit < $digit;
        $arabic += ($last_digit = $digit);
    }
    $arabic;
}

sub Roman($) {
    my $arg = shift;
    my %roman2arabic = qw(I 1 V 5 X 10 L 50 C 100 D 500 M 1000);
    my %roman_digit = qw(1 IV 10 XL 100 CD 1000 MMMMMM);
    my @figure = reverse sort keys %roman_digit;
    $roman_digit{$_} = [split(//, $roman_digit{$_}, 2)] foreach @figure;
    0 < $arg and $arg < 4000 or return undef;
    my $roman="";
    my $x;
    foreach (@figure) {
        my ($digit, $i, $v) = (int($arg / $_), @{$roman_digit{$_}});
        if (1 <= $digit and $digit <= 3) {
            $roman .= $i x $digit;
        } elsif ($digit == 4) {
            $roman .= "$i$v";
        } elsif ($digit == 5) {
            $roman .= $v;
        } elsif (6 <= $digit and $digit <= 8) {
            $roman .= $v . $i x ($digit - 5);
        } elsif ($digit == 9) {
            $roman .= "$i$x";
        }
        $arg -= $digit * $_;
        $x = $i;
    }
    $roman;
}

sub roman($) {
    lc Roman shift;
}

sub cmpChrs () {
    my $a2_c = lc($b);
    my $b2_c = lc($a);
    if($a2_c =~ /^\d+$/ && !($b2_c =~ /^\d+$/)) {
        return 1;
    }
    if($b2_c =~ /^\d+$/ && !($a2_c =~ /^\d+$/)) {
        return -1;
    }
    if($a2_c =~ /^[ivxym]+$/ && !($b2_c =~ /^[ivxym]+$/)) {
        return 1;
    }
    if($b2_c =~ /^[ivxym]+$/ && !($a2_c =~ /^[ivxym]+$/)) {
        return -1;
    }
    if($a2_c eq 'm' && ($b2_c eq 'y' || $b2_c eq 'x')) {
        return -1;
    }
    if($b2_c eq 'm' && ($a2_c eq 'y' || $a2_c eq 'x')) {
        return 1;
    }
    if($a2_c =~ /^[ivx]+$/ && $b2_c =~ /^[ivx]+$/) {
        $a2_c = "chr" . $a2_c;
        $b2_c = "chr" . $b2_c;
    }
   if($a2_c =~ /$b2_c/) {
	return -1;
    }
    if($b2_c =~ /$a2_c/) {
	return 1;
    }
    # dealing with roman numerals starts here

    if($a2_c =~ /chr([ivx]+)/ && $b2_c =~ /chr([ivx]+)/) {
	$a2_c =~ /chr([ivx]+)/;
	my $a2_roman = $1;
	$b2_c =~ /chr([ivx]+)/;
	my $b2_roman = $1;
	my $a2_arabic = arabic($a2_roman);
    	my $b2_arabic = arabic($b2_roman);
	if($a2_arabic > $b2_arabic) {
	    return -1;
	} 
	if($a2_arabic < $b2_arabic) {
	    return 1;
	}
	if($a2_arabic == $b2_arabic) {
	    my $tempa = $a2_c;
	    my $tempb = $b2_c;
	    $tempa =~ s/chr([ivx]+)//;
	    $tempb =~ s/chr([ivx]+)//;
	    my %temphash;
	    $temphash{$tempa}=1;
	    $temphash{$tempb}=1;
	    foreach my $tempkey (sort {cmpChrs($a,$b)} keys %temphash) {
		if($tempkey eq $tempa) {
		    return 1;
		} else {
		    return -1;
		}
	    }
	}
    }

    if($b2_c =~ /chr([ivx]+)/ && !($a2_c =~ /chr([a-z]+)/) && !($a2_c =~ /chr(\d+)/)) {
	return -1;
    }
    if($a2_c =~ /chr([ivx]+)/ && !($b2_c =~ /chr([a-z]+)/) && !($b2_c =~ /chr(\d+)/)) {
	return 1;
    }
    # roman numerals ends here
    if($a2_c =~ /chr(\d+)$/ && $b2_c =~ /chr.*_/) {
        return 1;
    }
    if($b2_c =~ /chr(\d+)$/ && $a2_c =~ /chr.*_/) {
        return -1;
    }
    if($a2_c =~ /chr([a-z])$/ && $b2_c =~ /chr.*_/) {
        return 1;
    }
    if($b2_c =~ /chr([a-z])$/ && $a2_c =~ /chr.*_/) {
        return -1;
    }
    if($a2_c =~ /chr(\d+)/) {
        my $numa = $1;
        if($b2_c =~ /chr(\d+)/) {
            my $numb = $1;
            if($numa < $numb) {return 1;}
	    if($numa > $numb) {return -1;}
	    if($numa == $numb) {
		my $tempa = $a2_c;
		my $tempb = $b2_c;
		$tempa =~ s/chr\d+//;
		$tempb =~ s/chr\d+//;
		my %temphash;
		$temphash{$tempa}=1;
		$temphash{$tempb}=1;
		foreach my $tempkey (sort {cmpChrs($a,$b)} keys %temphash) {
		    if($tempkey eq $tempa) {
			return 1;
		    } else {
			return -1;
		    }
		}
	    }
        } else {
            return 1;
        }
    }
    if($a2_c =~ /chrx(.*)/ && ($b2_c =~ /chr(y|m)$1/)) {
	return 1;
    }
    if($b2_c =~ /chrx(.*)/ && ($a2_c =~ /chr(y|m)$1/)) {
	return -1;
    }
    if($a2_c =~ /chry(.*)/ && ($b2_c =~ /chrm$1/)) {
	return 1;
    }
    if($b2_c =~ /chry(.*)/ && ($a2_c =~ /chrm$1/)) {
	return -1;
    }
    if($a2_c =~ /chr\d/ && !($b2_c =~ /chr[^\d]/)) {
	return 1;
    }
    if($b2_c =~ /chr\d/ && !($a2_c =~ /chr[^\d]/)) {
	return -1;
    }
    if($a2_c =~ /chr[^xy\d]/ && (($b2_c =~ /chrx/) || ($b2_c =~ /chry/))) {
        return -1;
    }
    if($b2_c =~ /chr[^xy\d]/ && (($a2_c =~ /chrx/) || ($a2_c =~ /chry/))) {
        return 1;
    }
    if($a2_c =~ /chr(\d+)/ && !($b2_c =~ /chr(\d+)/)) {
        return 1;
    }
    if($b2_c =~ /chr(\d+)/ && !($a2_c =~ /chr(\d+)/)) {
        return -1;
    }
    if($a2_c =~ /chr([a-z])/ && !($b2_c =~ /chr(\d+)/) && !($b2_c =~ /chr[a-z]+/)) {
        return 1;
    }
    if($b2_c =~ /chr([a-z])/ && !($a2_c =~ /chr(\d+)/) && !($a2_c =~ /chr[a-z]+/)) {
        return -1;
    }
    if($a2_c =~ /chr([a-z]+)/) {
        my $letter_a = $1;
        if($b2_c =~ /chr([a-z]+)/) {
            my $letter_b = $1;
            if($letter_a lt $letter_b) {return 1;}
	    if($letter_a gt $letter_b) {return -1;}
        } else {
            return -1;
        }
    }
    my $flag_c = 0;
    while($flag_c == 0) {
        $flag_c = 1;
        if($a2_c =~ /^([^\d]*)(\d+)/) {
            my $stem1_c = $1;
            my $num1_c = $2;
            if($b2_c =~ /^([^\d]*)(\d+)/) {
                my $stem2_c = $1;
                my $num2_c = $2;
                if($stem1_c eq $stem2_c && $num1_c < $num2_c) {
                    return 1;
                }
                if($stem1_c eq $stem2_c && $num1_c > $num2_c) {
                    return -1;
                }
                if($stem1_c eq $stem2_c && $num1_c == $num2_c) {
                    $a2_c =~ s/^$stem1_c$num1_c//;
                    $b2_c =~ s/^$stem2_c$num2_c//;
                    $flag_c = 0;
                }
            }
        }
    }
    if($a2_c le $b2_c) {
	return 1;
    }
    if($b2_c le $a2_c) {
	return -1;
    }


    return 1;
}

1;
