package DJob::DistribJobTasks::GSNAPTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;
use CBIL::TranscriptExpression::SplitBamUniqueNonUnique qw(splitBamUniqueNonUnique);

@ISA = (DJob::DistribJob::Task);
use strict;
# [name, default (or null if reqd), comment]
my @properties = 
    (
     ["mateA",   "",     "full path to reads file"],
     ["mateB",   "none",     "full path to paired reads file (optional)"],
     ["genomeDatabase",   "",     "full path to the genome database"],
     ["iitFile",   "none",     "full path to the iit file for splice sites"],
     ["gtfFile",   "none",     "full path to the gtf file (rRNAs removed)"],
     ["maskFile",   "none",     "full path to the gtf masked file (rRNAs removed); required for HTseq"],
     ["sraSampleIdQueryList", "none", "Comma delimited list of identifiers that can be used to retrieve SRS samples"],
     ["extraGsnapParams", "none", "GSNAP parameters other than default"],
     ["outputFileBasename", "results", "Base name for the results file"],
     ["nPaths",   "30",     "Limits the number of nonunique mappers printed to a max of [30]"],
     ["deleteIntermediateFiles", "true", "[true]|false: if true then deletes intermediate files to save space"],
     ["quantify", "true", "[true]|false: if true then runs HTSeq"],
     ["writeCovFiles", "true", "[true]|false: if true then runs bamutils"],
     ["isStrandSpecific", "false", "[true]|false"],
     ["quantifyJunctions", "true", "[true]|false: if true then runs gsnapSam2Junctions"],
     ["topLevelFastaFaiFile", "none", "required if writeCovFiles is turned on"],
     ["topLevelGeneFootprintFile", "none", "required if quantify is true"],
     ["hasKnownSpliceSites", "true", "if true gsnap will use the -s flag"]
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
    
    if($sidlist && $sidlist ne 'none'){ ##have a value and other than default
	my $mateA = $self->getProperty('mateA');
	my $mateB = $self->getProperty('mateB');
	
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
	    &runCmd("getFastqFromSra.pl --workingDir $inputDir --readsOne $mateA --readsTwo $mateB --sampleIdList '$sidlist'");
	}
	
	
    }
    else {
	my $mateA = $self->getProperty('mateA');
	my $mateB = $self->getProperty('mateB');
	my $basemateA = basename($mateA);
	my $basemateB = basename($mateB);
	if (! -e $mateB) {
	    $baseName = $basemateA;
	}
	else {
	    my @A= split "", $basemateA;
	    my @B= split "", $basemateB;
	    my $count = 0;
	    foreach my $element (@A) {
		if ($element eq $B[$count]) {
		    $baseName.=$element;
		    $count ++;
		}
		else {
		    last;
		}
	    }
	}
    }
    my $mateA = $self->getProperty('mateA');
    my $mateB = $self->getProperty('mateB');
    
######## SO I ONLY WANT TO RUN THIS FOR RNASEQ DATA
    
    
#### so if the files are fasta I want to skip this step. I think I will also add a warning that the files were fasta and so if not expecting microarray data then you should create fast
#q files and re run. 
    
    
#perl one liner to use perl -ne 'my $type; if ( /^>/ ) { $type="is fasta"} elsif ( /^@/ ) {$type="is fastq"} else {exit}; print "$type\n"
    my $type;
    if (-e $mateA) {
	open IN, "$mateA" or die "cant open $mateA";
	while (my $line =<IN>) {
	    if ( $line =~ /^>/ ) {
		$type = "fasta";
		last;
	    } 
	    elsif ( $line =~ /^@/ ) {
		$type = "fastq";
		last;
	    } 
	    else {
		print " ERROR: the first line doesnt start with > or @ its $line line ends and so cant determine if fasta or fastq\n";
		last;
		
	    }
	}
    }
    close IN;
   
    if (-e $mateB) {
	my $typeB;
	open  IN2, "$mateB" or die "cant open $mateB";
	while (my $line =<IN2>) {
	    if ( $line =~/^>/ ) {
		$typeB = "fasta";
	    }
            elsif ( $line =~ /^@/ ) {
		$typeB = "fastq";
	    }
            else {
                last;
            }
        }
	if ($type ne $typeB) {
	    die "mateA and mateB are not the same type of files one is fasta and one is fastq";
	}

    }	
	close IN2;

### However I want to check the fastq files and if they are created files (all quality score are I) then I want to force -phred33 for trimmomatic. 
    
    if ($type eq "fasta") {
	print "WARNING : fasta file detected ... if microarray data then please ignore warning otherwise please convert fasta to fastq for RNASeq data\n\n";
    }
    
    else {
#HERE I WANT TO DEAL WITH THE FASTQ QUAL ISSUE
	my $is_fake;
	if (-e $mateA) {
	    open IN3, "$mateA" or die "cant open $mateA";
	    my $count = 0;
	    while (my $line =<IN3>) {
		chomp $line;
		$count ++;
		if ($count == 4) {
		    $count = 0;
		    if ($line !~ /^I+$/) { #is anything but just Is
			$is_fake = 0;
			last;
		    }
		    else {
			$is_fake = 1;
		    }
		}
		else {
		    next;
		}
	    }
	}
	    close IN3;
	
	
	
	
	
	print "running FastQC on raw reads output files can be found in the main results folder \n";
	&runCmd("fastqc $mateA $mateB -o $inputDir");	    
	
	
	#want to do the trimming here and want to set the properties mateA and mateB here. this propery it already set for those with no sidlist
	if ($is_fake ==0) {
	    
	    if((-e "$mateA")&& (-e "$mateB") && ($mateB ne 'none')){
		print "running Paired End Trimmomatic to remove any adaptors if different chemistry than  TruSeq2 (as used in GAII machines) and TruSeq3 (as used by HiSeq and MiSeq machines) please supply custom adaptor fasta";
		&runCmd("java -jar \$eupath_dir/workflow-software/software/Trimmomatic/0.36/trimmomatic.jar PE -trimlog ${inputDir}/trimLog $mateA $mateB -baseout ${inputDir}/${baseName} ILLUMINACLIP:\$GUS_HOME/data/DJob/DistribJobTasks/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	    }
	    elsif((-e "$mateA")&&((! -e $mateB) || ($mateB eq 'none'))) {
		print "running Single End Trimmomatic to remove any adaptors if different chemistry than  TruSeq2 (as used in GAII machines) and TruSeq3 (as used by HiSeq and MiSeq machines) please supply custom adaptor fasta";
		&runCmd("java -jar \$eupath_dir/workflow-software/software/Trimmomatic/0.36/trimmomatic.jar SE -trimlog ${inputDir}/trimLog $mateA ${inputDir}/${baseName}_1P ILLUMINACLIP:\$GUS_HOME/data/DJob/DistribJobTasks/All_adaptors-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	    }
	    else {
		"ERROR: print reads files not found in $inputDir or not retrieved from SRA";
	    }
	} 
	elsif ($is_fake ==1) {
	    if((-e "$mateA")&& (-e "$mateB") && ($mateB ne 'none')){
		print "running Paired End Trimmomatic to remove any adaptors if different chemistry than  TruSeq2 (as used in GAII machines) and TruSeq3 (as used by HiSeq and MiSeq machines) please supply custom adaptor fasta";
		&runCmd("java -jar \$eupath_dir/workflow-software/software/Trimmomatic/0.36/trimmomatic.jar PE  -phred33 -trimlog ${inputDir}/trimLog $mateA $mateB -baseout ${inputDir}/${baseName} ILLUMINACLIP:\$GUS_HOME/data/DJob/DistribJobTasks/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	    }
	    elsif((-e "$mateA")&&((! -e $mateB) || ($mateB eq 'none'))) {
		print "running Single End Trimmomatic to remove any adaptors if different chemistry than  TruSeq2 (as used in GAII machines) and TruSeq3 (as used by HiSeq and MiSeq machines) please supply custom adaptor fasta";
		&runCmd("java -jar \$eupath_dir/workflow-software/software/Trimmomatic/0.36/trimmomatic.jar SE -phred33 -trimlog ${inputDir}/trimLog $mateA ${inputDir}/${baseName}_1P ILLUMINACLIP:\$GUS_HOME/data/DJob/DistribJobTasks/All_adaptors-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");
	    }
	    else {
		"ERROR: print reads files not found in $inputDir or not retrieved from SRA";
	    }
	}
	
	else { 
	    die "could not determine if FASTQ has experimental quality score or if we have created a fastq from a fasta\n";
	}
	
	
	
	
	
	my $trimmedA = $inputDir."/".$baseName."_1P";
	my $trimmedB = $inputDir."/".$baseName."_2P";
	if ((! -e $mateB || $mateB eq 'none')) {
	    $self->setProperty('mateB',"$mateB");
	}
	else {
	    $self->setProperty('mateB',"$trimmedB");
	}
	$self->setProperty('mateA',"$trimmedA");
	print "running FastQC on trimmed reads output files can be found in the main results folder\n";
	&runCmd("fastqc $trimmedA $trimmedB -o $inputDir");	    
    }
    
}

sub initNode {
    my ($self, $node, $inputDir) = @_;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $reads = $self->getProperty('mateA');
    my $paired = $self->getProperty('mateB');

    if (-e "$reads.gz"){
      print "unzipping $reads.gz\n";
      `gunzip $reads.gz`;
    }

    if (-e "$paired.gz"){
      print "unzipping $paired.gz\n";
      `gunzip $paired.gz`;
    }

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

    my $dashSParam;
    my $hasKnownSpliceSites = $self->getProperty("hasKnownSpliceSites");
    if($hasKnownSpliceSites && lc($hasKnownSpliceSites) eq 'true') {
      $dashSParam = "-s $iitFile";
    }
    my $cmd;
    if (-e $mateB) {
	 $cmd = "gsnap $extraGsnapParams --force-xs-dir -q $q  --quiet-if-excessive -N 1 $dashSParam -A sam -n $nPaths -D $databaseDirectory -d $databaseName  $mateA $mateB";
    }
    elsif (! -e $mateB) {
	 $cmd = "gsnap $extraGsnapParams --force-xs-dir -q $q  --quiet-if-excessive -N 1 $dashSParam -A sam -n $nPaths -D $databaseDirectory -d $databaseName  $mateA";
    }
    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $self->runCmdOnNode($node, "samtools view -Sb $nodeExecDir/subtask.output > $mainResultDir/${subTaskNum}_node.bam 2>>$nodeExecDir/subtask.stderr");

    return $node->getErr();
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;

  my $outputFileBasename = $self->getProperty("outputFileBasename");
    print "output file base is $outputFileBasename\n";
  my @bams = glob "$mainResultDir/*_node.bam";

  die "Did not find  bam files in $mainResultDir/*_node.bam" unless(scalar @bams > 0);
  opendir DIR, $inputDir;
  my @list = readdir DIR;
  closedir (DIR);
  foreach my $file (@list) {
      if ($file =~ m/.html$/) {
  print "moving fastqc files to results directory\n";
  &runCmd("mv $inputDir/$file  $mainResultDir/"); 
      }
  }

  if(scalar @bams > 1) {
    $self->runCmdOnNode($node, "samtools merge $mainResultDir/${outputFileBasename}.bam $mainResultDir/*_node.bam");
  }
  else {
    $self->runCmdOnNode($node, "cp $bams[0] $mainResultDir/${outputFileBasename}.bam");
  }

  # sort bams by location
  $self->runCmdOnNode($node, "samtools sort -o $mainResultDir/${outputFileBasename}_sorted.bam $mainResultDir/${outputFileBasename}.bam");

  $self->runCmdOnNode($node, "samtools view -bh -F 4 -f 8 $mainResultDir/${outputFileBasename}.bam > $mainResultDir/pair1.bam");
  $self->runCmdOnNode($node, "samtools view -bh -F 8 -f 4 $mainResultDir/${outputFileBasename}.bam > $mainResultDir/pair2.bam");
  $self->runCmdOnNode($node, "samtools view -b -F 12 $mainResultDir/${outputFileBasename}.bam > $mainResultDir/pairs.bam");

  $self->runCmdOnNode($node, "samtools merge $mainResultDir/trimmed.bam $mainResultDir/pair*");
  $self->runCmdOnNode($node, "samtools sort -n -o $mainResultDir/${outputFileBasename}_sortedByName.bam $mainResultDir/trimmed.bam");

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
    my $maskedFile = $self->getProperty("maskFile");
    my $topLevelGeneFootprintFile = $self->getProperty("topLevelGeneFootprintFile");

   # Cufflinks
    # if($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
    #     $self->runCmdOnNode($node, "cufflinks --no-effective-length-correction --compatible-hits-norm --library-type fr-firststrand -o $mainResultDir -G $maskedFile $mainResultDir/${outputFileBasename}_sorted.bam");
    #     rename "$mainResultDir/genes.fpkm_tracking", "$mainResultDir/genes.cuff.firststrand.fpkm_tracking";
    #     rename "$mainResultDir/isoforms.fpkm_tracking", "$mainResultDir/isoforms.cuff.firststrand.fpkm_tracking";
    #     $self->runCmdOnNode($node, "cufflinks --no-effective-length-correction --compatible-hits-norm --library-type fr-secondstrand -o $mainResultDir -G $maskedFile $mainResultDir/${outputFileBasename}_sorted.bam");
    #     rename "$mainResultDir/genes.fpkm_tracking", "$mainResultDir/genes.cuff.secondstrand.fpkm_tracking";
    #     rename "$mainResultDir/isoforms.fpkm_tracking", "$mainResultDir/isoforms.cuff.secondstrand.fpkm_tracking";
    # }
    # else {
    #     $self->runCmdOnNode($node, "cufflinks --no-effective-length-correction --compatible-hits-norm --library-type fr-unstranded -o $mainResultDir -G $maskedFile $mainResultDir/${outputFileBasename}_sorted.bam");
    #     rename "$mainResultDir/genes.fpkm_tracking", "$mainResultDir/genes.cuff.unstranded.fpkm_tracking";
    #     rename "$mainResultDir/isoforms.fpkm_tracking", "$mainResultDir/isoforms.cuff.unstranded.fpkm_tracking";
    # }
    
    # HTSeq
    my @modes = ('union');
#    my @modes = ('union', 'intersection-nonempty', 'intersection-strict');
    if ($isStrandSpecific && lc($isStrandSpecific) eq 'true') {

	for (my $i=0; $i<@modes; $i++) {
	    my $mode = $modes[$i];
	    $self->runCmdOnNode($node, "htseq-count --format=bam --order=name --stranded=reverse --type=exon --idattr=gene_id --mode=$mode $mainResultDir/${outputFileBasename}_sortedByName.bam $maskedFile > $mainResultDir/genes.htseq-$mode.firststrand.counts");
	    $self->runCmdOnNode($node, "htseq-count --format=bam --order=name --stranded=yes --type=exon --idattr=gene_id --mode=$mode $mainResultDir/${outputFileBasename}_sortedByName.bam $maskedFile > $mainResultDir/genes.htseq-$mode.secondstrand.counts");

	    $self->runCmdOnNode($node, "makeFpkmFromHtseqCounts.pl --geneFootprintFile $topLevelGeneFootprintFile --countFile $mainResultDir/genes.htseq-$mode.firststrand.counts --fpkmFile $mainResultDir/genes.htseq-$mode.firststrand.fpkm --antisenseCountFile $mainResultDir/genes.htseq-$mode.secondstrand.counts --antisenseFpkmFile $mainResultDir/genes.htseq-$mode.secondstrand.fpkm");
	}
    }
    else {
      for (my $i=0; $i<@modes; $i++) {
        my $mode = $modes[$i];
        $self->runCmdOnNode($node, "htseq-count --format=bam --order=name --stranded=no --type=exon --idattr=gene_id --mode=$mode $mainResultDir/${outputFileBasename}_sortedByName.bam $maskedFile > $mainResultDir/genes.htseq-$mode.unstranded.counts");
	$self->runCmdOnNode($node, "makeFpkmFromHtseqCounts.pl --geneFootprintFile $topLevelGeneFootprintFile --countFile $mainResultDir/genes.htseq-$mode.unstranded.counts --fpkmFile $mainResultDir/genes.htseq-$mode.unstranded.fpkm");
      }
    }
  }

  # Junctions
  if($quantifyJunctions && lc($quantifyJunctions eq 'true')) {
    $self->runCmdOnNode($node, "gsnapSam2Junctions.pl  --is_bam  --input_file $mainResultDir/${outputFileBasename}_sorted.bam --output_file $mainResultDir/junctions.tab");
  }

  # COVERAGE PLOTS
  if($writeCovFiles && lc($writeCovFiles) eq 'true') {


#    my $topLevelFastaFaiFile = $self->getProperty("topLevelFastaFaiFile");
#    unless(-e $topLevelFastaFaiFile) {
#      die "Top Level Genome fa.fai File $topLevelFastaFaiFile does not exist";
#    }

      $self->runCmdOnNode($node, "samtools index $mainResultDir/${outputFileBasename}_sorted.bam");

    my $mateB = $self->getProperty('mateB');

    my $isPairedEnd = 1;
      if((lc($mateB) eq 'none') || (! -e $mateB)){
	  $isPairedEnd = 0;
      }
    my $strandSpecific = 0;
    if ($isStrandSpecific && lc($isStrandSpecific) eq 'true') {
      $strandSpecific = 1;
    }
 

    $self->runCmdOnNode($node, "gsnapSplitBam.pl --mainResultDir $mainResultDir --strandSpecific $strandSpecific --isPairedEnd $isPairedEnd --bamFile $mainResultDir/${outputFileBasename}_sorted.bam");
      
#    my $splitExpDir = splitBamUniqueNonUnique($mainResultDir, $strandSpecific, $isPairedEnd, "$mainResultDir/${outputFileBasename}_sorted.bam");
  }

  return 1;
}

1;
