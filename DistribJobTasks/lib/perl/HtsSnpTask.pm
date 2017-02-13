package DJob::DistribJobTasks::HtsSnpTask;

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
	["fastaFile", "", "full path of fastaFile for genome"],
	["mateA", "none", "full path to file of reads"],
	["mateB", "none", "full path to file of paired ends reads"],
	["sraSampleIdQueryList", "none", "Comma delimited list of identifiers that can be used to retrieve SRS samples"],
	["outputPrefix", "result", "prefix of output files"],
	["bwaIndex", "none", "full path of the bwa indices .. likely same as fasta file"],
	["bowtieIndex", "none", "full path of the bowtie indices .. likely same as fasta file"],
	["isColorspace", "false", "input sequence reads are in SOLiD colorspace.  Quality files must be exactly matename.qual"],
	["varscan", "/genomics/eupath/workflow-software/software/VarScan/2.2.10/VarScan.jar", "full path to the varscan jar file"],
	["gatk", "/genomics/eupath/workflow-software/software/gatk/1.5.31/GenomeAnalysisTK.jar", "full path to the GATK jar file"],
	["bowtie2", "default", "full path to the bowtie2 bin dir"],
	["strain", "", "strain to be put into the GFF file"],
	["consPercentCutoff", "60", "minimum allele percent for calling consensus base"],
	["snpPercentCutoff", "20", "minimum allele percent for calling SNPs"],
	["editDistance", "0.04", "mismatch in bwa"],
        ["snpsOnly", "false", "if true then doesn't compute consensus or indels"],
	["deleteIntermediateFiles", "true", "[true]|false: if true then deletes intermediate files to save space"],
        ["noTrimming", "false", "true|[false]: if true then no adaptor or quality trimming is carried out"]
);

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;
    ##need to download fastq from sra if sample ids passed in.
    my $baseName;
    my $sidlist = $self->getProperty('sraSampleIdQueryList');
    my $isColorspace= $self->getProperty('isColorspace');
    my $noTrimming = $self->getProperty('noTrimming');
    if($sidlist && $sidlist ne 'none'){ ##have a value and other than default
	my $mateA = $self->getProperty('mateA');
	my $mateB = $self->getProperty('mateB');
	if(!$mateA || $mateA eq 'none'){
	    $mateA = $isColorspace eq 'true' ? "$inputDir/reads_1.csfasta" : "$inputDir/reads_1.fastq";
	    $self->setProperty('mateA',"$mateA");
	    $mateB = $isColorspace eq 'true' ? "$inputDir/reads_2.csfasta" : "$inputDir/reads_2.fastq";
	    $self->setProperty('mateB',"$mateB");
	    $baseName = "reads";
	}
	if(-e "$mateA"){
	    print "reads file $mateA already present so not retrieving from SRA\n";
	}
	else{  ##need to retrieve here
	    print "retrieving reads from SRA for '$sidlist'\n";
	    my $sraCmd = "getDataFromSra.pl --workingDir $inputDir --readsOne $mateA --readsTwo $mateB --sampleIdList '$sidlist'";
	    if($isColorspace eq 'true'){
		$sraCmd .= " --isColorspace";
	    }
	    &runCmd($sraCmd);
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

    
###### MAKE SURE THE REST DOESNT HAPPEN IF ITS COLOURSPACE!!!!!!!!!! 
    if ($isColorspace eq 'false') {  
	print "running FastQC on raw reads output files can be found in the main results folder \n";
	&runCmd("fastqc $mateA $mateB -o $inputDir");
	if ($noTrimming eq 'true') {
	}
	else {
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

    my $fastaFile = $self->getProperty ("fastaFile");
    my $mateA = $self->getProperty ("mateA");
    my $mateB = $self->getProperty ("mateB");
    my $outputPrefix = $self->getProperty ("outputPrefix");
    my $bwaIndex = $self->getProperty ("bwaIndex");
    my $bowtieIndex = $self->getProperty ("bowtieIndex");
    my $gatk = $self->getProperty ("gatk");
    my $varscan = $self->getProperty ("varscan");
    my $strain = $self->getProperty ("strain");
    my $consPercentCutoff = $self->getProperty ("consPercentCutoff");
    my $snpPercentCutoff = $self->getProperty ("snpPercentCutoff");
    my $editDistance = $self->getProperty ("editDistance");
    my $snpsOnly = $self->getProperty ("snpsOnly");
    my $wDir = "/eupath/data/EuPathDB/workflows/TestFlow/jane/SNP_tests/PNG17/new/master/mainresult";
    my $bowtie2 = $self->getProperty ("bowtie2");

    $consPercentCutoff = $snpPercentCutoff if $snpPercentCutoff > 60;

    if ($fastaFile !~ /\.fa/ || $fastaFile !~ /\.fasta/) {
        my $tempFile = $fastaFile;
        $tempFile =~ s/\.\w+$/\.fa/;
        `ln -s $fastaFile $tempFile` unless (-e $tempFile);
        $fastaFile = $tempFile;
    }
    
    my $cmd = "runHTS_SNPs.pl --fastaFile $fastaFile --mateA $mateA".(-e "$mateB" ? " --mateB $mateB" : "");
    $cmd .= " --outputPrefix $outputPrefix --varscan $varscan --bowtieIndex $bowtieIndex";
    $cmd .= " --gatk $gatk --bowtie2 $bowtie2";
    if($self->getProperty('isColorspace') eq 'true'){
      $cmd .= " --isColorspace";
    }
    $cmd .= " --strain $strain --consPercentCutoff $consPercentCutoff --snpPercentCutoff $snpPercentCutoff";
    $cmd .= " --editDistance $editDistance".($snpsOnly eq 'false' ? "" : " --snpsOnly");
    $cmd .= " --workingDir $wDir" . ($self->getProperty('deleteIntermediateFiles') eq 'true' ? " --deleteIntermediateFiles" : "");
      
#    print "Returning command: $cmd\n";
#    exit(0);  ##for testing
    return $cmd;
}

##cleanup materDir here and remove extra files that don't want to transfer back to compute node
sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
}

sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;
  my $sidlist = $self->getProperty('sraSampleIdQueryList');
  if($sidlist && $sidlist ne 'none' && $self->getProperty('deleteIntermediateFiles') eq 'true'){ ##have a value and other than default so reads were retrieved from sra
    my $mateA = $self->getProperty('mateA');
    my $mateB = $self->getProperty('mateB');
    unlink($mateA) if -e "$mateA";
    unlink($mateB) if -e "$mateB";
  }
}

1;
