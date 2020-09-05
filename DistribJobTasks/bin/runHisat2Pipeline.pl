#!@perl@
use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use Getopt::Long;
use CBIL::Util::Utils;

my($mateA,$mateB,$sampleName,$delIntFiles,$hisat2Index,$hisat2,$extraHisatParams, $ppn, $maxIntronLen, $quantify, $maskedFile, $topLevelGeneFootprintFile, $junctions, $isStranded, $writeCov);
my $workingDir = ".";
#TODO check and update these
&GetOptions( 
            "mateA|ma=s"=> \$mateA,
            "mateB|mb=s" => \$mateB,
            "hisat2=s" => \$hisat2,
            "hisat2Index|x=s" => \$hisat2Index,
            "sampleName|s=s" => \$sampleName,
            "workingDir|w=s" => \$workingDir,
            "extraHisatParams|e:s" => \$extraHisatParams,
            "deleteIntermediateFiles!" => \$delIntFiles,
            "ppn=s" => \$ppn,
            "maxIntronLen=s" => \$maxIntronLen,
            "quantify!" => \$quantify,
            "maskedFile=s" => \$maskedFile,
            "topLevelGeneFootprintFile=s" => \$topLevelGeneFootprintFile,
            "junctions!" => \$junctions,
            "isStranded!" => \$isStranded,
            "writeCov!" => \$writeCov
            );

die "MateA file $mateA not found\n".&getParams() unless -e "$mateA";
die "MateB file $mateB not found\n".&getParams() if ($mateB && !-e "$mateB");
die "Hisat2 index must be specified\n".&getParams() unless (-e "$hisat2Index.1.ht2"); 
die "You must provide a sample name\n".&getParams() unless $sampleName;
##should add in usage
$hisat2 = $hisat2 eq 'default' ? 'hisat2' : $hisat2;  ##if not specified then hisat2 must be in path.

#change hisat2 params back to correct form
$extraHisatParams =~ s/_/-/g;
$extraHisatParams =~ s/#/ /g;


open(L,">>$workingDir/runHisat2Mapping.log");

select L;
$| = 1;

print L "runHisat2Mapping.pl run starting ".&getDate()."\n\n";

print L "Determining input file type ".&getDate()."\n\n";
my $mateAType = &getFileType($mateA);
my $mateBType;

if ($mateB) {
    $mateBType = &getFileType($mateB);
    if ($mateAType ne $mateBType) {
        die "Mate A $mateA and mate B $mateB files are not of the same type.\n"
    }
}
print L "Input file type is $mateAType ".&getDate()."\n\n";

if ($mateB) {
    print L "Running in paired end mode".&getDate()."\n\n";
} else {
    print L "Mate B not specified. Running in single end mode.".&getDate()."\n\n";
}

print L "Running FASTQC on raw reads".&getDate()."\n\n";

my ($mateAEncoding, $mateBEncoding);
$mateAEncoding = &runFastQC($mateA, $workingDir);

if ($mateB) {
    $mateBEncoding = &runFastQC($mateB, $workingDir);
    if ($mateAEncoding ne $mateBEncoding) {
        die "ERROR: The two read files use different phred encoding\n";
    }
}

print L "Files are using $mateAEncoding for quality encoding".&getDate()."\n\n";

print L "Running Trimmomatic".&getDate()."\n\n";
&runCmd("mkdir -p  $workingDir/trimmedReads");

if ($mateB) {
    print L "Running Trimmomatic in paired-end mode. Trimmomatic will remove adaptors from TruSeq2 and TruSeq3 compatible kits.".&getDate()."\n\n";
    &runCmd("java -jar \$eupath_dir/workflow-software/software/Trimmomatic/0.36/trimmomatic.jar PE -trimlog ${workingDir}/trimLog $mateA $mateB -$mateAEncoding -baseout ${workingDir}/trimmedReads/${sampleName} ILLUMINACLIP:\$GUS_HOME/data/DJob/DistribJobTasks/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20");
} else {
    print L "Running Trimmomatic in single-end mode. Trimmomatic will remove adaptors from TruSeq2 and TruSeq3 compatible kits.".&getDate()."\n\n";
    &runCmd("java -jar \$eupath_dir/workflow-software/software/Trimmomatic/0.36/trimmomatic.jar SE -trimlog ${workingDir}/trimLog -$mateAEncoding $mateA ${workingDir}/trimmedReads/${sampleName}_1P ILLUMINACLIP:\$GUS_HOME/data/DJob/DistribJobTasks/All_adaptors-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20");
}

my $trimmedA = "$workingDir/trimmedReads/$sampleName"."_1P";
my $trimmedB = "$workingDir/trimmedReads/$sampleName"."_2P";

die "Trimmed read file $trimmedA cannot be located.  Please check trimming was successful.\n\n" unless (-e $trimmedA);

if ($mateB) {
    die "Trimmed read file $trimmedB cannot be located.  Please check trimming was successful\n\n" unless (-e $trimmedB);
}

print L "Re-running FASTQC on trimmed reads".&getDate()."\n\n!";
&runFastQC($trimmedA, $workingDir);

if ($mateB) {
    &runFastQC($trimmedB, $workingDir);
}


print L "Aligning with Hisat2".&getDate()."\n\n";

my $type;
if ($mateAType eq "fastq") {
    print L "Aligning using fastq input".&getDate()."\n\n";
    $type = "-q";
} elsif ($mateAType eq "fasta") {
    print L "Aligning using fasta input".&getDate()."\n\n";
    $type = "-f";
}

my $cmd;
if ($mateB) {
    print L "Aligning in paired end mode\n\n".&getDate()."";
    $cmd = "$hisat2 $extraHisatParams -p $ppn -x $hisat2Index $type --max-introlen $maxIntronLen -1 $trimmedA -2 $trimmedB | samtools view -bS - | samtools sort -T $workingDir/$sampleName - > $workingDir/${sampleName}_sorted.bam";
} else {
    print L "Aligning in single end mode\n\n".&getDate()."";
    $cmd = "$hisat2 $extraHisatParams -p $ppn -x $hisat2Index $type --max-intronlen $maxIntronLen -U $trimmedA | samtools view -bS - | samtools sort -T $workingDir/$sampleName - > $workingDir/${sampleName}_sorted.bam";
}

print L "Alignment command: $cmd".&getDate()."\n\n";
&runCmd($cmd);

if ($quantify) {
    print L "Preparing for quantification".&getDate()."\n\n";
    &runCmd("samtools view -bh -F 4 -f 8 $workingDir/${sampleName}_sorted.bam > $workingDir/pair1.bam");
    &runCmd("samtools view -bh -F 8 -f 4 $workingDir/${sampleName}_sorted.bam > $workingDir/pair2.bam");
    &runCmd("samtools view -b -F 12 $workingDir/${sampleName}_sorted.bam > $workingDir/pairs.bam");
    &runCmd("samtools merge $workingDir/trimmed.bam $workingDir/pair*");

    &runCmd("samtools sort -n $workingDir/trimmed.bam > $workingDir/${sampleName}_sortedByName.bam");

    print L "Starting read quantification".&getDate()."\n\n";

    die "Masked gtf file $maskedFile cannot be opened for reading\n\n$!\n" unless defined($maskedFile);
    die "Gene footprint file $topLevelGeneFootprintFile cannot be opened for reading\n\n$!\n" unless defined($topLevelGeneFootprintFile);

    if ($isStranded) {
        print L "Quantifying in stranded mode".&getDate()."\n\n";
        print L "Counting unique reads".&getDate()."\n\n";

        &runCmd("htseq-count -a 0 -f bam -s reverse -t exon -i gene_id $workingDir/${sampleName}_sortedByName.bam $maskedFile > $workingDir/genes.htseq-union.firststrand.counts");
        &runCmd("htseq-count -a 0 -f bam -s yes -t exon -i gene_id $workingDir/${sampleName}_sortedByName.bam $maskedFile > $workingDir/genes.htseq-union.secondstrand.counts");

        print L "Counting nonunique reads".&getDate()."\n\n";
        &runCmd("htseq-count -a 0 -f bam --nonunique all -s reverse -t exon -i gene_id $workingDir/${sampleName}_sortedByName.bam $maskedFile > $workingDir/genes.htseq-union.firststrand.nonunique.counts");
        &runCmd("htseq-count -a 0 -f bam --nonunique all -s yes -t exon -i gene_id $workingDir/${sampleName}_sortedByName.bam $maskedFile > $workingDir/genes.htseq-union.secondstrand.nonunique.counts");

        print L "Calculating TPM".&getDate()."\n\n";
        &runCmd("makeTpmFromHtseqCountsDJob.pl --geneFootprintFile $topLevelGeneFootprintFile --senseUniqueCountFile $workingDir/genes.htseq-union.firststrand.counts --senseNUCountFile $workingDir/genes.htseq-union.firststrand.nonunique.counts --senseUniqueTpmFile $workingDir/genes.htseq-union.firststrand.tpm --senseNUTpmFile $workingDir/genes.htseq-union.firststrand.nonunique.tpm --antisenseUniqueCountFile $workingDir/genes.htseq-union.secondstrand.counts --antisenseNUCountFile $workingDir/genes.htseq-union.secondstrand.nonunique.counts --antisenseUniqueTpmFile $workingDir/genes.htseq-union.secondstrand.tpm --antisenseNUTpmFile $workingDir/genes.htseq-union.secondstrand.nonunique.tpm");

    } else { 
        print L "Quantifying in unstranded mode".&getDate()."\n\n";
        print L "Counting unique reads".&getDate()."\n\n";
        &runCmd("htseq-count -a 0 -f bam -s no -t exon -i gene_id $workingDir/${sampleName}_sortedByName.bam $maskedFile > $workingDir/genes.htseq-union.unstranded.counts");

        print L "Counting nonunique reads".&getDate()."\n\n";
        &runCmd("htseq-count -a 0 -f bam --nonunique all -s no -t exon -i gene_id $workingDir/${sampleName}_sortedByName.bam $maskedFile > $workingDir/genes.htseq-union.unstranded.nonunique.counts");

        print L "Calculating TPM".&getDate()."\n\n";
        &runCmd("makeTpmFromHtseqCountsDJob.pl --geneFootprintFile $topLevelGeneFootprintFile --senseUniqueCountFile $workingDir/genes.htseq-union.unstranded.counts --senseNUCountFile $workingDir/genes.htseq-union.unstranded.nonunique.counts --senseUniqueTpmFile $workingDir/genes.htseq-union.unstranded.tpm --senseNUTpmFile $workingDir/genes.htseq-union.unstranded.nonunique.tpm");
    }
}

if ($junctions) {
    print L "Quantifying junctions".&getDate()."\n\n";
    &runCmd("gsnapSam2Junctions.pl --is_bam --input_file $workingDir/${sampleName}_sorted.bam --output_file $workingDir/junctions.tab");
}

if ($writeCov) {
    print L "Writing coverage files".&getDate()."\n\n";
 
    my $isPairedEnd = (-e $mateB) ? 1 : 0;
    my $strandSpecific = $isStranded ? 1 : 0;  

    &runCmd("samtools index $workingDir/${sampleName}_sorted.bam");
    &runCmd("gsnapSplitBam.pl --mainResultDir $workingDir --strandSpecific $strandSpecific --isPairedEnd $isPairedEnd --bamFile $workingDir/${sampleName}_sorted.bam");  
}


sub getFileType {
    my $file = shift;
    my $type;
    open IN, "$file" or die "Cannot open $file\n";
    while (my $line = <IN>) {
        if ($line =~ /^>/) { 
            $type = "fasta";
            last;
        }
        elsif ($line =~ /^@/) {
            $type = "fastq";
            last;
        }
        else { 
            last;
        }
    }
    close IN;
    die "Could not determine if file $file is in fasta or fastq format\n" unless defined ($type);
    return $type;
}

sub runFastQC {
    my ($file, $workingDir) = @_;
    my ($fastQcFolder) = ( $file =~ m/([^\/]+)$/ );
    $fastQcFolder =~ s/\.fastq$//;
    $fastQcFolder =~ s/\.fq$//;
    $fastQcFolder .= "_fastqc";
    &runCmd("fastqc $file -o $workingDir --extract");
    my $encoding = &phred("$workingDir/$fastQcFolder/fastqc_data.txt");
    return $encoding;
}

sub phred {
    (my $fastqcFile) = @_;
    my %encoding = (
    "sanger" => "phred33",
    "illumina 1.3" => "phred64",
    "illumina 1.4" => "phred64",
    "illumina 1.5" => "phred64",
    "illumina 1.6" => "phred64",
    "illumina 1.7" => "phred64",
    "illumina 1.8" => "phred33",
    "illumina 1.9" => "phred33",
    "illumina 2"   => "phred33",
    "illumina 3"   => "phred33",
    "solexa" => "phred64"
    );

    my $phred;
    open (FH, $fastqcFile) or die "Cannot open $fastqcFile to determine phred encoding: $!\n";
    while (<FH>) {
        my $line = $_;
        if ($line =~ /encoding/i) {
            foreach my $format (keys %encoding) {
                if ($line =~ /$format/i) {
                    if (! defined $phred) {
                        $phred = $encoding{$format};
                    } elsif ($phred ne $encoding{$format}) {
                        $phred = "Error: more than one encoding type";
                    }
                }
            }
            if (! defined $phred) {
                $phred = "Error: format not recognized on encoding line";
            }
        }
    }
    if (! defined $phred) {
        $phred = "Error: encoding line not found in file";
    }
    close(FH);
    if ($phred =~ /error/i) {
        die "ERROR: Could not determine phred encoding: $phred\n\n";
    }
    return $phred;
} 
            
if($delIntFiles){
  print L "deleting extra files\n";
  &runCmd("rm -r $workingDir/trimmedReads");

  &runCmd("rm $workingDir/*pair*bam $workingDir/trimmed.bam $workingDir/${sampleName}_sortedByName.bam");  
}

close L;

sub getParams {
  return &getDate().": runHisat2Pipeline.pl ... parameter values:\n\tmateA=$mateA\n\tmateB=$mateB\n\tsampleName=$sampleName\n\thisat2Index=$hisat2Index\n\textraHisatParams=$extraHisatParams\n\tmaxIntronLen=$maxIntronLen\n\tmaskedfile=$maskedFile\n\ttopLevelGeneFootprintFile=$topLevelGeneFootprintFile\n\tworkingDir=$workingDir\n\n";

}

sub getDate {
  my $date = `date`;
  chomp $date;
  return $date;
}
