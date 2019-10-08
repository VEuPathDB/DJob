#!@perl@
use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use warnings;
use Getopt::Long;
use CBIL::Util::Utils;
use File::Basename;

my ($genomicSeqsFile, $bamFile, $gtfFile, $geneFootprintFile, $samtoolsIndex, $sampleName, $window, $snpsClusterDir);
my $workingDir = ".";

#get args
&GetOptions(
            "genomicSeqsFile|g=s"=> \$genomicSeqsFile,
            "bamFile|b=s"=> \$bamFile,
            "gtfFile=s"=>\$gtfFile,
            "geneFootprintFile=s" => \$geneFootprintFile,
            "samtoolsIndex|s=s"=>\$samtoolsIndex,
            "workingDir|w=s"=>\$workingDir,
            "sampleName=s"=>\$sampleName,
            "window=s"=>\$window,
            "snpsClusterDir=s"=>\$snpsClusterDir,
            );
die "Genomic Sequence file not found\n".&getParams() unless -e $genomicSeqsFile;
die "Bam file not found\n".&getParams() unless -e $bamFile;
die "Gtf file not found\n" .&getParams() unless -e $gtfFile;
die "Gene footprint file not found\n" .&getParams() unless -e $geneFootprintFile;
die "Samtools index not found\n" .&getParams() unless -e $samtoolsIndex;
die "you must provide a sample name\n".&getParams() unless $sampleName;

#open and select logging fh
open(L,">>$workingDir/runCNVTasks.log");

select L;
$| = 1;

print L "runCNVTasks.pl run starting ".&getDate()."\n\n";

my $out = "$workingDir/$sampleName.bed";
my $bamFileBase = basename($bamFile);


###sort bam file by name###
print L &getDate(). ": sorting bam file for htseq-count...\n";
my $cmd = "samtools sort -n $bamFile > $workingDir/namesort_$bamFileBase";
print L &getDate(). ": $cmd\n";
if (-e "$workingDir/namesort_$bamFileBase") {
    print L &getDate(). ": Sorting succeeded in previous run\n\n";
} else {
    &runCmd ($cmd);
}

###run htseq count###
print L &getDate(). ": running htseq-count...\n";
$cmd = "htseq-count -f bam -s no -t CDS -i gene_id -a 0 $workingDir/namesort_$bamFileBase $gtfFile > $workingDir/counts_$sampleName.txt";
print L &getDate(). ": $cmd\n";
&runCmd($cmd);
unless (-e "$workingDir/counts_$sampleName.txt") {
    die "Counts file $workingDir/counts_$sampleName.txt was not successfully created\n";
}

###calculate FPKMS###
print L &getDate(). ": calculating fpkms...\n";
$cmd = "makeFpkmFromHtseqCounts.pl --geneFootprintFile $geneFootprintFile --countFile $workingDir/counts_$sampleName.txt --fpkmFile $workingDir/fpkm_$sampleName.txt";
print L &getDate(). ": $cmd\n";
&runCmd($cmd);
unless (-e "$workingDir/fpkm_$sampleName.txt") {
    die "FPKM file $workingDir/fpkm_$sampleName.txt was not successfully created\n";
}

###make binned coverage file###
print L &getDate(). ": creating binned coverage file\n";
my $bedFile = &createBed($samtoolsIndex, $window, $workingDir, $sampleName);
unless (-e $bedFile) {
    die "Bed file $bedFile was not successfully created\n";
}
my $coverage = &getCoverage($bedFile, $bamFile, $out, $workingDir);

###tidy up###
close L;
&cleanUp($snpsClusterDir, $bedFile, $workingDir);

sub createBed {
    my ($index, $winLen, $workingDir, $sampleName) = @_;
    my $bedfile = "$workingDir/$sampleName"."_$winLen.bed";
    open (OUT, ">$bedfile") or die "Cannot write to temporary file $bedfile\n$!\n";
    open (IN, "$index") or die "Cannot open samtools index $index for reading\n$!\n";
    while (<IN>) {
        my ($chr, $length, $cumulative, $lineLength, $lineBLength) = split(/\t/, $_);
        die "Chromosome and length are not defined for line $. in $index. [$chr, $length]\n" unless(defined($chr) && defined($length));
        for (my $i=1; $i+$winLen<=$length; $i+=$winLen){
            if ($i+$winLen == $length) {
                printf OUT "%s\t%d\t%d\n", $chr, $i, $length;
            } else {
                printf OUT "%s\t%d\t%d\n", $chr, $i, $i+$winLen-1;
            }
        }
        printf OUT "%s\t%d\t%d\n", $chr, $length-($length % $winLen)+1, $length unless ($length % $winLen <= 1);
    }
    close(IN);
    close(OUT);
    return $bedfile
}

sub getCoverage {
    my ($bed, $bam, $out, $workingDir) = @_;
    my $genomeFile = &getGenomeFile($bam, $workingDir);
    my @coverageBed = split(/\n/, &runCmd("bedtools coverage -counts -sorted -g $genomeFile -a $bed -b $bam"));
    my $totalMapped = &runCmd("samtools view -c -F 4 $bam");
    open (OUT, ">$out") or die "Cannot write output file\n$!\n";
    foreach (@coverageBed) {
        my ($chr, $start, $end, $mapped, $numNonZero, $lenB, $propNonZero) = split(/\t/, $_);
        die "Chromosome, start and end coordinates or number of mapped reads are not defined\n" unless (defined($chr) && defined($start) && defined($end) && defined($mapped));
        $mapped = $mapped/$totalMapped;
        printf OUT "%s\t%d\t%d\t%g\n", $chr, $start, $end, $mapped;
    }
    close OUT;
}

sub getGenomeFile {
    my ($bam, $workingDir) = @_;
    open (G, ">$workingDir/genome.txt") or die "Cannot open genome file $workingDir/genome.txt for writing\n";
    my @header = split(/\n/, &runCmd("samtools view -H $bam"));
    foreach my $line (@header) {
        if ($line =~ m/\@SQ\tSN:/) {
            $line =~ s/\@SQ\tSN://;
            $line =~ s/\tLN:/\t/;
            print G "$line\n";
        }
    }
    close G;
    return "$workingDir/genome.txt";
}
sub getParams {
  return &getDate().": runCNVTasks.pl ... parameter values:\n\tgenomicSeqsFile=$genomicSeqsFile\n\tbamFile=$bamFile\n\tgtfFile=$gtfFile\n\tsampleName=$sampleName\n\tworkingDir=$workingDir\n\n";
}

sub cleanUp {
    my ($snpsClusterDir, $bedFile, $workingDir) = @_;
    &runCmd("/bin/rm $bedFile");
    &runCmd("/bin/rm $workingDir/genome.txt");
    &runCmd("/bin/rm -r $snpsClusterDir");
    &runCmd("/bin/rm $workingDir/namesort_$bamFileBase");
}

sub getDate {
  my $date = `date`;
  chomp $date;
  return $date;
}
