#!@perl@
use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use warnings;
use Getopt::Long;
use CBIL::Util::Utils;

my ($genomicSeqsFile, $bamFile, $gtfFile, $samtoolsIndex, $sampleName, $window, $snpsClusterDir);
my $workingDir = ".";

#get args
&GetOptions(
            "genomicSeqsFile|g=s"=> \$genomicSeqsFile,
            "bamFile|b=s"=> \$bamFile,
            "gtfFile=s"=>\$gtfFile,
            "samtoolsIndex|s=s"=>\$samtoolsIndex,
            "workingDir|w=s"=>\$workingDir,
            "sampleName=s"=>\$sampleName,
            "window=s"=>\$window,
            "snpsClusterDir=s"=>\$snpsClusterDir,
            );
die "Genomic Sequence file not found\n".&getParams() unless -e $genomicSeqsFile;
die "Bam file not found\n".&getParams() unless -e $bamFile;
die "Gtf file not found\n" .&getParams() unless -e $gtfFile;
die "Samtools index not found\n" .&getParams() unless -e $samtoolsIndex;
die "you must provide a sample name\n".&getParams() unless $sampleName;

#open and select logging fh
open(L,">>$workingDir/runCNVTasks.log");

select L;
$| = 1;

print L "runCNVTasks.pl run starting ".&getDate()."\n\n";

my $out = "$workingDir/$sampleName.bed";


###run cufflinks###
print L &getDate(). ": running Cufflinks...\n";
my $cmd = "cufflinks -u -N -p 4 -b $genomicSeqsFile -G $gtfFile -o $workingDir/Cufflinks $bamFile";
print L &getDate(). ": $cmd\n";
if (-e "$workingDir/Cufflinks/genes.fpkm_tracking") {
    print L &getDate(). "Cufflinks command succeeded in previous run\n\n";
} else {
    &runCmd ($cmd); print L "\n";
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
}

sub getDate {
  my $date = `date`;
  chomp $date;
  return $date;
}
