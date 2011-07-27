#!@perl@

use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use Getopt::Long;
use CBIL::Util::Utils;

open(LOG,">identifySNPsFromBamFile.log");

#print STDERR "args:\n"; foreach my $a (@ARGV){ print STDERR " '$a'\n"; }

##generates the following files
## $fileBase.snps  ## called SNPs
## $fileBase.cons  ## consensus sequence
## $fileBase.bam   ## sorted bam file
## $filebase.bam.bai  ## index for sorted bam file

## files that can be deleted
## $samFile  ##can be regenerated from sorted bam file and is very large
    ## should probably be deleted only if this entire script runs properly so outside this script
## $fileBase.pileup

my $debug = 0;

$| = 1;

my ($genomeFastaFile, $perlScriptsDir, $varscan, $samtools, $bamFile);

&GetOptions("genomeFastaFile=s" => \$genomeFastaFile, 
            "perlScriptsDir=s" => \$perlScriptsDir, 
            "varScanJarFile|vsjf=s" => \$varscan,
            "samtoolsPath|sp=s" => \$samtools,
            "bamFile|sf=s" => \$bamFile,
            );

die "requires --genomeFastaFile --varScanJarFile --samtoolsPath --bamFile\n" unless -e $genomeFastaFile && -e $varscan && -e $samtools && -e $bamFile;

my $fileBase = $bamFile;
$fileBase =~ s/\.bam$//;

print LOG "starting: ".`date`."\n";

##samtools to generate sorted bam file
my $cmd = "$samtools view -u $bamFile | $samtools sort - $fileBase";
&runCmd($cmd);
## check for error and return proper error code if failed
print LOG "finished BAM file generation ".`date` . "$cmd\n\n";

##samtools to generate pileup file
$cmd = "$samtools pileup -f $genomeFastaFile $fileBase.bam > $fileBase.pileup";
&runCmd($cmd);
print LOG "finished pileup file generation ".`date` . "$cmd\n\n";

##run varscan to identify SNPs
$cmd = "java -jar $varscan pileup2snp $fileBase.pileup --p-value .01 --min-coverage 5 > $fileBase.snps";
&runCmd($cmd);
print LOG "finished varscan run to call SNPs ".`date` . "$cmd\n\n";

## samtools to index bam file
$cmd = "$samtools index $fileBase.bam";
&runCmd($cmd);
print LOG "finished indexing BAM file ".`date` . "$cmd\n\n";

## samtools to index genome fasta file? .. should probably be done separately as once / genome.

## VarScan to generate consensus file so can use to generate combined snps
$cmd = "java -jar $varscan pileup2cns $fileBase.pileup > $fileBase.cons";
&runCmd($cmd);
print LOG "finished varscan run to create consensus file ".`date` . "$cmd\n\n";

## now cleanup by deleting the pileup file
$cmd = "/bin/rm $fileBase.pileup";
&runCmd($cmd);
print LOG "deleted pileup file ".`date` . "$cmd\n\n";


## the SAM file can also be deleted but should probably be done outside this script as this is one of the inputs to this script.
