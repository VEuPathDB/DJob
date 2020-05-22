#!@perl@
use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use CBIL::TranscriptExpression::CalculationsForTPM qw(doTPMCalculation);

my ($verbose,$geneFootprintFile,$senseUniqueCountFile,$senseUniqueTpmFile,$antisenseUniqueCountFile,$antisenseUniqueTpmFile, $senseNUCountFile, $senseNUTpmFile, $antisenseNUCountFile, $antisenseNUTpmFile);
&GetOptions("verbose!"=>\$verbose,
            "geneFootprintFile=s"=> \$geneFootprintFile,
            "senseUniqueCountFile=s" => \$senseUniqueCountFile,
            "senseUniqueTpmFile=s" => \$senseUniqueTpmFile,
            "antisenseUniqueCountFile=s" => \$antisenseUniqueCountFile,
            "antisenseUniqueTpmFile=s" => \$antisenseUniqueTpmFile,
            "senseNUCountFile=s" => \$senseNUCountFile,
            "senseNUTpmFile=s" => \$senseNUTpmFile,
            "antisenseNUCountFile=s" => \$antisenseNUCountFile,
            "antisenseNUTpmFile=s" => \$antisenseNUTpmFile
    ); 

if(!$geneFootprintFile || !$senseUniqueCountFile || !$senseUniqueTpmFile || !$senseNUCountFile || !$senseNUTpmFile){
	die "usage: makeTpmFromhtseqCountsDJob.pl --geneFootprintFile <geneFootprintFile> --senseUniqueCountFile <input file for unique counts (sense if strand-specific)> --senseUniqueTpmFile <output file for unique counts (sense if strand-specific)> --senseNUCountFile <input file for non-unique counts (sense if strand specific)> --senseNUTpmFile <output file for nonunique counts (sense if strand specific)> [--antisenseUniqueCountFile <input file for unique counts (antisense if strand-specific)> --antisenseUniqueTpmFile <output file for unique counts(antisense if strand-specific)> --antisenseNUCountFile <input file for nonunique counts (antisense if strand specific)> --antisenseNUTpmFile <output file for nonunique counts (antisense if strand specific)>]\n";
}

my $geneLengths;
open(IN, "<$geneFootprintFile");
my $line = <IN>;
while ($line=<IN>) {
    chomp($line);
    my ($project, $gene, $length, @rest) = split(/\t/, $line);
    $geneLengths->{$gene} = $length;
}
close(IN);

&doTPMCalculation ($geneLengths, $senseUniqueCountFile, $senseNUCountFile, $antisenseUniqueCountFile, $antisenseNUCountFile, $senseUniqueTpmFile, $senseNUTpmFile, $antisenseUniqueTpmFile, $antisenseNUTpmFile);

