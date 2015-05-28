# Take as input count files (two if strand-specific data, else one) from HTSeq and a gene footprint file and outputs fpkm files (two if strand-specific data, else one). The denominator in the fpkm is the product of the length of the gene footprint and the total counts (sense+antisense total if strand-specific).

use strict;
use warnings;
use Getopt::Long;

my ($verbose,$geneFootprintFile,$countFile,$fpkmFile,$antisenseCountFile,$antisenseFpkmFile);
&GetOptions("verbose!"=>\$verbose,
            "geneFootprintFile=s"=> \$geneFootprintFile,
            "countFile=s" => \$countFile,
            "fpkmFile=s" => \$fpkmFile,
            "antisenseCountFile=s" => \$antisenseCountFile,
            "antisenseFpkmFile=s" => \$antisenseFpkmFile,
    ); 

if(!$geneFootprintFile || !$countFile || !$fpkmFile){
	die "usage: makeFpkmFromhtseqCounts.pl --geneFootprintFile <geneFootprintFile> --countFile <input file (sense if strand-specific)> --fpkmFile <output file (sense if strand-specific)> [--antisenseCountFile <input file (antisense if strand-specific)> --antisenseFpkmFile <output file (antisense if strand-specific)>]\n";
}

my $isStrandSpecific = 0;
if ($antisenseCountFile) {
    $isStrandSpecific = 1;
}

my %specialCounters = ('__no_feature'=>1, '__ambiguous'=>1, '__too_low_aQual'=>1, '__not_aligned'=>1, '__alignment_not_unique'=>1);

my %geneLengths;
open(IN, "<$geneFootprintFile");
my $line = <IN>;
while ($line=<IN>) {
    chomp($line);
    my ($project, $gene, $length, @rest) = split(/\t/, $line);
    $geneLengths{$gene} = $length;
}
close(IN);

my $counts = 0;
open(IN, "<$countFile");
while (my $line=<IN>) {
    chomp($line);
    my ($id, $count) = split(/\t/, $line);
    if ($specialCounters{$id}) {
	next;
    }
    $counts += $count;
}
close(IN);

if ($isStrandSpecific) {
    open(IN, "<$antisenseCountFile");
    while (my $line=<IN>) {
	chomp($line);
	my ($id, $count) = split(/\t/, $line);
	if ($specialCounters{$id}) {
	    next;
	}
	$counts += $count;
    }
    close(IN);
}

print STDERR ("Total counts: $counts\n") if $verbose;

open(IN, "<$countFile");
open(OUT, ">$fpkmFile");

while (my $line=<IN>) {
    chomp($line);
    my ($id, $count) = split(/\t/, $line);
    if ($specialCounters{$id}) {
	next;
    }
    my $fpkm = $count*10**9/($counts*$geneLengths{$id});
    print OUT ("$id\t$fpkm\n");
}
close(IN);
close(OUT);

if ($isStrandSpecific) {
    open(IN, "<$antisenseCountFile");
    open(OUT, ">$antisenseFpkmFile");
    
    while (my $line=<IN>) {
	chomp($line);
	my ($id, $count) = split(/\t/, $line);
	if ($specialCounters{$id}) {
	    next;
	}
	my $fpkm = $count*10**9/($counts*$geneLengths{$id});
	print OUT ("$id\t$fpkm\n");
    }
    close(IN);
    close(OUT);
}
