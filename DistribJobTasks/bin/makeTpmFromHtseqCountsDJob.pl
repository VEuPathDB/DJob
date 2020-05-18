#!@perl@

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my ($verbose,$geneFootprintFile,$countFile,$tpmFile,$antisenseCountFile,$antisenseTpmFile);
&GetOptions("verbose!"=>\$verbose,
            "geneFootprintFile=s"=> \$geneFootprintFile,
            "countFile=s" => \$countFile,
            "tpmFile=s" => \$tpmFile,
            "antisenseCountFile=s" => \$antisenseCountFile,
            "antisenseTpmFile=s" => \$antisenseTpmFile,
    ); 

if(!$geneFootprintFile || !$countFile || !$tpmFile){
	die "usage: makeTpmFromhtseqCountsDJob.pl --geneFootprintFile <geneFootprintFile> --countFile <input file (sense if strand-specific)> --tpmFile <output file (sense if strand-specific)> [--antisenseCountFile <input file (antisense if strand-specific)> --antisenseTpmFile <output file (antisense if strand-specific)>]\n";
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

&doTPMCalculation ($geneLengths, $countFile, $antisenseCountFile, $tpmFile, $antisenseTpmFile);


sub _calcRPK {
    my %specialCounters = ('__no_feature'=>1, '__ambiguous'=>1, '__too_low_aQual'=>1, '__not_aligned'=>1, '__alignment_not_unique'=>1);
    my ($geneLengths, $countFile, $rpkSum) = @_;
    my $rpkHash;
    open (IN, "<$countFile") or die "Cannot open file $countFile. Please check and try again\n$!\n";
    while (<IN>) {
        my ($geneId, $count) = split /\t/, $_;
        chomp($count);
        next if ($specialCounters{$geneId});
        my $geneLength = $geneLengths->{$geneId}/1000;
        my $rpk = $count/$geneLength;
        $rpkSum += $rpk;
        $rpkHash->{$geneId} = $rpk;
    }
    close IN;
    return ($rpkSum, $rpkHash);
}

sub _calcTPM {
    my ($rpkHash, $rpkSum) = @_;
    $rpkSum = $rpkSum/1000000;
    my $tpmHash;
    while (my($geneId, $rpk) = each %{$rpkHash}) {
        my $tpm =  $rpk/$rpkSum;
        $tpmHash->{$geneId} = $tpm;
    }
    return $tpmHash;
}

sub _writeTPM {
    my ($tpmFile, $tpmHash) = @_;
    open (OUT, ">$tpmFile") or die "Cannot  open TPM file $tpmFile for writing. Please check and try again.\n$!\n";
    while (my($geneId, $tpm) = each %{$tpmHash}) {
        print OUT ("$geneId\t$tpm\n");
    }
    close OUT;
}

sub doTPMCalculation {
    my ($geneLengths, $countFile, $antisenseCountFile, $tpmFile, $antisenseTpmFile) = @_;
    my ($rpkSum, $senseRpkHash) = &_calcRPK($geneLengths, $countFile, 0);
    if ($antisenseCountFile) {
        if ($antisenseTpmFile) {
            ($rpkSum, my $antisenseRpkHash) = &_calcRPK($geneLengths, $antisenseCountFile, $rpkSum);
            my $antisenseTpmHash = &_calcTPM($antisenseRpkHash, $rpkSum);
            &_writeTPM($antisenseTpmFile, $antisenseTpmHash);
        } else {
            die "An antisense count file $antisenseCountFile has been provided, but no antisense TPM file has been specified for writing output. Please add this using the --antisenseTpmFile flag.\n";
        }
    }

    my $tpmHash = &_calcTPM($senseRpkHash, $rpkSum);
    &_writeTPM($tpmFile, $tpmHash);
}
