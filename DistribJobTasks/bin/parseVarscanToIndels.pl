#!@perl@

## parses varscan pileup2indel output to generate a gff file for loading into the database

use strict;
use Getopt::Long;
use lib "$ENV{GUS_HOME}/lib/perl";

my $file; 
my $strain;
my $indelOutput = '';
my $percentCutoff = 60;
my $coverage = 5;

&GetOptions("file|f=s" => \$file, 
            "percentCutoff|pc=i"=> \$percentCutoff,
            "coverage|c=i"=> \$coverage,
            "indelOutput|io=s"=> \$indelOutput,
            "strain|s=s"=> \$strain,
            );

if (! -e $file || !$strain || !$indelOutput ){
die &getUsage();
}

my $pc = $percentCutoff > 1 ? $percentCutoff / 100 : $percentCutoff;

open(G,">$indelOutput");

open(F, "$file") || die "unable to open file $file\n";
my $ct = 0;
while(<F>){
  next if /^Chrom\s+Position/;
  chomp;
  $ct++;
  print STDERR "  Processing $ct\n" if $ct % 10000 == 0;
  my @t = split("\t",$_);
  chop $t[6];
  if ($t[5] / ($t[4] + $t[5]) > $pc && $t[6] > 20){ 
    &generateGffIndel(\@t);
  }
}

close F;
close G;

sub generateGffIndel {
  my($l) = @_;
  my ($t,$allele) = ($l->[3] =~ /\/(.)(\w+)/);
  print STDERR "ERROR: unable to determine allele from ".join("\t",@$l)."\n" unless $allele;
  my $type = $t eq '+' ? 'insertion' : 'deletion';
  my $id = "NGS_indel.".$l->[0] .".".$l->[1];
  my $cov = $l->[4] + $l->[5];
  return unless $cov >= $coverage;  ##minimum coverage must be met.
  my $perc = (int($l->[5] / $cov * 1000) / 10) . '%';
  my $pvalue = $l->[11];
  my $qual = $l->[10];
  print G "$l->[0]\t$type\t$type\t$l->[1]\t$l->[1]\t.\t+\t.\tID $id; Allele \"$strain:$allele:$cov:$perc:$pvalue:$qual\"\n";

}

sub getUsage {
return <<endOfUsage;
parseVarscanToIndels.pl usage:

  parseVarscanToIndels.pl --file|f <varscan file> --strain <strain for indels>  --percentCutoff|pc <frequency percent cutoff [60]> --pvalueCutoff|pvc <pvalue cutoff [0.01]> --indelOutut|io <output file for indels in GFF format>
endOfUsage
}
