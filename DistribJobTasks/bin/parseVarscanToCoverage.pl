#!/usr/bin/perl

use strict;
use Getopt::Long;
use lib "$ENV{GUS_HOME}/lib/perl";

my $file; 
my $out;
my $percentCutoff = 60;
my $coverageCutoff = 5;

&GetOptions("file|f=s" => \$file, 
            "percentCutoff|pc=i"=> \$percentCutoff,
            "coverageCutoff|cc=i"=> \$coverageCutoff,
            "outputFile|o=s"=> \$out,
            );

if (! -e $file || !$out){
die &getUsage();
}

if($file =~ /\.gz$/) {
  open(F, "zcat $file|") || die "unable to open file $file for reading: $!";
}
else {
  open(F, "$file") || die "unable to open file $file for reading: $!";
}


open(O,"|sort -k 1,1 -k 2,2n > $out") or die "Cannot open file for writing: $!";


my ($spanStart, $prevSeq, $prevLoc);

my @stats;
while(<F>){
  next if /^Chrom\s+Position/;
  chomp;

  my @a = split(/\t/, $_);

  chop $a[6]; # chop off the % sign
  my $varPercent = 100 - $a[6];
  my $coverage = $a[4];

  my $hasCoverage = $a[2] eq $a[3] && $varPercent >= $percentCutoff && $coverage >= $coverageCutoff;

  my $isSameSequence = $prevSeq ? $prevSeq eq $a[0]  : 1;

  # start span (doesn't matter what the sequence is 
  if(!$spanStart && $hasCoverage) {
    $spanStart = $a[1];
    push @stats, "$coverage:$varPercent";
  }
  # end span (sequence transition)
  elsif($spanStart && $hasCoverage && (!$isSameSequence || $prevLoc + 1 != $a[1])) {
    print O "$prevSeq\t$spanStart\t$prevLoc\t" . join(',', @stats) . "\n";
    $spanStart = $a[1];;
    @stats = ();
    push @stats, "$coverage:$varPercent";
  }
  # end span (no coverge )
  elsif($spanStart && !$hasCoverage) {
    print O "$prevSeq\t$spanStart\t$prevLoc\t" . join(',', @stats) . "\n";
    $spanStart = undef;
    @stats = ();
  }
  elsif($hasCoverage) { 
    push @stats, "$coverage:$varPercent";
  }
  else {}

  $prevSeq = $a[0];
  $prevLoc = $a[1];
}

if($spanStart) {
  print O "$prevSeq\t$spanStart\t$prevLoc\t" . join(',', @stats) . "\n";
}

close F;
close O;



sub getUsage {
return <<endOfUsage;
parseVarscanToCoverage.pl usage:

  parseVarscanToCoverage.pl --file|f <varscan file> --percentCutoff|pc <frequency percent cutoff [60]> --outputFile|o <output File for coverage> 
endOfUsage
}
