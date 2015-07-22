#!/usr/bin/perl


use strict;


## fuzztran -sequence mmtv.seq -pattern 'DYXX[LI]X(6,12)YXX[LI]' -mismatch 0 -frame 6 -outfile mmtv.fuzzpro

my $pattern = 'DXXYXX[LI]X(6,20)YXX[LI]';
#my $pattern = 'DYXX[LI]X(6,20)YXX[LI]';

my $input = shift;
open(IN, "$input");

open(O,">runFuzzTran.out");
foreach my $file (<IN>){
  chomp $file;
  
  if(!-e $file){
#    print STDERR "Unable to find file $file\n";
    next;
  }
  
  print STDERR "running fuzztran on $file\n";

  my $cmd = "/genomics/share/pkg/emboss-2.8.0/bin/fuzztran -sequence $file -pattern '$pattern' -mismatch 0 -frame 6 -outfile fuzzpro.out"; 
  print STDERR "$cmd\n";
  my $err = `$cmd`;
  
  my @loc;
  my $start = 0;
  my $hits = "";
  my $name = "";
  open(F,"fuzzpro.out");
  while(<F>){
    if(/Sequence:\s(\S+)/){
      $name = $1;
      $start = 1;
    }
    next unless $start;
    if(/^\s*(\d+)\s+(\d+).*\s(\S+)\s*$/){
      my($start,$end,$match) = ($1,$2,$3);
      next if $match =~ /(\*|X)/;
      push(@loc,[$start,$end,$match]);
    }
    $hits .= $_;
  }
  close F;
  
  if(scalar(@loc) > 0){
    print O "$hits\n\n";
    open(F,"$file");
    my $sequence = 0;
    while(<F>){
      next if /^\>/;
      chomp $_;
      $sequence .= $_;
    }
    foreach my $l (@loc){
      my $start = $l->[0] - 200;
      $start = $start < 0 ? 0 : $start;
      my $end = $l->[1] + 200;
      $end = $end > length($sequence) ? length($sequence) : $end;
      print O "\>$name\_$l->[0]-$l->[1] $file $start - $end $l->[2]\n".substr($sequence,$start,$end - $start)."\n\n";
    }
  }
}
close O;
close IN;
