#!/usr/bin/perl

## parses varscan output to generate a fasta file for loading into the database

use strict;
use Getopt::Long;
use lib "$ENV{GUS_HOME}/lib/perl";
use CBIL::Bio::SequenceUtils;

my $file; 
my $strain;
my $referenceFasta;
my $faOutput = '';
my $percentCutoff = 60;
my $pvalueCutoff = .01;

&GetOptions("file|f=s" => \$file, 
            "percentCutoff|pc=i"=> \$percentCutoff,
            "pvalueCutoff|pvc=s"=> \$pvalueCutoff,
            "fastaOutput|fo=s"=> \$faOutput,
            "strain|s=s"=> \$strain,
            "referenceFasta|rf=s"=> \$referenceFasta,
            );

if (! -e $file || !$strain || !$faOutput || !-e $referenceFasta){
die &getUsage();
}

open(O,">$faOutput");

##first open fasta file and get lengths of all sequences so can fill in Ns at the end.
open(A, "$referenceFasta") || die &getUsage;
my $refLen = 0;
my $sid;
my %seqLen;
while(<A>){
  if(/^\>(\S+)/){
    $seqLen{$sid} = $refLen if $sid;
    $refLen = 0;
    $sid = $1;
  }else{
    $_ =~ s/\s//g;
    $refLen += length($_);
  }
}
print STDERR "Have lengths for ".scalar(keys%seqLen)." sequences\n";


my $pos = 0;
my $cons = "";
my $inDeletion = 0;
my $id;
my %ctGaps;
open(F, "$file") || die "unable to open file $file\n";
while(<F>){
  next if /^Chrom\s+Position/;
  chomp;
  print STDERR "  Processing $id: $pos\n" if $pos % 100000 == 0;

  my @tmp = split("\t",$_);
  print STDERR "ERROR: $pos = $tmp[1]: $_\n" if $pos == $tmp[1];
  if($inDeletion > 0){  ##in a deletion
    if($id eq $tmp[0] && $seqLen{$id} > $pos){
#      print STDERR " In deletion $tmp[1]\n";
      $pos += 1;
      $inDeletion -= 1;
      $cons .= '-';
      next;
    }
  }
  if($id ne $tmp[0]){  ##starting new sequence
    if($id){
      if($pos < $seqLen{$id}){
        $cons .= &fillInNs($seqLen{$id} - $pos);
      }
      print O CBIL::Bio::SequenceUtils::makeFastaFormattedSequence("$id.$strain",$cons);  
#      print STDERR "$id has $ctGaps{$id} gaps\n";
    }
    $cons = "";
    $pos = 0;
    $id = $tmp[0];
    $inDeletion = 0;
    print STDERR "Processing $id\n";
  }
  if($pos + 1 < $tmp[1]){  ##need to fill in some Ns.
    $cons .= &fillInNs($tmp[1] - ($pos + 1));
    $ctGaps{$id}++;
#    print STDERR "$id: $pos ($tmp[1]) - filled in ".($tmp[1] - $pos - 1)." Ns\n";
  }
  if(length($tmp[3]) > 1){  ##deal with the indel here.  need to put something in the file but also create gff entry for this indel.
    ## note that if deletion then should just add - for each base deleted.
    my $b = $tmp[-1];
#    my($a,$b) = split("\/",$tmp[3]);
    if($b =~ /^-/){ ##deletion
      print STDERR "Deletion: $_\n";
      $inDeletion = length($b) - 1;
    }elsif($b =~ /^\+/){ ##insertion
      print STDERR "Insertion: $_\n";
    }else{
      print STDERR "\nERROR: not dealing with putative indel properly ... appending '$tmp[2]'\n$_\n\n";
    }
    $cons .= $tmp[2];
  }else{
    $cons .= $tmp[3];
  }
  $pos = $tmp[1];
}

##need to do the last one here!!
if($id){
  if($pos < $seqLen{$id}){
    $cons .= &fillInNs($seqLen{$id} - $pos);
  }
  print O CBIL::Bio::SequenceUtils::makeFastaFormattedSequence("$id.$strain",$cons);  
}

close F;
close O;

sub fillInNs {
  my($num) = @_;
  my $ret;
  for(my $a = 0;$a < $num; $a++){
    $ret .= 'N';
  }
  return $ret;
}

sub getUsage {
return <<endOfUsage;
parseVarscanToGFF.pl usage:

  parseVarscanToConsensus.pl --file|f <varscan file> --strain <strain for snps>  --referenceFasta <fasta file of reference sequence> --percentCutoff|pc <frequency percent cutoff [60]> --pvalueCutoff|pvc <pvalue cutoff [0.01]> --fastaOutput|fo <output File for fasta consensus> --indelOutut|io <output file for insertions in GFF format>
endOfUsage
}
