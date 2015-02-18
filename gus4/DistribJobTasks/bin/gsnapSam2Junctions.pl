#!/usr/bin/perl

# Modified from sam2junctions.pl by Gregory R Grant (University of Pennsylvania, 2010)

use strict;

use Getopt::Long;

my ($minsize, $uniqueSamFile, $nuSamFile, $outputFile);

&GetOptions("min_size=i"=> \$minsize,
            "unique_sam=s" => \$uniqueSamFile,
            "nu_sam=s" => \$nuSamFile,
            "output_file=s" => \$outputFile,
    );

$minsize = 0 unless($minsize);;

unless(-e $uniqueSamFile && -e $nuSamFile) {
  die "usage:  gsnapSam2Junctions.pl --unique_sam=s --nu_sam=s [--min_size=i]\n";
}

my %junctions;


foreach("$uniqueSamFile|unique", "$nuSamFile|nu") {

  my ($samFile, $mapperType) = split(/\|/, $_);

  open(INFILE, $samFile) or die "\nError: Cannot open '$samFile' for reading\n\n";

  while(my $line = <INFILE>) {
    chomp($line);

    my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @tags) = split(/\t/, $line);

    next unless($cigar =~ /N/);
#    $sam_chr =~ s/:[^:]*$//;

    my $xs; 
    foreach my $tag (@tags) {
      next unless($tag =~ /XS:A:([+-])/);
      $xs = $1;
    }

    die "XS:A tag not found for read $qname" unless($xs);

    my $running_offset = 0;
    while($cigar =~ s/^(\d+)([^\d])//) {
      my $num = $1;
      my $type = $2;
      
	if($type eq 'N') {
          my$start = $pos + $running_offset - 1;
          my $end = $start + $num + 1;
          my $junction = "$rname:$start-$end";

          if($num >= $minsize) {
            $junctions{$junction}->{$xs}->{$mapperType}++;
          }
	}

      if($type eq 'N' || $type eq 'M' || $type eq 'D') {
        $running_offset = $running_offset + $num;
      }
    }
  }
  close(INFILE);
}


open(OUT, "> $outputFile") or die "Cannot open file $outputFile for writing: $!";


print OUT "Junction\tStrand\tUnique\tNU\n"; 

foreach my $junction (keys %junctions) {
  foreach my $xs (keys %{$junctions{$junction}}) {

    my $uniqueCount = $junctions{$junction}->{$xs}->{unique};
    my $nuCount = $junctions{$junction}->{$xs}->{nu};

    $uniqueCount = 0 unless($uniqueCount);
    $nuCount = 0 unless($nuCount);

    print OUT "$junction\t$xs\t$uniqueCount\t$nuCount\n";
  }
}

close OUT;

1;
