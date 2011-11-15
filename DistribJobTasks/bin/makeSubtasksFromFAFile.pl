#!@perl@

use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use Getopt::Long;

my $debug = 0;

$| = 1;

my ($fileName,$stdin,$directory,$subtaskSize);
my $outStem = "subTask.fa";

&GetOptions("fileName=s" => \$fileName, ##just one fasta file as input
            "outStem=s" => \$outStem,
            "stdin!" => \$stdin,        ##take fasta on stdin
            "subtaskSize=i" => \$subtaskSize,        ##take fasta on stdin
            "directory=s" => \$directory  ##directory for output
            );

die "invalid directory '$directory'" unless -d "$directory";
die "you must specify either --stdin or --fileName <inputfile>" unless ($stdin || -e $fileName);
die "you must specify --subtaskSize <size>" unless $subtaskSize;

my $ct = 0;
my $stnum = 1;
open(OUT,">$directory/$outStem.$stnum");
if($stdin){
  while(<STDIN>){
    if(/^\>seq/){
      $ct++;
      if($ct != 1 && $ct % $subtaskSize == 1){
        close(OUT);
        $stnum++;
        open(OUT,">$directory/$outStem.$stnum");
      }
    }
    print OUT $_;
  }
}else{
  open(F,"$fileName");
    while(<F>){
    if(/^\>/){
      if($ct != 1 && $ct++ % $subtaskSize == 1){
        close(OUT);
        $stnum++;
        open(OUT,">$directory/$outStem.$stnum");
      }
      $ct++;
    }
    print OUT $_;
  }
}
close(OUT);

##print total number of sequences
open(N,">$directory/sequenceCount");
print N "$ct\n";
close(N);
