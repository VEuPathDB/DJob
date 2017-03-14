#!@perl@

use strict;
use Getopt::Long;

my($numGood,$annotFile,$pwmFile,$sequenceFile,$primaFile,$cutoff);

&GetOptions("pwmFile|p=s" => \$pwmFile, 
            "annotFile|a=s"=> \$annotFile,
            "sequenceFile|s=s" => \$sequenceFile,
            "cutoff|c=f" => \$cutoff,
            "numGood|n=i" => \$numGood,
            "primaFile|f=s" => \$primaFile,  ## in excel output from getGoodPrimaHits.pl -e
            );

$pwmFile = "/genomics/share/pkg/bio/Expander-PRIMA/PRIMA/$pwmFile" if $pwmFile =~ /^transfac\d\.dat/ || $pwmFile eq 'jaspar.mat' || $pwmFile eq 'jasparAndTransfac4.dat';

die "usage: runGoodPrimaHits.pl -p <pwmFile> -a <annotation file> -s <sequenceFile> -f <prima output file> --c <pValue cutoff> --n <number of prima hits>\n" unless(-e "$pwmFile" && -e "$sequenceFile" && -e "$primaFile");

##first get goodPrima hits
my $goodHitsCmd = "~brunkb/bin/getGoodPrimaHits.pl --e ".($cutoff ? "--c $cutoff " : "").($numGood ? "--n $numGood " : "")."$primaFile > goodHits.tab";
system($goodHitsCmd);

my $findElementsCmd = "~brunkb/bin/findMultipleTfElement.pl --s $sequenceFile --p $pwmFile --f goodHits.tab -a $annotFile > goodSites.tab"; 
system($findElementsCmd);

my $coocCmd = "~brunkb/bin/computePrimaOccurence.pl -f goodSites.tab -c 1.5 > goodHitsCooccurence.tab";
system($coocCmd);
