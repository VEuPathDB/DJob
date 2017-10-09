#!@perl@

use strict; 
use lib "$ENV{GUS_HOME}/lib/perl";
use CBIL::Util::Sra;
use Cwd;
use Getopt::Long;

my($doNotGetFastq,$workingDir,$readsOne,$readsTwo,$sampleIdList,$hasPairedEnds);

&GetOptions("doNotGetFastq!" => \$doNotGetFastq,
            "workingDir=s" => \$workingDir,
            "readsOne=s" => \$readsOne,
            "readsTwo=s" => \$readsTwo,
            "sampleIdList=s" => \$sampleIdList,
	    "pairs=s" => \$hasPairedEnds,
           );

my $cwd = getcwd();
chdir($workingDir) if $workingDir;

my @tmp;
foreach my $s (split(/,\s*/,$sampleIdList)){
  push(@tmp,$s);
}
&getFastqForSampleIds(\@tmp, "$readsOne", "$readsTwo", $doNotGetFastq, $hasPairedEnds);

chdir($cwd) if $workingDir;
