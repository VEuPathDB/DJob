#!@perl@

use strict; 
use lib "$ENV{GUS_HOME}/lib/perl";
use CBIL::Util::Sra;
use Cwd;
use Getopt::Long;

my($doNotGetFastq,$workingDir,$readsOne,$readsTwo,$sampleIdList, $isColorspace, $hasPairedEnds);

&GetOptions("doNotGetFastq!" => \$doNotGetFastq,
            "workingDir=s" => \$workingDir,
            "readsOne=s" => \$readsOne,
            "readsTwo=s" => \$readsTwo,
            "sampleIdList=s" => \$sampleIdList,
            "isColorspace!" => \$isColorspace,
	    "pairs=s" =>\$hasPairedEnds,
            );

my $cwd = getcwd();
chdir($workingDir) if $workingDir;

my @tmp;
foreach my $s (split(/,\s*/,$sampleIdList)){
  push(@tmp,$s);
}

if($isColorspace){
    &getCsForSampleIds(\@tmp, "$readsOne", "$readsTwo", $doNotGetFastq);
    
}else{
    &getFastqForSampleIds(\@tmp, "$readsOne", "$readsTwo", $doNotGetFastq, $hasPairedEnds);
}

chdir($cwd) if $workingDir;
