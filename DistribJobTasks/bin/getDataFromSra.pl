#!@perl@

use strict; 
use lib "$ENV{GUS_HOME}/lib/perl";
use CBIL::Util::Sra;
use Cwd;
use Getopt::Long;

my($doNotGetFastq,$workingDir,$readsOne,$readsTwo,$sampleIdList, $isColorspace);

&GetOptions("doNotGetFastq!" => \$doNotGetFastq,
            "workingDir=s" => \$workingDir,
            "readsOne=s" => \$readsOne,
            "readsTwo=s" => \$readsTwo,
            "sampleIdList=s" => \$sampleIdList,
            "isColorspace!" => \$isColorspace,
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
    &getFastqForSampleIds(\@tmp, "$readsOne", "$readsTwo", $doNotGetFastq);
}

chdir($cwd) if $workingDir;
