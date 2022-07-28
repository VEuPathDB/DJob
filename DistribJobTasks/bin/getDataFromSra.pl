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
if ($workingDir){
  chdir($workingDir) or die "$!: could not change directory to $workingDir";
}
printf("$0 downloading fastqs to %s\n", getcwd());
my @tmp;
foreach my $s (split(/,\s*/,$sampleIdList)){
  push(@tmp,$s);
}

if($isColorspace){
    &getCsForSampleIds(\@tmp, "$readsOne", "$readsTwo", $doNotGetFastq);
    
}else{
    if ($workingDir =~ m/FungiDB|AmoebaDB|CryptoDB|GiardiaDB|PiroplasmaDB|TrichDB|VectorBase|HostDB|MicrosporidiaDB|PlasmoDB|TriTrypDB|ToxoDB/){
     &getFastqForSraRunId("$sampleIdList", $hasPairedEnds);
    }else{
     &getFastqForSampleIds(\@tmp, "$readsOne", "$readsTwo", $doNotGetFastq, $hasPairedEnds);
    }
}

