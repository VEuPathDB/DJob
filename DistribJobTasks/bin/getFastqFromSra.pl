#!@perl@

use strict; 
use lib "$ENV{GUS_HOME}/lib/perl";
use CBIL::Util::Sra;
use Cwd;
use Getopt::Long;

my($doNotGetFastq,$workingDir,$readsOne,$readsTwo,$sampleIdList,$hasPairedEnds,$studyId,$deflineVars);

&GetOptions("doNotGetFastq!" => \$doNotGetFastq,
            "workingDir=s" => \$workingDir,
            "readsOne=s" => \$readsOne,
            "readsTwo=s" => \$readsTwo,
            "sampleIdList=s" => \$sampleIdList,
	    "pairs=s" => \$hasPairedEnds,
            "studyId=s" => \$studyId,
            "deflineVars=s" => \$deflineVars,
           );

if ($workingDir){
  chdir($workingDir) or die "$!: cannot change dir to $workingDir";
}
printf("$0 downloading fastqs to %s\n", getcwd());

if (defined $sampleIdList) {
  my @tmp;
  foreach my $s (split(/,\s*/,$sampleIdList)){
    push(@tmp,$s);
  }
  CBIL::Util::Sra::getFastqForSampleIds(\@tmp, "$readsOne", "$readsTwo", $doNotGetFastq, $hasPairedEnds);
} elsif (defined $studyId) {
  CBIL::Util::Sra::getFastqForStudyId($studyId, $hasPairedEnds, $doNotGetFastq, $deflineVars);
}
