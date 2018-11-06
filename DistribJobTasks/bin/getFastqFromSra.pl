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

my $cwd = getcwd();
chdir($workingDir) if $workingDir;

if (defined $sampleIdList) {
  my @tmp;
  foreach my $s (split(/,\s*/,$sampleIdList)){
    push(@tmp,$s);
  }
  &getFastqForSampleIds(\@tmp, "$readsOne", "$readsTwo", $doNotGetFastq, $hasPairedEnds);
} elsif (defined $studyId) {
  print STDERR "calling CBIL::sra.pm";
  &getFastqForStudyId($studyId, $hasPairedEnds, $doNotGetFastq, $deflineVars);
}

chdir($cwd) if $workingDir;
