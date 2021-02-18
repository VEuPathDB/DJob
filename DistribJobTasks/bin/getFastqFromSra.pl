#!@perl@

use strict; 
use lib "$ENV{GUS_HOME}/lib/perl";
use CBIL::Util::Sra;
use Cwd;
use Getopt::Long;

my($doNotGetFastq,$workingDir,$readsOne,$readsTwo,$sampleIdList,$hasPairedEnds,$studyId,$deflineVars, $sampleAndRunIdsPath);

&GetOptions("doNotGetFastq!" => \$doNotGetFastq,
            "workingDir=s" => \$workingDir,
            "readsOne=s" => \$readsOne,
            "readsTwo=s" => \$readsTwo,
            "sampleIdList=s" => \$sampleIdList,
	    "pairs=s" => \$hasPairedEnds,
            "studyId=s" => \$studyId,
            "deflineVars=s" => \$deflineVars,
            "sampleAndRunIdsPath=s" => \$sampleAndRunIdsPath,
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
} elsif (defined $sampleAndRunIdsPath){
  my @runs;
  open(my $fh, "<", $sampleAndRunIdsPath) or die "$!: $sampleAndRunIdsPath"; 
  while(<$fh>){
    chomp;
    my ($srs, $srr, $libraryLayout) = split "\t";
    die "Bad line: $_ - Need sample then run, at least" unless $srs && $srr;
    $libraryLayout //= $hasPairedEnds ? 'PAIRED' : 'SINGLE';
    push @runs, [$srs, $srr, $libraryLayout];
  }
  die "No runs in $sampleAndRunIdsPath" ? unless @runs;
  CBIL::Util::Sra::getFastqForStudyRuns(\@runs, $hasPairedEnds, $doNotGetFastq, $deflineVars);

}
