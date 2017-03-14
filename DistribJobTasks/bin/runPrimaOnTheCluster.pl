#!@perl@


##Need to set up the controller_liniac.prop and task.prop files then run liniacsubmit

use strict;
use Getopt::Long;
use File::Basename;

my($prefix,$numNodes,$pwm,$targetFile,$bgFile,$masterDir,$promoterFile);

&GetOptions("targetFile|t=s" => \$targetFile, 
            "bgFile|b=s" => \$bgFile,
            "pwm|p=s" => \$pwm,
            "masterDir|m=s" => \$masterDir,
            "promoterFile|s=s" => \$promoterFile,
            "prefix=s" => \$prefix,
            "numNodes|n=i" => \$numNodes,
            );

$pwm = "/genomics/share/pkg/bio/Expander-PRIMA/PRIMA/$pwm" if $pwm =~ /^transfac\d\.dat/ || $pwm eq 'jaspar.mat' || $pwm eq 'jasparAndTransfac4.dat';

if(!-e $targetFile){print STDERR "TargetFile $targetFile not found\n";}
if(!-e $bgFile){print STDERR "Background file $bgFile not found\n";}
if(!-e $pwm){print STDERR "pwm $pwm not found\n";}
if(!-e $promoterFile){print STDERR "promoterFile $promoterFile not found\n";}

die "usage: runPrimaOnTheCluster.pl --targetFile|t <targetlist> --bgFile|b <bg list> --pwm|p <file of pwms> --promoterFile|s <file of promoter sequences> --numNodes|n <Number of nodes to run on..[numPWM / 20] --masterDir|m <directory to work in> --prefix <seq prefix to use for tmp file [NM_]>\n" unless (-e $targetFile && -e $bgFile && -e $pwm && -e $promoterFile);

#print STDERR "number of nodes requested: $numNodes\n";
$masterDir = 'master' unless $masterDir;

$prefix = "NM_" unless $prefix;

my $dir = $ENV{PWD};

my $numPwm = `grep -c // $pwm`;
chomp $numPwm;
my $nodeNum = $numNodes; 
my $subtasksize = 10;
if($numNodes){
  ##need to set the subtasksize...want to run only once / processor
  $subtasksize = int($numPwm / ($numNodes * 2)) + 1;
}else{
  $nodeNum = int(($numPwm / 20)+1);
}

##want to determine time...with the mus promoters, the algorithm is #bases/setSize/13000.
open(B,"$bgFile");
my %bg;
while(<B>){
  chomp;
  $bg{$_} = 1;
}
close B;
open(S, "$promoterFile");
my $p = 0;
my $seq = "";
my $bases = 0;
while(<S>){
  if(/^\>Locus=(\d+)/){
    $p = $bg{$1} ? 1 : 0;
    $bases += length($seq);
    $seq = "";
  }else{
    $seq .= $_ if $p;
  }
}
$bases += length($seq);
close S;

my $time = int(($bases / 1300000) * $subtasksize );

#print STDERR "Working directory: $dir\n";

##first remove the leftovers...
my $rmCmd = "/bin/rm -r $masterDir 0*";
system($rmCmd);

##now create the master directory
mkdir($masterDir) || die "Unable to create directory $masterDir\n";
chdir($masterDir);

##first the controller_liniac.prop

open(C,">controller_liniac.prop");
print C "masterdir=$dir/$masterDir/result\n";
print C "inputdir=$dir/$masterDir\n";
print C "nodedir=/scratch/user/$ENV{USER}\n";
print C "slotspernode=2\n";
print C "subtasksize=$subtasksize\n";
print C "taskclass=DJob::DistribJobTasks::RunPRIMA\n";
print C "nodeclass=DJob::DistribJob::BprocNode\n";
print C "restart=no\n";
close C;

##now the task.prop
open(T,">task.prop");
$pwm = "$dir/$pwm" unless $pwm =~ /^\//;
print T "transfacFile=$pwm\n";
print T "outputFile=prima.out\n";
$promoterFile = "$dir/$promoterFile" unless $promoterFile =~ /^\//;
print T "promoterFile=$promoterFile\n";
$bgFile = "$dir/$bgFile" unless $bgFile =~ /^\//;
print T "bgFile=$bgFile\n";
$targetFile = "$dir/$targetFile" unless $targetFile =~ /^\//;
print T "targetFile=$targetFile\n";
close T;


##now kick off liniacsubmit!
##first determine how many nodes to reserve

my $cmd = "liniacsubmit $nodeNum $time $dir/$masterDir/controller_liniac.prop";
print STDERR "\n$cmd\n";
system($cmd);


##now want to wait and run the code to summarize the results...
##first need to get the sequences...
open(F,"$targetFile");
my %n;
my $p = 0;
while(<F>){
  chomp;
  $n{$_} = 1;
} 
close F;
my $targetSeq = "targets.seq";
open(O,">$targetSeq");
open(F,"$promoterFile");
while(<F>){
  if(/^\>Locus=(\d+)/){
    $p = $n{$1} ? 1 : 0;
    print O ">$prefix$1\n" if $p;
  }else{
    print O if $p;
  }
}
close O;
close F;

##now wait until the run is complete and then run summary
#while(1){
#  sleep 300;
#  if(-e "$masterDir/result/mainresult/prima.out"){
#    my $fin = `grep -c Hyper $masterDir/result/mainresult/prima.out`;
#    if($fin >= ($numPwm - 2)){
#      system("~brunkb/bin/runGoodPrimaHits.pl --p $pwm --s $targetSeq --f $masterDir/result/mainresult/prima.out --n 20");
#      last;
#    }
#  }
#}
