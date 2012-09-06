#!/usr/bin/perl

use strict;
use Getopt::Long;
use CBIL::Util::Utils;

my($fastaFile,$mateA,$mateB,$bwaIndex,$strain,$delIntFiles,$bowtieIndex,$isColorspace);
my $out = "result";
my $varscan = "/genomics/eupath/eupath-tmp/software/VarScan/2.2.10/VarScan.jar";
my $gatk = "/genomics/eupath/eupath-tmp/software/gatk/1.5.31/GenomeAnalysisTK.jar";
my $consPercentCutoff = 60; ## use to generate consensus
my $snpPercentCutoff = 20; ## use for snps and indels
my $editDistance = 0.04;
my $workingDir = ".";

&GetOptions("fastaFile|f=s" => \$fastaFile, 
            "mateA|ma=s"=> \$mateA,
            "mateB|mb=s" => \$mateB,
            "outputPrefix|o=s" => \$out,
            "bwaIndex|b=s" => \$bwaIndex,
            "bowtieIndex|x=s" => \$bowtieIndex,
            "varscan|v=s" => \$varscan,
            "gatk|g=s" => \$gatk,
            "strain|s=s" => \$strain,
            "consPercentCutoff|cpc=s" => \$consPercentCutoff,
            "snpPercentCutoff|spc=s" => \$snpPercentCutoff,
            "editDistance|ed=s" => \$editDistance,
            "workingDir|w=s" => \$workingDir,
            "deleteIntermediateFiles!" => \$delIntFiles,
            "isColorspace!" => \$isColorspace,
            );

die "varscan jar file not found\n".&getParams() unless -e "$varscan";
die "mateA file not found\n".&getParams() unless -e "$mateA";
die "mateB file not found\n".&getParams() if ($mateB && !-e "$mateB");
die "fasta file not found\n".&getParams() unless -e "$fastaFile";
die "either bwa or bowtie2 indices must be specified\n".&getParams() unless (-e "$bwaIndex.amb" || -e "$bowtieIndex.1.bt2" || ($isColorspace && -e "$bowtieIndex.1.ebwt")); 
die "you must provide a strain\n".&getParams() unless $strain;
##should add in usage

$ENV{_JAVA_OPTIONS}="-Xmx2g";

$snpPercentCutoff = $consPercentCutoff unless $snpPercentCutoff;

open(L,">>$workingDir/runHTS_SNPs.log");

select L;
$| = 1;

print L "runHTS_SNPs.pl run starting ".&getDate()."\n\n";

## indices ..  bwa index -p <prefix for indices> <fastafile>  ##optional -c flag if colorspace

## do we want to check to see if bwa indices are present and if not, create?
## NO, would cause race condition.  workflow should create ahead of time so there for all datasets

## do want to make is so that can be restarted and picks up where left off.  look for output file presence of subsequent cmd and don't run if present.

my $tmpOut = $out . "_tmp";

my $cmd;

##aligning with Bowtie2
if( -e "$bowtieIndex.1.bt2"){  
  ##NOTE: need to remove path to my install once mark puts into place
  $cmd = "(~brunkb/software/bowtie/bowtie2-2.0.0-beta7/bowtie2 --end-to-end --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' -x $bowtieIndex -1 $mateA ".(-e "$mateB" ? "-2 $mateB " : "")."-S $workingDir/$tmpOut.sam) >& $workingDir/$tmpOut.bowtie.log";
  print L &getDate().": $cmd\n";
  if(-e "$workingDir/complete" || -e "$workingDir/$tmpOut.bam"){ print L "  succeeded in previous run\n\n";
  }else{ &runCmd($cmd); print L "\n"; }

### aligning with bowtie 1 if colorspace
}elsif($isColorspace && -e "$bowtieIndex.1.ebwt"){  
  die "if isColorspace=true you must provide qual files that are named exactly like the reads files but with .qual appended to the read file name" unless -e "$mateA.qual";
  if(-e "$mateB"){  ## pairedEnd
    $cmd = "(bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' $bowtieIndex -1 $mateA --Q1 $mateA.qual -2 $mateB --Q2 $mateB.qual > $workingDir/$tmpOut.sam) >& $workingDir/$tmpOut.bowtie.log";
  }else{  ##single end
    $cmd = "(bowtie -f -C -a -S -n 3 --best --strata --sam-RG 'ID:EuP' --sam-RG 'SM:TU114' --sam-RG 'PL:Illumina' $bowtieIndex $mateA -Q $mateA.qual > $workingDir/$tmpOut.sam) >& $workingDir/$tmpOut.bowtie.log";
  }
  print L &getDate().": $cmd\n";
  if(-e "$workingDir/complete" || -e "$workingDir/$tmpOut.bam"){ print L "  succeeded in previous run\n\n";
  }else{ &runCmd($cmd); print L "\n"; }

##aligning with BWA
}elsif(-e "$bwaIndex.amb"){  
  $cmd = "(bwa aln -t 4 -n $editDistance $bwaIndex $mateA > $workingDir/$tmpOut.mate1.sai) >& $workingDir/$tmpOut.bwa_aln_mate1.err";
  print L &getDate().": $cmd\n";
  if(-e "$workingDir/complete" || -e "$workingDir/$tmpOut.sam"){ print L "  succeeded in previous run\n\n";
   }else{ &runCmd($cmd); print L "\n"; }
  
  if(-e "$mateB"){
    $cmd = "(bwa aln -t 4 -n $editDistance $bwaIndex $mateB > $workingDir/$tmpOut.mate2.sai) >& $workingDir/$tmpOut.bwa_aln_mate2.err";
    print L &getDate().": $cmd\n";
    if(-e "$workingDir/complete" || -e "$workingDir/$tmpOut.sam"){ print L "  succeeded in previous run\n\n";
   }else{ &runCmd($cmd); print L "\n"; }
    
    $cmd = "(bwa sampe -r '\@RG\tID:EuP\tSM:$strain\tPL:Illumina' $bwaIndex $workingDir/$tmpOut.mate1.sai $workingDir/$tmpOut.mate2.sai $mateA $mateB > $workingDir/$tmpOut.sam) >& $workingDir/$tmpOut.bwa_sampe.err";
    print L &getDate().": $cmd\n";
    if(-e "$workingDir/complete" || -e "$workingDir/$tmpOut.bam"){ print L "  succeeded in previous run\n\n";
    }else{ &runCmd($cmd); print L "\n"; }
  }else{
    print L &getDate().": Aligning in single end mode only\n";
    
    $cmd = "(bwa samse -r '\@RG\tID:EuP\tSM:$strain\tPL:Illumina' $bwaIndex $workingDir/$tmpOut.mate1.sai $mateA > $workingDir/$tmpOut.sam) >& $workingDir/$tmpOut.bwa_samse.err";
    print L &getDate().": $cmd\n";
    if(-e "$workingDir/complete" || -e "$workingDir/$tmpOut.bam"){ print L "  succeeded in previous run\n\n";
    }else{ &runCmd($cmd); print L "\n"; }
  }
}

$cmd = "samtools faidx $fastaFile";
print L &getDate().": $cmd\n";
if(-e "$fastaFile.fai"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print L "\n"; }

$cmd = "(samtools view -t $fastaFile.fai -uS $workingDir/$tmpOut.sam | samtools sort - $workingDir/$tmpOut) >& $workingDir/$tmpOut.samtools_view.err";
print L &getDate().": $cmd\n";
if(-e "$workingDir/complete" || -e "$workingDir/$tmpOut.bam.bai"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print L "\n"; }

$cmd = "samtools index $workingDir/$tmpOut.bam";
print L &getDate().": $cmd\n";
if(-e "$workingDir/complete" || -e "$workingDir/forIndelRealigner.intervals"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print L "\n"; }

$cmd = "java -Xmx2g -jar $gatk -I $workingDir/$tmpOut.bam -R $fastaFile -T RealignerTargetCreator -o $workingDir/forIndelRealigner.intervals >& $workingDir/realignerTargetCreator.err";
print L &getDate().": $cmd\n";
if(-e "$workingDir/complete" || -e "$workingDir/$out.bam"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print L "\n"; }

$cmd = "java -Xmx2g -jar $gatk -I $workingDir/$tmpOut.bam -R $fastaFile -T IndelRealigner -targetIntervals $workingDir/forIndelRealigner.intervals -o $workingDir/$out.bam >& $workingDir/indelRealigner.err";
print L &getDate().": $cmd\n";
if(-e "$workingDir/complete" || -e "$workingDir/$out.pileup"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print L "\n"; }

## I wonder if would be good to sort again here?

$cmd = "(samtools mpileup -f $fastaFile -B $workingDir/$out.bam > $workingDir/$out.pileup) &> $workingDir/$out.pileup.err";
print L &getDate().": $cmd\n";
if(-e "$workingDir/complete" || -e "$workingDir/$out.varscan.snps"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print L "\n"; }

my $pc = $consPercentCutoff / 100;
my $mpc = $snpPercentCutoff / 100;

$cmd = "(java -Xmx2g -jar $varscan pileup2snp $workingDir/$out.pileup --p-value 0.01 --min-coverage 5 --min-var-freq $mpc > $workingDir/$out.varscan.snps ) >& $workingDir/$out.varscan_snps.err";
print L &getDate().": $cmd\n";
if(-e "$workingDir/complete" || -e "$workingDir/$out.SNPs.gff"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print L "\n"; }

$cmd = "parseVarscanToGFF.pl --f $workingDir/$out.varscan.snps --strain $strain --pc $snpPercentCutoff --o $workingDir/$out.SNPs.gff >& $workingDir/$out.parseVarscan.err";
print L &getDate().": $cmd\n";
if(-e "$workingDir/complete" || -e "$workingDir/$out.varscan.cons"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print L "\n"; }

$cmd = "(java -Xmx2g -jar $varscan pileup2cns $workingDir/$out.pileup --p-value 0.01 --min-coverage 5 --min-var-freq $pc > $workingDir/$out.varscan.cons ) >& $workingDir/$out.varscan_cons.err";
print L &getDate().": $cmd\n\n";
if(-e "$workingDir/complete" || -e "$workingDir/$out.consensus.fa"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print "\n"; }

##now parse to generate the consensus fasta file and gff file of inserts.
$cmd = "parseVarscanToConsensus.pl --file $workingDir/$out.varscan.cons --strain $strain --referenceFasta $fastaFile --fastaOutput $workingDir/$out.consensus.fa --indelOutput $workingDir/$out.insertions.gff --percentCutoff $consPercentCutoff >& parseToConsensus.stderr";
print L &getDate().": $cmd\n\n";
if(-e "$workingDir/complete"){ print L "  succeeded in previous run\n\n";
}else{ &runCmd($cmd); print "\n"; }

print L &getDate().": run COMPLETE\n\n";

&runCmd("echo complete > $workingDir/complete");

##should cleanup unneeded files before exiting so don't transfer too much back
## can delete all $tmpOut* and .err files for starters.

if($delIntFiles){
  print L "deleting extra files\n";
  system("/bin/rm $workingDir/$tmpOut.*");
  system("/bin/rm $workingDir/*.err");
  unlink("$workingDir/result.pileup");
  unlink("$workingDir/result.varscan.snps");
  unlink("$workingDir/result.varscan.cons");
  unlink("$workingDir/forIndelRealigner.intervals");
}

close L;

sub getParams {
  return &getDate().": runBWA_HTS.pl ... parameter values:\n\tfastaFile=$fastaFile\n\tbwaIndex=$bwaIndex\n\t\tOR\n\tbowtieIndex=$bowtieIndex\n\tmateA=$mateA\n\tmateB=$mateB\n\toutputPrefix=$out\n\tstrain=$strain\n\tconsPercentCutoff=$consPercentCutoff\n\tsnpPercentCutoff=$snpPercentCutoff\n\teditDistance=$editDistance\n\tvarscan=$varscan\n\tgatk=$gatk\n\tworkingDir=$workingDir\n\n";
}

sub getDate {
  my $date = `date`;
  chomp $date;
  return $date;
}
