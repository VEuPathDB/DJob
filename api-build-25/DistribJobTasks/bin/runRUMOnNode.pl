#!@perl@

use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use Getopt::Long;
use CBIL::Util::Utils;

open(LOG,">runRumOnNode.log");

#print STDERR "args:\n"; foreach my $a (@ARGV){ print STDERR " '$a'\n"; }

my $debug = 0;

$| = 1;

my ($genomeBowtieIndex,$readsFile,$genomeFastaFile,$transcriptBowtieIndex,
$geneAnnotationFile,$qualFile,$bowtieExec,$blatExec,$mdustExec,$perlScriptsDir,$fastBlat,
$limitNU,$minBlatIdentity,$countMismatches,$alignToTranscriptome,$inputType,$pairedEnd,$createSAMFile,$numInsertions,$subtaskNumber,$mainResultDir,$readlength,$minoverlap,$minlength,$match_length_cutoff,$strandSpecific);

&GetOptions("readsFile=s" => \$readsFile,
            "qualFile=s" => \$qualFile, 
            "readlength=s" => \$readlength, 
            "genomeFastaFile=s" => \$genomeFastaFile, 
            "genomeBowtieIndex=s" => \$genomeBowtieIndex, 
            "transcriptBowtieIndex=s" => \$transcriptBowtieIndex, 
            "geneAnnotationFile=s" => \$geneAnnotationFile, 
            "bowtieExec=s" => \$bowtieExec, 
            "blatExec=s" => \$blatExec, 
            "mdustExec=s" => \$mdustExec, 
            "perlScriptsDir=s" => \$perlScriptsDir, 
            "limitNU=i" => \$limitNU, 
            "minlength=s" => \$minlength, 
            "pairedEnd=s" => \$pairedEnd,
            "strandSpecific=s" => \$strandSpecific,
            "minBlatIdentity=i" => \$minBlatIdentity, 
            "numInsertions=i" => \$numInsertions, 
            "createSAMFile=i" => \$createSAMFile, 
            "mainResultDir=s" => \$mainResultDir, 
            "subtaskNumber=s" => \$subtaskNumber,
            "countMismatches=i" => \$countMismatches,
            "matchLengthCutoff=i" => \$match_length_cutoff 
            );


print LOG "starting: ".`date`; 

my $cmd = "$bowtieExec -a --best --strata -f $genomeBowtieIndex $readsFile -v 3 --suppress 6,7,8 -p 1 > genomeBowtie.out";
&runCmd($cmd);
## check for error and return proper error code if failed
print LOG "finished first bowtie run: ".`date` . "$cmd\n\n";

$cmd = "perl $perlScriptsDir/make_GU_and_GNU.pl genomeBowtie.out btGenomeU btGenomeNU $pairedEnd";
&runCmd($cmd);
print LOG "finished parsing genome bowtie run: ".`date` . "$cmd\n\n";

if($transcriptBowtieIndex){
  $cmd = "$bowtieExec -a --best --strata -f $transcriptBowtieIndex $readsFile -v 3 --suppress 6,7,8 -p 1 > transcriptBowtie.out";
  &runCmd($cmd);
  print LOG "finished second bowtie run: ".`date` . "$cmd\n\n";

  $cmd = "perl $perlScriptsDir/make_TU_and_TNU.pl transcriptBowtie.out $geneAnnotationFile btTranscriptU btTranscriptNU $pairedEnd";
  &runCmd($cmd);
  print LOG "finished parsing transcriptome bowtie run: ".`date` . "$cmd\n\n";
  
  if($minlength eq 'false') {
      $cmd = "perl $perlScriptsDir/merge_GU_and_TU.pl btGenomeU btTranscriptU btGenomeNU btTranscriptNU bowtieUnique combinedNU $pairedEnd -readlength $readlength";
  } else {
      $cmd = "perl $perlScriptsDir/merge_GU_and_TU.pl btGenomeU btTranscriptU btGenomeNU btTranscriptNU bowtieUnique combinedNU $pairedEnd -readlength $readlength -minoverlap $minlength";
  }
  &runCmd($cmd);
  print LOG "finished merging TU and GU: ".`date` . "$cmd\n\n";

  $cmd = "perl $perlScriptsDir/merge_GNU_and_TNU_and_CNU.pl btGenomeNU btTranscriptNU combinedNU bowtieNU";
  &runCmd($cmd);
  print LOG "finished merging GNU, TNU and CNU: ",`date` . "$cmd\n\n";
}else{
  &runCmd("mv btGenomeU bowtieUnique");
  &runCmd("mv btGenomeNU bowtieNU");
}

$cmd = "perl $perlScriptsDir/make_unmapped_file.pl $readsFile bowtieUnique bowtieNU blatInput.fa $pairedEnd";
&runCmd($cmd);
print LOG "finished making R: ".`date` . "$cmd\n\n";
$cmd = "$blatExec $genomeFastaFile blatInput.fa blat.out -minIdentity=$minBlatIdentity -minScore=20 -stepSize=5";
&runCmd($cmd);
print LOG "finished first BLAT run: ".`date` . "$cmd\n\n";

$cmd = "$mdustExec blatInput.fa > mdust.out";
&runCmd($cmd);
print LOG "finished running mdust on R: ".`date` . "$cmd\n\n";

if($minlength eq 'false') {
    $cmd = "perl $perlScriptsDir/parse_blat_out.pl blatInput.fa blat.out mdust.out blatUnique blatNU".($transcriptBowtieIndex ? "" : " -dna")." -num_insertions_allowed $numInsertions";
} else {
    $cmd = "perl $perlScriptsDir/parse_blat_out.pl blatInput.fa blat.out mdust.out blatUnique blatNU -match_length_cutoff $minlength ".($transcriptBowtieIndex ? "" : " -dna")." -num_insertions_allowed $numInsertions";
}
&runCmd($cmd);
print LOG "finished parsing first BLAT run: ".`date` . "$cmd\n\n";

if($minlength eq 'false') {
    $cmd = "perl $perlScriptsDir/merge_Bowtie_and_Blat.pl bowtieUnique blatUnique bowtieNU blatNU merge_Unique_temp merge_NU_temp $pairedEnd -readlength $readlength";
} else {
    $cmd = "perl $perlScriptsDir/merge_Bowtie_and_Blat.pl bowtieUnique blatUnique bowtieNU blatNU merge_Unique_temp merge_NU_temp $pairedEnd -readlength $readlength -minoverlap $minlength";
}
&runCmd($cmd);
print LOG "finished merging Bowtie and Blat: ".`date` . "$cmd\n\n";

if($minlength eq 'false') {
    $cmd = "perl $perlScriptsDir/RUM_finalcleanup.pl merge_Unique_temp merge_NU_temp merge_Unique_temp2 merge_NU_temp2 $genomeFastaFile "." $mainResultDir/sam_header.$subtaskNumber ".($countMismatches ? " -countmismatches" : "");
} else {
    $cmd = "perl $perlScriptsDir/RUM_finalcleanup.pl merge_Unique_temp merge_NU_temp merge_Unique_temp2 merge_NU_temp2 $genomeFastaFile "." $mainResultDir/sam_header.$subtaskNumber ".($countMismatches ? " -countmismatches" : "") . " -match_length_cutoff $minlength";
}
&runCmd($cmd);
print LOG "finished cleaning up final results: ".`date` . "$cmd\n\n";

$cmd = "perl $perlScriptsDir/sort_RUM_by_id.pl merge_NU_temp2 merge_NU_idsorted";
&runCmd($cmd);

$cmd = "perl $perlScriptsDir/removedups.pl merge_NU_idsorted merge_NU_temp3 merge_Unique_temp2";
&runCmd($cmd);

if($limitNU > 0) {
    $cmd = "perl $perlScriptsDir/limit_NU.pl merge_NU_temp3 $limitNU > $mainResultDir/RUM_NU.$subtaskNumber";
} else {
    $cmd = "mv merge_NU_temp3 $mainResultDir/RUM_NU.$subtaskNumber";
}
&runCmd($cmd);

$cmd = "perl $perlScriptsDir/sort_RUM_by_id.pl merge_Unique_temp2 $mainResultDir/RUM_Unique.$subtaskNumber";
&runCmd($cmd);
print LOG "finished sorting final results: ".`date` . "$cmd\n\n";

##Note am writing these final results directly back to the server ...
$cmd = "perl $perlScriptsDir/limit_NU.pl merge_NU_temp2 $limitNU > merge_NU_temp3";
print LOG "finished limiting NU: ".`date` . "$cmd\n\n";
&runCmd($cmd);

if($createSAMFile){
$cmd = "perl $perlScriptsDir/rum2sam.pl $mainResultDir/RUM_Unique.$subtaskNumber $mainResultDir/RUM_NU.$subtaskNumber $readsFile $qualFile $mainResultDir/RUM_sam.$subtaskNumber";
  &runCmd($cmd);
  print LOG "finished creating SAM file: ".`date` . "$cmd\n\n";
}

$cmd = "perl $perlScriptsDir/sort_RUM_by_location.pl $mainResultDir/RUM_Unique.$subtaskNumber $mainResultDir/RUM_Unique.sorted.$subtaskNumber > $mainResultDir/chr_counts_u.$subtaskNumber";
&runCmd($cmd);
$cmd = "perl $perlScriptsDir/sort_RUM_by_location.pl $mainResultDir/RUM_NU.$subtaskNumber $mainResultDir/RUM_NU.sorted.$subtaskNumber > $mainResultDir/chr_counts_nu.$subtaskNumber";
&runCmd($cmd);

if($transcriptBowtieIndex) {
    if($strandSpecific eq 'true') {
	$cmd = "perl $perlScriptsDir/rum2quantifications.pl $geneAnnotationFile $mainResultDir/RUM_Unique.sorted.$subtaskNumber $mainResultDir/RUM_NU.sorted.$subtaskNumber $mainResultDir/quant.ps.$subtaskNumber -countsonly -strand p";
	&runCmd($cmd);
	$cmd = "perl $perlScriptsDir/rum2quantifications.pl $geneAnnotationFile $mainResultDir/RUM_Unique.sorted.$subtaskNumber $mainResultDir/RUM_NU.sorted.$subtaskNumber $mainResultDir/quant.ms.$subtaskNumber -countsonly -strand m";
	&runCmd($cmd);
	$cmd = "perl $perlScriptsDir/rum2quantifications.pl $geneAnnotationFile $mainResultDir/RUM_Unique.sorted.$subtaskNumber $mainResultDir/RUM_NU.sorted.$subtaskNumber $mainResultDir/quant.pa.$subtaskNumber -countsonly -strand p -anti";
	&runCmd($cmd);
	$cmd = "perl $perlScriptsDir/rum2quantifications.pl $geneAnnotationFile $mainResultDir/RUM_Unique.sorted.$subtaskNumber $mainResultDir/RUM_NU.sorted.$subtaskNumber $mainResultDir/quant.ma.$subtaskNumber -countsonly -strand m -anti";
	&runCmd($cmd);
    } else {
	$cmd = "perl $perlScriptsDir/rum2quantifications.pl $geneAnnotationFile $mainResultDir/RUM_Unique.sorted.$subtaskNumber $mainResultDir/RUM_NU.sorted.$subtaskNumber $mainResultDir/quant.$subtaskNumber -countsonly";
	&runCmd($cmd);
    }
}

print LOG "runRUMOnNode.pl complete: ".`date`;

close LOG;

