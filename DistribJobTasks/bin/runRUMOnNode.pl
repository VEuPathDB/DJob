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
$limitNU,$minBlatIdentity,$countMismatches,$alignToTranscriptome,$inputType,$pairedEnd,$createSAMFile,$numInsertions,$subtaskNumber,$mainResultDir);

&GetOptions("readsFile=s" => \$readsFile, ##??just one fasta file as input??
            "qualFile=s" => \$qualFile, 
            "genomeFastaFile=s" => \$genomeFastaFile, 
            "genomeBowtieIndex=s" => \$genomeBowtieIndex, 
            "transcriptBowtieIndex=s" => \$transcriptBowtieIndex, 
            "geneAnnotationFile=s" => \$geneAnnotationFile, 
            "bowtieExec=s" => \$bowtieExec, 
            "blatExec=s" => \$blatExec, 
            "mdustExec=s" => \$mdustExec, 
            "perlScriptsDir=s" => \$perlScriptsDir, 
            "limitNU=i" => \$limitNU, 
            "pairedEnd=s" => \$pairedEnd,
            "minBlatIdentity=i" => \$minBlatIdentity, 
            "numInsertions=i" => \$numInsertions, 
            "createSAMFile=i" => \$createSAMFile, 
            "mainResultDir=s" => \$mainResultDir, 
            "subtaskNumber=s" => \$subtaskNumber, 
            "countMismatches=i" => \$countMismatches 
            );


print LOG "starting: ".`date`; 

&runCmd("$bowtieExec -a --best --strata -f $genomeBowtieIndex $readsFile -v 3 --suppress 6,7,8 -p 1 > genomeBowtie.out");
## check for error and return proper error code if failed
print LOG "finished first bowtie run: ".`date`;

&runCmd("perl $perlScriptsDir/make_GU_and_GNU.pl genomeBowtie.out btGenomeU btGenomeNU $pairedEnd");
print LOG "finished parsing genome bowtie run: ".`date`;

if($transcriptBowtieIndex){
  &runCmd("$bowtieExec -a --best --strata -f $transcriptBowtieIndex $readsFile -v 3 --suppress 6,7,8 -p 1 > transcriptBowtie.out");
  print LOG "finished second bowtie run: ".`date`;

  &runCmd("perl $perlScriptsDir/make_TU_and_TNU.pl transcriptBowtie.out $geneAnnotationFile btTranscriptU btTranscriptNU $pairedEnd");
  print LOG "finished parsing transcriptome bowtie run: ".`date`;
  
  &runCmd("perl $perlScriptsDir/merge_GU_and_TU.pl btGenomeU btTranscriptU btGenomeNU btTranscriptNU bowtieUnique combinedNU $pairedEnd");
  print LOG "finished merging TU and GU: ".`date`;

  &runCmd("perl $perlScriptsDir/merge_GNU_and_TNU_and_CNU.pl btGenomeNU btTranscriptNU combinedNU bowtieNU");
  print LOG "finished merging GNU, TNU and CNU: ",`date`;
}else{
  &runCmd("mv btGenomeU bowtieUnique");
  &runCmd("mv btGenomeNU bowtieNU");
}

&runCmd("perl $perlScriptsDir/make_unmapped_file.pl $readsFile bowtieUnique bowtieNU blatInput.fa $pairedEnd");
print LOG "finished making R: ".`date`;

&runCmd("$blatExec $genomeFastaFile blatInput.fa blat.out -minIdentity=$minBlatIdentity");
print LOG "finished first BLAT run: ".`date`;

&runCmd("$mdustExec blatInput.fa > mdust.out");
print LOG "finished running mdust on R: ".`date`;

&runCmd("perl $perlScriptsDir/parse_blat_out.pl blatInput.fa blat.out mdust.out blatUnique blatNU".($transcriptBowtieIndex ? "" : " -dna")." -num_insertions_allowed $numInsertions");
print LOG "finished parsing first BLAT run: ".`date`;

&runCmd("perl $perlScriptsDir/merge_Bowtie_and_Blat.pl bowtieUnique blatUnique bowtieNU blatNU merge_Unique_temp merge_NU_temp $pairedEnd");
print LOG "finished merging Bowtie and Blat: ".`date`;

&runCmd("perl $perlScriptsDir/RUM_finalcleanup.pl merge_Unique_temp merge_NU_temp merge_Unique_temp2 merge_NU_temp2 $genomeFastaFile -faok".($countMismatches ? " -countmismatches" : ""));
print LOG "finished cleaning up final results: ".`date`;

&runCmd("perl $perlScriptsDir/sort_RUM.pl merge_Unique_temp2 $mainResultDir/RUM_Unique.$subtaskNumber");
##Note am writing these final results directly back to the server ...
&runCmd("perl $perlScriptsDir/limit_NU.pl merge_NU_temp2 $limitNU > RUM_NU_temp3");
&runCmd("perl $perlScriptsDir/sort_RUM.pl RUM_NU_temp3 $mainResultDir/RUM_NU.$subtaskNumber");
print LOG "finished sorting final results: ".`date`;

if($createSAMFile){
  &runCmd("perl $perlScriptsDir/rum2sam.pl RUM_Unique RUM_NU $readsFile $qualFile $mainResultDir/RUM_sam.$subtaskNumber");
  print LOG "finished creating SAM file: ".`date`;
}

print LOG "runRUMOnNode.pl complete: ".`date`;

close LOG;

