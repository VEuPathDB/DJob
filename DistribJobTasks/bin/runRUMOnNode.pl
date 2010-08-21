#!@perl@


BEGIN {
  unshift(@INC,"/genomics/share/lib/perl");
}

use lib "$ENV{GUS_HOME}/lib/perl";
use CBIL::Bio::Blast::BlastAnal;
use strict;
use Getopt::Long;
use File::Basename;

open(LOG,">runRumOnNode.log");

#print STDERR "args:\n"; foreach my $a (@ARGV){ print STDERR " '$a'\n"; }

my $debug = 0;

$| = 1;

my ($genomeBowtieIndex,$readsFile,$genomeFastaFile,$transcriptBowtieIndex,
$geneAnnotationFile,$qualFile,$bowtieExec,$blatExec,$mdustExec,$perlScriptsDir,$fastBlat,
$limitNU,$minBlatIdentity,$countMismatches,$alignToTranscriptome,$inputType,$pairedEnd,$createSAMFile);

&GetOptions("readsFile=s" => \$readsFile, ##??just one fasta file as input??
            "qualFile=s" => \$qualFile, 
            "genomeFastaFile=s" => \$genomeFastaFile, 
            "genomeBowtieIndex=s" => \$genomeBowtieIndex, 
            "transcriptBowtieIndex=s" => \$transcriptBowtieIndex, 
            "geneAnnotationFile=s" => \$geneAnnotationFile, 
            "bowtieExec=s" => \$bowtieBinDir, 
            "blatExec=s" => \$blatExec, 
            "mdustExec=s" => \$mdustExec, 
            "perlScriptsDir=s" => \$perlScriptsDir, 
            "limitNU=i" => \$limitNU, 
            "pairedEnd=i" => \$pairedEnd,
            "minBlatIdentity=i" => \$minBlatIdentity, 
            "numInsertions=i" => \$numInsertions, 
            "createSAMFile=i" => \$createSAMFile, 
            "countMismatches=i" => \$countMismatches 
            );


print LOG "starting: ".`date`; 

system("$bowtieExec -a --best --strata -f $genomeBowtieIndex $readsFile -v 3 --suppress 6,7,8 -p 1 > genomeBowtie.out")
## check for error and return proper error code if failed
&handleError($?,"ERROR running genome bowtie");
print LOG "finished first bowtie run: ".`date`;

system("perl $perlScriptsDir/make_GU_and_GNU.pl genomeBowtie.out btGenomeU btGenomeNU $pairedEnd");
&handleError($?,"ERROR parsing genome bowtie");
print LOG "finished parsing genome bowtie run: ".`date`;

if($transcriptBowtieIndex){
  system("$bowtieExec -a --best --strata -f $transcriptBowtieIndex $readsFile -v 3 --suppress 6,7,8 -p 1 > transcriptBowtie.out");
  &handleError($?,"ERROR running transcript bowtie");
  print LOG "finished second bowtie run: ".`date`;

  system("perl $perlScriptsDir/make_TU_and_TNU.pl transcriptBowtie.out $geneAnnotationFile btTranscriptU btTranscriptNU $pairedEnd");
  &handleError($?,"ERROR parsing transcriptome bowtie run");
  print LOG "finished parsing transcriptome bowtie run: ".`date`;
  
  system("perl $perlScriptsDir/merge_GU_and_TU.pl btGenomeU btTranscriptU btGenomeNU btTranscriptNU bowtieUnique combinedNU $pairedEnd");
  &handleError($?,"ERROR merging TU and GU");
  print LOG "finished merging TU and GU: ".`date`;

  system("perl $perlScriptsDir/merge_GNU_and_TNU_and_CNU.pl btGenomeNU btTranscriptNU combinedNU bowtieNU");
  &handleError($?,"ERROR merging GNU, TNU and CNU");
  print LOG "finished merging GNU, TNU and CNU: ",`date`;
}else{
  system("mv btGenomeU bowtieUnique");
  system("mv btGenomeNU bowtieNU");
}

system("perl $perlScriptsDir/make_unmapped_file.pl $readsFile bowtieUnique bowtieNU blatInput.fa $pairedEnd");
&handleError($?,"ERROR making R");
print LOG "finished making R: ".`date`;

system("$blatExec $genomeFastaFile blatInput.fa blat.out -minIdentity=$minBlatIdentity");
&handleError($?,"ERROR running first BLAT");
print LOG "finished first BLAT run: ".`date`;

system("$mdustExec blatInput.fa > mdust.out");
&handleError($?,"ERROR running mdust on R");
print LOG "finished running mdust on R: ".`date`;

system("perl $perlScriptsDir/parse_blat_out.pl blatInput.fa blat.out mdust.out blatUnique blatNU".($transcriptBowtieIndex ? "" : " -dna")." -num_insertions_allowed $numInsertions");
&handleError($?,"ERROR parsing first BLAT run");
print LOG "finished parsing first BLAT run: ".`date`;

system("perl $perlScriptsDir/merge_Bowtie_and_Blat.pl bowtieUnique blatUnique bowtieNU blatNU merge_Unique_temp merge_NU_temp $pairedEnd");
&handleError($?,"ERROR merging Bowtie and Blat");
print LOG "finished merging Bowtie and Blat: ".`date`;

system("perl $perlScriptsDir/RUM_finalcleanup.pl merge_Unique_temp merge_NU_temp merge_Unique_temp2 merge_NU_temp2 $genomeFastaFile".($countMismatches ? " -countmismatches" : ""));
&handleError($?,"ERROR cleaning up final results");
print LOG "finished cleaning up final results: ".`date`;

system("perl $perlScriptsDir/sort_RUM.pl RUM_Unique_temp2 RUM_Unique");
&handleError($?,"ERROR with first sort_RUM.pl run");
system("perl $perlScriptsDir/limit_NU.pl RUM_NU_temp2 $limitNU > RUM_NU_temp3");
&handleError($?,"ERROR with limit_NU.pl script");
system("perl $perlScriptsDir/sort_RUM.pl RUM_NU_temp3 RUM_NU");
&handleError($?,"ERROR with second sort_RUM.pl run");
print LOG "finished sorting final results: ".`date`;

if($createSAMFile){
  system("$perlScriptsDir/rum2sam.pl RUM_Unique RUM_NU $readsFile $qualFile RUM_sam");
  &handleError($?,"ERROR creating SAM file");
  print LOG "finished creating SAM file: ".`date`;
}

print LOG "runRUMOnNode.pl complete: ".`date`;

close LOG;

sub handleError {
  my($err,$msg) = @_;
  if($err){
    print LOG "$msg\nError Code: $err\n";
    exit($err);
  }
}
