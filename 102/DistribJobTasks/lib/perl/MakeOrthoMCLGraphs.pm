package DJob::DistribJobTasks::MakeOrthoMCLGraphs;

## Brian Brunk (ApiDB)
## NOTE: the subtaskSize must be = 1!
## configFile must be created .... will not be created in the script

use DJob::DistribJob::Task;
use CBIL::Util::Utils;
use File::Basename;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["orthomclExe",    "",     "directory containing the orthoMCL executables"],
 ["genomesFile",    "",     "full path to the genomes file"],
 ["bpoFile",   "",   "full path to the genomes file"],
 ["dataDir",   "",   "directory where in_paralog files and output files are found/written"],
 );
my @commands;

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    ##check to make certain the subtasksize = 1
    die "the subtasksize MUST = 1\n" unless $self->{subTaskSize} == 1;
    return $self;
}

##if able, check to see if the bpo.idx file exists and create if not .... for now create ahead of time
sub initServer {
    my ($self, $inputDir) = @_;
    my $bpoFile = $self->getProperty('bpoFile');
    my $fn = basename($bpoFile);
    my $idxFile = $fn;
    $idxFile =~ s/\./_/g;
    $idxFile = "$idxFile.idx";
    my $dir = dirname($bpoFile);
    die "you must index the bpoFile before running this task\n" unless (-e "$dir/$idxFile");
    $self->{bpoIdxFile} = "$dir/$idxFile";
}

## copy the bpo and bpo.idx files to the nodes
## also copy the orthomcl.conf file to the node as this won't change
sub initNode {
    my ($self, $node, $inputDir) = @_;
    my $nodeDir = $node->getDir();
    my $bpoFile = $self->getProperty('bpoFile');
    $node->runCmd("cp $bpoFile $nodeDir");
    $node->runCmd("cp $self->{bpoIdxFile} $nodeDir");
}

# note that this is called before initServer so need to create file here ...
## generate the file of commands (genomes to be compared)... just one pair / line space delimited
## unless clever, also generate the files that will be copied to the nodes for the genomes if needed
sub getInputSetSize {
  my ($self, $inputDir) = @_;
  ##first prepare the input files
  my $taskDir = "$inputDir/taskFiles";
  if(!-d "$taskDir"){
    mkdir("$taskDir");
  }
  my $genomeFile = $self->getProperty('genomesFile');
  my $numGenomes = `wc -l $genomeFile`;
  chomp $numGenomes;
  die "you must provide at least two genomes in the genome file\n" unless $numGenomes >= 2;
  if(-e "$inputDir/subtask.list" && `wc -l $inputDir/subtask.list` == ($numGenomes * ($numGenomes - 1))/2){
    print STDERR "Files and subtask.list already created\n";
  }else{  #need to create the files ...
    print STDERR "Creating files and subtask.list\n";
    open(F,"$genomeFile");
    my @gg = <F>;
    close F;
    open(O,">$inputDir/subtask.list") || die "unable to open $inputDir/subtask.list for writing\n";
    for(my $a=0;$a < scalar(@gg) - 1;$a++){
      my($one) = ($gg[$a] =~ /^(\w+)/);
      for(my $b = $a+1;$b < scalar(@gg);$b++){
        my($two) = ($gg[$b] =~ /^(\w+)/);
        my $fn = "$one"."_"."$two.gg";
        print O "$fn\n";
        open(G,">$taskDir/$fn") || die "unable to open '$taskDir/$fn for writing\n";
        print G "$gg[$a]$gg[$b]";
        close G;
      }
    }
    close O;
  }
  
  open(F, "$inputDir/subtask.list");
  while(<F>){
    next if /^\s*\#/;
    chomp;
    push(@commands,$_);
  }
  $self->{commands} =  \@commands;
  return scalar(@{$self->{commands}});
}

#copy the correct file to the node ... will have already made it in the getInputSetSize method
sub initSubTask {
  my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir,$subTask) = @_;

  my $file = $commands[$start];
  system("echo $file > $subTaskDir/subtask.file"); ##so have record on server
  ##also need to create config file and copy to the node
  open(F,">$subTaskDir/orthomcl.conf");
  print F $self->getConfigString($nodeSlotDir,$node->getDir());
  close F;
  $node->runCmd("cp $inputDir/taskFiles/$file $nodeSlotDir/subtask.gg");
  $node->runCmd("cp $subTaskDir/orthomcl.conf $nodeSlotDir/");
}

sub makeSubTaskCommand { 
  my ($self, $node, $inputDir, $nodeExecDir) = @_;
  
  my $orthomcl = $self->getProperty('orthomclExe');
  my $cmd = "$orthomcl $nodeExecDir/orthomcl.conf";
  return $cmd;
}

##nothing to be done here as writing back onto server directly
sub integrateSubTaskResults {
  my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
  return 1;
}

sub getConfigString {
  my($self,$slotDir,$nodeDir) = @_;
  my $bpoFile = basename($self->getProperty('bpoFile'));
  my $dataDir = $self->getProperty('dataDir');
  return "
#DIR is your working dir, bpo and gg files should be put within this dir
DIR=$slotDir
GG_FILE=$slotDir/subtask.gg
BPO_FILE=$nodeDir/$bpoFile
#BPO_IDX_MODE is the mode to index your BPO file. Options: all or taxon
BPO_IDX_MODE=all

#PVALUE_CUTOFF=1e-5
#PIDENT_CUTOFF=0
PMATCH_CUTOFF=50

#MAX_WEIGHT=316
MAX_WEIGHT=181

#MCL=/Users/fengchen/mcl-02-063/shmcl/mcl
#MCL=/files/cbil2/orthoMCL/program/orthomcl_engine/mcl-02-063/shmcl/mcl
MCL=/home/praveenc/orthomcl_engine_new/mcl-02-063/shmcl/mcl
INFLATION=1.5

# here is where the in_paralog files and results are stored ....
FORMER_RUN_SUBDIR=$dataDir
";

}
1;
