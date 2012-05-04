package DJob::DistribJobTasks::BlastSimilarityTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["blastVendor",   "wu",   "(wu | ncbi) [wu]"],
 ["blastBinDir",   "",   "eg, /genomics/share/bin"],
 ["dbFilePath",      "",     "full path to database file"],
 ["inputFilePath",   "",     "full path to input file"],
 ["dbType",          "",     "p or n (not nec. if cdd run)"],
 ["pValCutoff",      "1e-5",  "[1e-5]"],
 ["lengthCutoff",    "10",    "[10]"],
 ["percentCutoff",   "20",    "[20]"],
 ["blastProgram",    "",     "rpsblast if cdd | any wu-blast"],
 ["regex",           "'(\\S+)'",     "regex for id on defline after the >"],
 ["blastParamsFile", "",    "file holding blast params -relative to inputdir"],
 ["saveGoodBlastFiles",   "no",    "If yes then blast results that meet parse are saved"],
 ["blastFileDirPath",   "$ENV{HOME}/blastFiles",    "Must specify a directory to save blast files into if saveGoodBlastFiles=yes [$ENV{HOME}/blastFiles]"],
 ["doNotExitOnBlastFailure", "no", "if 'yes' then prints error in output file rather than causing subtask to fail"],
 ["copyDbToNode", "no", "(yes | [no]) if 'yes' then copies the blast indices to the local nodeDir on node ... may be faster in some contexts"],
 ["printSimSeqsFile", "no", "print output in format for sqlldr for similarsequences"]
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;

    ##deal with saving blast files if desired
    if($self->getProperty("saveGoodBlastFiles") =~ /yes/i){
      ##check to make certain directory exists and if doesn't, create it
      my $fileDir = $self->getProperty("blastFileDirPath");
      print STDERR "Saving good blast files to $fileDir\n";
      if(!-d "$fileDir"){
        mkdir($fileDir) || die "directory '$fileDir' for storing blast files can not be created\n";
      }
    }


    my $blastBin = $self->getProperty("blastBinDir");
    my $dbFilePath = $self->getProperty("dbFilePath");
    my $dbType = $self->getProperty("dbType");
    my $blastProgram = $self->getProperty("blastProgram");
    my $blastVendor = $self->getProperty("blastVendor");

    if (-e "$dbFilePath.gz") {
	&runCmd("gunzip $dbFilePath.gz");
    }

    die "blastBinDir $blastBin doesn't exist" unless -e $blastBin;
    die "dbFilePath $dbFilePath doesn't exist" unless -e "$dbFilePath";

    # run if we don't have indexed files or they are older than seq file
    if ($blastProgram eq 'rpsblast') {  
      my $cwd = &getcwd();
      chdir(dirname($dbFilePath));
      my @ls = `ls -rt $dbFilePath.mn $dbFilePath.rps`;
      map { chomp } @ls;
      if (scalar(@ls) != 2 || $ls[0] ne "$dbFilePath.mn") {
        &runCmd("$blastBin/copymat -r T -P $dbFilePath");
      }
      
      @ls = `ls -rt $dbFilePath $dbFilePath.p*`;
      map { chomp } @ls;
      if (scalar(@ls) != 6 || $ls[0] ne $dbFilePath) {
        &runCmd("$blastBin/formatdb -o T -i $dbFilePath");
      }
      chdir $cwd;
      
   } elsif($blastVendor eq "ncbi"){  ##need to format as per ncbi formatdb ...
     my @ls = `ls -rt $dbFilePath $dbFilePath.$dbType*`;
     map { chomp } @ls;
     if (scalar(@ls) < 4 || $ls[0] ne $dbFilePath) {
       &runCmd("$blastBin/formatdb -i $dbFilePath -p ".($dbType eq 'p' ? 'T' : 'F'));
     }
   } else {
     my @ls = `ls -rt $dbFilePath $dbFilePath.x$dbType*`;
     map { chomp } @ls;
     if (scalar(@ls) != 4 || $ls[0] ne $dbFilePath) {
       &runCmd("cat $dbFilePath | perl -pe 'unless (/^>/){s/J/X/g;}' > ${dbFilePath}.replaceJ");
       &runCmd("rm -fr $dbFilePath");
       &runCmd("cat ${dbFilePath}.replaceJ | perl -pe 'unless (/^>/){s/O/X/g;}' > $dbFilePath");
       &runCmd("rm -fr ${dbFilePath}.replaceJ");
       &runCmd("$blastBin/xdformat -$dbType $dbFilePath");
     }
   } 
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

    return if $self->getProperty("copyDbToNode") eq 'no';

    my $blastProgram = $self->getProperty("blastProgram");
    my $dbFilePath = $self->getProperty("dbFilePath");
    my $nodeDir = $node->getDir();

    $node->runCmd("cp $dbFilePath.* $nodeDir");
    if ($blastProgram eq 'rpsblast') {  
	$node->runCmd("cp $dbFilePath $nodeDir");
    }
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->getProperty("inputFilePath");

    if (-e "$fastaFileName.gz") {
	&runCmd("gunzip $fastaFileName.gz");
    }

    $self->{fastaFile} = CBIL::Bio::FastaFileSequential->new($fastaFileName);
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;
#    print STDERR "start=$start, end=$end\n";

    if(!$subTask->getRedoSubtask()){
      my $blastParamsFile = $self->getProperty("blastParamsFile");
      &runCmd("cp $inputDir/$blastParamsFile $serverSubTaskDir");
      $self->{fastaFile}->writeSeqsToFile($start, $end, "$serverSubTaskDir/seqsubset.fsa");
    }
    $node->runCmd("touch $serverSubTaskDir/seqsubset.fsa.touch",1);
    $node->runCmd("/bin/rm $serverSubTaskDir/seqsubset.fsa.touch",1);

    $node->runCmd("cp -r $serverSubTaskDir/* $nodeExecDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $blastBin = $self->getProperty("blastBinDir");
    my $blastVendor = $self->getProperty("blastVendor");
    my $lengthCutoff = $self->getProperty("lengthCutoff");
    my $pValCutoff = $self->getProperty("pValCutoff");
    my $percentCutoff = $self->getProperty("percentCutoff");
    my $blastProgram = $self->getProperty("blastProgram");
    my $regex = $self->getProperty("regex");
    my $blastParamsFile = $self->getProperty("blastParamsFile");
    my $dbFilePath = $self->getProperty("dbFilePath");
    my $saveGood = $self->getProperty("saveGoodBlastFiles");
    my $blastFilePath = $self->getProperty("blastFileDirPath");
    my $doNotExitOnBlastFailure = $self->getProperty("doNotExitOnBlastFailure");
    my $dbFile = $self->getProperty("copyDbToNode") eq 'no' ? $dbFilePath : $node->getDir() . "/" . basename($dbFilePath);


    my $cmd =  "blastSimilarity  --blastBinDir $blastBin --database $dbFile --seqFile $nodeExecDir/seqsubset.fsa --lengthCutoff $lengthCutoff --pValCutoff $pValCutoff --percentCutoff $percentCutoff --blastProgram $blastProgram --blastVendor $blastVendor --regex $regex --blastParamsFile $nodeExecDir/$blastParamsFile".($saveGood =~ /yes/i ? " --saveGoodBlastFiles --blastFileDir $blastFilePath" : "").($doNotExitOnBlastFailure =~ /yes/i ? " --doNotExitOnBlastFailure" : "");

    my $printSimSeqsFile = $self->getProperty("printSimSeqsFile");
    $cmd .= " --printSimSeqsFile" if $printSimSeqsFile eq 'yes';

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    $node->runCmd("cat $nodeExecDir/blastSimilarity.out >> $mainResultDir/blastSimilarity.out");
    return 1 if $node->getErr();
    $node->runCmd("cat $nodeExecDir/blastSimilarity.log >> $mainResultDir/blastSimilarity.log");
}


## concatenate files here so are in order
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;

}

1;
