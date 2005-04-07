package DJob::DistribJobTasks::BlastSimilarityTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFile;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["blastBinDir",   "",   "eg, /genomics/share/bin"],
 ["dbFilePath",      "",     "full path to database file"],
 ["inputFilePath",   "",     "full path to input file"],
 ["dbType",          "",     "p or n (not nec. if cdd run)"],
 ["pValCutoff",      "1e-5",  ""],
 ["lengthCutoff",    "10",    ""],
 ["percentCutoff",   "20",    ""],
 ["blastProgram",    "",     "rpsblast if cdd | any wu-blast"],
 ["regex",           "'(\\S+)'",     "regex for id on defline after the >"],
 ["blastParamsFile", "",    "file holding blast params -relative to inputdir"],
 ["saveGoodBlastFiles",   "no",    "If yes then blast results that meet parse are saved"],
 ["blastFileDirPath",   "$ENV{HOME}/blastFiles",    "Must specify a directory to save blast files into if saveGoodBlastFiles=yes [$ENV{HOME}/blastFiles]"],
 ["doNotExitOnBlastFailure", "no", "if 'yes' then prints error in output file rather than causing subtask to fail"],
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

    if (-e "$dbFilePath.gz") {
	&runCmd("gunzip $dbFilePath.gz");
    }

    die "blastBinDir $blastBin doesn't exist" unless -e $blastBin;

    # run if we don't have indexed files or they are older than seq file
    if ($blastProgram eq 'rpsblast') {  
	die "dbFilePath $dbFilePath doesn't exist" unless -e "$dbFilePath";
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

    } else {
	die "dbFilePath $dbFilePath doesn't exist" unless -e $dbFilePath;
	my @ls = `ls -rt $dbFilePath $dbFilePath.x$dbType*`;
	map { chomp } @ls;
	if (scalar(@ls) != 4 || $ls[0] ne $dbFilePath) {
	    &runCmd("$blastBin/xdformat -$dbType $dbFilePath");
	}
    }
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

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

    print "Creating index for $fastaFileName (may take a while)\n";
    $self->{fastaFile} = CBIL::Bio::FastaFile->new($fastaFileName);
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir) = @_;

    my $blastParamsFile = $self->getProperty("blastParamsFile");
    &runCmd("cp $inputDir/$blastParamsFile $serverSubTaskDir");
    $self->{fastaFile}->writeSeqsToFile($start, $end, "$serverSubTaskDir/seqsubset.fsa");

    $node->runCmd("cp -r $serverSubTaskDir/* $nodeExecDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $blastBin = $self->getProperty("blastBinDir");
    my $lengthCutoff = $self->getProperty("lengthCutoff");
    my $pValCutoff = $self->getProperty("pValCutoff");
    my $percentCutoff = $self->getProperty("percentCutoff");
    my $blastProgram = $blastBin . "/" . $self->getProperty("blastProgram");
    my $regex = $self->getProperty("regex");
    my $blastParamsFile = $self->getProperty("blastParamsFile");
    my $dbFilePath = $self->getProperty("dbFilePath");
    my $saveGood = $self->getProperty("saveGoodBlastFiles");
    my $blastFilePath = $self->getProperty("blastFileDirPath");
    my $doNotExitOnBlastFailure = $self->getProperty("doNotExitOnBlastFailure");

    my $dbFile = $node->getDir() . "/" . basename($dbFilePath);

    my $cmd =  "blastSimilarity  --blastBinDir $blastBin --database $dbFile --seqFile $nodeExecDir/seqsubset.fsa --lengthCutoff $lengthCutoff --pValCutoff $pValCutoff --percentCutoff $percentCutoff --blastProgram $blastProgram --regex $regex --blastParamsFile $nodeExecDir/$blastParamsFile".($saveGood =~ /yes/i ? " --saveGoodBlastFiles --blastFileDir $blastFilePath" : "").($doNotExitOnBlastFailure =~ /yes/i ? " -doNotExitOnBlastFailure" : "");

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    $node->runCmd("cat $nodeExecDir/blastSimilarity.out >> $mainResultDir/blastSimilarity.out");
    $node->runCmd("cat $nodeExecDir/blastSimilarity.log >> $mainResultDir/blastSimilarity.log");
}
1;
