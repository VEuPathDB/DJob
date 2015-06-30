package DJob::DistribJobTasks::MultiBlastSimilarityTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use Cwd;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

# Similar to BlastSimilarityTask, but processes a set of fasta files that need to be blasted against themselves.  These are provided in a tar file.
# Only supports NCBI blast

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["blastVendor",   "wu",   "(wu | ncbi) [wu]"],
 ["blastBinDir",   "default",   "eg, /genomics/share/bin"],
 ["fastasTarPath",   "",     "full path to tar file containing fasta files"],
 ["dbType",          "",     "p or n (not nec. if cdd run)"],
 ["pValCutoff",      "1e-5",  "[1e-5]"],
 ["lengthCutoff",    "10",    "[10]"],
 ["percentCutoff",   "20",    "[20]"],
 ["blastProgram",    "",     "rpsblast if cdd | any wu-blast"],
 ["regex",           "'(\\S+)'",     "regex for id on defline after the >"],
 ["blastParamsFile", "",    "file holding blast params -relative to inputdir"],
 ["doNotExitOnBlastFailure", "no", "if 'yes' then prints error in output file rather than causing subtask to fail"],
 ["printSimSeqsFile", "no", "print output in format for sqlldr for similarsequences"]
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;

    my $blastBin = $self->getProperty("blastBinDir");
    my $fastasTarPath = $self->getProperty("fastasTarPath");
    my $dbType = $self->getProperty("dbType");
    my $blastProgram = $self->getProperty("blastProgram");
    my $blastVendor = $self->getProperty("blastVendor");

    die "blastBinDir $blastBin doesn't exist" unless ( $blastBin eq 'default' ||  -e $blastBin);
    die "fastasTarPath '$fastasTarPath' doesn't exist" unless -d "$fastasTarPath";
}

sub initNode {
    my ($self, $node, $inputDir) = @_;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $tarredDir = $self->getProperty("fastasTarPath");
    $tarredDir =~ /(.*)\.tar\.gz/ || die "property fastasTarPath must be a .tar.gz file";
    my $fastaDir = $1;

    &runCmd("tar -xzf $tarredDir");

    opendir(DIR, $fastaDir) || die "can't open inputDir '$fastaDir'\n";
    my @files = readdir(DIR);
    $self->{fastaFiles} = \@files;
    closedir(DIR);
    return scalar(@files);
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir, $subTask) = @_;

    die "initSubTask:  end must equal start.  start=$start end=$end\n" unless $start = $end;

    my $dbFilePath = $self->{fastaFiles}->[$start];
    my $blastBin = $self->getProperty("blastBinDir");
    my $dbType = $self->getProperty("dbType");
    &runCmd(($blastBin eq 'default' ? "" : "$blastBin/") . "formatdb -i $dbFilePath -p ".($dbType eq 'p' ? 'T' : 'F'));
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
    my $doNotExitOnBlastFailure = $self->getProperty("doNotExitOnBlastFailure");


    my $cmd =  "blastSimilarity  --blastBinDir $blastBin --database $dbFilePath --seqFile $dbFilePath --lengthCutoff $lengthCutoff --pValCutoff $pValCutoff --percentCutoff $percentCutoff --blastProgram $blastProgram --blastVendor $blastVendor --regex $regex --blastParamsFile $nodeExecDir/$blastParamsFile".($doNotExitOnBlastFailure =~ /yes/i ? " --doNotExitOnBlastFailure" : "");

    my $printSimSeqsFile = $self->getProperty("printSimSeqsFile");
    $cmd .= " --printSimSeqsFile" if $printSimSeqsFile eq 'yes';

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    $node->runCmd("cat $nodeExecDir/blastSimilarity.out >> $mainResultDir/blastSimilarity.out");
    return 1 if $node->getErr();
    $node->runCmd("cat $nodeExecDir/blastSimilarity.log >> $mainResultDir/blastSimilarity.log");
    return 0;
}


## concatenate files here so are in order
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;
##  print "BlastSimilarity:  Cleaning up server method called using node: ".$node->getNum()."\n" if $node;
}

1;
