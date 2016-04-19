package DJob::DistribJobTasks::BlastMatrixTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["blastBinDir",   "default",   "eg, /genomics/share/pkg/bio/ncbi-blast/latest"],
 ["dbFilePath",    "",   "full path to database file"],
 ["inputFilePath", "",   "full path to input file"],
 ["lengthCutoff",  "40",   ""],
 ["pValCutoff",    "1e-5", ""],
 ["percentCutoff", "92",   ""],
 ["endSlop",       "15",   ""],
 ["maskRepeats",   "n",  "y or n"],

 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}


# called once 
sub initServer {
    my ($self, $inputDir) = @_;

    my $blastBin = $self->getProperty("blastBinDir");
    my $dbFilePath = $self->getProperty("dbFilePath");

    if (-e "$dbFilePath.gz") {
	&runCmd("gunzip $dbFilePath.gz");
    }

    die "blastBinDir $blastBin doesn't exist" unless ( $blastBin eq 'default' ||  -e $blastBin);
    die "dbFilePath $dbFilePath doesn't exist" unless -e $dbFilePath;

    my @ls = `ls -rt $dbFilePath $dbFilePath.xn*`;
    map { chomp } @ls;  
    if (scalar(@ls) != 4 || $ls[0] ne $dbFilePath) {
	&runCmd(($blastBin eq 'default' ? "" : "$blastBin/") ."xdformat -n $dbFilePath");
    }
}

sub initNode {
    my ($self, $node, $inputDir) = @_;
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->getProperty("inputFilePath");

    if (-e "$fastaFileName.gz") {
	&runCmd("gunzip $fastaFileName.gz");
    }

    print "Counting sequences in $fastaFileName\n";
    $self->{fastaFile} = CBIL::Bio::FastaFileSequential->new($fastaFileName);
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir,$subTask) = @_;

    if(!$subTask->getRedoSubtask()){
      $self->{fastaFile}->writeSeqsToFile($start, $end, 
					"$subTaskDir/seqsubset.fsa");
    }

    $self->runCmdOnNode("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $blastBin = $self->getProperty("blastBinDir");
    my $lengthCutoff = $self->getProperty("lengthCutoff");
    my $pValCutoff = $self->getProperty("pValCutoff");
    my $percentCutoff = $self->getProperty("percentCutoff");
    my $endSlop = $self->getProperty("endSlop");
    my $maskRepeats = $self->getProperty("maskRepeats") eq "y"? "--maskRepeats" : "";
    my $dbFilePath = $self->getProperty("dbFilePath");

    my $cmd = "blastMatrix --blastBinDir $blastBin --db $dbFilePath --seqFile $nodeExecDir/seqsubset.fsa --lengthCutoff $lengthCutoff --pValCutoff $pValCutoff --percentCutoff $percentCutoff --endSlop $endSlop $maskRepeats";
   
    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $self->runCmdOnNode("cat $nodeExecDir/blastMatrix.out >> $mainResultDir/blastMatrix.out");
}
1;
