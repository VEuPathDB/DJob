package DJob::DistribJobTasks::ArrayOligoSelectorTask;

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
 ["dbFilePath",      "",     "full path to database file"],
 ["dbType",          "genomic",     "database type (rna|genomic)"],
 ["inputFilePath",   "",     "full path to input file"],
 ["alignmentAlg",    "blat",    "algorithm for comparing to genome (blast|blat|gfclient)"],
 ["maskByLowerCase", "no", "ignore sequence that has been masked and is now lower case when picking oligos (no|yes)"],
 ["length",   "70",    "length of oligos"],
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
# I think simplest here if nothing is done....arrayoligoselector creates the blast indices itself
# if necessary the first time is run...need to check if recreates each subsequent time, if so
# then do something to decrease the time for this...

# NOTE: the getInputSize method is called before the initServer method!!

sub initServer {
    my ($self) = @_;
    return 1;
}

sub initNode {
    my ($self, $node, $inputDir) = @_;
}

## Note that this method is also doing something that could be done in initServer..ie
## creating the index files for the input database...inherited from Steve
## NOTE:  
sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->{props}->getProp("inputFilePath");

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
      $self->{fastaFile}->writeSeqsToFile($start, $end, "$subTaskDir/seqsubset.fsa");
    }

    $self->runCmdOnNode($node, "cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $binPath = '/genomics/share/pkg/bio/ArrayOligoSelector/latest';
    my $script = $self->{props}->getProp("dbType") =~ /^r/i ? 'Pick70_script1' : 'Pick70_script1_contig';
    my $dbFile = $self->{props}->getProp("dbFilePath");
    my $alignAlg = $self->{props}->getProp("alignmentAlg");
    my $length = $self->{props}->getProp("length");
    my $mask = $self->getProperty("maskByLowerCase");

    my $cmd = "bash $binPath/$script $nodeExecDir/seqsubset.fsa $dbFile $length $mask".($script =~ /contig/ ? " $alignAlg" : "");
#    print STDERR "command:\n$cmd\n\n";

    return $cmd;

}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $self->runCmdOnNode($node, "cat $nodeExecDir/output* >> $mainResultDir/output0");
}
1;
