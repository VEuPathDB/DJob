package DJob::DistribJobTasks::GenomeAlignWithGfClientTask;

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
 ["gaBinPath",   "",   "eg, /genomics/share/bin/blat"],
 ["targetDirPath",      "",     "full path to directory containing *.nib files"],
 ["queryPath",   "",     "full path to input file"],
 ["nodePort", "", "port used on port for gfServer and gfClient"],
 ["queryType", "dna", "type of query ([dna]|prot)"],
 ["maxIntron", "", "the maximum length allowed for gaps that correspond to introns"]
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;
    my $gaBinPath = $self->getProperty("gaBinPath");
    die "gaBinPath $gaBinPath doesn't exist" unless -e $gaBinPath;
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

    my $targetDirPath = $self->getProperty("targetDirPath");
    my $nodeDir = $node->getDir();

    my $gaBinPath = $self->getProperty("gaBinPath"); 
   
    my $port = $self->getProperty("nodePort");
    my $queryType = $self->getProperty("queryType");

    $node->runCmd("startGfServer --binPath $gaBinPath --nodePort $port --targetDir $targetDirPath".($queryType eq 'prot' ? " --trans" : ""));

}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->getProperty("queryPath");

    if (-e "$fastaFileName.gz") {
	&runCmd("gunzip $fastaFileName.gz");
    }

    print "Creating index for $fastaFileName (may take a while)\n";
    $self->{fastaFile} = CBIL::Bio::FastaFileSequential->new($fastaFileName);

    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir) = @_;

    $self->{fastaFile}->writeSeqsToFile($start, $end, "$subTaskDir/seqsubset.fsa");

    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand {
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $gaBinPath = $self->getProperty("gaBinPath"); #the path of the gfClient script
    my $targetPath = $self->getProperty("targetDirPath"); #path of the dir with the .nib files
    my $port = $self->getProperty("nodePort");
    my $nodeDir = $node->getDir();
    my $maxIntron = $self->getProperty("maxIntron");
    my $queryType = $self->getProperty("queryType");
    my $dbType = $queryType eq 'dna' ? 'dna' : 'dnax';

    my $cmd = "${gaBinPath}/gfClient -nohead -maxIntron=$maxIntron -t=$dbType -q=$queryType -dots=10 localhost $port $targetPath seqsubset.fsa out.psl";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $node->runCmd("cat $nodeExecDir/out.psl >> $mainResultDir/out.psl");
}

sub cleanUpNode {
    my ($self,$node) = @_;

    my $port = $self->getProperty('nodePort');

    $node->runCmd($self->getProperty('gaBinPath')."/gfServer stop localhost $port");
}


1;
