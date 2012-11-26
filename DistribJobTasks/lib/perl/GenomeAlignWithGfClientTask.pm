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
 ["gaBinPath",   "default",   "eg, /genomics/share/bin/blat"],
 ["targetDirPath",      "",     "full path to directory containing *.2bit file"],
 ["queryPath",   "",     "full path to input file"],
 ["queryType", "dna", "type of query ([dna]|prot)"],
 ["blatParams", "none", "additional params to be passed to blat"],
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
    die "gaBinPath $gaBinPath doesn't exist" unless ($gaBinPath eq 'default' || -e $gaBinPath);
    my $targetDirPath = $self->getProperty("targetDirPath");
    die "targetDirPath $targetDirPath doesn't exist" unless -e $targetDirPath;
    my $twoBitFile = "$targetDirPath/genomicSeqs.2bit";
    die "There is no 2bit file in $targetDirPath" unless -e  $twoBitFile;
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

    my $targetDirPath = $self->getProperty("targetDirPath");
    my $nodeDir = $node->getDir();

    my $gaBinPath = $self->getProperty("gaBinPath"); 
   
    my $queryType = $self->getProperty("queryType");


    my $cmd = "startGfServer --binPath $gaBinPath --nodePort $node->{gfport} --targetDir $targetDirPath".($queryType eq 'prot' ? " --trans" : "");

#    print STDERR "$cmd\n";
    print $node->runCmd($cmd)."\n";
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
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir,$subTask) = @_;

    if(!$subTask->getRedoSubtask()){
      $self->{fastaFile}->writeSeqsToFile($start, $end, "$subTaskDir/seqsubset.fsa");
    }
    $node->runCmd("touch $subTaskDir/seqsubset.fsa.touch",1);
    $node->runCmd("/bin/rm $subTaskDir/seqsubset.fsa.touch",1);

    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand {
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $gaBinPath = $self->getProperty("gaBinPath"); #the path of the gfClient script
    my $targetPath = $self->getProperty("targetDirPath"); #path of the dir with the .2bit file
    my $port = $node->{gfport};
    my $nodeDir = $node->getDir();
    my $maxIntron = $self->getProperty("maxIntron");
    my $queryType = $self->getProperty("queryType");
    my $dbType = $queryType eq 'dna' ? 'dna' : 'dnax';
    my $tmpBlatParams = $self->getProperty('blatParams');
    my $blatParams = $tmpBlatParams eq 'none' ? "" : $tmpBlatParams;

    my $cmd = ($gaBinPath eq 'default' ? "" : "$gaBinPath/")."gfClient -nohead -maxIntron=$maxIntron -t=$dbType -q=$queryType -dots=10 $blatParams localhost $port $targetPath seqsubset.fsa out.psl";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $node->runCmd("cat $nodeExecDir/out.psl >> $mainResultDir/out.psl");
    return 1 if $node->getErr();
}

sub cleanUpNode {
    my ($self,$node) = @_;

    my $port = $node->{gfport};
    my $gaBinPath = $self->getProperty('gaBinPath');

    $node->runCmd(($gaBinPath eq 'default' ? "" : "$gaBinPath/")."gfServer stop localhost $port");
}


1;
