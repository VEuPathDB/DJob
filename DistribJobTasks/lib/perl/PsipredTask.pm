package DJob::DistribJobTasks::PsipredTask;

use CBIL::Util::Utils;
use CBIL::Bio::FastaFile;

use File::Basename;

use DJob::DistribJob::Task;

@ISA = (DJob::DistribJob::Task);

use strict;

my @properties = (["psipredDir", "", "full path to the psipred dir"],
                  [ "dbFilePath", "", "subject file path"],
                  [ "inputFilePath", "", "query file path"],
		  );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

sub initServer {
    my ($self, $inputDir) = @_;

    return(1);
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

    my $dbFilePath = $self->{props}->getProp("dbFilePath");
    my $nodeDir = $node->getDir();

    my $nrfilt = $dbFilePath . "*";

    $node->runCmd("cp $nrfilt $nodeDir/");
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
    my ($self, $start, $end, $node, $inputDir, $subTaskDir, $nodeSlotDir) = @_;

    $self->{fastaFile}->writeSeqsToFile($start, $end, "$subTaskDir/seqsubset.fsa");

    open(FILE, "$subTaskDir/seqsubset.fsa") || die "Cannot open file $subTaskDir/seqsubset.fsa for reading: $!";

    my $line = <FILE>;
    close(FILE);

    chomp($line);
    my ($newName) = $line =~ /\>([a-zA-Z0-9_\.]+) /;
    $newName = $newName . ".fsa";

    $node->runCmd("cp -r $subTaskDir/seqsubset.fsa $nodeSlotDir/$newName");
}


sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $runpsipred = $self->{props}->getProp("psipredDir") . "/runpsipred";
    my $nrFilt = $self->{props}->getProp("dbFilePath");

    my $dbFile = $node->getDir() . "/" . basename($nrFilt);

    my $cmd = "$runpsipred $nodeExecDir/*.fsa $dbFile ";
    print STDERR "command:\n$cmd\n\n";

    return $cmd;
}


sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $node->runCmd("mv $nodeExecDir/*.ss2 $mainResultDir/");
}


1;
