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
		  [ "ncbiBinDir", "", "fullpath to the ncbi bin dir"],
		  );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

sub initServer {
    my ($self, $inputDir) = @_;

    my $ncbiBinDir = $self->getProperty("ncbiBinDir");
    my $dbFilePath = $self->getProperty("dbFilePath");

    my @ls = `ls -rt $dbFilePath.p*`;
    unless (scalar(@ls) > 0) {
	print "Formatting database (this may take a while)\n";
	&runCmd("$ncbiBinDir/formatdb -i $dbFilePath -p T");
    }

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
    $self->{basename}->{$node->getNum()}->{$nodeSlotDir} = $newName; #store the specific base filename 
    $newName = $newName . ".fsa";

    $node->runCmd("cp -r $subTaskDir/seqsubset.fsa $nodeSlotDir/$newName");
}


sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $runpsipred = $self->{props}->getProp("psipredDir") . "/runpsipred";
    my $nrFilt = $self->{props}->getProp("dbFilePath");

    my $dbFile = $node->getDir() . "/" . basename($nrFilt);
    my $inputFile = $self->{basename}->{$node->getNum()}->{$nodeExecDir} . ".fsa";

    my $cmd = "$runpsipred $nodeExecDir/$inputFile $dbFile ";
    print STDERR "command:\n$cmd\n\n";

    return $cmd;
}


sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    my $outputFile = $self->{basename}->{$node->getNum()}->{$nodeExecDir} . ".ss2";

    ##if move fails want to return 1 so that the task will fail the subtask 
    ##check to see that outputfile exists using ls
    return 1 unless $node->runCmd("ls $nodeExecDir/$outputFile",1);
    $node->runCmd("mv $nodeExecDir/$outputFile $mainResultDir/");
}


1;
