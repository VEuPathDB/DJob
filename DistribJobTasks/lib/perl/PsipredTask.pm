package DJob::DistribJobTasks::PsipredTask;

use CBIL::Util::Utils;
use CBIL::Bio::FastaFileSequential;
use CBIL::Util::Utils;

use File::Basename;

use DJob::DistribJob::Task;

@ISA = (DJob::DistribJob::Task);

use strict;

my @properties = (["psipredDir", "default", "full path to the psipred dir"],
                  [ "dbFilePath", "", "subject file path"],
                  [ "inputFilePath", "", "query file path"],
		  [ "ncbiBinDir", "", "fullpath to the ncbi bin dir"],
		  );

sub new {
use CBIL::Util::Utils;
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
      $self->{fastaFile}->writeSeqsToFile($start, $end, "$subTaskDir/seqsubset.fsa");
    }

    open(FILE, "$subTaskDir/seqsubset.fsa") || die "Cannot open file $subTaskDir/seqsubset.fsa for reading: $!";
    
    my $line = <FILE>;
    close(FILE);

    chomp($line);
    my ($newName) = $line =~ /\>([a-zA-Z0-9_\.]+) /;

    $newName = $newName . ".fsa";
    $self->runCmdOnNode("touch $subTaskDir/seqsubset.fsa");

    $self->runCmdOnNode("cp -r $subTaskDir/seqsubset.fsa $nodeSlotDir/$newName");
}


sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $dir = $self->{props}->getProp("psipredDir");
    $dir = $dir eq 'default' ? getProgramDir("runpsipred_single") : $dir;
    
    my $runpsipred =  "$dir/runpsipred_single";
    my $nrFilt = $self->{props}->getProp("dbFilePath");

    my $dbFile = $nrFilt  . "/" . basename($nrFilt);

    my $cmd = "$runpsipred $nodeExecDir/*.fsa $dbFile ";
    print STDERR "command:\n$cmd\n\n";

    return $cmd;
}


sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    ##check to see that outputfile exists using ls
    if ($self->runCmdOnNode("ls $nodeExecDir/*.ss2",1;) {

    	$self->runCmdOnNode("cp $nodeExecDir/*.ss2 $mainResultDir/");
	}
}


1;
