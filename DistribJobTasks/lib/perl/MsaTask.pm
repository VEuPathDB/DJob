package DJob::DistribJobTasks::MsaTask;

use CBIL::Util::Utils;
use DJob::DistribJob::Task;


@ISA = (DJob::DistribJob::Task);

use strict;


my @properties = 
    (
     ["muscleBinDir", "", "full path to muscle program"],
     ["inputFileDir", "", "directory containing tarballs of files, each file contains sequences in fasta format, to be aligned"]
     );


sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}


sub initServer {
    my ($self, $inputDir) = @_;

    my $muscleBinDir = $self->getProperty("muscleBinDir");

    my $inputFileDir = $self->getProperty("inputFileDir");


    die "muscleBinDir $muscleBinDir doesn't exist" unless -e $muscleBinDir;
    die "inputFileDir $inputFileDir doesn't exist" unless -e $inputFileDir;

}


sub initNode {
    my ($self, $node, $inputDir) = @_;

}


sub getInputSetSize {
  my ($self, $inputDir) = @_;

  my $inputFileDir = $self->getProperty("inputFileDir");

  my @fileArr;

  opendir(DIR, $inputFileDir) || die "Can't open directory $inputFileDir";

  while (defined (my $file = readdir (DIR))) {
    next if ($file eq "." || $file eq "..");
    push(@fileArr,$file);
  }

  my $count = @fileArr;

  $self->{fileArray} = @fileArr;

  return $count;
}


sub initSubTask {
    my ($self, $start, $stop, $node, $inputDir, $serverSubTaskDir,$nodeExecDir) = @_;


    my @files = @{$self->{fileArray}}[$start..$stop];

    foreach my $file (@files) {
      $node->runCmd("cp $file $serverSubTaskDir");

    }
    $node->runCmd("cp -r $serverSubTaskDir/* $nodeExecDir");

}


sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $muscleDir = $self->getProperty("muscleBinDir");

    my $cmd = "runMuscle --muscleDir $muscleDir --inputFileDir $nodeExecDir";

    return $cmd;

}

sub runMuscle {
  my ($self, $nodeExecDir,$muscleDir) = @_;

  opendir(DIR, $nodeExecDir) || die "Can't open directory 'nodeExecDir'";

  my @tarBall;

  while (defined (my $file = readdir (DIR))) {
    next if ($file eq "." || $file eq "..");
    system("$muscleDir/muscle -in $nodeExecDir/$file -out $nodeExecDir/${file}.msa -clw -log $nodeExecDir/${file}.log");

    push(@tarBall,$nodeExecDir/${file}.msa); 
  }

  my $list = join (' ',@tarBall);

  chdir $nodeExecDir;

#how can I get the subtask number?

  `tar -zcf $outDir/tarrBall_${tarNum}.tar.gz $list`;
}


# cleanUpNode is an optional method that is called when the node has completed
# to allow the user to stop processes that they may have started on the node.
# the files and directory structure on the nodes are already cleaned up.
sub cleanUpNode {
  my($self,$node) = @_;
}

# cleanUpServer is an optional method that is called by the controller after
# all nodes have completed.  This allows users to run some additional analysis
# on the server that may clean up or further analyze the combined  results of
# the run.
sub cleanUpServer {
  my($self, $inputDir, $mainResultDir) = @_;
  return 1;
}

1;
