package DJob::DistribJobTasks::GenomeAlignTask;

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
 ["targetListPath",      "",     "full path to file list target files"],
 ["queryPath",   "",     "full path to input file"],
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;
    my $gaBinPath = $self->{props}->getProp("gaBinPath");
    die "gaBinPath $gaBinPath doesn't exist" unless -e $gaBinPath;
}

sub initNode {
    my ($self, $node, $inputDir) = @_;

}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->{props}->getProp("queryPath");

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

  my $gaBinPath = $self->{props}->getProp("gaBinPath");
  my $targetListPath = $self->{props}->getProp("targetListPath");
  my $paramsPath = $inputDir . '/params.prop';

  my $cmd =  "blatSearch --blatBinPath $gaBinPath --targetListPath $targetListPath --seqPath $nodeExecDir/seqsubset.fsa --paramsPath $paramsPath";

  return $cmd;
}


sub integrateSubTaskResults {
  my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

  my @files = split /\n/, $self->runCmdOnNode($node, "ls -1 $nodeExecDir");

  # expect output file name to be like blat-chr5.psl
  my @seqs = grep(/\-.+\./i, @files);
  foreach my $seq (@seqs) {
    chomp $seq;
    my $outdir = "$mainResultDir/per-seq";
    my $outfile = "$outdir/$seq";
    $self->runCmdOnNode($node, "mkdir -p $outdir") unless -d $outdir;

    # workaround to odd bpsh problem:
    # cat >> $file fails if $file does not exist,
    # even if we "touch $file" first.
    if (-e $outfile) {
	$self->runCmdOnNode($node, "cat $nodeExecDir/$seq >> $outfile");
    } else {
	$self->runCmdOnNode($node, "cp $nodeExecDir/$seq $outfile");
    }

 }
}
1;
