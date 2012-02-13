# Task for TMHMM
# Center for Biological Sequence Analysis
# http://www.cbs.dtu.dk/services/TMHMM/TMHMM2.0b.guide.php

package DJob::DistribJobTasks::TmhmmTask;

use DJob::DistribJob::Task;
use CBIL::Bio::FastaFileSequential;
use File::Basename;
use CBIL::Util::Utils;

@ISA = (DJob::DistribJob::Task);

use strict;

# [name, default (or null if reqd), comment]
my @properties = 
(
 ["tmhmmProgram",  "",  "eg: /usr/local/TMHMM/bin/tmhmm"],
 ["inputFilePath", "",  "full path to input file"],
 ["tmhmmOptions",  "",  ""],
 );

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

sub initServer {
    my ($self, $inputDir) = @_;

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

    $node->runCmd("cp -r $subTaskDir/* $nodeSlotDir");
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

    my $tmhmmProgram = $self->getProperty("tmhmmProgram");
    my $tmhmmOptions = $self->getProperty("tmhmmOptions");
    my $inputFilePath = $self->getProperty("inputFilePath");
    my $cmd = "runTMHMM_distribjob -binPath $tmhmmProgram -workdir=$nodeExecDir $tmhmmOptions -seqFile $nodeExecDir/seqsubset.fsa -outFile tmhmm.out";

    return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;

    $node->runCmd("cat $nodeExecDir/tmhmm.out >> $mainResultDir/tmhmm.out");
    $node->runCmd("cat $nodeExecDir/subtask.output >> $mainResultDir/task.log");
}
1;
