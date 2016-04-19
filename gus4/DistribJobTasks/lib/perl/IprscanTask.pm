package DJob::DistribJobTasks::IprscanTask;

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
	["seqfile", "", "Input sequence file (Protein or Nucleotide)"],
	["outputfile", "", "Output result file name"],
	["seqtype", "p", "Sequence type: n for DNA/RNA, p for protein"],
	["trlen", "false", "Transcript length threshold"],
	["appl", "", "Comma-separated names of applications. Supported applications: blastprodom coils "
					 . "fprintscan gene3d hmmpanther hmmpfam hmmpir hmmsmart hmmtigr profilescan "
					 . "scanregexp seg superfamily"],
	["crc", "false", "Checks interpro compatible databases for data consistency. Set this to false "
					. " (default) to allow updating individual databases without updating interproscan"],
	["email", $ENV{USER} . "\@pcbi.upenn.edu", "Submitter's email address to send job status"],
	["taxo", "false", "Activate the Taxonomy filter for abbreviated taxonomy"],
	["datadir", "false", "Full path to the Interproscan data directory"]
	
	# Notes:
	#-------
	# Always set output format to 'xml', because this Task handles only xml
	#
	# Always use goterms and iprlookup, because we want simple and consistent xml parsing when
	# combining the results.
	#
);

sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties);
    return $self;
}

# called once 
sub initServer {
    my ($self, $inputDir) = @_;

}

sub initNode {
    my ($self, $node, $inputDir) = @_;
	
}

sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $fastaFileName = $self->getProperty("seqfile");

    if (-e "$fastaFileName.gz") {
	&runCmd("gunzip $fastaFileName.gz");
    }

    print "Counting sequences in $fastaFileName\n";
    $self->{fastaFile} = CBIL::Bio::FastaFileSequential->new($fastaFileName);
    return $self->{fastaFile}->getCount();
}

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;

    if(!$subTask->getRedoSubtask()){
      $self->{fastaFile}->writeSeqsToFile($start, $end, "$serverSubTaskDir/seqsubset.fsa");
    }
    $self->runCmdOnNode("touch $serverSubTaskDir/seqsubset.fsa.touch",1);
    $self->runCmdOnNode("/bin/rm $serverSubTaskDir/seqsubset.fsa.touch",1);

    $self->runCmdOnNode("cp -r $serverSubTaskDir/seqsubset.fsa $nodeExecDir");

	# print "Created subtask fasta file at: $nodeExecDir/seqsubset.fsa\n";
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

	my @cmd_params;

	# Required parameters
	push @cmd_params, "-i $nodeExecDir/seqsubset.fsa";
#	push @cmd_params, "-email " . $self->getProperty ("email");
	push @cmd_params, "-email dontcare\@dontcare.nowhere.com";
	push @cmd_params, "-o $nodeExecDir/" . $self->getProperty ("outputfile");

	# Task defaults (which override default iprscan configuration in case of a conflict)
	push @cmd_params, "-format xml";
	push @cmd_params, "-iprlookup";
	push @cmd_params, "-goterms";
	push @cmd_params, "-verbose";
    

	# Optional parameters (Use these values only if specified; otherwise use iprscan configuration)
	my $appl = $self->getProperty ("appl");
	if ($appl) {
		my $appl_param = "";
		foreach my $app (split (/,\s*/, $appl)) {
			$appl_param .= "-appl $app ";
		}
		push @cmd_params, $appl_param;
	}
	
	my $seqtype = $self->getProperty ("seqtype");
	$seqtype and push @cmd_params, "-seqtype $seqtype";
	
	my $trlen = $self->getProperty ("trlen");
	($trlen && $trlen !~ /false/i) and push @cmd_params, "-trlen $trlen";

	my $crc = $self->getProperty ("crc");
	(!$crc || $crc =~ /false/i) and push @cmd_params, "-nocrc";
	
	my $cmd = "$ENV{IPRSCAN_HOME}/bin/iprscan -cli " . join (" ", @cmd_params);
	# print "Returning command: $cmd\n";
	return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    my $outputfile = $self->getProperty ("outputfile");
    my $output_part_file = $outputfile;
    $output_part_file =~ s/\.xml$//i;
    $output_part_file .= "_" . $subTaskNum . ".xml";
    $self->runCmdOnNode ("cp $nodeExecDir/$outputfile $mainResultDir/$output_part_file");
}

sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;
}

1;
