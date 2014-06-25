package DJob::DistribJobTasks::Iprscan5Task;

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
	["appl", "", "Comma-separated names of applications. Supported applications: ProDom coils "
					 . "PIRSF PfamA SMART TIGRFAM PrositeProfiles SuperFamily Gene3d "
					 . "scanregexp seg"],
	["email", $ENV{USER} . "\@pcbi.upenn.edu", "Submitter's email address to send job status"],
	["taxo", "false", "Activate the Taxonomy filter for abbreviated taxonomy"],
        ["copyDbToNode", "no", "(yes | [no]) if 'yes' then copies the blast indices to the local nodeDir on node ... may be faster in some contexts"],
	["datadir", "false", "Full path to the Interproscan data directory"]
	
	# Notes:
	#-------
	# Always set output format to 'TSV'
	#
	# Always use goterms and iprlookup
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
	
    return if $self->getProperty("copyDbToNode") eq 'no';

	# This determines what application uses which files. If adding a new application,
	# make sure this hash is populated appropriately. Look at the individual [app].conf
	# files in IPRSCAN_HOMe/conf directory to figure out the files used by each app.
	# Make sure the index files are copied over too.
	my %appl_files = ( prodom=>'prodom/2006.1/prodom.*',
						pirsf=>'pirsf/2.80/pirsf.dat;pirsf/2.80/sf_*;pirsf/2.80/sf.*',
						pfama=>'pfam/26.0/Pfam*',
						smart=>'smart/6.2/smart.*',
						tigrfam=>'tigrfam/10.1/TIGRFAMs_*',
						prints=>'prints/41.1/prints.*;prints/41.1/FingerPRINTSparser.db',
						prositeprofiles=>'prosite/20.80/prosite*',
						superfamily=>'superfamily/1.75/superfamily.*',
						tmhmm=>''
						);
	
	# Copying "data" files to the nodes avoids expensive network I/O
	# Create entire iprscan_home because that's what the individual applications expect
	# Copy over the data files, but the rest of the dir can be symlinked.
	my $nodeIprscanHome = $node->getDir() . "/iprscan_home";
	$node->runCmd ("mkdir $nodeIprscanHome");
	$node->runCmd ("ln -s $ENV{IPRSCAN_HOME}/bin $nodeIprscanHome/bin");
	$node->runCmd ("ln -s $ENV{IPRSCAN_HOME}/lib $nodeIprscanHome/lib");
	$node->runCmd ("mkdir $nodeIprscanHome/tmp");
	$node->runCmd ("mkdir $nodeIprscanHome/data");
	
	my $appls = $self->getProperty("appl");
	if ($appls) {
		#copy only the data files that are need for the given applications
		foreach my $appl ( split (/,\s*/, $appls)) {
			foreach my $appl_db_file (split (/;/, $appl_files{$appl})) {
			        my($filename, $directories, $suffix) = fileparse($appl_db_file);
			        $node->runCmd ("mkdir $nodeIprscanHome/data/$directories");
				$node->runCmd ("cp $ENV{IPRSCAN_HOME}/data/$appl_db_file $nodeIprscanHome/data/$directories");
			}
		}
	} else {
		# the default is to run all the configured iprscan applications 
		# (see IPRSCAN_HOME/conf/iprscan.conf)
		# Simpler to just copy over all the data files instead of parsing iprscan.conf and picking
		# only the configured applications.

		foreach my $appl (values %appl_files) {
			foreach my $appl_db_file (split (/;/, $appl_files{$appl})) {
                                my($filename, $directories, $suffix) = fileparse($appl_db_file);
			        $node->runCmd ("mkdir $nodeIprscanHome/data/$directories");
				$node->runCmd ("cp $ENV{IPRSCAN_HOME}/data/$appl_db_file $nodeIprscanHome/data/$directories");
			}
		}
	}
	
	# All the applications read their data from $IPRSCAN_HOME/data
	$node->runCmd ("export IPRSCAN_HOME=$nodeIprscanHome");
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
    $node->runCmd("touch $serverSubTaskDir/seqsubset.fsa.touch",1);
    $node->runCmd("/bin/rm $serverSubTaskDir/seqsubset.fsa.touch",1);

    $node->runCmd("cp -r $serverSubTaskDir/seqsubset.fsa $nodeExecDir");

	# print "Created subtask fasta file at: $nodeExecDir/seqsubset.fsa\n";
}

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;

	my @cmd_params;

	# Required parameters
	push @cmd_params, "-i $nodeExecDir/seqsubset.fsa";
#	push @cmd_params, "-email " . $self->getProperty ("email");
#	push @cmd_params, "-email dontcare\@dontcare.nowhere.com";
	push @cmd_params, "-o $nodeExecDir/" . $self->getProperty ("outputfile");

	# Task defaults (which override default iprscan configuration in case of a conflict)
	push @cmd_params, "-f TSV";
	push @cmd_params, "-iprlookup";
	push @cmd_params, "-goterms";
#	push @cmd_params, "-verbose";
    

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

	#my $crc = $self->getProperty ("crc");
	#(!$crc || $crc =~ /false/i) and push @cmd_params, "-nocrc";
	
	my $cmd = "$ENV{IPRSCAN_HOME}/interproscan.sh " . join (" ", @cmd_params);
	print "Returning command: $cmd\n";
	return $cmd;
}

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
    my $outputfile = $self->getProperty ("outputfile");
    my $output_part_file = $outputfile;
    $output_part_file =~ s/\.tsv$//i;
    $output_part_file .= "_" . $subTaskNum . ".tsv";
    $node->runCmd ("cp $nodeExecDir/$outputfile $mainResultDir/$output_part_file");
    return 1 if $node->getErr();
}

sub cleanUpServer {
  my($self, $inputDir, $mainResultDir, $node) = @_;
}

1;
