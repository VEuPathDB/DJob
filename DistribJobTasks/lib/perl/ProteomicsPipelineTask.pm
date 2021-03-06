package DJob::DistribJobTasks::ProteomicsPipelineTask; 

use CBIL::Util::Utils;
use Data::Dumper;
use File::Basename;

use DJob::DistribJob::Task; 
@ISA = (DJob::DistribJob::Task); 

use strict; 

# Declare the properties that will be passed to your task in the task.prop
# properties file.  This file is the only way for the user to pass parameters 
# to your task.  
#
# The task.prop file contains one property per line, where a property is of 
# the following form (lines beginning with # are comments):
#   name=value
#
# The user will place the file in the input directory (which is specified by 
# the inputdir property of the controller property file). The controller 
# parses the file, and validates the properties based on your declaration.
# It throws an error if the file includes invalid property names, or 
# fails to include required properties.
#
# Use @properties to declare the properties your task needs.  Each element in 
# @properties is an array containing:
#    [name_of_property, default_value, description_of_property]
# If the default value is "" then the property is required.  
#
# The properties are made available to your task in a CBIL::Util::PropertySet 
# object in $self->{props}.  Access a property by calling the getProperty() method:
#   $self->getProperty('name_of_property');
#
my @properties =
    (
     ["configDir", "$ENV{GUS_HOME}/config/ProteoAnnotator", "directory containing config files for proteoAnnotator"],
     ["mgfFile", "", "directory containing the mgf file"],
     ["searchConfigFile", "/inputFiles/searchCriteria.txt", "name of the search criteria file to be used" ],
     ["databaseConfigFile", "/inputFiles/databaseCriteria.txt", "name of the database criteria file to be used" ],
     ["mgfSplitterSizeLimit", "1000000", "number of lines to use for each split mgf file"],
     );
# Construct a new instance of your task.  This method should be copied 
# verbatim into your task.
sub new {
    my $self = &DJob::DistribJob::Task::new(@_, \@properties,1);

    ##would be good here to check to see that all files exist .. die if not ...
    die("configDir does not exist\n:$!") unless (-d $self->getProperty("configDir"));
    my $configDir =  $self->getProperty("configDir");
    die("mgfFile does not exist\n") unless (-e $self->getProperty("mgfFile"));
    my $searchConfigFile =$configDir.$self->getProperty("searchConfigFile");
    my $databaseConfigFile =$configDir.$self->getProperty("databaseConfigFile");
    die("searchConfigFileName does not exist in the folder \n :$!") unless (-e $searchConfigFile);
    die("databaseConfigFileName does not exist in the folder \n :$!") unless (-e $databaseConfigFile); 
    return $self;
}

# Initialize the server.  This method is called once by the controller, before
# any processing on the nodes begins.  Use this method to do any
# pre-processing (such as indexing) on data required by your task.  Your task
# may then access this data directly on the server, or, possibly more
# efficiently, on the nodes' local disks if you copy this data to the nodes 
# during node initialization.
# 
# There are only two correct ways for you to locate resources that need
# initializing.  You may have declared a property that contains the location
# (eg, matrixfile) or, you may require the user to place the resource in the
# input directory.
#
# param inputDir The input directory where you may expect the user to place resources.
#

sub initServer {
    my ($self, $inputDir) = @_;
    ##here we should split up the mgf file into as many subparts 
#    unless(-d "$inputDir/subtasks") {
#	mkdir("$inputDir/subtasks") || die "unable to make dir $inputDir/subtasks" ;    
#    }
#check if Successful (make a status file that can be checked)
    ##split file
    my $mgfFile = $self->getProperty("mgfFile");
    my $limit = $self->getProperty("mgfSplitterSizeLimit");
    &runCmd("perl gus_home/bin/mgfSplitter.pl --inputMgfFile $mgfFile --outputDir $inputDir/subtasks --limit $limit") unless -e "$inputDir/subtasks/success.txt"; 
    my @files = glob("$inputDir/subtasks/*.mgf");
    $self->{"files"} = \@files;
#    print STDERR "My File List is ".join(',',@files)."\n";
    $self->{"inputSetSize"} = scalar(@files);
    $self->{size}=$self->getInputSetSize();
#    print STDERR $self->getInputSetSize();
    if(system("echo $self->{'inputSetSize'} >$inputDir/subtasks/success.txt") !=0)
      {
        die "unable to create file $inputDir/subtasks/success.txt";
      }
}
 
# Initialize the local disk on a node.  This method is called once per node
# by the controller, before any processing on the nodes begins.  
# Use this method to copy data to the node from the server, and to initialize
# that data as needed.  If the user runs more that one subtask per node,
# each of those subtasks will share this data.  It must not be written to
# by the subtasks.
#
# There are only two correct ways for you to locate resources to copy onto
# the node.  You may have declared a property that contains the location
# (eg, matrixfile) or, you may require the user to place the resource in the
# input directory.
#
# There is only one valid location on the node where you can place data:  the
# "node directory" provided by $node->getDir()
#
# The only way to access the node directory is by executing commands on the
# node itself.  To do this, use $self->runCmdOnNode($node, "my command").  
#
# param node The node object.
# param inputDir The input directory where you may expect the user to place resources.
#

sub initNode {
    my ($self, $node, $inputDir) = @_;
    #do nothing
}
 
# Provide the size of the input set.  The controller needs to know this so
# it knows when the task is complete.
#
# param inputDir The input directory where you may find the input file.
#

sub getInputSetSize {
    my ($self, $inputDir) = @_;
    return $self->{"inputSetSize"};
}

# Initialize subtask-specific directories on the server and/or the node for use
# when the subtask runs.  This method is called by the controller once for each
# subtask just before the subtask is run.
# 
# There are two directories available to the subtask to find its input.
# $serverSubTaskDir is on the server.  $nodeSubTaskDir is on the node.  Place 
# in one or both of these directories the input needed by the subtask. Use 
# $nodeSubTaskDir if having the input locally on the node will speed the 
# application.
#
# If you use $nodeSubTaskDir, it may be useful to first put the input 
# into $serverSubTaskDir, and then copy it to $nodeSubTaskDir.  This way you 
# have a mirror of $nodeSubTaskDir on the server. This makes it easy to
# peruse the subtask's input on the server without needing to go to the node.
# (The server subtask dirs are found in the master/running/subtask_xxx dir)
#
# The controller creates the directory $serverSubTaskDir/result, which is
# a convenient place to run from and to store results. (If you copy the 
# contents of $serverSubTaskDir to $nodeSubTaskDir this directory will be
# included and you can use it as your nodeRunDir.  See the runSubTask method 
# below).
#
# $start and $end are indices into the input set. They are  provided by the 
# controller to specify the subset of the input used by the current subtask.  
# Use these to access that subset from the input file found in $inputDir 
# and write the subset to a file in $serverSubTaskDir and/or $nodeSubTaskDir. 
# This file will be the input to the subtask.
#
# param start The index of the first element in the subset (starting at 0).
# param end The index of the first element beyond the subset (ending at sz).
# param node The Node object on which the subtask will run.
# param inputDir The input dir on the server where the task's input is found.
# param serverSubTaskDir The subtask specific input dir on the server.
# param nodeSubTaskDir The subtask specific input dir on the node.
# 

sub initSubTask {
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, $nodeExecDir,$subTask) = @_;
     ##check here to make certain that subtask size = 1
    die "You must set the subtask size to 1\n" unless $start == $end-1;
    
    my $mgfFile = $self->{"files"}->[$start];
    my $configFileDir =  $self->getProperty("configDir");
    $self->runCmdOnNode($node, "mkdir $nodeExecDir/mgfFiles/");
    my $cmd = "cp $mgfFile $nodeExecDir/mgfFiles/subtaskMgf.mgf";
#    print STDERR "initSubTask : my command is :".$cmd."\n";
    $self->runCmdOnNode($node, "$cmd");
    my $cmdLine = "find $configFileDir/ -maxdepth 1 -exec ln -s {} $nodeExecDir \\;";
    $self->runCmdOnNode($node, "$cmdLine");
} 

# Actually run the subtask by issuing a command on a node. This method
# is called once by the controller for each subtask, after the subtask
# is initialized.
#
# The command you run must:
#   - get its input from $inputDir, $serverSubTaskDir or $nodeSubTaskDir
#   - place all results, logs and temporary files in the directory that it is
#     run from.
#   - return non-zero exit status on failure.
#
# To run the subtask, you just need to return the command in the makeSubTaskCommand method.  
#
# NOTE:  All STDOUT will be redirected into a file called subtask.output
#        All STDERR will be redirected into a file called subtask.stderr
#        Users can NOT redirect either STDERR or STDOUT to files of their choice
#
# param $node The Node object on which the subtask will run.
# param $inputDir The input dir on the server where the task's input is found.
# param nodeExecDir The subtask specific dir on the node where the command will be run.
# 

sub makeSubTaskCommand { 
    my ($self, $node, $inputDir, $nodeExecDir) = @_;    
    my $configFileDir =  $self->getProperty("configDir");
    my $searchConfigFile =$configFileDir.$self->getProperty("searchConfigFile");
    my $databaseConfigFile =$configFileDir.$self->getProperty("databaseConfigFile");
    my $mgfFile = "$nodeExecDir/mgfFiles/subtaskMgf.mgf";
    my $fileSuccess = $self->runCmdOnNode($node, "ls -la $nodeExecDir/mgfFiles/subtaskMgf.mgf");
    $self->runCmdOnNode($node, "source /gpfs/fs121/h/jcade/setenv");
    if($fileSuccess)
    {
#	print STDERR "my ls :". $fileSuccess."\n";
    }
    my $jarFile =  "$ENV{GUS_HOME}/lib/java/proteoannotator.jar";
   
    $self->runCmdOnNode($node, "rm $nodeExecDir/output -r"); 
    
    my $cmd = "java -jar $jarFile single_mode -searchInput $searchConfigFile -databaseInput $databaseConfigFile -inputMgf $mgfFile -outputResultDir $nodeExecDir/output/";

    return $cmd; 
} 

# Merge subTask results from $nodeExecDir into the main result in
# $mainResultDir as needed.  This method is called once by the controller
# for each subtask after it completes successfully.
#
# Each subtask produces results in $nodeExecDir.   These
# results must be placed or merged into the main result, which is stored
# in $mainResultDir.
#more 
# param subTaskNum The index of this subtask.  Use this if you want to store
#                  subtask results by subtask number.
# param node the node obect that this subtask was run on
# param nodeExecDir subtask specific dir on node where command was executed...result files here
# param mainResultDir The directory in which the main result is stored.
#

sub integrateSubTaskResults {
    my ($self, $subTaskNum, $node, $nodeExecDir, $mainResultDir) = @_;
     ##could check to make sure that the number of files is correct, otherwise throw an error.
    ##need to go back and look at how to get this to copy into failures ... simply return non-zero in this method
    my @files=$self->runCmdOnNode($node, "ls $nodeExecDir/output/dirsubtaskMgf/");
   print STDERR Dumper @files;
    print STDERR "My subtask number = $subTaskNum\n";
    foreach my $file (@files){	
	chomp $file;
	my $cmdLine = undef;
	my $uniqueFileName = undef;
	print STDERR Dumper (\@files);
	if ( $file=~m/^Summary(_\d+)?\.txt/ || $file=~m/^FinalOutput(_\d+)?\.txt/ || $file=~m/^FinalOutput_Verbose(_\d+)?\.txt/)
	{
	    $uniqueFileName=$file;
	    print STDERR "Orginal File Name : $uniqueFileName \n\n ";
	    $uniqueFileName=~s/_\d+/_$subTaskNum/;
	    print STDERR "New File Name : $uniqueFileName \n\n ";
	    $cmdLine = "cp $nodeExecDir/output/dirsubtaskMgf/$file $mainResultDir/$uniqueFileName";
	}
	$self->runCmdOnNode($node, "$cmdLine") if $cmdLine;
	
    }

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
  my($self, $inputDir, $mainResultDir,$node) = @_;
  my $jarFile =  "$ENV{GUS_HOME}/lib/java/proteoannotator.jar";
  my $mainFiles = $self->runCmdOnNode($node, "ls $mainResultDir");
  print STDERR $mainFiles;
  $self->runCmdOnNode($node, "java -jar $jarFile create_summary -resultDir $mainResultDir -summaryFile $mainResultDir/WholeSummary.txt");
  
  #generate summary
  return 1;
} 1;
