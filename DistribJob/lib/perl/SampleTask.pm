package DistribJob::SampleTask;

use Common::Utils;
use DistribJob::Task;

@ISA = (DistribJob::Task);

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
# The properties are made available to your task in a Common::PropertySet 
# object in $self{props}.  Access a property by calling the getProp() method:
#   $self{props}->getProp('name_of_property');
#
my @properties = 
    (
     ["msg", "", "a message to print to the output"]
     );

# Construct a new instance of your task.  This method should be copied 
# verbatim into your task.
sub new {
    my $self = &DistribJob::Task::new(@_, \@properties);
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
# node itself.  To do this, use $node->runCmd("my command").  For example,
# to copy a file from the inputdir to the node dir:
#    my $nodeDir = $node->getDir();
#    $node->runCmd("cp $inputdir/myfile $nodeDir");
#
# param node The node object.
# param inputDir The input directory where you may expect the user to place resources.
#
sub initNode {
    my ($self, $node, $inputDir) = @_;

}

# Provide the size of the input set.  The controller needs to know this so
# it knows when the task is complete.
#
# param inputDir The input directory where you may find the input file.
#
sub getInputSetSize {
    my ($self, $inputDir) = @_;

    my $count;
    open(F, "$inputDir/inputset");
    while (<F>) {
	$count++ if (/\w+/)
    }
    close(F);
    return $count;
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
    my ($self, $start, $end, $node, $inputDir, $serverSubTaskDir, 
	$nodeSubTaskDir) = @_;

    &runCmd("cp $inputDir/task.prop $serverSubTaskDir");

    my $index = 0;
    open(IN, "$inputDir/inputset");
    open(OUT, ">$serverSubTaskDir/inputset");
    while (<IN>) {
	if ($index >= $start && $index < $end) {
	    print OUT;
	}
	$index++;
    }
    close(IN);
    close(OUT);

    $node->runCmd("cp -r $serverSubTaskDir/* $nodeSubTaskDir");
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
# To run the subtask, use the $node->execSubTask() method.  It has three
# arguments ($nodeRunDir, $serverResultDir, $cmd): 
#   
#  1) $nodeRunDir: the directory on the node from which to run the command. 
#     The controller changes into that dir before running the command.
#  2) $serverResultDir: the directory on the server to copy the contents of
#     $nodeRunDir when the command is complete.
#  3) $cmd: the command to run
#
# After the command completes successfully, the controller copies the
# contents of $nodeRunDir to $serverResultDir.
#
# param node The Node object on which the subtask will run.
# param inputDir The input dir on the server where the task's input is found.
# param serverSubTaskDir The subtask specific input dir on the server.
# param nodeSubTaskDir The subtask specific input dir on the node.
# 
sub runSubTask { 
    my ($self, $node, $inputDir, $serverSubTaskDir, $nodeSubTaskDir) = @_;

    my $msg = $self->{props}->getProp("msg");

    my $cmd = "samplecmd $nodeSubTaskDir/inputset $msg";

    $node->execSubTask("$nodeSubTaskDir/result", "$serverSubTaskDir/result", $cmd);
}

# Merge subTask results from $subTaskResultDir into the main result in
# $mainResultDir as needed.  This method is called once by the controller
# for each subtask after it completes successfully.
#
# Each subtask produces results in $subTaskResultDir (which was the second
# argument to $node->execSubTask() called above in runSubTask()).  These
# results must be placed or merged into the main result, which is stored
# in $mainResultDir.
#
# The controller deletes $subTaskResultDir after this method is called.
#
# param subTaskNum The index of this subtask.  Use this if you want to store
#                  subtask results by subtask number.
# param subTaskResultDir The directory in which all files produced by the
#                        subtask's command are available.
# param mainResultDir The directory in which the main result is stored.
#
sub integrateSubTaskResults {
    my ($self, $subTaskNum, $subTaskResultDir, $mainResultDir) = @_;

    &runCmd("cat $subTaskResultDir/answer >> $mainResultDir/answer");
}

1;
