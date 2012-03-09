#!/bin/sh
# The shell script to run the Proteomics Pipeline
# Set the PIPELINE_HOME in .profile file

echo "Starting the Proteomics Pipeline For A Single MGF.."
echo "Proteomics Pipeline Path = $PIPELINE_HOME"

PIPELINE_PATH="$PIPELINE_HOME"

INPUT_FILE_TEMPLATE="$1"
INPUT_MGF="$2"
OUTPUT_DIR="$3"
PARSER_INPUT_FILE="$4"

echo "$INPUT_MGF_DIR"

# Change to the Pipeline directory
cd "$PIPELINE_PATH"

echo "$PWD"

java -jar pipelineForCluster.jar $INPUT_FILE_TEMPLATE $INPUT_MGF $OUTPUT_DIR = $PARSER_INPUT_FILE  =

