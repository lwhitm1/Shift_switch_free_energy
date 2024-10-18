#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <parent_directory> <xvg_file> <log_file_prefix>"
    exit 1
fi

# Assign command-line arguments to variables
PARENT_DIR="$1"
XVG_FILE="$2"
LOG_FILE_PREFIX="$3"

# Set the path to the semi_analysis.py script
SCRIPT_PATH="/gpfs/alpine1/scratch/liwh2139/more_lambdas/new_cutoffs/85_9/semi_analysis.py"

# Loop through each subdirectory in the parent directory
for SUBDIR in "$PARENT_DIR"/*; do
    if [ -d "$SUBDIR" ]; then
        echo "Running analysis script in: $SUBDIR"
        cd "$SUBDIR" || continue
        # Run the analysis command
        python "$SCRIPT_PATH" --temp 300 --xvg_file "$XVG_FILE" > "${LOG_FILE_PREFIX}_$(basename "$SUBDIR")_analysis.log" 2>&1
        cd - > /dev/null
    fi
done

