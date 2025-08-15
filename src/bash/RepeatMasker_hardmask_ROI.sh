 #!/bin/bash

#Author: lahumada
#Date: 2022

# Function: repeatmasker_run
# Description: Hardmasks FASTA files using RepeatMasker and organizes the output
#
# This function runs RepeatMasker on all FASTA files in a specified input directory,
# using given parameters. It then organizes the output by moving masked files and
# annotation files to separate subdirectories in the output directory:
#   - masked_files: Contains all .masked files produced by RepeatMasker
#   - repeat_annotations: Contains all .out.gff files produced by RepeatMasker
#
# Usage: repeatmasker_run <input_ROI_raw> <output_ROI_hardmasked> <species_arg_repeatmasker> [threads]
#
# Parameters:
#   $1 (input_ROI_raw): Directory containing input FASTA files to be processed
#   $2 (output_ROI_hardmasked): Directory where RepeatMasker output will be stored
#   $3 (species_arg_repeatmasker): Species argument for RepeatMasker (e.g., "primates")
#   $4 (threads): Optional. Number of parallel jobs to run
#
# Example:
#   repeatmasker_run "/path/to/input_ROI" "/path/to/output_ROI_hardmasked" "primates" 4
#
# Notes: 
# - Ensure that RepeatMasker is installed and accessible in the current environment
# - RepeatMasker is run with options: -s -norna -nolow -gff -u


repeatmasker_run(){
    local input_ROI_raw="$1"
    local output_ROI_hardmasked="$2"
    local species_arg_repeatmasker="$3"
    local threads="${4:-false}"
    
    # Subdirectories to move the masked ROI and the annotation of repeats
    local roi_hardmasked_folder="$output_ROI_hardmasked/roi_hardmasked"
    local repeat_annotations_folder="$output_ROI_hardmasked/repeat_annotations"

    # Create directories if they don't exist
    mkdir -p "$roi_hardmasked_folder"
    mkdir -p "$repeat_annotations_folder"
    
    
    files=($input_ROI_raw/*.fasta)
    if [ ! -e "${files[0]}" ]; then
    	echo "No FASTA files found in $input_ROI_raw"
    	exit 1
    fi
    
    cd $output_ROI_hardmasked
    
    for file in "${files[@]}"; do
        echo "Processing $file"
        
        # Prepare the RepeatMasker command
        local cmd="RepeatMasker -s -norna -nolow -species \"$species_arg_repeatmasker\" -dir \"$output_ROI_hardmasked\" -gff -u \"$file\""
        
        # Add threads option if specified
        if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            cmd+=" -pa $threads"
        fi
        
        # Run RepeatMasker
        eval $cmd
    done
    
    # Move all .masked files to masked_folder
    for file in "$output_ROI_hardmasked"/*.masked; do
        mv "$file" "$roi_hardmasked_folder"
    done

    # Move all .out.gff files to annotations_folder
    for file in "$output_ROI_hardmasked"/*.out.gff; do
        mv "$file" "$repeat_annotations_folder"
    done
    
    echo "MASKED_FOLDER=$roi_hardmasked_folder"
    echo "ANNOTATIONS_FOLDER=$repeat_annotations_folder"
    
    cd
}

# Export the function so it can be called from outside the script
export -f repeatmasker_run
