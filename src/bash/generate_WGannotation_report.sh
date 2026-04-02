 #!/bin/bash

#Author: lahumada
#Date: 2026


# Function: generate_WGannotation_report
# Description: Generates a report of whole-genome minimap2 annotation results and cleans up empty files
# This function scans minimap2 output directories to identify species folders containing successful # alignments and generates a summary report of minimap2 whole-genome annotation results
#           - Lists species directories with >1 non-empty PAF files (successful alignments)
#           - Shows relative filenames of successful alignment files
#           - Cleans up empty output files from minimap2 runs

# Usage: generate_WGannotation_report <input_minimap2_dir> <output_report_dir>
#
# Parameters:
#   $1 (input_minimap2_dir): Directory containing species subdirectories with minimap2 PAF output files
#   $2 (output_report_dir): Directory where the report file will be saved (creates report_WGannotation.txt)

# Example:
#   generate_WGannotation_report "/path/to/minimap2/raw" "/path/to/output"
#
# Notes:
# - Uses POSIX sh for maximum compatibility (no bash arrays)
#
# Example report output:
#   Pan_paniscus (2 files)
#   CM003395.1_genome_NKG2_cDNA_minimap2.paf
#   CM003389.1_genome_NKG2_cDNA_minimap2.paf
#   ...

generate_WGannotation_report(){
    local input_minimap2_dir="$1"
    local output_report_dir="$2"

    cd "$input_minimap2_dir"
    
    # Find files with matches found by minimap2 (>0 bytes) and generate report
    find . -mindepth 1 -maxdepth 1 -type d -links 2 -exec sh -c '
  	dir="$1"
  	basename="${dir##*/}"
  	files=$(find "$dir" -maxdepth 1 -type f -size +0c 2>/dev/null)
  	count=$(echo "$files" | grep -c .)
  	if [ "$count" -gt 1 ]; then
    		echo "$basename ($count files)"
    		echo "$files" | sed "s|$dir/||g"
    		echo
  	fi
    ' sh {} \; > "$output_report_dir/report_WGannotation.txt"

    # Delete empty files = 0 bytes
    find . -type f -size 0c -delete
    
    cd -
}
    
# Export the function so it can be called from outside the script
export -f generate_WGannotation_report
