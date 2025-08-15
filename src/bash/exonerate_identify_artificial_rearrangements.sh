#!/bin/bash

#Author: lahumada
#Date: 2022

# Function: exonerate_run_identify_rearrangements
# Description: Uses exonerate to identify potential artificial rearrangements by comparing 
# exon-by-exon between each predicted gene and the exons of the reference genes
#
# Usage: exonerate_run_identify_rearrangements <subject_file> <query_file> <output_dir> <name_genome> [threads]
#
# Parameters:
#   $1 (subject_file): File with FASTA exon sequences of reference
#   $2 (query_file): File with FASTA exon sequences of predictions
#   $3 (output_dir): Directory where exonerate output files will be stored
#   $4 (name_genome): Name of the genome in which the genes were annotated (e.g. "Macaca mulatta")
#   $5 (threads): Optional. Number of CPU threads to use for exonerate
#
# Example:
#   exonerate_run_identify_rearrangements "/path/to/reference_sequences.fasta" "/path/to/predicted_sequences.fasta" "/path/to/exonerate_results" "Macaca_mulatta" 4
#
# Notes:
# - Ensure that exonerate is installed and accessible in the current environment.
# - The exonerate output is formatted using --ryo to produce tab-separated values: query ID, target ID, and percentage identity.


exonerate_run_identify_rearrangements(){
	local subject_file="$1"
	local query_file="$2"
	local output_dir="$3"
    	local name_genome="$4"
    	local threads="${5:-false}"
	
	# Prepare the exonerate command
 	local cmd="exonerate --model affine:local \"$query_file\" \"$subject_file\" --showvulgar no --showalignment no --showcigar no  --showtargetgff no --ryo '%qi\t%ti\t%pi\n'"

	# Add threads option if specified
        if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
	    cmd+=" -c $threads"
	fi
	
	# Make output file path
	local outfile="$output_dir/Exonerate_identify_artificial_rearrangements_${name_genome}"
	
	# Run exonerate and post-process output
	eval "$cmd" | sed '/^Command line:/d; /^Hostname:/d; /^-- completed exonerate analysis/d' > "$outfile"
	
	echo "$outfile" 
}

# Export the function so it can be called from outside the script
export -f exonerate_run_identify_rearrangements
