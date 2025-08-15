 #!/bin/bash

#Author: lahumada
#Date: 2022


# Function: gmap_build_run
# Description: Processes a list of genomes or fasta files using gmap_build
#
# This function reads a list with the genome names from a file and runs gmap_build on each input FASTA file. 
# gmap_build generates a genomic map database for each input FASTA file.
#
# This script will be used in two steps of the pipeline:
#   1. Find flanking genes in input FASTA genome files
#   2. Find query genes to annotate in input hardmasked FASTA ROI files
#
# Usage: process_genomes <genome_list> <gmap_build_dir> <input_genome_fasta> <type_of_fasta> [threads]
#
# Parameters:
#   $1 (genome_list): Path to a TXT file containing a list of genome names, one per line
#   $2 (gmap_build_dir): Directory where gmap_build will store its output
#   $3 (input_fasta): Directory containing the input FASTA files
#   $4 (type_of_fasta): Specifies the type of FASTA file as input. Possible values:
#       - "genome"   
#       - "ROI"
#   $5 (threads): Optional. Number of threads to use for gmap_build. If not specified or set to "false", 
#                 gmap_build will use its default thread settings.
#
# Examples:
#   gmap_build_run "/path/to/genome_list.txt" "/path/to/gmap_build_output" "/path/to/fasta_files" "genome"
#   gmap_build_run "/path/to/genome_list.txt" "/path/to/gmap_build_output" "/path/to/fasta_files" "genome" 4
#
# Note: Ensure that gmap_build (gmap) is installed and accessible in the current environment


gmap_build_run() {
    local genome_list="$1"
    local gmap_build_dir="$2"
    local input_fasta="$3"
    local type_of_fasta="$4"
    local threads="${5:-false}"

    while read genome; do
        echo "$genome"

        # Dynamically find the matching file for $genome
        fasta_file=$(compgen -G "$input_fasta/${genome}*")
        
        # Prepare the gmap_build command
        local cmd="gmap_build -d ${genome} -D $gmap_build_dir $fasta_file"
        
        # Add threads option if specified
        if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            cmd+=" -t $threads"
        fi
        
        # Run gmap_build
        eval $cmd
        
    done < "$genome_list"
}

# Export the function so it can be called from outside the script
export -f gmap_build_run
