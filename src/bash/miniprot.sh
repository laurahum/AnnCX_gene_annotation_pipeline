#!/bin/bash

#Author: lahumada
#Date: 2026

# Function: miniprot_run  
# Description: Runs miniprot on a list of genomes or fasta files against a query sequence

# This script will be used to:
# 1. Find query genes to annotate in input hardmasked FASTA ROI files
#   - query: protein
# 2. Find query genes to annotate in input soft or hardmasked FASTA whole genome files
#   - query: protein

# Usage: miniprot_run <input_target_fasta_dir> <input_query_fasta> <output_dir> <genome_list> <gene_to_annotate> <max_intron> <type_of_fasta> [threads]
#
# Parameters:
#   $1 (input_target_fasta_dir): Directory containing target FASTA files
#   $2 (input_query_fasta): Query FASTA file (protein)
#   $3 (output_dir): Directory for GFF output files
#   $4 (genome_list): Path to TXT file with genome names, one per line
#   $5 (gene_to_annotate): gene name (e.g. "NKG2")
#   $6 (max_intron): Sets maximum intron size in basepairs  (default = 7000)
#   $7 (type_of_fasta): Specifies the type of FASTA file as input. Possible values:
#       - "genome"   
#       - "ROI"
#   $8 (threads): Optional threads

# Example:
#   miniprot_run "/path/ROI_hardmasked" "/path/cDNA.fasta" "/path/output" "/path/genome_list.txt" "NKG2" 4000 "ROI" 8
#
# Notes:
# - Ensure that miniprot is installed and accessible in the current environment
# - miniprot is run with the following key options:
#   - -G: Maximum intron length
#   - -t: Number of threads (optional)
# - Outputs GFF format
#
# Output files:
#   ${genome}_${type_of_fasta}_${gene_to_annotate}_prot_miniprot.gff3



miniprot_run(){
    local input_target_fasta_dir="$1"
    local input_query_fasta="$2"
    local output_dir="$3"
    local genome_list="$4"
    local gene_to_annotate="$5"
    local max_intron="$6"
    local type_of_fasta="$7"
    local threads="${8:-false}"

    while read genome; do
        echo "$genome"
        
        # Dynamically find the matching file for $genome
        target_fasta=$(compgen -G "$input_target_fasta_dir/${genome}*")
        
        # Prepare miniprot command
        local cmd="miniprot -G $max_intron --gff \"$target_fasta\" \"$input_query_fasta\""
                
        # Add threads option if specified
        if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            cmd+=" -t $threads"
        fi
                          
        # Run miniprot → GFF output
        eval $cmd > "$output_dir/${genome}_${type_of_fasta}_${gene_to_annotate}_prot_miniprot.gff3"
        echo "Saved: ${genome}_${type_of_fasta}_${gene_to_annotate}_prot_miniprot.gff3"
        
    done < "$genome_list"
}

# Export the function so it can be called from outside the script
export -f miniprot_run

