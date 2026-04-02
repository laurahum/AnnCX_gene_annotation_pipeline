#!/bin/bash

#Author: lahumada
#Date: 2026

# Function: miniprot_index
# Description: Makes miniprot index on a list of genomes or fasta files

# This script will be used to:
# Make index for input soft or hardmasked FASTA whole genome files

# Usage: miniprot_index <input_target_fasta_dir> <output_dir> <genome_list> <gene_to_annotate> <type_of_fasta> [threads]
#
# Parameters:
#   $1 (input_target_fasta_dir): Directory containing target FASTA files
#   $2 (output_dir): Directory for index MPI output files
#   $3 (genome_list): Path to TXT file with genome names, one per line
#   $4 (gene_to_annotate): gene name (e.g. "NKG2")
#   $5 (type_of_fasta): Specifies the type of FASTA file as input. Possible values:
#       - "genome"   
#   $6 (threads): Optional threads

# Example:
#   miniprot_index "/path/ROI_hardmasked" "/path/output" "/path/genome_list.txt" "NKG2" "genome" 8
#
# Notes:
# - Ensure that miniprot is installed and accessible in the current environment
# - miniprot is run with the following key options:
#  	 -t: Number of threads (optional)
# - Outputs MPI format
#
# Output files:
#   ${genome}_${type_of_fasta}_${gene_to_annotate}_prot_miniprot.gff3


miniprot_index(){
    local input_target_fasta_dir="$1"
    local output_dir="$2"
    local genome_list="$3"
    local gene_to_annotate="$4"
    local type_of_fasta="$5"
    local threads="${6:-false}"

    cd "$output_dir"
    while read genome; do
        echo "$genome"
        
        # Dynamically find the matching file for $genome
        target_fasta=$(compgen -G "$input_target_fasta_dir/${genome}*")
        
        # Prepare miniprot command
        local cmd="miniprot -d ${genome}_${type_of_fasta}_${gene_to_annotate}_miniprot.mpi \"$target_fasta\""
                
        # Add threads option if specified
        if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            cmd+=" -t $threads"
        fi
                          
        # Run miniprot → MPI output
        eval $cmd
        echo "Saved: ${genome}_${type_of_fasta}_${gene_to_annotate}_miniprot.mpi"
        
    done < "$genome_list"
    cd
}

# Export the function so it can be called from outside the script
export -f miniprot_index
