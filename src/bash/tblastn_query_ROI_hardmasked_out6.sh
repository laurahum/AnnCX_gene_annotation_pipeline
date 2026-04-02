 #!/bin/bash

#Author: lahumada
#Date: 2022


# Function: tblastn_run
# Description: Performs TBLASTN searches against hardmasked ROI makeblastdb databases
#
# This function runs TBLASTN for each genome listed in a single contig list file.
# It uses pre-built BLAST databases of hardmasked ROI sequences, generated using
# makeblastdb, and a specific query file containing protein sequences.
#
# Usage: blastn_run <makeblastdb_dir> <query_dir> <output_dir> <single_contig_list> [threads] <gene_to_annotate>
#
# Parameters:
#   $1 (makeblastdb_dir): Directory containing the pre-built BLAST databases
#   $2 (query_dir): Directory containing the query FASTA file (cDNA)
#   $3 (output_dir): Directory where BLAST output files will be stored
#   $4 (single_contig_list): TXT file containing a list of single-contig genomes
#   $5 (gene_to_annotate): Name of the genes of interest to be annotated, used to label output files (Example: "NKG2")  
#   $6 (threads): Optional. Number of CPU threads to use for tblastn
#
# Example:
#   tblastn_run "/path/to/blastdb" "/path/to/protein_queries" "/path/to/tblastn_results" "/path/to/single_contig_list.txt" "NKG2" 4
#
# Notes: 
# - Ensure that TBLASTN is installed and accessible in the current environment
# - TBLASTN output format is set to outfmt 6 with specific fields:
# "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs"


tblastn_run(){
    local makeblastdb_dir="$1"
    local query_dir="$2"
    local output_dir="$3"
    local single_contig_list="$4"
    local gene_to_annotate="$5"
    local threads="${6:-false}"

    while read genome; do
        echo "$genome"
    
        # Prepare the tblastn command
        local cmd="tblastn -db \"$makeblastdb_dir/${genome}\" -query \"$query_dir\" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs\""
        
        # Add threads option if specified
        if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            cmd+=" -num_threads $threads"
        fi
        
        # Run tblastn
        eval $cmd > "$output_dir/${genome}_ROI_${gene_to_annotate}_protein_tblastn_out6"
    
    done < "$single_contig_list"
}

# Export the function so it can be called from outside the script
export -f tblastn_run
