#!/bin/bash

#Author: lahumada
#Date: 2022

# Function: blastn_run_identify_gene
# Description: This post-processing feature uses BLASTN to search the sequences of reference vs 
# predicted genes to identify genes found by the genomic annotation process.
#
# Usage: blastn_run_identify_gene <subject_file> <query_file> <output_dir> <type_fasta_seq_query> <type_fasta_seq_subject> <name_genome> [gapopen] [gapextend]
#
# Parameters:
#   $1 (subject_file): File with FASTA sequences of reference
#   $2 (query_file): File with FASTA sequences of predictions
#   $3 (output_dir): Directory where BLAST output files will be stored
#   $4 (type_fasta_seq_query): Type of FASTA query sequence being used (e.g., "cDNA", "exon_all")
#   $5 (type_fasta_seq_subject): Type of FASTA subject sequence being used (e.g., "cDNA", "exon_all")
#   $6 (name_genome): Name of the genome in which the genes were annotated (e.g., "Macaca_mulatta")
#   $7 (gapopen): Int to run the -gapopen open_penalty argument in BLASTN (e.g., 1)
#   $8 (gapextend): Int to run the -gapextend extend_penalty argument in BLASTN (e.g., 1)
#
# Example:
#   blastn_run_identify_gene "/path/to/reference_sequences.fasta" "/path/to/predicted_sequences.fasta" "/path/to/blastn_results" "cDNA" "cDNA" "Macaca_mulatta" 1 1
#
# Notes: 
# - Ensure that BLASTN is installed and accessible in the current environment.
# - BLASTN output format is set to outfmt 6 with specific fields:
#   "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
# - The function uses a gap opening penalty of 1 and a gap extension penalty of 1.



blastn_run_identify_gene(){
    local subject_file="$1"
    local query_file="$2"
    local output_dir="$3"
    local type_fasta_seq_query="$4"
    local type_fasta_seq_subject="$5"
    local name_genome="$6"
    local gapopen="$7"
    local gapextend="$8"
    
    # Prepare the blastn command
    local cmd="blastn -subject \"$subject_file\" -query \"$query_file\" -gapopen \"$gapopen\" -gapextend \"$gapextend\" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\""
    
    # Make output file path
    local outfile="$output_dir/Identify_pred2ref_${name_genome}_${type_fasta_seq_query}vs${type_fasta_seq_subject}"
    
    # Run blastn
    eval $cmd > "$outfile"
    
    echo "$outfile"    
}

# Export the function so it can be called from outside the script
export -f blastn_run_identify_gene
