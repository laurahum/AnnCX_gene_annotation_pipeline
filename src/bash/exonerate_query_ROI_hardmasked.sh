 #!/bin/bash

#Author: lahumada
#Date: 2022


# Function: exonerate_run
# Description: Performs gene annotation using Exonerate on hardmasked ROI sequences
#
# This function runs Exonerate for each genome listed in a single contig list file.
# It uses transcript sequences as queries to annotate genes in the corresponding 
# hardmasked ROI sequences. The results are output in GFF format.
#
# Usage: exonerate_run <ROI_hardmasked_dir> <query_dir> <output_dir> <single_contig_list> <gene_to_annotate> [max_intron] [threads]
#
# Parameters:
#   $1 (ROI_hardmasked_dir): Directory containing hardmasked ROI FASTA files
#   $2 (query_dir): Directory containing the query transcript FASTA file(s)
#   $3 (output_dir): Directory where Exonerate output files will be stored
#   $4 (single_contig_list): File containing a list of single-contig genomes
#   $5 (gene_to_annotate): Name of the gene being annotated, used to label output files (Example: "NKG2")
#   $6 (max_intron): Optional. Sets maximum intron size in basepairs  (default = 7000)
#   $7 (threads): Optional. Number of CPU cores to use for parallel processing (default: 1)
#
# Example:
#   exonerate_run "/path/to/hardmasked_ROI" "/path/to/transcript_queries" "/path/to/exonerate_results" "/path/to/single_contig_list.txt" "NKG2" 7000 4
#
# Notes:
# - Ensure that Exonerate is installed and accessible in the current environment
# - Exonerate is run with the following key options:
#   - --model cdna2genome: Aligns transcript to genomic DNA
#   - --refine region: Refines alignments in specific regions
#   - --bestn 1: Reports only the best alignment
#   - --percent 80: Requires 80% sequence identity


exonerate_run(){
    local ROI_hardmasked_dir="$1"
    local query_dir="$2"
    local output_dir="$3"
    local single_contig_list="$4"
    local gene_to_annotate="$5"
    local max_intron="$6"
    local threads="${7:-false}"
    
    while read genome; do
        echo "$genome"
        
        # Dynamically find the matching file for $genome
        fasta_file=$(compgen -G "$ROI_hardmasked_dir/${genome}*")
    
        # Prepare the exonerate command
        local cmd="exonerate \
        --refine region \
        --maxintron \"$max_intron\" \
        --bestn 1 \
        --percent 80 \
        --model cdna2genome \
        --fastasuffix .fasta \
        --querytype dna \
        --targettype dna \
        --query \"$query_dir\" \
        --target \"$fasta_file\" \
        --showvulgar no \
        --showalignment no \
        --showcigar no \
        --showtargetgff yes"
        
        # Add threads option if specified
        if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            cmd+=" -c $threads"
        fi
        
        # Run exonerate
        eval $cmd > "$output_dir/${genome}_ROI_${gene_to_annotate}_transcript_exonerate.gff"
    
    done < "$single_contig_list"
}


# Export the function so it can be called from outside the script
export -f exonerate_run
