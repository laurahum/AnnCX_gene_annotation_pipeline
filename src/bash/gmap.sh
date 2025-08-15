 #!/bin/bash

#Author: lahumada
#Date: 2022


# Function: gmap_run
# Description: Runs GMAP on a list of genomes or fasta files against a query sequence (gene, cDNA or exon)
# The genomes must have been processed beforehand using gmap_build

# This function reads a list with the genome names from a file and runs gmap on each input FASTA file 
# against a specified query sequence. The output is in GFF3 format.
#
# This script will be used in two steps of the pipeline:
#   1. Find flanking genes in input FASTA genome files
#       - query: gene or cDNA
#   2. Find query genes to annotate in input hardmasked FASTA ROI files
#       - query: cDNA or exons
#
# Usage: gmap_run <genome_list> <gmap_build_dir> <input_query_fasta> <output_dir> <gene_to_annotate> <type_of_query> <type_of_fasta> [threads]
#
# Parameters:
#   $1 (genome_list): Path to a TXT file containing a list of genome names, one per line
#   $2 (gmap_build_dir): Directory containing the GMAP databases for each genome
#   $3 (input_query_fasta): Directory containing the query FASTA file
#   $4 (output_dir): Directory where the output GFF3 files will be saved
#   $5 (gene_to_annotate): Name of the gene type to be annotated, used to label output files
#      The values can be either "flanking" for the flanking genes or 
#      the name of the genes of interest to be annotated (Example: "NKG2")         
#   $6 (type_of_query): Specifies the type of query sequence. Possible values:
#      - "flanking": Full gene sequence -> Recommended for queries using flanking genes
#      - "cDNA": For queries using cDNA sequences of genes to be annotated
#      - "exon": For queries using exon sequences of genes to be annotated
#   $7 (type_of_fasta): Specifies the type of FASTA file as input. Possible values:
#       - "genome"   
#       - "ROI"
#   $8 (threads): Optional. Number of threads to use for gmap. If not specified or set to "false", 
#                 gmap will use its default thread settings.
#
# Examples:
#   gmap_run "/path/to/genome_list.txt" "/path/to/gmap_databases" "/path/to/query_fasta" "/path/to/output" "NKG2" "cDNA" "ROI"
#   gmap_run "/path/to/genome_list.txt" "/path/to/gmap_databases" "/path/to/query_fasta" "/path/to/output" "NKG2" "cDNA" "ROI" 4
#
# Note: Ensure that GMAP is installed and accessible in the current environment


gmap_run(){
    local genome_list="$1"
    local gmap_build_dir="$2"
    local input_query_fasta="$3"
    local output_dir="$4"
    local gene_to_annotate="$5"
    local type_of_query="$6"
    local type_of_fasta="$7"
    local threads="${8:-false}"

    while read genome; do
            echo "$genome"
            
            # Prepare the gmap command
            local cmd="gmap -d ${genome} -D $gmap_build_dir -A $input_query_fasta --format=gff3_gene"
            
            # Add threads option if specified
            if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
                cmd+=" -t $threads"
            fi
            
            # Run gmap
            eval $cmd > "$output_dir/${genome}_${gene_to_annotate}_${type_of_query}_${type_of_fasta}_gmap.gff3"

    done < "$genome_list"
}

# Export the function so it can be called from outside the script
export -f gmap_run
