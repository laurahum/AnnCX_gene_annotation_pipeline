 #!/bin/bash

#Author: lahumada
#Date: 2022


# Function: gmap_run
# Description: Runs GMAP on a list of genomes or fasta files against a query sequence (gene, transcript or exon)
# The genomes must have been processed beforehand using gmap_build

# This function reads a list with the genome names from a file and runs gmap on each input FASTA file 
# against a specified query sequence. The output is in GFF3 format.
#
# This script will be used for:
# This script will be used in two steps of the pipeline:
#   1. Find flanking genes in input FASTA genome files
#       - query: gene or transcript
#   2. Find query genes to annotate in input hardmasked FASTA ROI files
#       - query: transcript or exons
# The flags max_intron_middle and max_intron_ends use the same input from the arguments that the
# user provides at the start stating the max length intron.
#
# Usage: gmap_run <genome_list> <gmap_build_dir> <input_query_fasta> <output_dir> <gene_to_annotate> <type_of_query> <max_intron_middle> <max_intron_ends> [threads] <type_of_fasta> 
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
#      - "transcript": For queries using transcript sequences of genes to be annotated
#      - "exon": For queries using exon sequences of genes to be annotated
#   $7 (max_intron_middle): Optional. Sets maximum middle intron size in basepairs  (default = 7000)
#   $8 (max_intron_ends): Optional. Sets maximum ends intron size in basepairs  (default = 7000)
#   $9 (threads): Optional. Number of threads to use for gmap. If not specified or set to "false", 
#                 gmap will use its default thread settings.
#   $10 (type_of_fasta): Specifies the type of FASTA file as input. Possible values:
#       - "genome"   
#       - "ROI"
#
# Examples:
#   gmap_run "/path/to/genome_list.txt" "/path/to/gmap_databases" "/path/to/query_fasta" "/path/to/output" "NKG2" "transcript" "6000" "6000" 4 "ROI"
#   gmap_run "/path/to/genome_list.txt" "/path/to/gmap_databases" "/path/to/query_fasta" "/path/to/output" "NKG2" "exon" "6000" "6000" 4 "ROI"
#
# Note: Ensure that GMAP is installed and accessible in the current environment
# - GMAP is run with the following key options:
#	--max-intronlength-middle
#	--max-intronlength-ends
#	--split-large-introns


gmap_run(){
    local genome_list="$1"
    local gmap_build_dir="$2"
    local input_query_fasta="$3"
    local output_dir="$4"
    local gene_to_annotate="$5"
    local type_of_query="$6"
    local max_intron_middle="$7"
    local max_intron_ends="$8"
    local threads="${9:-false}"
    local type_of_fasta="${10}"

    while read genome; do
            echo "$genome"
            
            # Run GMAP to extract ROI
            if [[ "$type_of_fasta" == "flanking" ]]; then
            	# Prepare the gmap command
            	local cmd="gmap -d ${genome} -D $gmap_build_dir -A $input_query_fasta --format=gff3_gene"
            
            	# Add threads option if specified
            	if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
                	cmd+=" -t $threads"
            	fi
            
            	# Run gmap
            	eval $cmd > "$output_dir/${genome}_${gene_to_annotate}_${type_of_query}_${type_of_fasta}_gmap.gff3"

           	
            else # Run GMAP to find query in ROI
            	# Prepare the gmap command
            	local cmd="gmap -d ${genome} -D $gmap_build_dir -A $input_query_fasta --format=gff3_gene --max-intronlength-middle ${max_intron_middle} --max-intronlength-ends ${max_intron_ends} --split-large-introns"
            
            	# Add threads option if specified
            	if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
                	cmd+=" -t $threads"
            	fi
            
            	# Run gmap
           	eval $cmd > "$output_dir/${genome}_${gene_to_annotate}_${type_of_query}_${type_of_fasta}_gmap.gff3"
            
            fi

    done < "$genome_list"
}

# Export the function so it can be called from outside the script
export -f gmap_run
