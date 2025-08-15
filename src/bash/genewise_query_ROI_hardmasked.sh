 #!/bin/bash

#Author: lahumada
#Date: 2022


# Function: genewise_run
# Description: Performs gene annotation using Genewise on hardmasked ROI sequences
#
# This function runs Genewise for each genome listed in a single contig list file.
# It uses protein sequences as queries to annotate genes in the corresponding 
# hardmasked ROI sequences. The results are output in GFF format.
#
# Usage: genewise_run <ROI_hardmasked_dir> <query_dir> <output_dir> <single_contig_list> <gene_to_annotate>
#
# Parameters:
#   $1 (ROI_hardmasked_dir): Directory containing hardmasked ROI FASTA files
#   $2 (query_dir): Directory containing the query protein database file
#   $3 (output_dir): Directory where Genewise output files will be stored
#   $4 (single_contig_list): File containing a list of single-contig genomes
#   $5 (gene_to_annotate): Name of the genes of interest to be annotated, used to label output files (Example: "NKG2") 
#
# Example:
#   genewise_run "/path/to/hardmasked_ROI" "/path/to/protein_queries" "/path/to/genewise_results" "/path/to/single_contig_list.txt" "NKG2"
#
# Notes:
# - Ensure that Genewise is installed and accessible in the current environment
# - Must set environment variable WISECONFIGDIR 
#   Example: (genewise_env) user@computer:~$ export WISECONFIGDIR=/path/to/anaconda3/envs/genewise_env/share/wise2/wisecfg
# - Genewise is run with options: -prodb for protein database and -dnas for DNA sequences


genewise_run(){
    local ROI_hardmasked_dir="$1"
    local query_dir="$2"
    local output_dir="$3"
    local single_contig_list="$4"
    local gene_to_annotate="$5"
  
    while read genome; do
        echo "$genome"
        
        # Dynamically find the matching file for $genome
        fasta_file=$(compgen -G "$ROI_hardmasked_dir/${genome}*.fasta.masked")
	
        # Run genewise
    	genewisedb -prodb "$query_dir" -dnas "$fasta_file" -gff > "$output_dir/${genome}_ROI_${gene_to_annotate}_prot_genewise.gff"
	
	
    done < "$single_contig_list"   
}

# Export the function so it can be called from outside the script
export -f genewise_run


