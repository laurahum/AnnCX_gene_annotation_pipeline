 #!/bin/bash

#Author: lahumada
#Date: 2022


# Function: augustus_run
# Description: Performs gene prediction using AUGUSTUS on raw ROI sequences
#
# This function runs AUGUSTUS for each genome listed in a single contig list file.
# It uses a protein profile and species-specific parameters to predict genes in the 
# corresponding unmasked ROI sequences. The results are output in GFF3 format.
#
# Usage: augustus_run <ROI_dir> <proteinprofile> <output_dir> <single_contig_list> <gene_to_annotate> <species_arg_augustus>
#
# Parameters:
#   $1 (ROI_dir): Directory containing unmasked ROI FASTA files
#   $2 (proteinprofile): Path to the protein profile file for AUGUSTUS
#   $3 (output_dir): Directory where AUGUSTUS output files will be stored
#   $4 (single_contig_list): File containing a list of single-contig genomes
#   $5 (gene_to_annotate): Name of the gene being annotated, used to label output files (Example: "NKG2")
#   $6 (species_arg_augustus): Species argument for AUGUSTUS (Example: "human", "mouse")
#
# Example:
#   augustus_run "/path/to/ROI" "/path/to/protein_profile.prf1" "/path/to/augustus_results" "/path/to/single_contig_list.txt" "NKG2" "human"
#
# Notes:
# - Ensure that AUGUSTUS is installed and accessible in the current environment
# - AUGUSTUS is run with the following key options:
#   - --proteinprofile: Uses the specified protein profile for gene prediction
#   - --softmasking=0: Disables soft masking
#   - --gff3=on: Outputs results in GFF3 format
#   - --species: Uses species-specific parameters for gene prediction


augustus_run(){
    local ROI_dir="$1"
    local proteinprofile="$2"
    local output_dir="$3"
    local single_contig_list="$4"
    local gene_to_annotate="$5"
    local species_arg_augustus="$6"
    
    while read genome; do
    	echo "$genome"
    	
    	# Dynamically find the matching file for $genome
        fasta_file=$(compgen -G "$ROI_dir/${genome}*")
    	
    	# Run Augustus
    	augustus "$fasta_file" \
    	--proteinprofile="$proteinprofile" \
    	--softmasking=0	\
    	--gff3=on  \
    	--species="$species_arg_augustus" > "$output_dir/${genome}_ROI_${gene_to_annotate}_protprof_augustus.gff3"
	
    done < "$single_contig_list"
}

# Export the function so it can be called from outside the script
export -f augustus_run
