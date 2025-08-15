 #!/bin/bash

#Author: lahumada
#Date: 2025


# Function: make_augustus_protprof
# Description: Generates a protein profile for AUGUSTUS
#
# This function takes and input FASTA file with several protein sequences and:
#   1. Makes a multiple sequence alignment with mafft
#   2. Makes a protein profile with msa2prfl.pl (provided by AUGUSTUS installation)
#
# Usage: make_augustus_protprof <input_prot_query> <output_prot_prof>
#
# Parameters:
#   $1 (input_prot_query): File containing several FASTA protein sequences that will be used as query to find genes
#   $2 (output_prot_prof): Directory to save a protein profile to be used when running AUGUSTUS
#
# Example:
#   make_autustus_protprof "/path/to/input_prot_query.fasta" "/path/to/output_prot_prof/"
#
# Notes:
# - Ensure that AUGUSTUS is installed and accessible in the current environment
# - Ensure that mafft is installed and accessible in the current environment


make_augustus_protprof(){
    local input_prot_query="$1"
    local output_prot_prof="$2"
    local gene_to_annotate="$3"
    
    mafft "$input_prot_query" | msa2prfl.pl - > "${output_prot_prof}/${gene_to_annotate}_protprof.prf1"
    
    echo "PROTPROF_FOLDER=${output_prot_prof}/${gene_to_annotate}_protprof.prf1"
}

# Export the function so it can be called from outside the script
export -f make_augustus_protprof
