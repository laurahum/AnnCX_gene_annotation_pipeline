 #!/bin/bash

#Author: lahumada
#Date: 2022


# Function: makeblastdb_run
# Description: Creates BLAST databases from hardmasked ROI sequences
#
# This function processes hardmasked ROI sequences for single-contig genome
# and creates BLAST databases using makeblastdb. It reads a list of
# genomes and processes each corresponding hardmasked ROI file.
#
# Usage: makeblastdb_run <ROI_hardmasked_dir> <output_dir> <single_contig_list>
#
# Parameters:
#   $1 (ROI_hardmasked_dir): Directory containing hardmasked ROI FASTA files
#   $2 (output_dir): Directory where BLAST databases will be stored
#   $3 (single_contig_list): TXT file containing a list of single-contig genomes
#
# Example:
#   makeblastdb_run "/path/to/hardmasked_ROI" "/path/to/blast_db_output" "/path/to/single_contig_list.txt"
#
# Notes: 
# - Ensure that makeblastdb is installed and accessible in the current environment
# - makeblastdb is run with options: -dbtype nucl -parse_seqids



makeblastdb_run(){
    local ROI_hardmasked_dir="$1"
    local output_dir="$2"
    local single_contig_list="$3"

    while read genome; do
        echo "$genome"
        
        # Dynamically find the matching file for $genome
        fasta_file=$(compgen -G "$ROI_hardmasked_dir/${genome}*.fasta.masked")
            
        # Run makeblastdb on hard masked ROI sequences for single contig genome    
    	makeblastdb -in "$fasta_file" -dbtype nucl -parse_seqids -title "${genome}" -out "$output_dir/${genome}"
	
    done < "$single_contig_list"
}

# Export the function so it can be called from outside the script
export -f makeblastdb_run
