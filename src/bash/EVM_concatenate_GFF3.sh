 #!/bin/bash

#Author: lahumada
#Date: 2024


# Function: EVM_concatenate_GFF3
# Description: Concatenates EVM-formatted-GFF3 files for EvidenceModeler (EVM) input
#
# This function processes each genome in a provided list, concatenating the files corresponding to:
#  - transcript_alignment (BLASTN, GMAP (cDNA, exons), Exonerate)
#  Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/transcript_alignment.gff3 
#  - protein_alignment (TBLASTN, GeneWise)
#  Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/protein_alignments.gff3
# These concatenated files, together with the gene_prediction files from AUGUSTUS, are used as input to run EVM.
#
# Usage: EVM_concatenate_GFF3 <input_dir_GeneWise> <input_dir_tblastn> <input_dir_blastn> <input_dir_gmap_cDNA> <input_dir_Exonerate> <input_dir_gmap_exons> <output_dir> <single_contig_list>
#
# Parameters:
#   $1 (input_dir_GeneWise): Directory containing formatted filtered formatted-for-evm GeneWise GFF3 files
#   $2 (input_dir_tblastn): Directory containing formatted filtered formatted-for-evm tblastn GFF3 files
#   $3 (input_dir_blastn): Directory containing formatted filtered formatted-for-evm blastn GFF3 files
#   $4 (input_dir_gmap_cDNA): Directory containing filtered formatted-for-evm GMAP cDNA GFF3 files
#   $5 (input_dir_Exonerate): Directory containing formatted filtered formatted-for-evm Exonerate GFF3 files
#   $6 (input_dir_gmap_exons): Directory containing filtered formatted-for-evm GMAP exon GFF3 files
#   $7 (output_dir): Directory where the concatenated GFF3 files will be stored
#   $8 (single_contig_list): File containing a list of genomes to process
#
# Example:
#   EVM_concatenate_GFF3 "/path/to/GeneWise" "/path/to/tblastn" "/path/to/blastn" "/path/to/gmap_cDNA" "/path/to/Exonerate" "/path/to/gmap_exons" "/path/to/output" "/path/to/genome_list.txt"
#
# - The function creates two output files per genome:
#   1. <genome>_protein_alignment.GFF3: Concatenation of GeneWise and tblastn results
#   2. <genome>_transcript_alignment.GFF3: Concatenation of blastn, gmap_cDNA, Exonerate, and gmap_exons results
#
# Internal Function:
# - process_genome: Handles the concatenation process for each individual genome


## Function to process each genome
process_genome() {
    local genome="$1"
    local input_dir_GeneWise="$2"
    local input_dir_tblastn="$3"
    local input_dir_blastn="$4"
    local input_dir_gmap_cDNA="$5"
    local input_dir_Exonerate="$6"
    local input_dir_gmap_exons="$7"
    local output_dir="$8"
    
    echo "Processing $genome..."

    # Protein alignment
    # Concatenate GeneWise, tblastn
    cat "$input_dir_GeneWise/${genome}"*_EVM.gff3 \
        "$input_dir_tblastn/${genome}"*_EVM.gff3 \
        > "$output_dir/${genome}_protein_alignment.GFF3"

    # Transcript alignment
    # Concatenate blastn, gmap_cDNA, gmap_exons and Exonerate
    cat "$input_dir_blastn/${genome}"*_EVM.gff3 \
        "$input_dir_gmap_cDNA/${genome}"*_EVM.gff3 \
        "$input_dir_Exonerate/${genome}"*_EVM.gff3 \
        "$input_dir_gmap_exons/${genome}"*_EVM.gff3 \
        > "$output_dir/${genome}_transcript_alignment.GFF3"

    echo "Completed processing $genome"
}


EVM_concatenate_GFF3(){
    local input_dir_GeneWise="$1"
    local input_dir_tblastn="$2"
    local input_dir_blastn="$3"
    local input_dir_gmap_cDNA="$4"
    local input_dir_Exonerate="$5"
    local input_dir_gmap_exons="$6"
    local output_dir="$7"
    local single_contig_list="$8"
    
    ## Main script
    # 1. Check if the genome file exists
    if [ ! -f "$single_contig_list" ]; then
        echo "Error: Species list file '$single_contig_list' not found."
        exit 1
    fi

    # 2. Process each genome
    while read genome; do
    
        echo "$genome"
    
        # Run function to concatenate the files for each genome
        process_genome "$genome" "$input_dir_GeneWise" "$input_dir_tblastn" "$input_dir_blastn" "$input_dir_gmap_cDNA" "$input_dir_Exonerate" "$input_dir_gmap_exons" "$output_dir"
    
    done < "$single_contig_list"

    echo "All genome processed successfully."
}

# Export the function so it can be called from outside the script
export -f EVM_concatenate_GFF3
