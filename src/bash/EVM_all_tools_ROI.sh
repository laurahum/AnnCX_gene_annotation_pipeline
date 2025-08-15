 #!/bin/bash

#Author: lahumada
#Date: 2024

# Function: EVM_run
# Description: Runs EvidenceModeler (EVM) on the output of AUGUSTUS, GMAP, BLAST, GeneWise and Exonerate
#
# This function executes EvidenceModeler for each genome listed in a provided file.
# It integrates various sources of evidence:
#   - gene_prediction: AUGUSTUS 
#   Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/gene_predictions.gff3
#   - transcript_alignment: BLASTN, GMAP (cDNA, exons), Exonerate 
#   Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/transcript_alignment.gff3 
#   - protein_alignment: TBLASTN, GeneWise
#   Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/protein_alignments.gff3
#
# Usage: EVM_run <ROI_hardmasked_dir> <input_Augustus_dir> <input_GFF3_dir> <output_dir> <output_dir_gff3> <EVM_weights> <single_contig_list> [threads]
#
# Parameters:
#   $1 (ROI_hardmasked_dir): Directory containing hardmasked ROI FASTA files
#   $2 (input_Augustus_dir): Directory containing AUGUSTUS gene predictions in formated-for-EVM GFF3 format
#   $3 (input_GFF3_dir): Directory containing both concatenated protein (GeneWise, tblastn) and transcript alignment (blastn, GMAP cDNA, GMAP exons, Exonerate in formated-for-EVM GFF3 files
#   $4 (output_dir): Directory where all EVM output will be stored
#   $5 (output_dir_gff3): Directory for storing EVM GFF3 output files
#   $6 (EVM_weights): File specifying weights for different evidence types in EVM
#   $7 (single_contig_list): File containing a list of genomes to process
#   $8 (threads): Optional. Number of CPUs to use for EVM processing (Default = 4)
#
# Example:
#   EVM_run "/path/to/masked_genomes" "/path/to/augustus" "/path/to/alignments" "/path/to/evm_output" "/path/to/gff3_output" "/path/to/weights.txt" "/path/to/genome_list.txt" 8
#
# Notes:
# - EVM must be installed and accessible
# - EVM is run with specific default parameters:
#   - segmentSize: 100000
#   - overlapSize: 10000
# - If threads is not specified, EVM will use its default CPU settings



EVM_run(){
    local ROI_hardmasked_dir="$1"
    local input_Augustus_dir="$2"
    local input_concatenated_dir="$3"
    local output_dir="$4"
    local output_dir_gff3="$5"
    local EVM_weights="$6"
    local single_contig_list="$7"
    local threads="${8:-false}"

    # Run EVM weights_2
    cd "$output_dir"
    while read genome; do
        echo "$genome"
        
        # Dynamically find the matching file for $genome
        input_fasta=$(compgen -G "$ROI_hardmasked_dir/${genome}*.fasta.masked")
        input_predictions=$(compgen -G "$input_Augustus_dir/${genome}*_EVM.gff3")
    
        # Prepare the EVM command
        local cmd="EVidenceModeler --sample_id \"$genome\" \
         --genome \"$input_fasta\" \
         --gene_predictions \"$input_predictions\" \
         --protein_alignments \"$input_concatenated_dir/${genome}_protein_alignment.GFF3\" \
         --transcript_alignments \"$input_concatenated_dir/${genome}_transcript_alignment.GFF3\" \
         --segmentSize 100000 \
         --overlapSize 10000 \
         --weights \"$EVM_weights\""
        
        # Add CPU option if specified
        if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            cmd+=" --CPU $threads"
        fi
       
        # Run EVM
        eval $cmd
        
    done < "$single_contig_list"

    # Move gff3 files to another directory
    for file in "$output_dir"/*.gff3; do
        cp "$file" "$output_dir_gff3"
    done

    cd
}

# Export the function so it can be called from outside the script
export -f EVM_run
