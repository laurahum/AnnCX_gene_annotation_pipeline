 #!/bin/bash

#Author: lahumada
#Date: 2024

# Function: EVM_run
# Description: Runs EvidenceModeler (EVM) on the output of AUGUSTUS, GMAP, BLAST, GeneWise, Exonerate, Minimap2, Miniprot
#
# This function executes EvidenceModeler for each genome listed in a provided file.
# It integrates various sources of evidence:
#   - gene_prediction: AUGUSTUS, GMAP (cDNA), minimap2 -> gene model, miniprot -> gene model
#   Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/gene_predictions.gff3
#   - transcript_alignment: BLASTN, GMAP (exons) -> exons, Exonerate, minimap2 -> exons
#   Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/transcript_alignment.gff3 
#   - protein_alignment: TBLASTN, GeneWise, GMAP (exons) -> CDS, miniprot -> CDS
#   Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/protein_alignments.gff3
#
# Usage: EVM_run <ROI_hardmasked_dir> <input_GFF3_dir> <output_dir> <output_dir_gff3> <EVM_weights> <single_contig_list> [threads]
#
# Parameters:
#   $1 (ROI_hardmasked_dir): Directory containing hardmasked ROI FASTA files
#   $2 (input_GFF3_dir): Directory containing concatenated protein, transcript alignment and gene_prediction files
#   $3 (output_dir): Directory where all EVM output will be stored
#   $4 (output_dir_gff3): Directory for storing EVM GFF3 output files
#   $5 (EVM_weights): File specifying weights for different evidence types in EVM
#   $6 (single_contig_list): File containing a list of genomes to process
#   $7 (threads): Optional. Number of CPUs to use for EVM processing (Default = 4)
#
# Example:
#   EVM_run "/path/to/masked_genomes" "/path/to/alignments" "/path/to/evm_output" "/path/to/gff3_output" "/path/to/weights.txt" "/path/to/genome_list.txt" 8
#
# Notes:
# - EVM must be installed and accessible
# - EVM is run with specific default parameters:
#   - segmentSize: 100000
#   - overlapSize: 10000
# - If threads is not specified, EVM will use its default CPU settings

EVM_run(){
    local ROI_hardmasked_dir="$1"
    local input_concatenated_dir="$2"
    local output_dir="$3"
    local output_dir_gff3="$4"
    local EVM_weights="$5"
    local single_contig_list="$6"
    local threads="${7:-false}"

    # Run EVM
    cd "$output_dir"
    while read genome; do
        echo "$genome"
        
        # Dynamically find the matching file for $genome
        input_fasta=$(compgen -G "$ROI_hardmasked_dir/${genome}*")
    
        # Prepare the EVM command
        local cmd="EVidenceModeler --sample_id \"$genome\" \
         --genome \"$input_fasta\" \
         --gene_predictions \"$input_concatenated_dir/${genome}_gene_predictions.GFF3\" \
         --protein_alignments \"$input_concatenated_dir/${genome}_protein_alignment.GFF3\" \
         --transcript_alignments \"$input_concatenated_dir/${genome}_transcript_alignment.GFF3\" \
         --segmentSize 100000 \
         --overlapSize 10000 \
         --search_long_introns 1 \
         --re_search_intergenic 1 \
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

    cd -
}

# Export the function so it can be called from outside the script
export -f EVM_run
