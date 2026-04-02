 #!/bin/bash

#Author: lahumada
#Date: 2024

# Function: EVM_concatenate_GFF3
# Description: Concatenates EVM formatted GFF3 files for EvidenceModeler (EVM) input.
#
# This function processes each genome in a provided list and creates three pergenome files:
#  protein_alignment (GeneWise, tblastn, GMAP CDS, miniprot CDS) #  Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/protein_alignments.gff3
#  transcript_alignment (BLASTN, Exonerate, GMAP exon, minimap2 exon) #  Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/transcript_alignment.gff3 
#  gene_predictions (AUGUSTUS, GMAP cDNA, minimap2 gene models, miniprot gene models) #  Reference: https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/gene_predictions.gff3
#
# Usage:
#   EVM_concatenate_GFF3 \
#     <input_dir_GeneWise> <input_dir_tblastn> <input_dir_blastn> \
#     <input_dir_gmap_cDNA> <input_dir_gmapExon_CDS> <input_dir_gmapExon_exons> \
#     <input_dir_Exonerate> <input_dir_miniprot_CDS> <input_dir_miniprot_model> \
#     <input_dir_augustus> <input_dir_minimap2_exons> <input_dir_minimap2_model> \
#     <output_dir> <single_contig_list>
#
# Parameters:
#   $1  (input_dir_GeneWise): Directory with GeneWise EVM‑formatted protein GFF3 files.
#   $2  (input_dir_tblastn):  Directory with tblastn EVM‑formatted protein GFF3 files.
#   $3  (input_dir_blastn):   Directory with blastn EVM‑formatted transcript GFF3 files.
#   $4  (input_dir_gmap_cDNA): Directory with GMAP cDNA EVM‑formatted gene model GFF3 files.
#   $5  (input_dir_gmapExon_CDS): Directory with GMAP exon‑based CDS EVM‑formatted protein GFF3 files.
#   $6  (input_dir_gmapExon_exons): Directory with GMAP exon‑based exon EVM‑formatted transcript GFF3 files.
#   $7  (input_dir_Exonerate): Directory with Exonerate EVM‑formatted transcript GFF3 files.
#   $8  (input_dir_miniprot_CDS): Directory with miniprot CDS EVM‑formatted protein GFF3 files.
#   $9  (input_dir_miniprot_model): Directory with miniprot gene‑model EVM‑formatted GFF3 files.
#   $10 (input_dir_augustus): Directory with AUGUSTUS EVM‑formatted gene prediction GFF3 files.
#   $11 (input_dir_minimap2_exons): Directory with minimap2 exon‑based EVM‑formatted transcript GFF3 files.
#   $12 (input_dir_minimap2_model): Directory with minimap2 gene‑model EVM‑formatted GFF3 files.
#   $13 (output_dir): Directory where concatenated GFF3 files will be written.
#   $14 (single_contig_list): Text file with one genome ID per line.
#
# Pergenome outputs:
#   1. <genome>_protein_alignment.GFF3
#      = GeneWise + tblastn + GMAP exon CDS + miniprot CDS
#   2. <genome>_transcript_alignment.GFF3
#      = blastn + Exonerate + GMAP exon exons + minimap2 exons
#   3. <genome>_gene_predictions.GFF3
#      = AUGUSTUS + GMAP cDNA + minimap2 models + miniprot models
#
# Internal Function:
#   process_genome: Handles the concatenation for a single genome.

## Function to process each genome
process_genome() {
    local genome="$1"
    local input_dir_GeneWise="$2"
    local input_dir_tblastn="$3"
    local input_dir_blastn="$4"
    local input_dir_gmapcDNA="$5"
    local input_dir_gmapExon_CDS="$6"
    local input_dir_gmapExon_exons="$7"
    local input_dir_Exonerate="$8"
    local input_dir_miniprot_CDS="$9"
    local input_dir_miniprot_model="${10}"   
    local input_dir_augustus="${11}"
    local input_dir_minimap2_exons="${12}"   
    local input_dir_minimap2_model="${13}"   
    local output_dir="${14}"
    
    echo "Processing $genome..."

    # Protein alignment: GeneWise + tblastn + GMAP exon CDS + miniprot CDS
    > "$output_dir/${genome}_protein_alignment.GFF3"  # Create empty file first
    [[ -d "$input_dir_GeneWise" && -n "$(ls -A "$input_dir_GeneWise/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_GeneWise/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_protein_alignment.GFF3"
    [[ -d "$input_dir_tblastn" && -n "$(ls -A "$input_dir_tblastn/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_tblastn/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_protein_alignment.GFF3"
    [[ -d "$input_dir_gmapExon_CDS" && -n "$(ls -A "$input_dir_gmapExon_CDS/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_gmapExon_CDS/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_protein_alignment.GFF3"
    [[ -d "$input_dir_miniprot_CDS" && -n "$(ls -A "$input_dir_miniprot_CDS/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_miniprot_CDS/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_protein_alignment.GFF3"

    # Transcript alignment: blastn + Exonerate + GMAP exon exons + minimap2 exons  
    > "$output_dir/${genome}_transcript_alignment.GFF3"
    [[ -d "$input_dir_blastn" && -n "$(ls -A "$input_dir_blastn/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_blastn/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_transcript_alignment.GFF3"
    [[ -d "$input_dir_Exonerate" && -n "$(ls -A "$input_dir_Exonerate/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_Exonerate/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_transcript_alignment.GFF3"
    [[ -d "$input_dir_gmapExon_exons" && -n "$(ls -A "$input_dir_gmapExon_exons/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_gmapExon_exons/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_transcript_alignment.GFF3"
    [[ -d "$input_dir_minimap2_exons" && -n "$(ls -A "$input_dir_minimap2_exons/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_minimap2_exons/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_transcript_alignment.GFF3"

    # Gene predictions: AUGUSTUS + GMAP cDNA + minimap2 models + miniprot models
    > "$output_dir/${genome}_gene_predictions.GFF3"
    [[ -d "$input_dir_augustus" && -n "$(ls -A "$input_dir_augustus/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_augustus/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_gene_predictions.GFF3"
    [[ -d "$input_dir_gmapcDNA" && -n "$(ls -A "$input_dir_gmapcDNA/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_gmapcDNA/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_gene_predictions.GFF3"
    [[ -d "$input_dir_minimap2_model" && -n "$(ls -A "$input_dir_minimap2_model/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_minimap2_model/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_gene_predictions.GFF3"
    [[ -d "$input_dir_miniprot_model" && -n "$(ls -A "$input_dir_miniprot_model/${genome}"*_EVM.gff3 2>/dev/null)" ]] && cat "$input_dir_miniprot_model/${genome}"*_EVM.gff3 >> "$output_dir/${genome}_gene_predictions.GFF3"
     
    echo "Completed processing $genome"
}



EVM_concatenate_GFF3(){
    local input_dir_GeneWise="$1"
    local input_dir_tblastn="$2"
    local input_dir_blastn="$3"
    local input_dir_gmapcDNA="$4"
    local input_dir_gmapExon_CDS="$5"
    local input_dir_gmapExon_exons="$6"
    local input_dir_Exonerate="$7"
    local input_dir_miniprot_CDS="$8"
    local input_dir_miniprot_model="$9"   
    local input_dir_augustus="${10}"
    local input_dir_minimap2_exons="${11}"   
    local input_dir_minimap2_model="${12}" 
    local output_dir="${13}"
    local single_contig_list="${14}"
    
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
        process_genome "$genome" "$input_dir_GeneWise" "$input_dir_tblastn" "$input_dir_blastn" "$input_dir_gmapcDNA" "$input_dir_gmapExon_CDS" "$input_dir_gmapExon_exons" "$input_dir_Exonerate" "$input_dir_miniprot_CDS" "$input_dir_miniprot_model" "$input_dir_augustus" "$input_dir_minimap2_exons" "$input_dir_minimap2_model" "$output_dir"
    
    done < "$single_contig_list"

    echo "All genome processed successfully."
}

# Export the function so it can be called from outside the script
export -f EVM_concatenate_GFF3
