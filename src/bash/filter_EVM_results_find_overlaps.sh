#!/bin/bash

#Author: lahumada
#Date: 2024

# Function: filter_EVM_results_find_overlaps
# Description: Filters EvidenceModeler (EVM) results by finding overlaps with other annotation tools

# This script is the second step (2/3) to filter EVM results.
# Process:
#   For each of the genome with the ROI in a single contig:
#   1. Sort the EVM results (only the gene features using the output of
#   'filter_EVM_results_get_only_genes.py) and the annotation results from running the tools in the pipeline
#   2. Uses bedtools to find the intersects between the EVM results and each of the tools 
#   3. Creates a GFF-like file 'combined_overlaps.txt' that displays the amount of overlaps per gene feature
#
# Usage: filter_EVM_results_find_overlaps <input_dir_EVM> <input_dir_augustus> <input_dir_blastn> <input_dir_tblastn> <input_dir_gmap_genes> <input_dir_gmap_exons> <input_dir_genewise> <input_dir_exonerate> <output_dir> <single_contig_list>
#
# Parameters:
#   $1  (input_dir_EVM): Directory containing EVM GFF3 result files after 1/3 filtering
#   $2  (input_dir_augustus): Directory containing AUGUSTUS formatted-for-EVM GFF3 result files
#   $3  (input_dir_blastn): Directory containing BLASTN formatted filtered formatted-for-EVM GFF3 result files
#   $4  (input_dir_tblastn): Directory containing TBLASTN formatted filtered formatted-for-EVM GFF3 result files
#   $5  (input_dir_gmap_genes): Directory containing GMAP cDNA filtered formatted-for-EVM GFF3 result files
#   $6  (input_dir_gmap_exons): Directory containing GMAP exons filtered formatted-for-EVM GFF3 result files
#   $7  (input_dir_genewise): Directory containing GeneWise formatted filtered formatted-for-EVM GFF3 result files
#   $8  (input_dir_exonerate): Directory containing Exonerate formatted filtered formatted-for-EVM GFF3 result files
#   $9  (output_dir): Directory where output files will be stored
#   $10 (single_contig_list): File containing a list of genome to process
#
# Example:
#   filter_EVM_results_find_overlaps "/path/to/evm" "/path/to/augustus" "/path/to/blastn" "/path/to/tblastn" "/path/to/gmap_genes" "/path/to/gmap_exons" "/path/to/genewise" "/path/to/exonerate" "/path/to/output" "/path/to/genome_list.txt"
#
# Notes:
# - Requires bedtools to be installed and accessible in the system PATH
# - Output includes sorted files for each tool and EVM, intersection results, and a combined overlaps file
# - The combined_overlaps.txt file is sorted by the number of overlaps in descending order

# Output files (per genome):
# - <genome>_EVM.sorted.gff3: Sorted EVM results
# - <genome>_tool[1-7].sorted.gff3: Sorted results for each tool
# - <genome>_result_tool[1-7]_overlap.gff3: Overlap results between EVM and each tool
# - <genome>_combined_overlaps.txt: Summary of overlaps for each gene feature


filter_EVM_results_find_overlaps(){
    local input_dir_EVM="$1"
    local input_dir_augustus="$2"
    local input_dir_blastn="$3"
    local input_dir_tblastn="$4"
    local input_dir_gmap_cDNA="$5"
    local input_dir_gmap_exons="$6"
    local input_dir_genewise="$7"
    local input_dir_exonerate="$8"
    local output_dir="$9"
    local single_contig_list="${10}"

    while read genome; do
        echo "Filter_EVM_results_find_overlaps for $genome"
        
    	# Find overlaps
    	
    	# 1. Sort the input files
    	EVM=($input_dir_EVM/${genome}*.gff3)
    	augustus=($input_dir_augustus/${genome}*)
    	blastn=($input_dir_blastn/${genome}*)
    	tblastn=($input_dir_tblastn/${genome}*)
    	gmap_cDNA=($input_dir_gmap_cDNA/${genome}*)
    	gmap_exons=($input_dir_gmap_exons/${genome}*)
    	genewise=($input_dir_genewise/${genome}*)
    	exonerate=($input_dir_exonerate/${genome}*)
	
    	sort -k1,1 -k4,4n $EVM > "$output_dir/${genome}_EVM.sorted.gff3"
    	sort -k1,1 -k4,4n $augustus > "$output_dir/${genome}_tool1.sorted.gff3"
    	sort -k1,1 -k4,4n $blastn > "$output_dir/${genome}_tool2.sorted.gff3"
    	sort -k1,1 -k4,4n $tblastn > "$output_dir/${genome}_tool3.sorted.gff3"
    	sort -k1,1 -k4,4n $gmap_cDNA > "$output_dir/${genome}_tool4.sorted.gff3"
    	sort -k1,1 -k4,4n $gmap_exons > "$output_dir/${genome}_tool5.sorted.gff3"
    	sort -k1,1 -k4,4n $genewise > "$output_dir/${genome}_tool6.sorted.gff3"
    	sort -k1,1 -k4,4n $exonerate > "$output_dir/${genome}_tool7.sorted.gff3"
    	
    	# Run bedtools to find the intersects
    	bedtools intersect -a "$output_dir/${genome}_EVM.sorted.gff3" -b "$output_dir/${genome}_tool1.sorted.gff3" -wa -u > "$output_dir/${genome}_result_tool1_overlap.gff3"
    	bedtools intersect -a "$output_dir/${genome}_EVM.sorted.gff3" -b "$output_dir/${genome}_tool2.sorted.gff3" -wa -u > "$output_dir/${genome}_result_tool2_overlap.gff3"
    	bedtools intersect -a "$output_dir/${genome}_EVM.sorted.gff3" -b "$output_dir/${genome}_tool3.sorted.gff3" -wa -u > "$output_dir/${genome}_result_tool3_overlap.gff3"
    	bedtools intersect -a "$output_dir/${genome}_EVM.sorted.gff3" -b "$output_dir/${genome}_tool4.sorted.gff3" -wa -u > "$output_dir/${genome}_result_tool4_overlap.gff3"
    	bedtools intersect -a "$output_dir/${genome}_EVM.sorted.gff3" -b "$output_dir/${genome}_tool5.sorted.gff3" -wa -u > "$output_dir/${genome}_result_tool5_overlap.gff3"
    	bedtools intersect -a "$output_dir/${genome}_EVM.sorted.gff3" -b "$output_dir/${genome}_tool6.sorted.gff3" -wa -u > "$output_dir/${genome}_result_tool6_overlap.gff3"
    	bedtools intersect -a "$output_dir/${genome}_EVM.sorted.gff3" -b "$output_dir/${genome}_tool7.sorted.gff3" -wa -u > "$output_dir/${genome}_result_tool7_overlap.gff3"

    	# Create a GFF-like file that states how many overlaps each gene feature in EVM results has
    	cat "$output_dir/${genome}_result_tool1_overlap.gff3" "$output_dir/${genome}_result_tool2_overlap.gff3" "$output_dir/${genome}_result_tool3_overlap.gff3" "$output_dir/${genome}_result_tool4_overlap.gff3" "$output_dir/${genome}_result_tool5_overlap.gff3" "$output_dir/${genome}_result_tool6_overlap.gff3" "$output_dir/${genome}_result_tool7_overlap.gff3" | sort | uniq -c | sort -k1,1nr > "$output_dir/${genome}_combined_overlaps.txt"
	
    done < "$single_contig_list"
}

# Export the function so it can be called from outside the script
export -f filter_EVM_results_find_overlaps
