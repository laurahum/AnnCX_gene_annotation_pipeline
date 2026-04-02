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
# Usage: filter_EVM_results_find_overlaps <input_dir_EVM> <input_dir_augustus> <input_dir_blastn> <input_dir_tblastn> <input_dir_gmap_genes> <input_dir_gmap_exons> <input_dir_genewise> <input_dir_exonerate> <input_dir_minimap2> <input_dir_miniprot> <output_dir> <single_contig_list>
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
#   $9  (input_dir_minimap2): Directory containing Minimap2 formatted filtered formatted-for-EVM GFF3 result files
#   $10 (input_dir_miniprot): Directory containing Miniprot formatted filtered formatted-for-EVM GFF3 result files
#   $11 (output_dir): Directory where output files will be stored
#   $12 (single_contig_list): File containing a list of genome to process
#
# Example:
#   filter_EVM_results_find_overlaps "/path/to/evm" "/path/to/augustus" "/path/to/blastn" "/path/to/tblastn" "/path/to/gmap_genes" "/path/to/gmap_exons" "/path/to/genewise" "/path/to/exonerate" "/path/to/minimap2" "/path/to/miniprot" "/path/to/output" "/path/to/genome_list.txt"
#
# Notes:
# - Requires bedtools to be installed and accessible in the system PATH
# - Output includes sorted files for each tool and EVM, intersection results, and a combined overlaps file
# - The combined_overlaps.txt file is sorted by the number of overlaps in descending order

# Output files (per genome):
# - <genome>_EVM.sorted.gff3: Sorted EVM results
# - <genome>_tool[1-8].sorted.gff3: Sorted results for each tool
# - <genome>_result_tool[1-8]_overlap.gff3: Overlap results between EVM and each tool
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
    local input_dir_minimap2="$9"
    local input_dir_miniprot="${10}"
    local output_dir="${11}"
    local single_contig_list="${12}"

    while read genome; do
        echo "Filter_EVM_results_find_overlaps for $genome"
        
        # 1. Sort EVM (always required)
        EVM_files=("$input_dir_EVM/${genome}"*.gff3)
        if [[ -f "${EVM_files[0]}" ]]; then
            sort -k1,1 -k4,4n "${EVM_files[0]}" > "$output_dir/${genome}_EVM.sorted.gff3"
        else
            echo "Warning: No EVM file found for $genome, skipping"
            continue
        fi

        # 2. Check and sort each tool's files (tool1=augustus, tool2=blastn, etc.)
        declare -A tool_dirs=(
            ["tool1"]="$input_dir_augustus"
            ["tool2"]="$input_dir_blastn" 
            ["tool3"]="$input_dir_tblastn"
            ["tool4"]="$input_dir_gmap_cDNA"
            ["tool5"]="$input_dir_gmap_exons"
            ["tool6"]="$input_dir_genewise"
            ["tool7"]="$input_dir_exonerate"
            ["tool8"]="$input_dir_minimap2"
            ["tool9"]="$input_dir_miniprot"
        )
        
        declare -a overlap_files=()
        
        for tool_num in "${!tool_dirs[@]}"; do
            tool_dir="${tool_dirs[$tool_num]}"
            tool_files=("$tool_dir/${genome}"*)
            
            if [[ -f "${tool_files[0]}" ]]; then
                sort -k1,1 -k4,4n "${tool_files[0]}" > "$output_dir/${genome}_${tool_num}.sorted.gff3"
                
                # Run bedtools intersect
                bedtools intersect -a "$output_dir/${genome}_EVM.sorted.gff3" \
                                  -b "$output_dir/${genome}_${tool_num}.sorted.gff3" \
                                  -wa -u > "$output_dir/${genome}_result_${tool_num}_overlap.gff3"
                overlap_files+=("$output_dir/${genome}_result_${tool_num}_overlap.gff3")
            else
                echo "Warning: No files for $tool_num ($tool_dir) for $genome, skipping"
            fi
        done

        # 3. Combine overlap results (only existing ones)
        if [[ ${#overlap_files[@]} -gt 0 ]]; then
            cat "${overlap_files[@]}" | sort | uniq -c | sort -k1,1nr > "$output_dir/${genome}_combined_overlaps.txt"
        else
            echo "No overlap files to combine for $genome" > "$output_dir/${genome}_combined_overlaps.txt"
        fi
        
    done < "$single_contig_list"
}

# Export the function so it can be called from outside the script
export -f filter_EVM_results_find_overlaps
