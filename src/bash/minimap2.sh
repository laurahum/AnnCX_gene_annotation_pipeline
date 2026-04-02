#!/bin/bash

#Author: lahumada
#Date: 2026

# Function: minimap2_run  
# Description: Runs minimap2 on a list of genomes or fasta files against a query sequence
# Uses minimap2 with splice-aware settings for gene annotation
# This script will be used to
# - (Step 4) Find query genes to annotate in input hardmasked FASTA ROI files
#       - query: transcript
# - (Step 0) Find query genes to annotate in input unmasked or softmasked FASTA genome files
#	- query: transcript

# Usage: minimap2_run <input_target_fasta_dir> <input_query_fasta> <output_dir> <genome_list> <gene_to_annotate> <max_intron> <type_of_query> <type_of_fasta> [threads]
#
# Parameters:
#   $1 (input_target_fasta_dir): Directory containing target FASTA files (ROI hardmasked)
#   $2 (input_query_fasta): Query FASTA file (transcript or CDS)
#   $3 (output_dir): Directory for PAF output files
#   $4 (genome_list): Path to TXT file with genome names, one per line
#   $5 (gene_to_annotate): gene name (e.g. "NKG2")
#   $6 (max_intron): Optional. Sets maximum intron size in basepairs  (default = 7000)
#   $7 (type_of_query): Specifies the type of query sequence. Possible values:
#      - "flanking": Full gene sequence -> Recommended for queries using flanking genes
#      - "transcript": For queries using transcript sequences of genes to be annotated
#   $8 (type_of_fasta): Specifies the type of target FASTA sequence. Possible values:
#      - "ROI": Search the extracted ROI
#      - "genome": Seach whole genome
#   $9 (threads): Optional threads

# Example:
#   minimap2_run "/path/ROI_hardmasked" "/path/transcript.fasta" "/path/output" "/path/genome_list.txt" "NKG2" 4000 "transcript" "ROI" 8
#
# Notes:
# - Ensure that minimap2 is installed and accessible in the current environment
# - minimap2 is run with the following key options:
#   - -G: Maximum intron length
#   - -x splice:hq: Splice-aware alignment for high-quality RNA/DNA mapping
#   - --end-bonus=10: Bonus score for end-to-end full-length alignments
#   - -uf: Forward-reverse orientation only of the query
#   - --cs, -MD, -c: generate informative CIGAR to construct gene models in GFF3 in postprocessing
#   - -t: Number of threads (optional)
# - Outputs PAF format
#
# Output files:
#   ${genome}_${type_of_fasta}_${gene_to_annotate}_${type_of_query}_minimap2.paf



minimap2_run(){
    local input_target_fasta_dir="$1"
    local input_query_fasta="$2"
    local output_dir="$3"
    local genome_list="$4"
    local gene_to_annotate="$5"
    local max_intron="$6"
    local type_of_query="$7"
    local type_of_fasta="$8"
    local threads="${9:-false}"
        
        # Logic to annotate extracted ROI
        if [[ "$type_of_fasta" == "ROI" ]]; then
        	while read genome; do
        		echo "$genome"
        	
        		# Dynamically find the matching file for $genome
	        	target_fasta=$(compgen -G "$input_target_fasta_dir/${genome}*")
        	
        		# Prepare minimap2 command
      		  	local cmd="minimap2 -G $max_intron -x splice:hq --end-bonus=10 -uf --cs -MD -c \"$target_fasta\" \"$input_query_fasta\""
                
      		  	# Add threads option if specified
      		  	if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            		cmd+=" -t $threads"
      		  	fi
                          
      		  	# Run minimap2 → PAF output
      		  	eval $cmd > "$output_dir/${genome}_${type_of_fasta}_${gene_to_annotate}_${type_of_query}_minimap2.paf"
      		  	echo "Saved: ${genome}_${type_of_fasta}_${gene_to_annotate}_${type_of_query}_minimap2.paf"
        
    		done < "$genome_list"
    	
    	# Logic to annotate whole genome
        else
        	while read genome; do
    		    	echo "$genome"
        	
        		# INPUT directories and files
        		# Directory with each chromosome/scaffold for $genome
        		target_subdir="$input_target_fasta_dir/$genome"
        	
        		# Get each FASTA file within $target_subdir
        		target_files=($(compgen -G "$target_subdir/*"))
        	
        		# OUTPUT directories
        		output_subdir="$output_dir/$genome"
        		mkdir $output_subdir
        		
        		for target_fasta in "${target_files[@]}"; do
        			# Prepare minimap2 command
        			local cmd="minimap2 -G $max_intron -x splice:hq --end-bonus=10 -uf --cs -MD -c \"$target_fasta\" \"$input_query_fasta\""
                
        			# Add threads option if specified
        			if [[ "$threads" != "false" && "$threads" =~ ^[0-9]+$ ]]; then
            				cmd+=" -t $threads"
        			fi
                          
                          	file_name="${target_fasta##*/}"
                          	file_name_no_ext="${file_name%.*}"
                          	
        			# Run minimap2 → PAF output
        			eval $cmd > "$output_subdir/${file_name_no_ext}_${type_of_fasta}_${gene_to_annotate}_${type_of_query}_minimap2.paf"
        			echo "Saved: ${file_name_no_ext}_${type_of_fasta}_${gene_to_annotate}_${type_of_query}_minimap2.paf"
		
			done
		        
    		done < "$genome_list"
        fi
}

# Export the function so it can be called from outside the script
export -f minimap2_run

