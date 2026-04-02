 #!/bin/bash

#Author: Charlotte Hoeltermann
#'The source for the awk code that turns blast output format 6 to GFF3 was taken 
#from https://github.com/raymondkiu/blastoutput2gff (Accessed in Juli 2021). 
#It was modified to fit the specific requirements of this context, and to allow 
#conversion of blast hits on the negative strand'
#Date: 2021
#Modified: lahumada



# Change the format for the files produced by BLASTN or TBLASTN for all genomes
# Functions: 
#   - format_blastn_out6_to_gff
#   - format_tblastn_out6_to_gff
# Description: Converts BLAST output (outfmt 6) to GFF format and performs additional formatting

# These functions process BLAST output files (blastn and tblastn) for each genome listed in a single contig list file.
# They convert the BLAST output to GFF format and perform additional formatting to:
#   - set the strand orientation 
#   - swap the start and end position if the start is > end position
#
# Usage: 
# format_blastn_out6_to_gff <input_dir_out6> <output_dir_gff> <output_dir_gff_formatted> <single_contig_list> <gene_to_annotate>
# format_tblastn_out6_to_gff <input_dir_out6> <output_dir_gff> <output_dir_gff_formatted> <single_contig_list> <gene_to_annotate>

# Parameters:
#   $1 (input_dir_out6): Directory containing BLAST output files in outfmt 6
#   $2 (output_dir_gff): Directory for initial GFF conversion output
#   $3 (output_dir_gff_formatted): Directory for final formatted GFF files
#   $4 (single_contig_list): File containing a list of single-contig genomes
#   $5 (gene_to_annotate): Name of the gene being annotated, used in output file names (Example: "NKG2")

# Example:
#   format_blastn_out6_to_gff "/path/to/blast_out6" "/path/to/gff_output" "/path/to/formatted_gff" "/path/to/single_contig_list.txt" "NKG2"



format_blastn_out6_to_gff(){
    local input_dir_out6="$1"
    local output_dir_gff="$2"
    local output_dir_gff_formatted="$3"
    local single_contig_list="$4"
    local gene_to_annotate="$5"
    
    # Process output files from script: blastn_query_ROI_hardmasked_out6.sh
    # 1. Turn blast outfmt 6 to .gff for viewing in Artemis:
    while read genome; do
        echo "$genome"
        echo "Format blastn step 1/2"
        
        # Dynamically find the matching file for $genome
        input_file=$(compgen -G "$input_dir_out6/${genome}*")
	
        # 1. awk to format the order of the tab-delimited fields
        # 2. sed to convert the words minus or plus for the strand orientation to - or +
        awk '{print $2"\tblast\tBLASTCDS\t"$9"\t"$10"\t"$11"\t"$13"\t.\tID=Gene"$1";Name="$2";pident="$3";length="$4";bitscore="$12";qcovs="$14";gapopen="$6}' "$input_file" | sed "s/minus/-/g" | sed "s/plus/+/g" > "$output_dir_gff/${genome}_ROI_${gene_to_annotate}_transcript_blastn.gff"
    
    done < "$single_contig_list"  
    
    # 2. Format gff files: Swap the start and end position if the start is > end position
    while read genome; do
        echo "$genome"
        echo "Format blastn step 2/2"
        
        while read line; do
        	# get content of line 7 (strand information):
            ln=$(echo $line | awk -F ' ' '{ print $7 }')
            # If "-" -> start > stop -> invert order
            if [ "$ln" = "-" ]; then
                echo $line | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9}' >> "$output_dir_gff_formatted/${genome}_ROI_${gene_to_annotate}_transcript_blastn_formatted.gff"
            # else "+" -> start < stop -> keep order
            else
                echo $line | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' >> "$output_dir_gff_formatted/${genome}_ROI_${gene_to_annotate}_transcript_blastn_formatted.gff"
            fi
        done < "$output_dir_gff/${genome}_ROI_${gene_to_annotate}_transcript_blastn.gff"

    done < "$single_contig_list" 

}


format_tblastn_out6_to_gff(){
    local input_dir_out6="$1"
    local output_dir_gff="$2"
    local output_dir_gff_formatted="$3"
    local single_contig_list="$4"
    local gene_to_annotate="$5"

    # Process output files from script: tblastn_query_ROI_hardmasked_out6.sh
    # 1. Turn blast outfmt 6 to .gff for viewing in Artemis:
    while read genome; do
        echo "$genome"
	echo "Format tblastn step 1/2" 
	        
        # Dynamically find the matching file for $genome
        input_file=$(compgen -G "$input_dir_out6/${genome}*")
           
        # 1. awk to format the order to the tab-delimited fields
        # The output of tblastn did not produce plus or minus, the column is N/A instead
        awk '{print $2"\tblast\tBLASTCDS\t"$9"\t"$10"\t"$11"\t"$13"\t.\tID=Gene"$1";Name="$2";pident="$3";length="$4";bitscore="$12";qcovs="$14";gapopen="$6}' "$input_file" > "$output_dir_gff/${genome}_ROI_${gene_to_annotate}_protein_tblastn.gff"

    done < "$single_contig_list"
    
    # 2. Format gff files: Swap the start and end position if the start is > end position
    while read genome; do
        echo "$genome"
	echo "Format tblastn step 2/2" 
	
        while read line; do
            # get start and stop
        	sta=$(echo $line | awk -F ' ' '{ print $4 }')
        	sto=$(echo $line | awk -F ' ' '{ print $5 }')
        	# If start > stop -> invert order -> name column 7 "-"
        	if [ "$sta" -gt "$sto" ]; then
        	    echo $line | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t-\t"$8"\t"$9}' >> "$output_dir_gff_formatted/${genome}_ROI_${gene_to_annotate}_protein_tblastn_formatted.gff"
        	# If start < end -> keep order
        	else
        	    echo $line | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t+\t"$8"\t"$9}' >> "$output_dir_gff_formatted/${genome}_ROI_${gene_to_annotate}_protein_tblastn_formatted.gff"
        	fi
        done < "$output_dir_gff/${genome}_ROI_${gene_to_annotate}_protein_tblastn.gff"   

    done < "$single_contig_list"
}


# Export the functions so they can be called from outside the script
export -f format_blastn_out6_to_gff
export -f format_tblastn_out6_to_gff
