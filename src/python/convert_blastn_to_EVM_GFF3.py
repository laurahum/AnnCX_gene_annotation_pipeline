#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 01:16:27 2024

@author: lahumada
"""

# This script converts formatted filtered GFF output from BLASTN (cDNA) to a format compatible with EVM
#
# 1. Convert features
# The output from BLASTN contains for features:
#     BLASTCDS
# from those features EVM (transcript alignment) accepts:
#     cDNA_match
# all the BLASTCDS features in BLASTN formatted filtered GFF file must be converted to cDNA_match features
#
# 2. Modify attributes column
# The attributes column of the BLASTN formatted filtered GFF files must be modified according to the example provided by EVM
# Give the entries associated to the same query output alignment, in order of entry, a counter in the attributes section:
# ex. ID=match.blastn.1 for the BLASTCDS that are associated to the first query match entry.
# The score field in the output will have the bitscore instead of the evalue
# as it is more informative and can be compared
#
# The following example file provided by EVM documentation was used as a guide to make the formatting changes
# https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/transcript_alignments.gff3
#
# Input: BLASTN formatted filtered GFF output (with features: BLASTCDS) 
# Output: BLASTN formated filtered EVM-compatible GFF3 file (with features: cDNA_match)

# Main function: convert_blastn_to_EVM_all_files(input_dir, output_dir, source)
# Note: source argument must be:
#   'blastn'



import re
import os

# Function to convert tblastn to EVM
def convert_blastn_to_EVM(input_file, output_file, source):
    """
    Convert BLASTN formatted-filtered-GFF output to GFF3 format compatible with EVM.
    BLASTCDS -> cDNA_match
        - input_file (str): Path to the input BLASTN output file
        - output_file (str): Path to the output GFF3 file compatible with EVM
        - source (str): Name of the tool used to generate the annotation output. Possible values:
            - blastn: run BLASTN with cDNA as query
    """
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        
        protein_counter = {}
        
        for line in infile:
            # Skip empty lines and comments
            if line.strip() and not line.startswith('#'):
                fields = line.strip().split('\t')
                
                # Process only BLASTCDS features
                if fields[2].upper() == 'BLASTCDS':
                    
                    # Extract necessary information for each field
                    seqid = fields[0]
                    start = fields[3]
                    end = fields[4]
                    score = fields[5]  # This score is the e-value
                    strand = fields[6]
                    phase = '.'  # EVM examples contain no phase
                    
                    # Extract gene ID and other attributes
                    attributes = fields[8]
                    id_match = re.search(r'ID=Gene([^;]+)', attributes)
                    pident_match = re.search(r'pident=([\d.]+)', attributes)
                    length_match = re.search(r'length=([\d]+)', attributes)
                    bitscore_match = re.search(r'bitscore=([\d.]+)', attributes)
                    
                    if id_match and pident_match and length_match and bitscore_match:
                        protein_name = id_match.group(1)
                        pident = pident_match.group(1)
                        length = length_match.group(1)
                        bitscore = bitscore_match.group(1)  # The output will have the field[5] value as bitscore instead of evalue
                        
                        # Assign or get the match number for this protein
                        if protein_name not in protein_counter:
                            protein_counter[protein_name] = len(protein_counter) + 1
                        match_number = protein_counter[protein_name]
                        
                        # Construct new ID
                        new_id = f"match.{source}.{match_number}"
                        
                        # Construct new attributes
                        new_attributes = f"ID={new_id};Target={protein_name};pident={pident};length={length};bitscore={bitscore}"
                        
                        # Construct and Write the reformatted line
                        new_line = f"{seqid}\t{source}\tcDNA_match\t{start}\t{end}\t{bitscore}\t{strand}\t{phase}\t{new_attributes}\n"
                        outfile.write(new_line)

# Main
def convert_blastn_to_EVM_all_files(input_dir, output_dir, source):
    '''
    Loop the convert_blastn_to_EVM function over all the formatted filtered BLASTN result annotation files for each genome
        - input_dir: Path to formatted-filtered-GFF BLASTN result annotation files
        - output_dir: Path to formatted-filtered-EVM-GFF BLASTN result annotation files
        - source: Name of the tool used to generate the annotation output. Possible values:
            - blastn: run BLASTN with cDNA as query
    '''
    # Loop the main function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename
        output_name = output_name.rsplit('.', 1)[0] + '_EVM.gff3'
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files
        convert_blastn_to_EVM(Filepath_input, Filepath_output, source)
        print("Saved_" + output_name + '\n')
