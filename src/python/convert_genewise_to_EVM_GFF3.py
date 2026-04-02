#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 00:40:08 2024

@author: lahumada
"""


# This script converts formatted filtered GFF output from GeneWise (protein) to a format compatible with EVM
#
# 1. Convert features
# The output from GeneWise contains for features:
#     cds
# from those features EVM (protein alignment) accepts:
#     nucleotide_to_protein_match
# all the cds features in GeneWise formatted filtered GFF file must be converted to nucleotide_to_protein_match features
#
# 2. Modify attributes column
# The attributes column of the GeneWise formatted filtered GFF files must be modified according to the example provided by EVM
# Give the entries associated to the same query output alignment, in order of entry, a counter in the attributes section:
# ex. ID=match.genewise.1 for the GeneWise that are associated to the first query match entry.
# The EVM example for genewise has no field for score
#
# The following example file provided by EVM documentation was used as a guide to make the formatting changes
# https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/protein_alignments.gff3
#
# Input: GeneWise formatted filtered GFF output (with features: cds)
# Output: GeneWise formatted filtered EVM-compatible GFF3 file (with features: nucleotide_to_protein_match)

# Main function: convert_genewise_to_EVM_all_files(input_dir, output_dir, source)
# Note: source argument must be:
#   'genewise'


import re
import os

# Function to convert genewise to EVM
def convert_genewise_to_EVM(input_file, output_file, source):
    """
    Convert GeneWise formatted-filtered-GFF output to GFF3 format compatible with EVM.
    cds -> nucleotide_to_protein_match
        - input_file (str): Path to the input GeneWise output file
        - output_file (str): Path to the output GFF3 file compatible with EVM
        - source (str): Name of the tool used to generate the annotation output. Possible values:
            - genewise: run GeneWise with protein as query
    """    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        
        # None = last gene model finished writing
        current_gene = None
        
        for line in infile:
            # Skip empty lines and comments
            if line.strip() and not line.startswith('#'):
                fields = line.strip().split('\t')
                
                # Process only CDS features
                if fields[2].lower() == 'cds':
                    
                    # Extract necessary information for each field
                    seqid = fields[0]
                    start = fields[3]
                    end = fields[4]
                    score = '.'  # EVM example for genewise has no field for score (website in header)
                    strand = fields[6]
                    phase = '.'  # EVM examples contain no phase
                    
                    # Extract gene number and target from attributes
                    attributes = fields[8]
                    gene_match = re.search(rf'{re.escape(seqid)}-genewise-prediction-(\d+\.\d+)', attributes)
                    target_match = re.search(r'Target=([^;]+)', attributes)
                    
                    if gene_match and target_match:
                        gene_number = gene_match.group(1)
                        target = target_match.group(1)
                        
                        # If this is a new gene, update the current_gene
                        if gene_number != current_gene:
                            current_gene = gene_number
                        
                        # Construct new ID
                        new_id = f"match.{source}.{gene_number}"
                        
                        # Construct new attributes
                        new_attributes = f"ID={new_id};Target={target}"
                        
                        # Construct and Write the reformatted line
                        new_line = f"{seqid}\t{source}\tnucleotide_to_protein_match\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{new_attributes}\n"
                        outfile.write(new_line)

# Main
def convert_genewise_to_EVM_all_files(input_dir, output_dir, source):
    '''
    Loop the convert_genewise_to_EVM function over all the formatted filtered GeneWise result annotation files for each genome
        - input_dir: Path to formatted-filtered-GFF GeneWise result annotation files
        - output_dir: Path to formatted-filtered-EVM-GFF GeneWise result annotation files
        - source: Name of the tool used to generate the annotation output. Possible values:
            - genewise: run GeneWise with protein as query
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
        convert_genewise_to_EVM(Filepath_input, Filepath_output, source)
        print("Saved_" + output_name)

