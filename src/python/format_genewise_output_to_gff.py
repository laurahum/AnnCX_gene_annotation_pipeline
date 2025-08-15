#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 23:48:03 2022

@author: lahumada
"""

### This script will edit the output files of running genewise 
### to make sure that start<end and to filter
### out the 'intron' and 'match' hits, leaving only the 'cds' hits.
### Removes any line that is not an annotation entry
### Makes sure that the start and end features are always start<end 
### Adds bitscore from the match hits on its associated cds features
### Adds a gene counter in the attributes 
### Adds the protein query used to generate the match to the attributes

## Main function: format_genewise_output_to_gff_all_files(input_dir, output_dir)

import os



# Define a function that does all the formatting
def format_genewise_output_to_gff(Filepath_input, Filepath_output):
    ''' 
    Format GeneWise raw output format to a GFF format:
        1. Reads the input file line by line.
        2. Search for lines starting with '>Results for' which correspond to the start of the annotation features found by GeneWise for a specific sequence used as query
            - Keeps track of gene entries (type 'match') by incrementing a gene counter.
            - Excludes entries of type 'intron' and 'match' since they interfere with Artemis visualization
            - Reformats CDS entries by adding gene number, target protein, and match score to attributes.
            - Ensures start position is always smaller than end position.
        3. Writes the GFF formatted entries to the output file.
    '''
    # Open the input file to be re-formated as read
    with open(Filepath_input, 'r') as file:
        lines = file.readlines()

    # Create a list to store the formatted output lines
    output_list = []
    gene_counter = 0
    current_gene_number = None
    current_protein = None
    current_score = None

    # Go through each line of the input file
    for line in lines:
        line = line.strip()

        # Check for protein information
        # This corresponds to the start of the annotation features found for a specific sequence used as query
        if line.startswith('>Results for'):
            current_protein = line.split()[2].split('(')[0]
            continue

        # Only process lines that have 9 elements
        if len(line.split('\t')) == 9:
            scaffold, source, type, start, end, score, strand, phase, attributes = line.split('\t')
            
            # Check if this is a new gene entry
            if type == 'match':
                gene_counter += 1
                current_gene_number = gene_counter
                current_score = score
            
            # Exclude the hits that are labeled as intron and match since they interfere with Artemis visualization
            if type not in ['intron', 'match']:
                # Update the attributes of CDS entries
                if type == 'cds':
                    attributes = f"{attributes}.{current_gene_number};Target={current_protein};match_score={current_score}"
                
                # Ensure start is always smaller than end
                if int(start) > int(end):
                    start, end = end, start
                
                # Create the formatted line
                formatted_line = f"{scaffold}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}"
                output_list.append(formatted_line)

    # Save the list of formatted lines to the output file
    with open(Filepath_output, 'w') as output_file:
        for element in output_list:
            output_file.write(element + '\n')
            

def format_genewise_output_to_gff_all_files(input_dir, output_dir):
    '''Loop the format_genewise_output_to_gff function over all the raw GeneWise result annotation files for each genome
        - input_dir: Path to raw GeneWise result annotation files
        - output_dir: Path to formatted-GFF GeneWise result annotation files
    '''
    #### Loop the function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename
        output_name = output_name.rsplit('.', 1)[0]    
        output_name = output_name + '_FORMATTED.gff'
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files in order to filter genes and save
        # output files:
        format_genewise_output_to_gff(Filepath_input, Filepath_output)
        print("Saved_" + output_name + '\n')
    