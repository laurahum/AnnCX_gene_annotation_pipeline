#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 23:24:02 2022

@author: lahumada
"""

### This script will edit the output files of running exonerate
### to get rid of anything that is not
### an annotation entry. Substitute all no annotation lines with '##'

## Main function: format_exonerate_output_to_gff_all_files (input_dir, output_dir)

import os



# Function that does all the formatting and saving to output file
def format_exonerate_output_to_gff (Filepath_input, Filepath_output):
    '''
       Format Exonerate raw output format to a GFF format:
        1. Reads the input file line by line.
        2. Gets rid of all the lines that are not an annotation entry
        3. Writes the GFF formatted entries to the output file.
    '''
    # Open the input file to be formated as read
    file=open(Filepath_input,'r')

    # Create a list of lists with each of the lines for the output file
    output_list =  []

    # Go through each line of the input file
    for line in file:
        # remove the enter from each line
        line=line.rstrip('\n')

        # Only the lines that have 9 elements
        if len(line.split('\t')) == 9: #make sure it's data line and not header
            # Define each element of the line what they are based on the column order of a gff file
            scaffold,source,type,start,end,score,strand,phase,attributes=line.split('\t')
            # All the gff elements as a list
            line_annotation=("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (scaffold,source,type,start,end,score,strand,phase,attributes))
            # Append to the general list (output_list)
            output_list.append(line_annotation)
        
    # Save the list of lists to an output file
    with open (Filepath_output, 'a+') as output_file:
        for element in output_list:
            output_file.write(element + '\n')
            
            
# Loop over all input files
def format_exonerate_output_to_gff_all_files (input_dir, output_dir):
    '''Loop the format_exonerate_output_to_gff function over all the raw Exonerate result annotation files for each genome
        - input_dir: Path to raw Exonerate result annotation files
        - output_dir: Path to formatted-GFF Exonerate result annotation files
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
        format_exonerate_output_to_gff(Filepath_input, Filepath_output)
        print("Saved_" + output_name)
