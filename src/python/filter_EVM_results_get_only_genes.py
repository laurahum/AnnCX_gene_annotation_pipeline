#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 20:43:30 2024

@author: lahumada
"""

## This script is the first step (1/3) to filter EVM results
## 1. Takes as input the raw EVM result files 
## 2. Writes a GFF3 file with only the gene entries in the EVM file

## Main function: filter_EVM_results_get_only_genes(input_dir, output_dir)


import os


# Function to figure out if it's an annotation line for a gene feature
def is_gene_line(line):
    ''' 
    Determine whether a given annotation line corresponds to a gene feature in a GFF3 file
    by checking if the provided line contains at least 3 tab-separated fields and if the
    thrid field (feature type) is a 'gene'.
        - line (str): Single line from a GFF3 file 
    
    Returns (bool): True if the line represents a gene feature, False otherwise
    '''
    fields = line.strip().split('\t')
    return len(fields) >= 3 and fields[2] == 'gene'


# Write to a new file only the gene features
def filter_genes(input_file, output_file):
    '''
    Extract and write only the gene feature lines from an input GFF3 file to an output file by
    reading the input file line by line and using the is_gene_line function to identify gene feature lines
        - input_file (str): Path to the input GFF3 files
        - output_file (str): Path where the filtered output will be saved
    '''
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if is_gene_line(line):
                outfile.write(line)
    print(f"Gene lines have been extracted to {output_file}")

# Main
def filter_EVM_results_get_only_genes(input_dir, output_dir):
    '''
    Process multiple EVM generated GFF3 files to extract gene features by iterating
    over the files in an input directory and apply the filter_genes function to each file
        - input_dir (str): Path to the directory containing EVM generated GFF3 files
        - output_dir (str): Path to the directory where the filtered fiels will be saved
    '''
    # Loop the main function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename
        output_name = output_name.rsplit('.', 1)[0] + '_genes_only.gff3'
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files
        filter_genes(Filepath_input, Filepath_output)
        print("Saved_" + output_name + '\n')
