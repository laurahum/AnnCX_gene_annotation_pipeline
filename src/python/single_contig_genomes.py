#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 20:00:00 2022

@author: lahumada
"""

## Considering only the list of genome in which both flanking genes are found:
# Use the GMAP results from the alignment of flanking regions to ROI and 
# divide the list of genome into two lists, one that contains genome where 
# the ROI is mapped to only one contig and genome in which the ROI is mapped 
# to two different contigs. Only the genome which have the ROI on a single 
# contig will be used for the downstream analysis. This is a limitation on how 
# well assembled is each genome.

# Main function: single_contig_genomes(gmap_flanking_dir, output_dir, input_genome_flanking, name_genes)

import os
import re
#from .read_gff3_to_df import read_gff3 # Custom function to read gff3

import pandas as pd
def read_gff3(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Keep lines that contain '###'
            if '###' in line:
                data.append(line)
            elif not line.startswith('#'): 
                # Skip comment lines
                fields = line.split('\t')
                # Only the fields that correspond to an annotation (have 9 fields)
                if len(fields) == 9:  
                    data.append(fields)
    
    # Create DataFrame only if there are data rows
    if not data:
        return None
    
    return pd.DataFrame(data)

# To find out if the ROI is mapped onto different chromosomes
# I must use the first column of the GFF3 file which has the chromosome information.

## Get the list of genome in which both flanking genes are found
def read_genome_flanking(input_genome_flanking):
    '''Read the TXT file with the list of genomes in which both flanking genes are found'''
    with open(input_genome_flanking, 'r') as genome_file:
        genome_list = [line.strip() for line in genome_file.readlines()]  # Remove '\n' from each line
    return genome_list


## Function that fills in the two lists of single and non-single contig genome
def single_or_not(filepath, filename, Single_contig, Non_single_contig, input_genome_flanking, name_genes):
    '''
    Check whether the ROI delimited by the flanking genes is contained 
    in a single contig or not and fill in the lists corresponding to either situation
    	- file_path = Path to a GFF3 file created by GMAP to annotate the flanking genes
        - file_name = Name of the GFF3 created by GMAP to annotate the flanking genes
        - Single_contig = list of genome in which the ROI is in a single contig
        - Non_single_contig = list of genome in which the ROI is not in a single contig
        - input_genome_flanking = list of genome in which both flanking genes were found by GMAP
        - name_genes = name of the gene type to be annotated
    '''

    # Read the gmap file
    gmap_gff3_file = read_gff3(filepath)
    
    # Drop nan values in rows
    gmap_gff3_file = gmap_gff3_file.dropna()
    
    ## Read list of genome with both flanking genes found:
    list_of_genome = read_genome_flanking(input_genome_flanking)    

    # Get the genome name from the GFF3 file name created by GMAP
    pattern = rf'^(.*?)_{name_genes}_flanking'
    match = re.search(pattern, filename)
    if match:
        genome = match.group(1)

    # Get only the column with the chromosome names
    chromosome_column_slice = gmap_gff3_file.iloc[:, 0]

    # Transform pandas.core.series.Series into a list
    chromosome_column_list = chromosome_column_slice.tolist()

    # Check if all the values in the list are the same using all()
    # If True = ROI is in a single contig
    # If False = ROI is in two different contig
    result = all(element == chromosome_column_list[0] for element in chromosome_column_list)
    
    # Only append if the genome is in the list of genome in which both flanking genes were found
    if genome in list_of_genome:
        if (result):
            Single_contig.append(genome)
        else:
            Non_single_contig.append(genome)


# Main function
def single_contig_genomes(gmap_flanking_dir, output_dir, input_genome_flanking, name_genes):
    '''
    Find in the GMAP results from the alignment of flanking genes to the genome FASTA files
    which genome contain the ROI in a single contig and which don't and write the
    name of the genome into two TXT files one for each situation.
    	- gmap_flanking_dir: Directory to GMAP GFF3 annotation files with the flanking genes
    	- output_dir: Directory to store the TXT files produced by this script:
    		- Single_contig_list.txt: genomes with the flanking genes in a single contig/scaffold/chromosome
    		- Non_single_contig_list.txt: genomes with the flanking genes in different contig/scaffold/chromosome
    	- input_genome_flanking: TXT file with genomes in which both flanking genes are found
    	- name_genes: name of the gene type to be annotated
    
    Return (result_filename_1): Path to the files containing the single contig genome
    '''        
    ## Make two empty lists: 
    # For the genome that have the ROI in the same contig:
    Single_contig = list()
    # For the genome that have the ROI in a different contig:
    Non_single_contig = list()
    
    ## Loop the function over the gff3 files created by gmap (genome vs ROI flanking regions)
    for filename in os.listdir(gmap_flanking_dir):     
        filepath = os.path.join(gmap_flanking_dir, filename)
        single_or_not(filepath, filename, Single_contig, Non_single_contig, input_genome_flanking, name_genes)

    # Order both lists alphabetically
    Single_contig.sort()
    Non_single_contig.sort()    
    
    # Single_contig
    result_filepath_1 = os.path.join(output_dir, 'Single_contig_list.txt')
    with open(result_filepath_1, 'w') as file:
    	file.writelines(f"{element}\n" for element in Single_contig)

    # Non_single_contig
    result_filepath_2 = os.path.join(output_dir, 'Non_single_contig_list.txt')
    with open(result_filepath_2, 'w') as file:
    	file.writelines(f"{element}\n" for element in Non_single_contig)

    return result_filepath_1
