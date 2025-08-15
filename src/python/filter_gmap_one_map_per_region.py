#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 20:06:33 2024

@author: lahumada
"""


## This script will filter the GMAP results of aligning query cDNA/exons vs ROI hardmasked FASTA

## Since GMAP produces gene model annotations, each gene model (labeled as 'gene' in the GFF3 files) is considered
## for filtering and refered to as 'genes' (which includes its associated features: mRNA, exon, CDS) throughout the script.
## GMAP aligns several gene model annotation entries on the same region, with overlapping coordinates and different score values, 
## and the objective of this filtering step is to leave only the best gene model annotation entries.

## Logic: If the ranges from two gene model annotation entries overlap, filter lower identity.
##  -> If both ranges have the same identity, filter one at random

## Main function: filter_gmap_one_map_per_region(input_dir, output_dir)

import os
import sys
import numpy as np
import re
from .ranges_overlap import find_ranges_overlap # Function to determine whether two ranges overlap
from .read_gff3_to_df import read_gff3 # Custom function to read gff3



# Function to know what genes are mapped to the same region
    # Check whether a gene range is contained in another gene range
    # Function that returns genes to be filtered out
def genes_same_region (gene_rows, list_gene_start, list_gene_end, list_gene_names, list_identities, pos_gene):
    """Returns index of genes to be filtered out.
    Compares each gene alignment range to the rest.
    If the gene range checked is within another gene range 
    the function filters the gene with lower identity.
    If both identity values are the same, filter
    one gene at random.

    Arguments:
        gene_rows = dataframe # only gene rows from input_file
        list_gene_start = list(int) # start of gene alignment
        list_gene_end = list(int) # end of gene alignment
        list_gene_names = list(str) # name of each gene
        list_identities = list(int) # identities of each gene
        pos_gene = list(int) # index of genes in gene_rows
        
    Prints the reason for the filtering of each annotation entry
    """

    list_gene_same_region = []
        
    for i in range(len(gene_rows)):
        # Range of gene alignment to be checked
        range_gene_1 = range(list_gene_start[i],list_gene_end[i]+1)
            
        for j in range(len(gene_rows)):
            if i != j: # Prevents comparing one gene to itself
            # Range of gene alignment to be checked
                range_gene_2  = range(list_gene_start[j],list_gene_end[j]+1)
                    
                # Check if gene [i] is within gene [j]
                if find_ranges_overlap(range_gene_1, range_gene_2):
                    identity_gene_1 = list_identities[i]
                    identity_gene_2 = list_identities[j]
                    
                    # If the identity of gene [i] is higher, filter gene [j]
                    if identity_gene_1 > identity_gene_2:
                        gene_pos_1 = pos_gene[j]
                        gene_name_1 = list_gene_names[j]
                        list_gene_same_region.append(gene_pos_1)
                        print("Add_gene_" + gene_name_1 + "_index_" + str(gene_pos_1))
                        
                    # Else, filter gene [i]
                    elif identity_gene_1 < identity_gene_2:
                        gene_pos_2 = pos_gene[i]
                        gene_name_2 = list_gene_names[i]
                        list_gene_same_region.append(gene_pos_2)
                        print("Add_gene_" + gene_name_2 + "_index_" + str(gene_pos_2))
       
                    else:
                        # If both identities are the same
                        # Create a new random number generator, with a set seed for reproducibility
                        rng = np.random.default_rng(seed=42)
                        # Randomly select an index
                        index_random = rng.choice([i,j], size=1)[0]
                        gene_pos_3 = pos_gene[index_random]
                        gene_name_3 = list_gene_names[index_random]
                        list_gene_same_region.append(gene_pos_3)
                        print("Add_gene_" + gene_name_3 + "_index_" + str(gene_pos_3))
        
       
    # Get list of unique genes to be filtered
    unique_genes = np.unique(list_gene_same_region)
    
    # Control: print what genes will be filtered
    for i in unique_genes:
        index_gene = list(pos_gene).index(i)
        name_gene = list_gene_names[index_gene]
        print("Filtered_gene_" + name_gene)
    
    return(unique_genes)


## Function to filter genes in the input gmap gff3 file (Filepath_gmap_gff3)
## and save output filtered files (Filepath_filtered_gff3):
def filter_genes_same_map_region (Filepath_input, Filepath_output):
    '''Finds which genes are mapped to the same region and filters the ones
    that have the lower identity.
    1. Load dataframe from Filepath_input
    2. Format data
        - Remove nan values
        - Convert columns with gene range (3 and 4) from float to int
        - Reset index
    3. Get gene alignment data
        - Gene name
        - Gene start
        - Gene end
        - Gene identity
    4. Split dataframe annotation entries by gene model (gene + mRNA + exons + CDS)
    5. gene_same_region() 
        - Returns index of genes to be filtered out
    6. Save the dataframes from the genes excluding the filtered genes
    '''
        
    # 1. Load input gff3 file to a dataframe
    input_file = read_gff3(Filepath_input)
   
    # If the annotation file read is empty, skip this genome
    if input_file is None:
        print (f"No data rows in {Filepath_input}. Skipping this file.")
        return False 

    # 2. Format the gff3 dataframe
    ##### - Drop nan values in rows
    input_file = input_file.dropna()
    
    ##### - Convert column 3 and 4 from float to int
    input_file.iloc[:,3] = input_file.iloc[:,3].astype(int)
    input_file.iloc[:,4] = input_file.iloc[:,4].astype(int)    
    
    ##### - Reset index
    input_file = input_file.reset_index()
    original_index = input_file.loc[:,'index']
    input_file = input_file.drop(columns=['index'])
    
    
    # 3. Get data from the gene alignment 
    ##### - Make a dataframe only with the rows containing the alignment for 'gene':
    # Get the indices for the rows containing the word 'gene'
    # Create a boolean mask to split later the dataframe
    mask_gene = input_file.iloc[:,2] == 'gene'
    # Get the rows that have 'gene' in the second column
    gene_rows = input_file[~mask_gene == False]
    # Index of the dataframe containing just the genes
    pos_gene = gene_rows.index


    ##### - Get the name and interval of each gene alignment
    list_gene_names = []
    list_gene_start = []
    list_gene_end = []
    
    info_column = gene_rows.iloc[:,-1]
    
    
    for i in range(len(gene_rows)):
        # Get start and end of the gene alignment
        start_align = int(gene_rows.iloc[i,3])
        end_align = int(gene_rows.iloc[i,4])
        
        list_gene_start.append(start_align)
        list_gene_end.append(end_align)
        
        # Get each gene name
        start_name = info_column[pos_gene[i]].find("=")
        end_name = info_column[pos_gene[i]].find(";")
        gene_name = info_column[pos_gene[i]][(start_name+1):end_name]
        
        list_gene_names.append(gene_name)
        
        
    ##### - Get the identity of each gene alignment
    # The identity is written within the mRNA row in the last column
    # Get the rows that have 'mRNA' in the second column
    mask_mRNA = input_file.iloc[:,2] == 'mRNA'
    pos_mRNA = np.flatnonzero(mask_mRNA)
    mRNA_rows = input_file.iloc[pos_mRNA]
    
    # Get the identiy
    list_identities = []
    
    for i in range(len(pos_mRNA)):
        last_column = mRNA_rows.iloc[i,-1]
        
        identity_pattern = r'identity=(\d+(?:\.\d+)?)'
        identity_match = re.search(identity_pattern, last_column)
        identity_value = float(identity_match.group(1)) if identity_match else None

        list_identities.append(identity_value)
    

    # 4. Split dataframe by gene alignment (gene + mRNA + exons + CDS)
    # For loop to fill a dictionary: each entry is the dataframe from each gene
    df_by_gene = {}
    
    for i in range(len(pos_gene)):
        # For the last element of the gene index list pos_gene:
        if pos_gene[i] == pos_gene[-1]:
            df_by_gene[i] = input_file.iloc[pos_gene[i]:len(input_file),:]
        # For the rest elements of the list:
        else:
            start = pos_gene[i]
            end = pos_gene[i+1]
            df_by_gene[i] = input_file.iloc[start:end,:]
   

    # Call function to return gene indexes to be filtered
    index_same_region = genes_same_region(gene_rows, list_gene_start, list_gene_end, list_gene_names, list_identities, pos_gene)
       

    # - Save dataframes of the rest of genes minus filtered genes
    with open(Filepath_output, "w") as new_file:
        
        # Iterate over the list of dictionaries for each gene alignment:            
        for i in range(len(df_by_gene)):
            # Get indexes of each gene dictionary in df_by_gene
            index_df = df_by_gene[i].index
            
            # Omit the gene indexes to be filtered:
            if (any(item in index_same_region for item in index_df) == False):
                # Save each dataframe using '.to_csv'
                df_each = df_by_gene[i].to_csv(header=False, index=False, sep='\t')
                new_file.write(df_each)
                new_file.write('### \n')

    return True


# Loop over all the input files
def filter_gmap_one_map_per_region(input_dir, output_dir):
    '''Loop the filter_genes_same_map_region function over all the raw GMAP result annotation files for each genome
        - input_dir: Path to raw-GFF3 GMAP result annotation files
        - output_dir: Path to filtered-GFF3 GMAP result annotation files
    '''
    # Loop the function 'filter_genes_same_map_region()' over the input files  
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        filtered_name = filename
        filtered_name = filtered_name.rsplit('.', 1)[0] 
        filtered_name = filtered_name + '_FILTERED.gff3'
        
        # Directory to each input and output file:
        Filepath_gmap_gff3 = os.path.join(input_dir, filename)
        Filepath_filtered_gff3 = os.path.join(output_dir, filtered_name)
        
        # Call function to process input files in order to filter genes and save
        # output files:
        if filter_genes_same_map_region (Filepath_gmap_gff3, Filepath_filtered_gff3):
            # The annotation file was correct and filtering was successful
            print(f"Saved {filtered_name}")
        else:
            # The annotation file was empty and it was skipped
            print(f"Skipped {filtered_name}")
        # print empty line for readability
        print()

