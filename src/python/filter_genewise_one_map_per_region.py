#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 19:41:55 2024

@author: lahumada
"""


## This script will filter the GeneWise results of aligning query protein vs ROI hardmasked FASTA

## Since GeneWise does not produce gene model annotations, each annotation entry in the GeneWise results is considered
## for filtering and refered to as 'genes' throughout the script. 
## After the previous formatting step, the annotation entries in GeneWise results correspond to the CDS features.
## GeneWise aligns several annotation entries on the same region, with overlapping coordinates and different score values, 
## and the objective of this filtering step is to leave only the best annotation entries.

## Logic: If two ranges overlap, filter lower bitscore.
##  -> If both ranges have the same bitscore, filter one at random

## Main function: filter_genewise_one_map_per_region(input_dir, output_dir)

import os
import re
import numpy as np
import pandas as pd
import sys
from .ranges_overlap import find_ranges_overlap # Function to determine whether two ranges overlap
from .read_gff3_to_df import read_gff3 # Custom function to read gff3




## Return index of genes to be filtered out
def genes_same_region(input_file, list_gene_start, list_gene_end, list_gene_names, list_bitscores):
    """Returns indices of CDS features to filter out based on bitscore.
    If ranges overlap, keep ONLY feature with highest bitscore.
    Ties broken randomly using fixed seed for reproducibility.
    
    Arguments:
    input_file = dataframe # full dataframe of CDS features
    list_gene_start = list(int) # start of each CDS alignment
    list_gene_end = list(int) # end of each CDS alignment  
    list_gene_names = list(str) # name of each CDS
    list_bitscores = list(float) # bitscores of each CDS
    
    Prints reason for filtering each annotation entry.
    """
    list_to_filter = []
    is_keeper = [True] * len(input_file)
    
    for i in range(len(input_file)):
        if not is_keeper[i]:  # Skip already filtered
            continue
            
        range_gene_1 = range(list_gene_start[i], list_gene_end[i]+1)
        bitscore_1 = list_bitscores[i]
        
        for j in range(len(input_file)):
            if i == j or not is_keeper[j]:  # Skip self and filtered
                continue
                
            range_gene_2 = range(list_gene_start[j], list_gene_end[j]+1)
            
            if find_ranges_overlap(range_gene_1, range_gene_2):
                bitscore_2 = list_bitscores[j]
                
                if bitscore_1 > bitscore_2:
                    # i beats j -> filter j
                    list_to_filter.append(j)
                    is_keeper[j] = False
                    print(f"1_Gene_{list_gene_names[j]}_index_{j}_overlaps_with_{list_gene_names[i]}_index_{i}_FILTERED")
                    
                elif bitscore_1 < bitscore_2:
                    # j beats i -> filter i and stop
                    list_to_filter.append(i)
                    is_keeper[i] = False
                    print(f"2_Gene_{list_gene_names[i]}_index_{i}_overlaps_with_{list_gene_names[j]}_index_{j}_FILTERED")
                    break
                    
                else:  # Tie
                    rng = np.random.default_rng(seed=42)
                    winner_idx = rng.choice([i, j], size=1)[0]
                    loser_idx = j if winner_idx == i else i
                    list_to_filter.append(loser_idx)
                    is_keeper[loser_idx] = False
                    winner_name = list_gene_names[winner_idx]
                    loser_name = list_gene_names[loser_idx]
                    print(f"5_Gene_{loser_name}_index_{loser_idx}_tie_with_{winner_name}_index_{winner_idx}")
                    if loser_idx == i:
                        break
    
    unique_genes = np.unique(list_to_filter)
    
    if len(unique_genes) == 0:
        print("A_No_genes_filtered")
    else:
        for i in unique_genes:
            print(f"B_Filtered_gene_{list_gene_names[i]}_index_{i}")
    
    return unique_genes.tolist()

   

## Function to filter genes in the input gff3 file (Filepath_input)
## and save output filtered files (Filepath_output):
def filter_genes_same_map_region (Filepath_input, Filepath_output):
    '''Finds which genes are mapped to the same region and filters the ones
    that have the lower bitscore.
    1. Load dataframe from Filepath_input
    2. Format data
        - Remove nan values
        - Convert columns with gene range (3 and 4) from float to int
    3. Get gene alignment data
        - Gene name
        - Gene start
        - Gene end
        - Gene bitscore
    4. genes_same_region() 
        - Returns genes to be filtered out
    5. Save the dataframes from the genes excluding the filtered genes
    '''
    
    # 1. Load input gff3 file to a dataframe
    input_file = read_gff3(Filepath_input)

    # If the annotation file read is empty, skip this genome
    if input_file is None:
        print (f"No data rows in {Filepath_input}. Skipping this file.")
        return False 

    # 2. Format the gff3 dataframe
    ##### - Fill nan values with 0
    input_file = input_file.fillna(0)
    
    ##### - Convert column 3 and 4 from float to int
    input_file[3] = pd.to_numeric(input_file[3], errors='coerce').astype('Int64')
    input_file[4] = pd.to_numeric(input_file[4], errors='coerce').astype('Int64')

    ##### - Get the name and interval of each gene alignment
    list_gene_names = []
    list_gene_start = []
    list_gene_end = []
    
    # info_column = input_file.iloc[:,-1]
    
    
    for i in range(len(input_file)):
        # Get start and end of the gene alignment
        start_align = int(input_file.iloc[i,3])
        end_align = int(input_file.iloc[i,4])
        
        list_gene_start.append(start_align)
        list_gene_end.append(end_align)
        
        # Get each gene name
        get_last_column = input_file.iloc[i,-1]
        
        match = re.search(r'Target=([^;]+)', get_last_column)
        
        if match:
    	    gene_name = match.group(1)
        else:
    	    gene_name = "unknown"
        
        list_gene_names.append(gene_name)
        
    
    ##### - Get the bitscore of each gene alignment
    list_bitscores = []
   
    for i in range(len(input_file)):
        get_last_column = input_file.iloc[i,-1]
        
        # Use regular expression to find the match_score value
        match = re.search(r'match_score=([\d.]+)', get_last_column)
        
        if match:
            bitscore_value = float(match.group(1))
            list_bitscores.append(bitscore_value)
        else:
            # None if no match_score is found
            list_bitscores.append(None)


    # 5. What genes are mapped to the same region?
    # Call function 'index_same_region' and return gene indexes to be filtered
    index_same_region = genes_same_region(input_file, list_gene_start, list_gene_end, list_gene_names, list_bitscores)
           
    # 6. Drop index of genes to be filtered from the original dataframe
    input_file_filtered = input_file.drop(index_same_region)
     
    # 7. Save filtered dataframe
    input_file_filtered.to_csv(Filepath_output, sep='\t', index=False, header=False)

    return True


# Loop over all the input files
def filter_genewise_one_map_per_region(input_dir, output_dir):
    '''Loop the filter_genes_same_map_region function over all the formatted GeneWise result annotation files for each genome
        - input_dir: Path to formatted-GFF GeneWise result annotation files
        - output_dir: Path to formatted-filtered-GFF GeneWise result annotation files
    '''
    # Loop the function 'filter_genes_same_map_region()' over the input files  
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        print(filename)
        filtered_name = filename
        filtered_name = filtered_name.rsplit('.', 1)[0] 
        filtered_name = filtered_name + '_FILTERED.gff'
    
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, filtered_name)
        
        # Call function to process input files in order to filter genes and save
        # output files:
        if filter_genes_same_map_region (Filepath_input, Filepath_output):
            # The annotation file was correct and filtering was successful
            print(f"Saved {filtered_name}")
        else:
            # The annotation file was empty and it was skipped
            print(f"Skipped {filtered_name}")
        # print empty line for readability
        print()
