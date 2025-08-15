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

## Logic: It two ranges overlap, filter lower bitscore.
##  -> If both ranges have the same bitscore, filter one at random

## Main function: filter_genewise_one_map_per_region(input_dir, output_dir)

import os
import re
import numpy as np
import sys
from .ranges_overlap import find_ranges_overlap # Function to determine whether two ranges overlap
from .read_gff3_to_df import read_gff3 # Custom function to read gff3




## Return index of genes to be filtered out
def genes_same_region (input_file, list_gene_start, list_gene_end, list_gene_names, list_bitscores):
    """Check whether a gene range is contained in another 
    gene range: Compares each gene alignment range to the rest. 
    If the gene range checked is within another gene range 
    the function filters the gene with lower bitscore.
    If both bitscore values are the same, filter
    one gene at random.
    
    Returns index of genes to be filtered out.
    
    Arguments:
    input_file = dataframe # full dataframe
    list_gene_start = list(int) # start of gene alignment
    list_gene_end = list(int) # end of gene alignment
    list_gene_names = list(str) # name of each gene
    list_bitscores = list(int) # bitscores of each gene
    
    Prints the reason for the filtering of each annotation entry
    """
    
    ###### Genes to be filtered
    # Two kind of genes to be filtered out:
    # 1. List for genes that overlap and have smaller bitscore
    list_gene_same_region = []
    
    # 2. List for genes that overlap and have the same bitscore but have lower bitscore
    list_gene_same_bitscore = []
    
    
    ###### Check the ranges:
    for i in range(len(input_file)):
        # Range of gene alignment to be checked
        range_gene_1 = range(list_gene_start[i],list_gene_end[i]+1)
        
        for j in range(len(input_file)):
            
            if i != j: # Prevents comparing one gene to itself
            # Range of gene alignment to be checked
                range_gene_2  = range(list_gene_start[j],list_gene_end[j]+1)
            
                # Check if gene [i] is within gene [j]
                if find_ranges_overlap(range_gene_1, range_gene_2):
                    
                    # Get bitscores for genes [i] and [j]
                    bitscore_gene_1 = list_bitscores[i]
                    bitscore_gene_2 = list_bitscores[j]
                    
                    # Get gene name for genes [i] and [j]
                    gene_name_1 = list_gene_names[i]
                    gene_name_2 = list_gene_names[j]
                
            
                    # 1. List for genes that overlap and have smaller bitscore
                    # If the bitscore of gene [i] is higher, filter gene [j]
                    if bitscore_gene_1 > bitscore_gene_2:
                        list_gene_same_region.append(j)
                        print("1_Gene_" + gene_name_2 + "_index_" + str(j) + "_overlaps_with_" + gene_name_1 + "_index_" + str(i) + "_FILTERED")
                        
                    elif bitscore_gene_1 < bitscore_gene_2:
                        list_gene_same_region.append(i)
                        print("2_Gene_" + gene_name_1 + "_index_" + str(i) + "_overlaps_with_" + gene_name_2 + "_index_" + str(j) + "_FILTERED")
                    
                    # 2. List for genes that overlap because the range is 
                    # exactly the same to other gene but have lower bitscore
                    else:
                        if ((any(i in sublist for sublist in list_gene_same_bitscore)) == False):
                            list_gene_same_bitscore.append([i,j])
                            print("5_Gene_" + gene_name_1 + "_index_" + str(i) + "_same_bitscore_as_" + gene_name_2 + "_index_" + str(j))

                            
                                                      
    ###### Return list of unique genes to be filtered    
     # No overlapping values found
    if (len(list_gene_same_region) == 0) and (len(list_gene_same_bitscore) == 0):
        unique_genes_1 = []
        
        print("A_No_genes_filtered")
        
        return(unique_genes_1)
   
    
    # Ony list_gene_same_region has values
    if (len(list_gene_same_region) > 0) and (len(list_gene_same_bitscore) == 0):
        unique_genes_1 = np.unique(list_gene_same_region)
        
        # Control: print what genes will be filtered
        for i in unique_genes_1:
            index_gene = i
            name_gene = list_gene_names[index_gene]
            print("B_Filtered_gene_" + name_gene + "_index_" + str(index_gene))
        
        return(unique_genes_1)
    
    # Both list_gene_same_region and list_gene_same_bitscore have values
    if (len(list_gene_same_region) > 0) and (len(list_gene_same_bitscore) > 0):
        unique_genes_1 = np.unique(list_gene_same_region)
        
        # Filter random element of each list in list_gene_same_bitscore
        filter_values = []
        for i in range(len(list_gene_same_bitscore)):
            # Create a new random number generator, with a set seed for reproducibility
            rng = np.random.default_rng(seed=42)
            # Randomly select an index
            index_random = rng.choice(list_gene_same_bitscore[i], size=1)[0]
            filter_values.append(index_random)
            
        unique_genes_3 = np.unique(filter_values)
        
        # Merge both lists
        gene_index_to_be_filtered = list(unique_genes_1) + list(unique_genes_3)
        unique_genes_filter = np.unique(gene_index_to_be_filtered)
        
        # Control: print what genes will be filtered
        for i in unique_genes_filter:
            index_gene = i
            name_gene = list_gene_names[index_gene]
            print("C_Filtered_gene_" + name_gene + "_index_" + str(index_gene))
        
        return(unique_genes_filter)
    
   

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
    input_file.iloc[:,3] = input_file.iloc[:,3].astype(int)
    input_file.iloc[:,4] = input_file.iloc[:,4].astype(int)

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
        gene_name = input_file.iloc[i,0]
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
