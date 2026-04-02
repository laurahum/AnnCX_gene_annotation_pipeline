#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 22:16:37 2024

@author: lahumada
"""

## This script takes as input augustus output files (GFF3 or GTF format) and filters them to keep
## only those genes that have a feature called protein_match, which are those
## that Augustus found with the help of the protein query that it was run with
## as protein profile.

## Main function: filter_augustus_protein_match(input_dir, output_dir)
## Returns:
        # (True) Genes protein_match found and saved to file
        # (False) No genes found -> Skipping this file


import os
import sys
import numpy as np
import pandas as pd
from .read_gff3_to_df import read_gff3

def find_indices_with_protein_match(df_by_gene):
    """
    Returns the indices of DataFrames in df_by_gene that contain a feature called 'protein_match'.
    Parameters:
    - df_by_gene: Dictionary where keys are indices and values are DataFrames.
    Returns:
    - A list of indices corresponding to DataFrames with 'protein_match'.
    """
    indices_with_protein_match = []

    for key, df in df_by_gene.items():
        # Check if 'protein_match' is in the feature column (column 2)
        for row in df.iloc[:,2]:
            if str(row).strip() == 'protein_match':
                indices_with_protein_match.append(key)
                break  # Found protein_match, no need to check rest of this gene

    return np.unique(indices_with_protein_match)

def filter_genes_same_map_region(input_file, output_file):
    '''Keeps only those genes that have a protein_match feature.
    Works with both GFF3 and GTF format Augustus output.
    
    1. Load dataframe from input_file
    2. Format data - handle GFF3/GTF (9 columns)
        - Remove nan values
        - Convert columns with gene range (3 and 4) from float to int
        - Reset index
    3. Find gene positions using 'gene' feature
    4. Split dataframe by gene model
    5. gene_keep() - Returns genes to be kept
    6. Save filtered genes to output_file
    
    Returns:
        (True) Genes protein_match found and saved to file
        (False) No genes found -> Skipping this file
    '''
     
    # 1. Load input gff3/gtf file to a dataframe
    gff3_file = read_gff3(input_file)
   
    # If the annotation file read is empty, skip this genome
    if gff3_file is None:
        print(f"No data rows in {os.path.basename(input_file)}. Skipping this file.")
        return False 

    # 2. Format the gff3/gtf dataframe
    # Drop nan values in rows
    gff3_file = gff3_file.dropna()
    
    ##### - Convert column 3 and 4 from float to int
    gff3_file[3] = pd.to_numeric(gff3_file[3], errors='coerce').astype('Int64') # start
    gff3_file[4] = pd.to_numeric(gff3_file[4], errors='coerce').astype('Int64') # end
    
    # Reset index
    gff3_file = gff3_file.reset_index(drop=True)
    
    # 3. Find gene positions
    # Get the rows that have 'gene' in the feature column (column 2, 0-indexed column 3)
    mask_gene = gff3_file.iloc[:,2] == 'gene'
    gene_rows = gff3_file[mask_gene]
    
    # If no genes found, skip this file
    if len(gene_rows) == 0:
        print(f"No gene features found in {os.path.basename(input_file)}. Skipping this file.")
        return False
    
    pos_gene = gene_rows.index.tolist()

    # 4. Split dataframe by gene alignment (gene + transcript + CDS + protein_match + etc)
    df_by_gene = {}
    
    for i in range(len(pos_gene)):
        # For the last gene in the list
        if pos_gene[i] == pos_gene[-1]:
            df_by_gene[i] = gff3_file.iloc[pos_gene[i]:]
        # For other genes
        else:
            start = pos_gene[i]
            end = pos_gene[i+1]
            df_by_gene[i] = gff3_file.iloc[start:end]

    # 5. Find genes with protein_match features
    index_protein_match = find_indices_with_protein_match(df_by_gene)
    
    print(f"Found {len(index_protein_match)}/{len(df_by_gene)} genes with protein_match features in {os.path.basename(input_file)}")

    
    # 6. Save only genes with protein_match features
    with open(output_file, "w") as new_file:
        for index in index_protein_match:
            df_each = df_by_gene[index].to_csv(header=False, index=False, sep='\t')
            new_file.write(df_each.rstrip())
            new_file.write('\n###\n')
    return True
    

def filter_augustus_protein_match(input_dir, output_dir):
    '''Loop the filter_genes_same_map_region function over all raw Augustus result annotation files.
    
    Args:
        input_dir: Path to raw GFF3/GTF Augustus result annotation files
        output_dir: Path to filtered GFF3 Augustus result annotation files
    
    Returns:
        (True) Genes protein_match found and saved to file
        (False) No genes found -> Skipping this file
    '''
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop over all input files
    for filename in os.listdir(input_dir):
        filtered_name = filename.rsplit('.', 1)[0] + '_FILTERED.gff3'
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, filtered_name)
            
        # Process file      
        result = filter_genes_same_map_region (Filepath_input, Filepath_output)
        if result == True:
            # The annotation file was correct and filtering was successful
            print(f"Saved {filtered_name}")
        else:
            # The annotation file was empty and it was skipped
            print(f"Skipped {filtered_name}")

        # print empty line for readability
        print()
