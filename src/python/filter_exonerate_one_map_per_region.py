#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 22:21:19 2022

@author: lahumada
"""

'''
This script will filter the Exonerate results of aligning query cDNA vs ROI hardmasked FASTA

Since Exonerate produces gene model annotations, each gene model (labeled as 'gene' in the GFF files) is considered
for filtering and refered to as 'genes' (which includes its associated features: exon, cds, intron, utr3, utr5, splice3, splice4, similarity)) 
throughout the script.
Exonerate aligns several gene model annotation entries on the same region, with overlapping coordinates and different score values, 
and the objective of this filtering step is to leave only the best gene model annotation entries.

Uses compound score: (score_norm × identity/100 × similarity/100) / (1 + normalized_avg_exon_indels)
where:
- score_norm = score internal exonerate / 10000.0
- normalized_avg_exon_indels = avg(indels_per_exon) / (1 + alignment_length)
- indels_per_exon = insertions + deletions

Logic: If gene models overlap, keep highest compound score. Ties broken randomly.

Main function: filter_exonerate_one_map_per_region(input_dir, output_dir)
'''


import os
import numpy as np
import pandas as pd
import re
from .ranges_overlap import find_ranges_overlap # Function to determine whether two ranges overlap
from .read_gff3_to_df import read_gff3 # Custom function to read gff3



# Function to find what genes are mapped to the same region
# Check whether a gene range is contained in another gene range
# Returns genes to be filtered out

def genes_same_region(gene_rows, list_gene_start, list_gene_end, list_gene_names, 
                     list_compound_scores, pos_gene):
    """Returns indices of gene models to filter out based on compound score.
    If gene ranges overlap, keep ONLY model with highest compound score.
    Ties broken randomly using fixed seed for reproducibility.
        Arguments:
        gene_rows = dataframe # only gene rows from input_file
        list_gene_start = list(int) # start of gene alignment
        list_gene_end = list(int) # end of gene alignment
        list_gene_names = list(str) # name of each gene
        list_compound_scores = list(float) # compound scores of each gene
        pos_gene = list(int) # index of genes in gene_rows
        
    Prints the reason for filtering each annotation entry.
    """    
    list_gene_same_region = []
    is_keeper = [True] * len(gene_rows)  # Track which genes are keepers
    
    for i in range(len(gene_rows)):
        if not is_keeper[i]:  # Skip if already marked for filtering
            continue
            
        range_gene_1 = range(list_gene_start[i], list_gene_end[i]+1)
        score_gene_1 = list_compound_scores[i]
        
        for j in range(len(gene_rows)):
            if i == j or not is_keeper[j]:  # Skip self and already-filtered genes
                continue
                
            range_gene_2 = range(list_gene_start[j], list_gene_end[j]+1)
            
            if find_ranges_overlap(range_gene_1, range_gene_2):
                score_gene_2 = list_compound_scores[j]
                
                if score_gene_1 > score_gene_2:
                    # Gene i beats gene j -> filter j
                    gene_pos_j = pos_gene[j]
                    gene_name_j = list_gene_names[j]
                    list_gene_same_region.append(gene_pos_j)
                    is_keeper[j] = False
                    print(f"Filter_gene_{gene_name_j}_index_{gene_pos_j}_score_{score_gene_2:.2f}_beaten_by_{list_gene_names[i]}")
                    
                elif score_gene_1 < score_gene_2:
                    # Gene j beats gene i -> filter i, break out
                    gene_pos_i = pos_gene[i]
                    gene_name_i = list_gene_names[i]
                    list_gene_same_region.append(gene_pos_i)
                    is_keeper[i] = False
                    print(f"Filter_gene_{gene_name_i}_index_{gene_pos_i}_score_{score_gene_1:.2f}_beaten_by_{list_gene_names[j]}")
                    break  # i is filtered, no need to check more opponents
                    
                else:  # Tie
                    rng = np.random.default_rng(seed=42)
                    winner_idx = rng.choice([i, j], size=1)[0]
                    loser_idx = j if winner_idx == i else i
                    gene_pos_loser = pos_gene[loser_idx]
                    gene_name_loser = list_gene_names[loser_idx]
                    list_gene_same_region.append(gene_pos_loser)
                    is_keeper[loser_idx] = False
                    print(f"Filter_gene_{gene_name_loser}_index_{gene_pos_loser}_tie_score_{score_gene_1:.2f}_loser_to_{list_gene_names[winner_idx]}")
                    if loser_idx == i:
                        break  # i lost tie, no need to check more opponents
    
    unique_genes = np.unique(list_gene_same_region)
    
    for i in unique_genes:
        index_gene = list(pos_gene).index(i)
        name_gene = list_gene_names[index_gene]
        print(f"Filtered_gene_{name_gene}_index_{i}")
    
    if len(unique_genes) == 0:
        print("No overlaps found - keeping all gene models")
    
    return unique_genes

    

def calculate_exonerate_compound_score(gene_row, exon_rows):
    """Exonerate compound score: (score_norm × identity × similarity) / (1 + normalized_avg_exon_indels)"""
    
    # Parse gene-level metrics
    score_norm = float(gene_row[5]) / 10000.0
    
    # Identity/similarity from attributes (column 9, 0-indexed = 8)
    attrs = str(gene_row[8])
    identity_match = re.search(r'identity\s+(\d+(?:\.\d+)?)', attrs)
    similarity_match = re.search(r'similarity\s+(\d+(?:\.\d+)?)', attrs)
    
    identity = float(identity_match.group(1))
    similarity = float(similarity_match.group(1))
    
    # Calculate avg exon indels (same logic)
    exon_indels = []
    for _, exon_row in exon_rows.iterrows():
        exon_attrs = str(exon_row[8])
        ins_match = re.search(r'insertions\s+(\d+)', exon_attrs)
        del_match = re.search(r'deletions\s+(\d+)', exon_attrs)
        
        insertions = float(ins_match.group(1))
        deletions = float(del_match.group(1))
        indels = insertions + deletions
        exon_indels.append(indels)
    
    avg_exon_indels = np.mean(exon_indels)
    alen = gene_row[4] - gene_row[3] + 1
    normalized_avg_exon_indels = avg_exon_indels / (1 + alen)
    
    compound_score = (score_norm * identity/100 * similarity/100) / (1 + normalized_avg_exon_indels)
    return compound_score


## Function to filter genes in the input gff3 file (Filepath_input)
## and save output filtered files (Filepath_output):
def filter_genes_same_map_region (Filepath_input, Filepath_output):
    '''Finds which genes are mapped to the same region and filters the ones
    that have the lower compound score.
    1. Load dataframe from Filepath_input
    2. Format data
        - Remove nan values
        - Convert columns with gene range (3 and 4) from float to int
    3. Get gene alignment data
        - Gene name
        - Gene start
        - Gene end
        - Gene compound score
    4. Split dataframe annotation entries by gene model (gene + exon + cds + intron + utr3 + utr5 + splice3 + splice4 + similarity)
    5. gene_same_region() 
        - Returns genes to be filtered out
    6. Save the dataframes from the genes excluding the filtered genes
    '''
    
    # 1. Load input GFF3
    input_file = read_gff3(Filepath_input)
    
    if input_file is None:
        print(f"No data rows in {Filepath_input}. Skipping this file.")
        return False

    # 2. Format data
    input_file = input_file.dropna()
    input_file[3] = pd.to_numeric(input_file[3], errors='coerce').astype('Int64')
    input_file[4] = pd.to_numeric(input_file[4], errors='coerce').astype('Int64')
    input_file = input_file.reset_index(drop=True)

    # 3. Get gene model data
    mask_gene = input_file.iloc[:,2] == 'gene'
    gene_rows = input_file[mask_gene]
    pos_gene = gene_rows.index.tolist()
    
    if len(gene_rows) == 0:
        print(f"No gene features in {os.path.basename(Filepath_input)}. Skipping.")
        return False

    # Extract gene names, ranges, and calculate compound scores
    list_gene_names = []
    list_gene_start = []
    list_gene_end = []
    list_compound_scores = []
    
    for i in range(len(gene_rows)):
        row = gene_rows.iloc[i]
        list_gene_start.append(int(row[3]))
        list_gene_end.append(int(row[4]))
        
        # Extract gene name from attributes (sequence name)
        attributes = str(row[8])
        name_match = re.search(r'sequence\s+[=|\s]\s*([^;\s]+)', attributes)
        gene_name = name_match.group(1) if name_match else f"gene_{i}"
        list_gene_names.append(gene_name)
        
        # Find corresponding exons for this gene model (slice by position)
        start_idx = pos_gene[i]
        end_idx = pos_gene[i+1] if i+1 < len(pos_gene) else len(input_file)
        gene_slice = input_file.iloc[start_idx:end_idx]
        exon_rows = gene_slice[gene_slice.iloc[:,2] == 'exon']
        
        # Calculate compound score
        compound_score = calculate_exonerate_compound_score(row, exon_rows)
        list_compound_scores.append(compound_score)

    print(f"Exonerate compound scores: min={min(list_compound_scores):.2f} max={max(list_compound_scores):.2f}")

    # 4. Split dataframe by gene model
    df_by_gene = {}
    for i in range(len(pos_gene)):
        if pos_gene[i] == pos_gene[-1]:
            df_by_gene[i] = input_file.iloc[pos_gene[i]:]
        else:
            df_by_gene[i] = input_file.iloc[pos_gene[i]:pos_gene[i+1]]

    # 5. Filter overlapping gene models
    index_same_region = genes_same_region(gene_rows, list_gene_start, list_gene_end, 
                                        list_gene_names, list_compound_scores, pos_gene)
    # 6. Save filtered gene models
    with open(Filepath_output, "w") as new_file:
        for i in range(len(df_by_gene)):
            index_df = df_by_gene[i].index
            if not any(item in index_same_region for item in index_df):
                df_each = df_by_gene[i].to_csv(header=False, index=False, sep='\t')
                new_file.write(df_each)
                new_file.write('###\n')
                
    return True


# Loop over all the input files
def filter_exonerate_one_map_per_region(input_dir, output_dir):
    '''Loop the filter_genes_same_map_region function over all the formatted Exonerate result annotation files for each genome
        - input_dir: Path to formatted-GFF Exonerate result annotation files
        - output_dir: Path to formatted-filtered-GFF Exonerate result annotation files
    '''
    ##### Loop the function 'filter_genes_same_map_region()' over the input files  
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
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
