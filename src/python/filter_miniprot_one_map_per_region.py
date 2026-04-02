#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 02:59:46 2026

@author: lahumada
"""

"""
This script will filter Miniprot gene model results.
Since Miniprot produces gene model annotations (mRNA → CDS → stop_codon), 
each gene model (labeled as 'mRNA' in the GFF3 files) is considered for filtering and referred 
to as 'genes' (which includes its associated features: CDS, stop_codon) throughout the script.

Miniprot aligns several gene model annotation entries on the same region, with overlapping 
coordinates and different score values, and the objective of this filtering step is to leave 
only the best gene model annotation entries using compound score: 
(mapq/60 × ms_norm × query_cov) / (1 + gaps_penalty)

where:
- ms_norm     = ms / (1 + qlen)            # Matching score normalized by query length
- query_cov   = nmatch / (1 + qlen)        # Query coverage  (match bases / query length)
- gaps_penalty= (alen - nmatch) / 100 # Normalized gaps/indels  (alen = align length)

Logic: If gene models overlap, keep highest compound score. Ties broken randomly.

Main function: filter_miniprot_one_map_per_region(input_dir, output_dir)
"""

import os
import sys
import numpy as np
import pandas as pd
import re
from .ranges_overlap import find_ranges_overlap # Function to determine whether two ranges overlap
from .read_gff3_to_df import read_gff3 # Custom function to read gff3


def genes_same_region(gene_rows, list_gene_start, list_gene_end, list_gene_names, 
                     list_compound_scores, pos_gene):
    """Returns indices of gene models to filter out based on compound score.
    If gene ranges overlap, keep ONLY model with highest compound score.
    Ties broken randomly using fixed seed for reproducibility.
    Arguments:
        gene_rows = dataframe # only mRNA rows from input_file
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


def calculate_miniprot_compound_score(query_len, match_bases, align_len, mapq, ms):
    """Calculate compound score: (mapq/60 × ms_norm × query_cov) / (1 + gaps_penalty)
    where:
    ms_norm     = ms / (1 + qlen)        # Matching score normalized by query length
    query_cov   = nmatch / (1 + qlen )   # Query coverage  
    gaps_penalty= (alen - nmatch) / 100  # Normalized gaps/indels
    """
    ms_norm = ms / (1 + query_len)
    query_cov = match_bases / (1 + query_len)
    gaps_penalty = (align_len - match_bases) / 100
    compound_score = (mapq/60 * ms_norm * query_cov) / (1 + gaps_penalty)
    return compound_score


def filter_genes_same_map_region(Filepath_input, Filepath_output):
    '''Filters Miniprot gene models by compound score when ranges overlap.
    1. Load and format GFF3
    2. Extract gene model data and calculate compound scores
    3. Split by gene model (mRNA+CDS+stop_codon)
    4. Filter overlapping models by compound score
    5. Save non-overlapping gene models
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

    # 3. Get gene model data (mRNA features for Miniprot)
    mask_gene = input_file.iloc[:,2] == 'mRNA'
    gene_rows = input_file[mask_gene]
    pos_gene = gene_rows.index.tolist()
    
    if len(gene_rows) == 0:
        print(f"No mRNA features in {os.path.basename(Filepath_input)}. Skipping.")
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
        
        # Extract gene name from ID attribute
        attributes = str(row[8])
        id_match = re.search(r'ID=([^;]+)', attributes)
        gene_name = id_match.group(1)
        list_gene_names.append(gene_name)
        
        # Extract compound score components from mRNA attributes
        query_len_match = re.search(r'qlen=([^;]+)', attributes)
        query_len = int(query_len_match.group(1))
        
        match_bases_match = re.search(r'nmatch=([^;]+)', attributes)
        match_bases = int(match_bases_match.group(1))
        
        align_len_match = re.search(r'alen=([^;]+)', attributes)
        align_len = int(align_len_match.group(1))
        
        mapq_match = re.search(r'mapq=([^;]+)', attributes)
        mapq = int(mapq_match.group(1))
        
        ms_match = re.search(r'ms=([^;]+)', attributes)
        ms = int(ms_match.group(1))
        
        # Calculate compound score
        compound_score = calculate_miniprot_compound_score(
            query_len, match_bases, align_len, mapq, ms
        )
        list_compound_scores.append(compound_score)

    print(f"Calculated compound scores: {list_compound_scores}")

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


def filter_miniprot_one_map_per_region(input_dir, output_dir):
    '''Loop over raw Miniprot GFF3 files and apply compound score filtering'''
    for filename in os.listdir(input_dir):
        filtered_name = filename.rsplit('.', 1)[0] + '_FILTERED.gff'
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, filtered_name)
        
        if filter_genes_same_map_region(Filepath_input, Filepath_output):
            print(f"Saved {filtered_name}")
        else:
            print(f"Skipped {filtered_name}")
        print()
