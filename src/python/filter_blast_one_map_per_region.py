#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 22:50:57 2022

@author: lahumada
"""

'''
This script will filter the BLAST results of aligning query prot/cDNA vs ROI hardmasked FASTA

Since BLAST does not produce gene model annotations, each annotation entry in the BLAST results is considered
for filtering and refered to as 'genes' throughout the script.
BLAST aligns several annotation entries on the same region, with overlapping coordinates and different score values, 
and the objective of this filtering step is to leave only the best annotation entries.

Filter BLAST outfmt 6 results converted to GFF3 using a compound score:

    compound_score = (pident/100 * qcovs/100 * bitscore_norm) / (1 + normalized_gapopen)

where:
    pident            = percentage identity (%)
    qcovs             = query coverage per subject (%)
    bitscore_norm     = bitscore / bitscore_scale
    normalized_gapopen = gapopen / alignment_length

Logic: for overlapping BLAST hits on the same target region, keep the hit with
the highest compound score. If there is a tie, pick one at random (fixed seed).

Main function: filter_blast_one_map_per_region(input_dir, output_dir)
'''


import os
import re
import numpy as np
import pandas as pd
from .ranges_overlap import find_ranges_overlap # Function to determine whether two ranges overlap
from .read_gff3_to_df import read_gff3 # Custom function to read gff3


def calculate_blast_compound_score(attributes,
                                   bitscore_scale=100.0):
    """
    Calculate BLAST compound score from the attributes column of a BLAST GFF3 entry.
    Expected attributes, e.g.: ID=...;Name=...;pident=90.233;length=215;bitscore=281;qcovs=97.5;gapopen=0

    compound_score = (pident * qcovs * bitscore_norm) / (1 + normalized_gapopen)

    Arguments
        - attributes (str): GFF3 attributes field for a BLAST feature.
    	- bitscore_scale (float): Divisor to normalize bitscore (typical bitscores ~100–500; 100 is reasonable).

    Returns (float): Compound score for this BLAST hit.
    """

    # Helper to pull float/int from attributes
    def get_attr_float(key):
        m = re.search(rf"{key}=([^;]+)", attributes)
        return float(m.group(1))

    def get_attr_int(key):
        m = re.search(rf"{key}=([^;]+)", attributes)
        return int(m.group(1))

    pident = get_attr_float("pident")      # % identity
    qcovs = get_attr_float("qcovs")        # % query coverage per subject
    bitscore = get_attr_float("bitscore")  # BLAST bitscore
    gapopen = get_attr_int("gapopen")      # number of gap openings
    alen = get_attr_int("length")      # alignment length

    bitscore_norm = bitscore / bitscore_scale if bitscore_scale > 0 else bitscore
    normalized_gapopen = gapopen / (1 + alen)

    compound_score = (pident/100 * qcovs/100 * bitscore_norm) / (1.0 + normalized_gapopen)

    return compound_score


def genes_same_region(blast_hits, list_start, list_end, list_names, list_scores):
    """
    Determine which BLAST hits to filter out when they overlap on the target.
    Keep ONLY highest compound score per overlap cluster. Ties broken randomly.
    """
    list_to_filter = []
    is_keeper = [True] * len(blast_hits)
    rng = np.random.default_rng(seed=42)
    
    for i in range(len(blast_hits)):
        if not is_keeper[i]:
            continue
            
        range_i = range(list_start[i], list_end[i] + 1)
        score_i = list_scores[i]
        
        for j in range(len(blast_hits)):
            if i == j or not is_keeper[j]:
                continue
                
            range_j = range(list_start[j], list_end[j] + 1)
            
            if find_ranges_overlap(range_i, range_j):
                score_j = list_scores[j]
                
                if score_i > score_j:
                    list_to_filter.append(j)
                    is_keeper[j] = False
                    print(f"BLAST_filter_{list_names[j]}_index_{j}_overlaps_{list_names[i]}_index_{i}_scores_{score_j:.3f}_vs_{score_i:.3f}")
                    
                elif score_i < score_j:
                    list_to_filter.append(i)
                    is_keeper[i] = False
                    print(f"BLAST_filter_{list_names[i]}_index_{i}_overlaps_{list_names[j]}_index_{j}_scores_{score_i:.3f}_vs_{score_j:.3f}")
                    break
                    
                else:  # Tie
                    winner_idx = rng.choice([i, j], size=1)[0]
                    loser_idx = j if winner_idx == i else i
                    list_to_filter.append(loser_idx)
                    is_keeper[loser_idx] = False
                    print(f"BLAST_filter_tie_{list_names[loser_idx]}_index_{loser_idx}_loser_to_{list_names[winner_idx]}_index_{winner_idx}_score_{score_i:.3f}")
                    if loser_idx == i:
                        break
    
    unique_to_filter = np.unique(list_to_filter)
    
    if len(unique_to_filter) == 0:
        print("No overlapping BLAST hits - keeping all")
    else:
        for i in unique_to_filter:
            print(f"A_Filtered_gene_{list_names[i]}_index_{i}")
    
    return unique_to_filter
    

def filter_genes_same_map_region(Filepath_input, Filepath_output):
    """
    Filter a single BLAST GFF3 file by compound score.

    Steps:
      1. Load GFF3 into a DataFrame.
      2. Select BLAST rows (e.g. type == 'BLASTCDS' or 'BLAST').
      3. Compute compound score for each BLAST feature.
      4. For overlapping hits on the same target, keep only the highest-scoring one.
      5. Write filtered GFF3 to Filepath_output.
    """

    gff_df = read_gff3(Filepath_input)

    if gff_df is None or gff_df.empty:
        print(f"No BLAST data in {os.path.basename(Filepath_input)}. Skipping this file.")
        return False

    # Ensure positions numeric
    gff_df = gff_df.dropna()
    gff_df[3] = pd.to_numeric(gff_df[3], errors='coerce').astype('Int64') # start
    gff_df[4] = pd.to_numeric(gff_df[4], errors='coerce').astype('Int64') # end

    gff_df = gff_df.reset_index(drop=True)

    # Filter to BLAST features (adjust types if needed)
    # Here assuming feature types in column 2 are 'BLASTCDS' / 'BLASTN' / 'BLAST'
    blast_mask = gff_df.iloc[:, 2].str.contains("BLAST", case=False, na=False)
    blast_hits = gff_df[blast_mask].reset_index(drop=True)

    if len(blast_hits) == 0:
        print(f"No BLAST features in {os.path.basename(Filepath_input)}. Skipping this file.")
        return False

    list_names = []
    list_start = []
    list_end = []
    list_scores = []

    for idx, row in blast_hits.iterrows():
        start = int(row[3])
        end = int(row[4])
        attributes = str(row[8])

        # Extract ID or Name for reporting
        m_id = re.search(r"ID=([^;]+)", attributes)
        name = m_id.group(1)

        score = calculate_blast_compound_score(attributes)

        list_names.append(name)
        list_start.append(start)
        list_end.append(end)
        list_scores.append(score)

    if len(list_scores) > 0:
        print(f"BLAST compound scores in {os.path.basename(Filepath_input)}: "
              f"min={min(list_scores):.3f} max={max(list_scores):.3f}")
    
        # Get original indices of blast_hits in gff_df
        blast_indices = blast_hits.index.tolist()
    
        # Filter using blast_hits indices (0-M)
        relative_indices_to_filter = genes_same_region(blast_hits, list_start, list_end, list_names, list_scores)
    
        # Map back to original gff_df indices
        index_same_region = [blast_indices[idx] for idx in relative_indices_to_filter]
    else:
        index_same_region = []

    # Drop index of genes to be filtered from the original dataframe
    input_file_filtered = gff_df.drop(index_same_region)
     
    # Save filtered dataframe
    input_file_filtered.to_csv(Filepath_output, sep='\t', index=False, header=False)
    
    return True


def filter_blast_one_map_per_region(input_dir, output_dir):
    """
    Loop over BLAST GFF3 files and apply compound-score-based filtering.
      - input_dir: path to BLAST GFF3 files
      - output_dir: path to filtered BLAST GFF3 files
    """
    for filename in os.listdir(input_dir):
        filtered_name = filename.rsplit('.', 1)[0] + '_FILTERED.gff3'
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, filtered_name)

        if filter_genes_same_map_region(Filepath_input, Filepath_output):
            # The annotation file was correct and filtering was successful
            print(f"Saved {filtered_name}")
        else:
            # The annotation file was empty and it was skipped
            print(f"Skipped {filtered_name}")
        # print empty line for readability
        print()
