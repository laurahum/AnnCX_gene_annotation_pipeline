#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 19:37:40 2024

@author: lahumada
"""

## This script is the third step (3/3) to filter EVM results
## 1. Take as input the EVM raw results and the combined_overlaps file generated 
## by 'filter_EVM_results_find_overlaps.sh'
## 2. Get the genes with a certain amount of overlaps from combined_overlaps file
## and then use their ids to grab the rest of the features from the EVM raw results
## 3. Write the filtered genes to a gff3

## Main function: filter_EVM_results_overlap_number_tools(input_dir_EVM_raw, input_dir_overlaps, ouput_dir, genome_list, overlaps_number)
## overlaps number set to default 4 in the main script and can be changed by the user

import re
import os
from collections import defaultdict


# Function to extract gene ID from the last field
# Entry is a list of strings that each correspond to a field of the annotation
def get_gene_id(entry):
    ''' 
    Extract gene ID from the attributes field (last column) of an annotation entry
    using regular expressions. The gene ID is a sequence of digits at the end of 
    the attributes string, following a period.
        - entry (list): List of strings, each corresponding to a field of an annotation entry
    
    Returns (str): Extracted gene ID
    '''
    attributes = entry[-1]
    match = re.search(r'\.(\d+)$', attributes)
    return match.group(1)


# Read gff file to list of entries
def read_gff_to_list_entries (input_file):
    ''' 
    Read GFF file line by line and return a list of annotation entries.
    Ignores lines that do not correspond to annotation entries (9 fields)
        input_file (str): Path to input GFF file
        
    Returns (list): List of lists where each inner list corresponds to an annotation entry
    '''
    list_of_entries = []

    with open(input_file, 'r') as infile:
    
        # Read line by line
        for line in infile:
            parts = line.strip().split('\t')
            print(parts)
        
            # Only lines that correspond to an annotation
            if len(parts) == 9:
                list_of_entries.append(parts)
       
    return list_of_entries
                

# Get the gene ids of the genes to keep
# Use the global variable for the number of overlaps so it's easy to change
def get_genes_id_to_keep (genes_input_file, overlaps_number):
    '''
    Identifies gene IDs to keep based on a specified overlap threshold.
    It reads the input file, extracts the overlap count for each gene and keeps
    the gene IDs where the overlap count is greater than or equal to the
    overlaps_number specified by the user.
        - gene_input_file (str): Path to the input TXT file (*_combined_overlaps.txt) 
        which is result from script filter_EVM_results_find_overlaps.sh in the second filtering step (2/3) 
        - overlaps_number (int): Number of tools to overlap for the filtering (Default = 4)
        
    Returns (list): List of integer gene IDs that meet the overlap threshold criteria
    '''
    list_of_entries = read_gff_to_list_entries(genes_input_file)
    
    # list of genes id to keep
    genes_to_keep = []
    
    for entry in list_of_entries:
        
        # Get the number of overlaps
        overlap_count = entry[0].strip()[0]
        print(overlap_count)
            
        # If the overlap count was 3, gene ids to keep
        if int(overlap_count) >= int(overlaps_number):
            gene_id = get_gene_id(entry)
            genes_to_keep.append(int(gene_id))
    
    return genes_to_keep


# Get the annotation lines from the EVM results corresponding to the gene ids to keep
def get_annotations_entry_to_keep (EVM_input_file, genes_input_file, overlaps_number):
    ''' 
    Get the annotation lines, for gene models, from the EVM raw results 
    corresponding to the gene IDs to keep. 
    It reads the EVM raw results, gets the list of gene IDs to keep and filters
    the EVM entries to only include those whose gene IDs are in the 'keep' list.
        - EVM_input_file (str): Path to the GFF3 EVM raw results file
        - genes_input_file (str): Path to the input TXT file (*_combined_overlaps.txt) 
        which is result from script filter_EVM_results_find_overlaps.sh in the second filtering step (2/3) 
        - overlaps_number (int): Number of tools to overlap for the filtering (Default = 4)
        
    Returns (list): List of annotation entries (each a list of strings) that correspond to gene models to keep
    '''
    list_of_entries_raw = read_gff_to_list_entries(EVM_input_file)
    gene_ids_keep = get_genes_id_to_keep (genes_input_file, overlaps_number)
    
    # list of annotation entries to keep
    list_of_entries_filtered = []
    
    for entry in list_of_entries_raw:
        entry_id = get_gene_id(entry)
        
        if int(entry_id) in gene_ids_keep:
            list_of_entries_filtered.append(entry)
    
    return list_of_entries_filtered


# Order the list of entries by annotation entry for gene -> mRNA -> exon -> CDS
def sort_gff_entries(list_of_entries):
    '''  
    Sorts GFF3 entries, for each gene model, as gene -> mRNA -> exon -> CDS
    It groups entries by gene ID, then sorts them based on the order of features specified
    before
        - list_of_entries (list): List of GFF3 entries where each entry is a list of strings
        
    Returns (list): Sorted list of GFF3 entries
    '''
    # Create a dictionary to group entries by gene ID
    gene_groups = defaultdict(lambda: defaultdict(list))

    # Group entries by gene and feature type
    for entry in list_of_entries:
        gene_id = get_gene_id(entry)
        feature_type = entry[2]
        gene_groups[gene_id][feature_type].append(entry)

    # Define the order of feature types
    feature_order = {'gene': 0, 'mRNA': 1, 'exon': 2, 'CDS': 3}

    # Sort and flatten the grouped entries
    sorted_entries = []
    for gene_id in sorted(gene_groups.keys()):
        gene_entries = gene_groups[gene_id]
        for feature_type in sorted(gene_entries.keys(), key=lambda x: feature_order.get(x, 4)):
            sorted_entries.extend(sorted(gene_entries[feature_type], key=lambda x: int(x[3])))  # Sort by start position

    return sorted_entries


# Write sorted list_of entries to a file:
def write_gff(list_of_entries, output_file, version):
    """
    Writes a list of GFF entries to a file.
        - list_of_entries (list): A list of GFF entries, where each entry is a list of strings.
        - output_file (str): The path to the output file.
    """   
    with open(output_file, 'w') as outfile:
        outfile.write (f'## Generated by AnnCX_gene_annotation_pipeline {version}' + '\n')
        
        for entry in list_of_entries:
            outfile.write('\t'.join(entry) + '\n')
        

# Processing of EVM annotation results
def process_EVM (EVM_input_file, gene_input_file, output_file, overlaps_number, version): 
    """
    Processes EVM annotation results:
    1. Filter entries: get_annotation_entry_to_keep
    2. Sort entries: sort_gff_entries
    3. Write filtered sorted entries: write_gff

        - EVM_input_file (str): Path to the EVM raw results file.
        - gene_input_file (str): Path to the file containing overlap info files for what genes to keep (*_combined_overlaps.txt) 
        - output_file (str): Path to the output EVM filtered results file.
        - overlaps_number (int): Number of tools to overlap for the filtering (Default = 4)
    """
    list_all_genome = get_annotations_entry_to_keep(EVM_input_file, gene_input_file, overlaps_number)
    sorted_list_all_genome = sort_gff_entries(list_all_genome)
    write_gff(sorted_list_all_genome, output_file, version)
    

# Find the input files that correspond to each genome in the list 
def find_input_files(genome_list, EVM_raw_folder, overlaps_folder):
    """
    Finds input files, EVM raw results files and overlap files, for each genome in the given list.
        - genome_list (list): A list of genome names.
        - EVM_raw_folder (str): Path to the folder containing EVM raw files.
        - overlaps_folder (str): Path to the folder containing overlap info files for what genes to keep (*_combined_overlaps.txt) 

    Returns (dict): Dictionary where keys are genome names and values are dictionaries
          containing paths to the corresponding EVM and overlaps files.
    """
    genome_files = {}
    
    for genome in genome_list:
        EVM_file = None
        overlaps_file = None
        
        # Search for reference file
        for filename in os.listdir(EVM_raw_folder):
            if genome in filename and filename.endswith('.gff3'):
                EVM_file = os.path.join(EVM_raw_folder, filename)
                break
        
        # Search for query file
        for filename in os.listdir(overlaps_folder):
            if genome in filename and filename.endswith('_combined_overlaps.txt'):
                overlaps_file = os.path.join(overlaps_folder, filename)
                break
        
        genome_files[genome] = {'EVM': EVM_file, 'overlaps': overlaps_file}
    
    return genome_files
    

# Main
def filter_EVM_results_overlap_number_tools(input_dir_EVM_raw, input_dir_overlaps, output_dir, genome_list, overlaps_number, version):
    """
    Processes EVM results for multiple genome based on overlap criteria.
    1. Reads a list of genome
    2. Finds the corresponding EVM raw input files,
    3. Filters the EVM results based on 'overlaps_number' criterion: process_EVM
    Keeps only those gene model annotation entries in EVM raw results that are found 
    by as many individual annotation tools as the overlaps_number identifies
    4. Writes the filtered results to output files.

        - input_dir_EVM_raw (str): Path to the directory containing raw EVM files.
        - input_dir_overlaps (str): Path to the directory containing overlap info files for what genes to keep (*_combined_overlaps.txt) 
        - genome_list (str): Path to a file containing a list of genome names.
        - overlaps_number (int): Number of tools to overlap for the filtering (Default = 4)
    """
    # Open the genome list
    with open(genome_list, 'r') as file:
        genome_list = file.read().splitlines()
        
    genome_dictionary = find_input_files(genome_list, input_dir_EVM_raw, input_dir_overlaps)


    # Loop the main function over the input files
    for genome, files in genome_dictionary.items():
        print(f"Filter_EVM_results_overlap_number_tools for: {genome}")  
        output_name = f"{genome}_consensus_filter_{str(overlaps_number)}.gff3"
        output_file = os.path.join(output_dir, output_name)
        
        process_EVM(files['EVM'], files['overlaps'], output_file, overlaps_number, version)
        print("Saved_" + output_name + '\n')
        
        
