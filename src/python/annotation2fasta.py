#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 23:02:30 2022

@author: lahumada
"""


## This script will read the annotation for each genome in a TXT list, 
## and extract the FASTA sequences for each gene annotated using pyfastx
##  for the following annotations:
    
## 1. gene
## 2. mRNA
## 3. CDS
##    3.1. each CDS sequence
##    3.2. all CDS joined: CDS_all ~ (cDNA)
## 4. exon
##    4.1. each exon sequence
##    4.2. all exon joined: exon_all
## 5. intron
##    5.1. each intron sequence inferred from the exon coordinates
## 6. protein
##    6.1. translated CDS_all

## Then the nucleotide or aminoacid sequences will be saved as FASTA files.
## One FASTA file per kind of annotation, each containing all the genome in
## the list.

## Main function: get_fasta_sequences_annotation(genome_file, gene_to_annotate, input_dir_ROI, input_dir_annotation, output_dir, genes_out = None)


import pandas as pd
import pyfastx
import os
import re
from Bio.Seq import Seq
#from .read_gff3_to_df import read_gff3 # Custom function to read gff3

def read_gff3(file_path):
    '''Function to read the gff3 file and return a pandas dataframe'''
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

# 2. List of genome, genes, types and levels
### List of types
types_annotation = ['gene', 'mRNA', 'CDS', 'exon', 'intron']
### List of levels
levels_annotation = ['start', 'end', 'strand']

### List of levels for List_of_fasta
# Make CDS_all and exon_all by concatenating strings once I get the nucleotide sequence
types_fasta = ['gene', 'mRNA', 'CDS', 'CDS_all', 'exon', 'exon_all', 'intron', 'protein'] 



# Make list of genome from a file
def get_list_genome(genome_file):
    """
    Read a file containing genome names and return a list of genome.
        - genome_file (str): Path to the file containing genome names, one per line.

    Returns (list): A list of genome names, with newline characters removed.
    """
    genome_txt = open(genome_file, 'r')
    genome_lines = genome_txt.readlines()
    genome_list = []

    for i in range(len(genome_lines)):
        genome_list.append(genome_lines[i].rstrip('\n'))
    
    return (genome_list)
    

# Make a list of dictionaries with the output directories
def create_output_dir(genome_name, gene_to_annotate, output_dir):
    """
    Create a list of dictionaries containing output file paths for each genome and FASTA type (from types_fasta).
        - genome_name (str): Name of the genome to create output paths for
        - gene_to_annotate (str): Name of the gene to be annotated (Example: "NKG2")
        - output_dir (str): Path to the directory where output files will be saved.

    Returns (list): List of dictionaries, each containing a FASTA type as key and its 
                    corresponding output file path as value.
                    
    Notes:
    - The function uses the global variable:
        - types_fasta (list of str): List of FASTA types to generate.
    """    
    ## Dictionary of output directories
    List_of_output_dir = list()

    for each_type_fasta in types_fasta:    
        file_name = [genome_name, gene_to_annotate, each_type_fasta, 'sequence.fasta']
        file_name = '_'.join(map(str, file_name))
    
        each_dict = {each_type_fasta:[os.path.join(output_dir, file_name)]}
        
        List_of_output_dir.append(each_dict)
    
    return (List_of_output_dir)


# Functions to populate data structures
#1. Populate List_of_annotation
def populate_annotation (df_by_gene, genes, List_of_annotation):
    """
    Populate the List_of_annotation structure with annotation information from GFF3 files.
    1. Extracts gene, mRNA, exon, and CDS information from each DataFrame in df_by_gene.
    2. Populates List_of_annotation with start, end, and strand information for each feature type.
    3. Calculates and adds intron information inferred from exon positions.

        - df_by_gene (list): List of pandas DataFrames, each containing annotation data for a gene model.
        - genes (list): List of gene names.
        - List_of_annotation (list): List of dictionaries to be populated with annotation data.
        
    Note:
    - It uses regular expressions to extract gene names from the 'attributes' field of the GFF3 data.
    """   
    ##### Populate the dictionaries in the data structure above:
    ## Fill in gene, mRNA, exon and CDS    
    for i in range(len(df_by_gene)):
        # Reset index of each dataframe
        df = df_by_gene[i].reset_index()
        original_index = df.loc[:,'index']
        df = df.drop(columns=['index'])
  
        for j in range(len(df)):
            # Make each row into a string and into a list
            each_row = df.iloc[j,:].to_string(header=False, index=False).split()
            # Create a variable for each element in the list
            scaffold,source,kind,start,end,score,strand,phase,attributes = each_row
            print(attributes)
        
            List_of_annotation[i][kind]['start'].append(int(start))
            List_of_annotation[i][kind]['end'].append(int(end))
            List_of_annotation[i][kind]['strand'].append(strand)  
            print(f"Extracted {genes[i]} for the {kind} annotation, start = {start}, end = {end}, strand = {strand}")
            
    # Fill introns
    for dic in range(len(List_of_annotation)):
        start_exon = sorted(List_of_annotation[dic]['exon']['start'])
        end_exon = sorted(List_of_annotation[dic]['exon']['end'])
        strand_exon = List_of_annotation[dic]['exon']['strand'][0]
    
        for i in range(len(end_exon)-1):
            List_of_annotation[dic]['intron']['start'].append(end_exon[i]+1)
            List_of_annotation[dic]['intron']['end'].append(start_exon[i+1]-1)
            List_of_annotation[dic]['intron']['strand'].append(strand_exon)
            print(f"Extracted {genes[dic]} for the intron annotation, start = {end_exon[i]+1}, end = {start_exon[i+1]-1}, strand = {strand_exon}")
               

# 2. Populate List_of_fasta
# Populate List_of_fasta using pyfastx      
def populate_fasta (List_of_annotation, List_of_fasta, ROI_fasta, ROI_entry_name, genes):
    """
    Populate the List_of_fasta structure with FASTA sequences (for various genetic elements: 
    CDS, CDS_all, exon, exon_all, intron, gene, mRNA) based on annotation data.
    The function handles strand-specific extractions and performs error checking for missing data.
    1. For CDS and exon:
       - Extracts individual sequences for each element -> CDS, exon
       - Combines all elements into a single sequence -> CDS_all, exon_all
       - Translates CDS_all -> protein
    2. For introns: Extracts individual sequences for each intron.
    3. For gene and mRNA: Extracts a single sequence for the entire element.

        - List_of_annotation (list): List of dictionaries containing annotation data for each gene.
        - List_of_fasta (list): List of dictionaries to be populated with FASTA sequences.
        - ROI_fasta (pysam.FastaFile): FASTA file object for the region of interest.
        - ROI_entry_name (str): Name of the entry sequence in the ROI FASTA file
        - genes (list): List of gene names.

    Note:
    - The function includes special handling for finding start codons in CDS sequences.
    - The function uses the global variable types_annotation (list): List of annotation types to process
    """
    for each_gene in range(len(genes)):
        for each_type in types_annotation:
        
            if each_type == 'CDS' or each_type == 'exon':
                start_list = List_of_annotation[each_gene][each_type]['start']
                end_list = List_of_annotation[each_gene][each_type]['end']
                strand_list = List_of_annotation[each_gene][each_type]['strand']
                
                if not start_list or not end_list or not strand_list:
                    print(f"Warning: Empty data for {genes[each_gene]} {each_type}. Skipping.")
                    continue
            
                each_type_tuple_list = []
            
                for i in range(len(start_list)):
                    start = start_list[i]
                    end = end_list[i]
                    strand = strand_list[i]
                    
                    fasta_seq = ROI_fasta.fetch(ROI_entry_name, (start,end), strand = strand)
                    
                    List_of_fasta[each_gene][each_type].append(fasta_seq)
                    print(f"Extracted {genes[each_gene]} for the fasta sequence of {each_type}_{i}")
            
                    each_type_tuple_list.append((start,end))

                each_type_tuple_list = sorted(each_type_tuple_list)

                fasta_seq = ROI_fasta.fetch(ROI_entry_name, each_type_tuple_list, strand = strand_list[0])   
            
                if each_type == 'CDS':
                    List_of_fasta[each_gene]['CDS_all'].append(fasta_seq)
                    print(f"Extracted {genes[each_gene]} for the fasta sequence of all {each_type}")
                    
                    # Translate to protein 
                    # Find the first occurrence of either ATG or GTG
                    start_codon_atg = fasta_seq.lower().find("atg")
                    start_codon_gtg = fasta_seq.lower().find("gtg")

                    start_codon = None
                    # Determine the first start codon position
                    if start_codon_atg == -1 and start_codon_gtg == -1:
                        print("No start codon found.")
                    elif start_codon_atg == -1:  # Only GTG found
                        start_codon = start_codon_gtg
                    elif start_codon_gtg == -1:  # Only ATG found
                        start_codon = start_codon_atg
                    else:  # Both found, take the first one
                        start_codon = min(start_codon_atg, start_codon_gtg)

                    # Slice the sequence from the identified start codon
                    sliced_sequence = fasta_seq[start_codon:]
                    # Translate
                    dna_seq = Seq(sliced_sequence)
                    prot_seq = dna_seq.translate(to_stop=True)
                    prot_seq = str(prot_seq)
                    List_of_fasta[each_gene]['protein'].append(prot_seq)
                    print(f"Translated {genes[each_gene]} for the fasta sequence of protein")
            
                if each_type == 'exon':
                    List_of_fasta[each_gene]['exon_all'].append(fasta_seq)
                    print(f"Extracted {genes[each_gene]} for the fasta sequence of all {each_type}")
                    
                    
            elif each_type == 'intron':
                start_list = List_of_annotation[each_gene][each_type]['start']
                end_list = List_of_annotation[each_gene][each_type]['end']
                strand_list = List_of_annotation[each_gene][each_type]['strand']
            
                if not start_list or not end_list or not strand_list:
                    print(f"Warning: Empty data for {genes[each_gene]} {each_type}. Skipping.")
                    continue
            
                for i in range(len(start_list)):
                    start = start_list[i]
                    end = end_list[i]
                    strand = strand_list[i]
               
                    fasta_seq = ROI_fasta.fetch(ROI_entry_name, (start,end), strand = strand)
                   
                    List_of_fasta[each_gene][each_type].append(fasta_seq)
                    print(f"Extracted {genes[each_gene]} for the fasta sequence of {each_type}_{i}")
                
            elif each_type == 'gene' or each_type == 'mRNA':
                start = List_of_annotation[each_gene][each_type]['start'][0]
                end = List_of_annotation[each_gene][each_type]['end'][0]
                strand = List_of_annotation[each_gene][each_type]['strand'][0]
                
                if not start or not end or not strand:
                    print(f"Warning: Empty data for {genes[each_gene]} {each_type}. Skipping.")
                    continue
                
                fasta_seq = ROI_fasta.fetch(ROI_entry_name, (start,end), strand = strand)
        
                List_of_fasta[each_gene][each_type].append(fasta_seq)
                print(f"Extracted {genes[each_gene]} for the fasta sequence of {each_type}")


# Function to initialize the output files
def initialize_output_files(List_of_output_dir):
    """
    Initialize empty output files for each file path in the provided list of directories
    to ensure that all output files exist and are empty before starting to write data to
    them in subsequent operations.
    1. Iterates through each dictionary in List_of_output_dir.
    2. For each dictionary, extracts the file path (which is the first and only
       element of the value list).
    3. Opens the file in write mode and immediately closes it, effectively creating
       an empty file or overwriting an existing file with an empty one.

        - List_of_output_dir (list): A list of dictionaries, where each dictionary contains
                                     a single key-value pair. key: file type, value: list containing a single file path.
    """
    for type_dict in List_of_output_dir:
        for file_path in type_dict.values():
            with open(file_path[0], 'w') as f:
                pass


# Write the fasta files
def write_fasta(each_genome_name, genes, List_of_fasta, List_of_output_dir):
    """
    Write FASTA sequences to output files for each genome, gene and sequence type.
    The function assumes that output files have been initialized by initialize_output_files
    1. Iterates through each gene and FASTA type.
    2. Determines the appropriate output file for each sequence type.
    3. Handles 'CDS', 'exon', and 'intron' types differently from other types:
       - For 'CDS', 'exon', and 'intron': Writes multiple sequences if available, each with a unique identifier.
       - For other types: Writes a single sequence.
    4. Formats the FASTA header using genome name, gene name, sequence type, and (if applicable) sequence number.
    5. Writes sequences in FASTA format, wrapping lines at 60 characters.
    6. Appends new sequences to existing files.

        - each_genome_name (str): Name of the genome being processed.
        - genes (list): List of gene names.
        - List_of_fasta (list): List of dictionaries containing FASTA sequences for each gene and type.
        - List_of_output_dir (list): List of dictionaries containing output file paths for each FASTA type.

    Note:
    - The function uses append mode ('a') for writing, allowing multiple entries in the same file.
    - The function uses the global variable:
        - types_fasta (list of str): List of FASTA types to generate.
    """
    for each_gene in genes:
        index_gene = genes.index(each_gene)
        for each_type_fasta in types_fasta:
            output_file = List_of_output_dir[types_fasta.index(each_type_fasta)][each_type_fasta][0]
            
            if each_type_fasta in ['CDS', 'exon', 'intron']:
                seq_list = List_of_fasta[index_gene][each_type_fasta]
                if not seq_list:
                    print(f"Warning: No data for {each_gene} {each_type_fasta}. Skipping.")
                    continue
 
                for i, seq_str in enumerate(seq_list, 1):
                    n = 60
                    seq_fasta_list = [seq_str[i:i+n] for i in range(0, len(seq_str), n)]
                    header_list = [each_genome_name, each_gene, each_type_fasta, i]
                    header_name = '>' + '_'.join(map(str, header_list))
                    
                    with open(output_file, 'a') as outfile:
                        outfile.write(header_name + '\n')
                        outfile.write('\n'.join(seq_fasta_list) + '\n')
                    
                    print(f"Written entry for {header_name}")
            else:
                if not List_of_fasta[index_gene][each_type_fasta]:
                    print(f"Warning: No data for {each_gene} {each_type_fasta}. Skipping.")
                    continue

                header_list = [each_genome_name, each_gene, each_type_fasta]
                header_name = '>' + '_'.join(map(str, header_list))
                seq_str = List_of_fasta[index_gene][each_type_fasta][0]
                n = 60
                seq_list = [seq_str[i:i+n] for i in range(0, len(seq_str), n)]
                
                with open(output_file, 'a') as outfile:
                    outfile.write(header_name + '\n')
                    outfile.write('\n'.join(seq_list) + '\n')
                
                print(f"Written entry for {header_name}")


# Main processing
def annotation_to_fasta (each_genome_name, annotation_path, List_of_output_dir, ROI_fasta, genes_out=None):
    """
    Process GFF3 annotation files and convert them to FASTA format.
    The function handles various annotation types (gene, mRNA, exon, CDS, intron) and generates 
    corresponding FASTA sequences. It uses helper functions (populate_annotation, populate_fasta, 
    write_fasta) to process the data and generate output files.
    1. Reads the GFF3 file into a dataframe using pandas.
    2. Filters out the GFF3 file (if genes_out is provided)
    3. Extracts gene information and creates a list of gene names.
    4. Splits the main dataframe by gene model alignment.
    5. Sorts exons and CDS by start position within each gene.
    6. Creates data structures (List_of_annotation and List_of_fasta) to store processed data.
    7. Calls helper functions to:
        a. populate_annotation: populate List_of_annotation with annotation data
        b. populate_fasta: populate List_of_fasta with FASTA sequences
        c. write_fasta: write output files

        - each_genome_name (str): Name of the genome being processed.
        - annotation_path (str): Path to the GFF3 annotation file.
        - List_of_output_dir (list): List of dictionaries containing output file paths for each FASTA type.
        - ROI_fasta (pysam.FastaFile): FASTA file object for the region of interest.
        - genes_out (list, optional, default = None): List of gene names to exclude from processing.

    Note:
        - The function uses the global variables:
        - types_annotation (list of str): List of annotation types to process.
        - types_fasta (list of str): List of FASTA types to generate.
        - levels_annotation (list of str): 'start', 'end', 'strand'
    """
    ##### Read files:
    ### gmap annotation with pandas
    gff3_file = read_gff3(annotation_path)
    
    # Skip if the gff3_file is empty or None
    if gff3_file is None or gff3_file.empty:
        print(f"Warning: GFF3 file is empty or not found for genome {each_genome_name}. Skipping...")
        return  # Exit the function early, skip this genome

    ### Name of fasta sequence
    ROI_entry_name = gff3_file[0][0]

    # Get rid of some genes (genes_out): genes KLRD1, KLRD2 and flanking genes
    if genes_out:
        for gene in genes_out:
            info_column = gff3_file.iloc[:,-1]
            gff3_file = gff3_file[~info_column.str.contains(gene)]
    else:
        print("No genes to process in genes_out")

    # Reset index in the main dataframe
    gff3_file = gff3_file.reset_index()
    original_index = gff3_file.loc[:,'index']
    gff3_file = gff3_file.drop(columns=['index'])
    
    ## Get index of the rows that are type = gene
    # bool
    mask_gene = gff3_file.iloc[:,2] == 'gene'
    # dataframe of only type = gene
    gene_df = gff3_file[~mask_gene == False]
    # indexes of genes
    pos_gene = gff3_file[mask_gene].index

    ### list of gene names
    genes = []

    # Fill in list of genes
    info_column = gene_df.iloc[:,-1]
    
    for i in range(len(gene_df)):
        name_match = re.search(r'ID=([^;]+)', info_column[pos_gene[i]])
        if name_match:
            name = name_match.group(1)
        else:
            print("Name not found in attributes")
        genes.append(name)
    
    #1. Split the main dataframe by gene alignment (gene + mRNA + exons + CDS)
    ### Split the main dataframe
    # For loop to fill a dictionary: each entry is the dataframe from each gene
    df_by_gene = {}
    
    for i in range(len(pos_gene)):
        # For the last element of the gene index list pos_gene:
        if pos_gene[i] == pos_gene[-1]:
            df_by_gene[i] = gff3_file.iloc[pos_gene[i]:len(gff3_file),:].copy()
        # For the rest elements of the list:
        else:
            start = pos_gene[i]
            end = pos_gene[i+1]
            df_by_gene[i] = gff3_file.iloc[start:end,:].copy()
    
    # Order the exons and CDS by start position
    for i in range(len(pos_gene)):
        # Convert start and end positions to int
        df_by_gene[i][3] = pd.to_numeric(df_by_gene[i][3], errors='coerce').astype('Int64')
        df_by_gene[i][4] = pd.to_numeric(df_by_gene[i][4], errors='coerce').astype('Int64')
        # Separate exons and CDS 
        exons = df_by_gene[i][df_by_gene[i][2] == 'exon']
        cds = df_by_gene[i][df_by_gene[i][2] == 'CDS']
        # Sort exons and CDS separately by start position (column 3)
        exons_sorted = exons.sort_values(by=3, ascending=False)
        cds_sorted = cds.sort_values(by=3, ascending=False)
        # Combine other rows, exons, and CDS
        other_rows = df_by_gene[i][~df_by_gene[i][2].isin(['exon', 'CDS'])]
        df_sorted = pd.concat([other_rows, cds_sorted]) 
        df_sorted = pd.concat([other_rows, exons_sorted, cds_sorted])
        # Substitute df for sorted df
        df_by_gene[i] = df_sorted

    #2. Create data structures: List_of_annotation and List_of_fasta
    ###### 2.1 List_of_annotation 
    ###### Create a data structure of nested dictionaries and lists to populate
    ###### with the data each level for each type and gene name
    List_of_annotation = list()

    for gene_name in genes: 
        each_dict = {}
        List_of_annotation.append(each_dict)   
    
    for i in range(len(genes)):
        for each_type in types_annotation:
            List_of_annotation[i][each_type] = {}
   
        for each_type in types_annotation:
            for each_level in levels_annotation:
                List_of_annotation[i][each_type][each_level] = []
            
    ###### 2.2 List_of_fasta
    ###### Create a data structure: list of dictionaries (one per gene) to populate
    ###### with the fasta sequences from pyfastx
    List_of_fasta = list()

    for gene_name in genes: 
        each_dict = {}
        List_of_fasta.append(each_dict)
    
    for i in range(len(genes)):
        for each_type in types_fasta:
            List_of_fasta[i][each_type] = [] 
        
    print("Edited annotation dataframe")

    ## 3. Run functions
    populate_annotation (df_by_gene, genes, List_of_annotation)
    populate_fasta (List_of_annotation, List_of_fasta, ROI_fasta, ROI_entry_name, genes)    
    write_fasta (each_genome_name, genes, List_of_fasta, List_of_output_dir)
 
    print("Processed annotation to fasta")


# Main function
def get_fasta_sequences_annotation(genome_file, gene_to_annotate, input_dir_ROI, input_dir_annotation, output_dir, genes_out = None):
    """
    Process multiple genome' GFF3 annotations and ROI FASTA files to generate FASTA sequences from the annotation coordinates.
    This function calls a series of helper functions:
    1. get_list_genome: Reads the list of genome from the genome_file.
    2. create_output_dir: Creates output directories for each FASTA type.
    3. initialize_output_files: Initializes output files.
    4. For each genome:
       a. Locate the corresponding ROI FASTA and GFF3 annotation files.
       b. Reads the ROI FASTA file using pyfastx.
       c. Calls annotation_to_fasta to process the annotation data and generate FASTA sequences.

        - genome_file (str): Path to a file containing a list of genome names.
        - gene_to_annotate (str): Name of the gene being annotated (Example: 'NKG2')
        - input_dir_ROI (str): Directory containing ROI FASTA files.
        - input_dir_annotation (str): Directory containing annotation files.
        - output_dir (str): Directory where FASTA output files will be saved.
        - genes_out (list, optional, default = None): List of gene names to exclude from processing.
    """
    print("Get FASTA sequences from the annotations")
    
    genome_list = get_list_genome(genome_file)

    ## Loop over the files and run functions
    fasta_patterns = ['.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa', '.frn']
    input_dir_ROI_list = [f for f in os.listdir(input_dir_ROI) 
                      if any(f.lower().endswith(p) for p in fasta_patterns) 
                      and not f.endswith('fasta.fxi')]
    input_dir_annotation_list = os.listdir(input_dir_annotation)

    # Loop over each genome annotation file
    for each_genome_name in genome_list:
        print(each_genome_name)
        
        # Create paths for the output FASTA files
        List_of_output_dir = create_output_dir(each_genome_name, gene_to_annotate, output_dir)
        initialize_output_files(List_of_output_dir)

        # This line returns the index in the list of the first item (therefore the [0] at the end) whose
        # value contains a certain string (name of the genome) 
        index_ROI = [idx for idx, s in enumerate(input_dir_ROI_list) if each_genome_name in s][0]
        index_annotation = [idx for idx, s in enumerate(input_dir_annotation_list) if each_genome_name in s][0]
   
        ##### Read and save ROI fasta sequence
        # Pyfastx: enter the data in List_of_annotation to get the fasta sequences
        ROI_file = input_dir_ROI_list[index_ROI]
        ROI_path = os.path.join(input_dir_ROI, ROI_file)
    
        # Read and save ROI fasta sequence
        # Pyfastx: enter the data in List_of_annotation to get the fasta sequences
        ROI_fasta = pyfastx.Fasta(ROI_path)

        annotation_name = input_dir_annotation_list[index_annotation]
        annotation_path = os.path.join(input_dir_annotation, annotation_name)
    
        annotation_to_fasta (each_genome_name, annotation_path, List_of_output_dir, ROI_fasta, genes_out)
                             
        # Remove fxi files created by pyfastx
        fxi_file = ROI_path + ".fxi"
        
        if os.path.exists(fxi_file):
            os.remove(fxi_file)
            print(f"Deleted index file: {fxi_file}")
        else:
            print(f"Index file not found: {fxi_file}")


genome_file='/home/lahumada/Desktop/AnnCX_gene_annotation_pipeline/examples/annotate2fasta/TXT/txt_genome.txt'
gene_to_annotate='annotation2fasta_test'
input_dir_ROI='/home/lahumada/Desktop/AnnCX_gene_annotation_pipeline/examples/genomic_sequences/genome'
input_dir_annotation='/home/lahumada/Desktop/AnnCX_gene_annotation_pipeline/examples/annotate2fasta/GFF3'
output_dir='/home/lahumada/Desktop/AnnCX_gene_annotation_pipeline/testing/output_revision/AnnCX_NKG2'
get_fasta_sequences_annotation(genome_file, gene_to_annotate, input_dir_ROI, input_dir_annotation, output_dir)
