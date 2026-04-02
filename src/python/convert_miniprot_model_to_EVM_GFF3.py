#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 19:13:47 2026

@author: lahumada
"""

'''
This script converts raw GFF3 output from miniprot to a GFF3 format compatible with EVM

Features to convert:
The output from miniprot contains annotations for features:
    mRNA, CDS, stop_codon
EVM (gene predictions) accepts:
    gene, mRNA, exon, CDS

1. mRNA coordinates -> gene coordinates  
2. Keep mRNA features as-is
3. CDS coordinates -> exon coordinates
4. Keep CDS features as-is
5. Skip stop_codon features

Modify attributes column to match EVM gene prediction format:
Generate hierarchical IDs: gene → mRNA → exon -> CDS with Parent relationships

Input: miniprot formatted filtered GFF3 (with features: mRNA, CDS, stop_codon)
Output: miniprot EVM-compatible GFF3 file (with features: gene, mRNA, exon, CDS)

Main function: convert_miniprot_model_to_EVM_all_files(input_dir, output_dir, source)
Note: source argument must be: 'miniprotModel'
'''

import os
import re

# Main function to convert the format
def convert_miniprot_to_EVM(input_file, output_file, source):
    """
    Convert miniprot raw GFF3 output format to GFF3 format compatible with EVM.
    mRNA, CDS, stop_codon -> gene, mRNA, exon, CDS (with EVM-standardized attributes)
        - input_file: Path to the raw GFF3 miniprot result annotation file
        - output_file: Path to the EVM GFF3 miniprot result annotation file
        - source (str): Name of the tool used to generate the annotation output. Possible values:
        	- miniprotModel
    """
     
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        
        # Counter for generating unique gene IDs       
        gene_id = 1
        # None = last gene model finished writing
        current_mRNA = None 
        
        # Lists to accumulate the CDS for each gene model
        cds_features = []

        for line in infile:
            line = line.strip()
            
            # Skip empty lines or comments
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            if len(fields) == 9:  # Annotation entry
                seqname, source_field, feature, start, end, score, strand, frame, attributes = fields
                
                if feature == 'mRNA':
                    # Write previous model in case it is still not written
                    if current_mRNA:
                        write_gene(outfile, current_mRNA, cds_features, gene_id, source)
                        gene_id += 1
                    # Start a new gene model
                    current_mRNA = fields
                    cds_features = []
                    
                elif feature == 'CDS':
                    cds_features.append(fields)
                    
                elif feature == 'stop_codon':
                    # Skip stop_codon features
                    pass

        # Write the last gene if exists
        if current_mRNA:
            write_gene(outfile, current_mRNA, cds_features, gene_id, source)


# Function to write a gene model
def write_gene(outfile, mRNA, cds_features, gene_id, source):
    """
    Write a single gene model in GFF3 format compatible with EVM.
        - outfile (file): Output file object
        - mRNA (list): mRNA feature fields (used for gene coordinates)
        - cds_features (list): List of CDS feature fields  
        - gene_id (int): Unique identifier for the gene
    """
    
    seqname, _, mRNA_feature, mRNA_start, mRNA_end, mRNA_score, mRNA_strand, mRNA_frame, mRNA_attributes = mRNA
    
    # Extract gene name from mRNA attributes
    gene_match = re.search(r'ID=([^;]+)', mRNA_attributes)
    gene_name = gene_match.group(1) if gene_match else f"gene_{gene_id}"
    
    # Write gene feature (using mRNA coordinates)
    outfile.write(f"{seqname}\t{source}\tgene\t{mRNA_start}\t{mRNA_end}\t{mRNA_score}\t{mRNA_strand}\t.\tID={seqname}.tPRED{gene_id:06d};Name={source} model {gene_name}\n")
    
    # Write mRNA feature 
    outfile.write(f"{seqname}\t{source}\tmRNA\t{mRNA_start}\t{mRNA_end}\t{mRNA_score}\t{mRNA_strand}\t.\tID={seqname}.m{gene_id:06d};Parent={seqname}.tPRED{gene_id:06d}\n")
    
    # Write exon features (using CDS coordinates)
    for i, cds in enumerate(cds_features, 1):
        seqname_c, source_c, feature_c, c_start, c_end, c_score, c_strand, c_frame, c_attributes = cds
        outfile.write(f"{seqname_c}\t{source}\texon\t{c_start}\t{c_end}\t{c_score}\t{c_strand}\t.\tID={seqname}.e{i:06d};Parent={seqname}.m{gene_id:06d}\n")
    
    # Write CDS features 
    for i, cds in enumerate(cds_features, 1):
        seqname_c, source_c, feature_c, c_start, c_end, c_score, c_strand, c_frame, c_attributes = cds
        outfile.write(f"{seqname_c}\t{source}\tCDS\t{c_start}\t{c_end}\t{c_score}\t{c_strand}\t{c_frame}\tID=cds_of_{seqname}.m{gene_id:06d};Parent={seqname}.m{gene_id:06d}\n")


# Main
def convert_miniprot_model_to_EVM_all_files(input_dir, output_dir, source):
    '''
    Loop the convert_miniprot_to_EVM function over all the raw miniprot result annotation files for each genome
        - input_dir: Path to raw-GFF3 miniprot result annotation files
        - output_dir: Path to EVM-GFF3 miniprot result annotation files
        - source: 'miniprotModel'
    '''
    # Loop the main function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename.rsplit('.', 1)[0] + '_EVM.gff3'
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files
        convert_miniprot_to_EVM(Filepath_input, Filepath_output, source)
        print("Saved_" + output_name)
        
