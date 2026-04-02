#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 19:13:47 2026

@author: lahumada
"""

'''
This script converts raw GFF3 output from minimap2 to a GFF3 format compatible with EVM

Features to convert:
The output from minimap2 contains annotations for features:
    gene, mRNA, exon
EVM (gene predictions) accepts:
    gene, mRNA, exon, CDS

1. Keep gene features as-is
2. Keep mRNA features as-is
3. Keep exon features as-is
4. Use exon features to create CDS features

Modify attributes column to match EVM gene prediction format:
Generate hierarchical IDs: gene → mRNA → exon -> CDS with Parent relationships

Input: minimap2 formatted filtered GFF3 (with features: gene, mRNA, exon)
Output: minimap2 EVM-compatible GFF3 file (with features: gene, mRNA, exon, CDS)

Main function: convert_minimap2_model_to_EVM_all_files(input_dir, output_dir, source)
Note: source must be 'minimap2Model'

'''

import os
import re

# Main function to convert the format
def convert_minimap2_to_EVM(input_file, output_file, source):
    """
    Convert minimap2 raw GFF3 output format to GFF3 format compatible with EVM.
    gene, mRNA, exon -> gene, mRNA, exon (with EVM-standardized attributes)
        - input_file: Path to the raw GFF3 minimap2 result annotation file
        - output_file: Path to the EVM GFF3 minimap2 result annotation file
        - source: must be 'minimap2Model'
    """
     
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        
        # Counter for generating unique gene IDs       
        gene_id = 1
        # None = last gene model finished writing
        current_gene = None 
        
        # Lists to accumulate the mRNA and exons for each gene model
        mRNAs = []
        exons = []

        for line in infile:
            line = line.strip()
            
            # Skip empty lines or comments
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            if len(fields) == 9:  # Annotation entry
                seqname, source_field, feature, start, end, score, strand, frame, attributes = fields
                
                if feature == 'gene':
                    # Write previous model in case it is still not written
                    if current_gene:
                        write_gene(outfile, current_gene, mRNAs, exons, gene_id, source)
                        gene_id += 1
                    # Start a new gene model
                    current_gene = fields
                    mRNAs = []
                    exons = []
                    
                elif feature == 'mRNA':
                    mRNAs.append(fields)
                    
                elif feature == 'exon':
                    exons.append(fields)

        # Write the last gene if exists
        if current_gene:
            write_gene(outfile, current_gene, mRNAs, exons, gene_id, source)


# Function to write a gene model
def write_gene(outfile, gene, mRNAs, exons, gene_id, source):
    """
    Write a single gene model in GFF3 format compatible with EVM.
        - outfile (file): Output file object
        - gene (list): Gene feature fields
        - mRNAs (list): List of mRNA feature fields
        - exons (list): List of exon feature fields
        - gene_id (int): Unique identifier for the gene
    """
    
    seqname, _, feature, start, end, score, strand, frame, attributes = gene
    
    # Extract gene name from attributes
    gene_match = re.search(r'Name=([^;]+)', attributes)
    gene_name = gene_match.group(1) if gene_match else f"gene_{gene_id}"
    
    # Write gene feature with EVM naming convention
    outfile.write(f"{seqname}\t{source}\tgene\t{start}\t{end}\t{score}\t{strand}\t.\tID={seqname}.tPRED{gene_id:06d};Name={source} model {gene_name}\n")
    
    # Write mRNA feature (use first mRNA)
    if mRNAs:
        seqname_m, source_m, feature_m, m_start, m_end, m_score, m_strand, m_frame, m_attributes = mRNAs[0]
        outfile.write(f"{seqname}\t{source}\tmRNA\t{m_start}\t{m_end}\t{m_score}\t{m_strand}\t.\tID={seqname}.m{gene_id:06d};Parent={seqname}.tPRED{gene_id:06d}\n")
    
    # Write exon features
    for i, exon in enumerate(exons, 1):
        seqname_e, source_e, feature_e, e_start, e_end, e_score, e_strand, e_frame, e_attributes = exon
        outfile.write(f"{seqname_e}\t{source}\texon\t{e_start}\t{e_end}\t{e_score}\t{e_strand}\t.\tID={seqname}.e{i:06d};Parent={seqname}.m{gene_id:06d}\n")
    
    # Write CDS features (same coordinates as exons)
    for i, exon in enumerate(exons, 1):
        seqname_e, source_e, feature_e, e_start, e_end, e_score, e_strand, e_frame, e_attributes = exon
        outfile.write(f"{seqname_e}\t{source}\tCDS\t{e_start}\t{e_end}\t{e_score}\t{e_strand}\t.\tID=cds_of_{seqname}.m{gene_id:06d};Parent={seqname}.m{gene_id:06d}\n")


# Main
def convert_minimap2_model_to_EVM_all_files(input_dir, output_dir, source):
    '''
    Loop the convert_minimap2_to_EVM function over all the raw minimap2 result annotation files for each genome
        - input_dir: Path to raw-GFF3 minimap2 result annotation files
        - output_dir: Path to EVM-GFF3 minimap2 result annotation files
        - source: 'minimap2Model'
    '''
    # Loop the main function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename.rsplit('.', 1)[0] + '_EVM.gff3'
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files
        convert_minimap2_to_EVM(Filepath_input, Filepath_output, source)
        print("Saved_" + output_name)
