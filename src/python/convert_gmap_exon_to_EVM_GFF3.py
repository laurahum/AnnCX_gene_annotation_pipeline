#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 22:23:31 2024

@author: lahumada
"""

'''
This script converts filtered GFF3 output from GMAP (exons) to a format compatible with EVM

1. Convert features
The output from GMAP contains for features:
    gene, mRNA, exon, CDS
from all those features EVM (transcript alignment) accepts:
    cDNA_match
the exon features in GMAP filtered GFF3 file must be converted to cDNA_match features

2. Modify attributes column
The attributes column of the GMAP filtered GFF3 files must be modified according to the example provided by EVM
Give the entries associated to the same gene model, in order of entry, a counter in the attributes section:
ex. ID=match.gmapcDNA.1 for the features that are associated to the first gene model entry.

The following example file provided by EVM documentation was used as a guide to make the formatting changes
https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/transcript_alignment.gff3

Input: GMAP filtered GFF3 output (with features: gene, mRNA, exon, CDS) 
Output: GMAP filtered EVM-compatible GFF3 file (with features: cDNA_match)

Main function: convert_gmap_exon_to_EVM_all_files(input_dir, output_dir, source)
Note: source argument must be: 'gmapExon'
'''

# import sys
import os
import re


# Main function to convert the format
def convert_gmap_to_EVM(input_file, output_file, source):
    """
    Convert GMAP filtered-GFF3 output to GFF3 format compatible with EVM.
    ONLY exon -> cDNA_match (gene provides grouping/ID)
        - input_file (str): Path to the input GMAP output file
        - output_file (str): Path to the output GFF3 file compatible with EVM
        - source (str): Name of the tool used to generate the annotation output. Possible values:
            - gmapExon: run GMAP with exons as query but keep only exon features
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Counter for generating unique gene IDs
        gene_id = 1  
        # Store the current gene being processed
        current_gene = None  
        # List to accumulate exon features only
        exons = []

        for line in infile:
            line = line.strip()
            
            # Skip empty lines or comments
            if not line or line.startswith('#'):
                continue
            
            # Check for end of gene model
            if line == "###":
                # If there is a previous gene, write its exons
                if current_gene and exons: 
                    write_gene_exons(outfile, current_gene, exons, gene_id, source)
                    gene_id += 1
                current_gene = None
                exons = []
                continue
                
            fields = line.split('\t')
            # Check if it is a valid annotation entry
            if len(fields) == 9:  
                seqname, source_field, feature, start, end, score, strand, frame, attributes = fields
                
                if feature == 'gene':  
                    # If there is a previous gene, write its exons
                    if current_gene and exons: 
                        write_gene_exons(outfile, current_gene, exons, gene_id, source)
                        gene_id += 1
                    # Start a new gene model
                    current_gene = fields  
                    # Reset exon list for new gene
                    exons = []
                    
                elif feature == 'mRNA':
                    # Explicitly skip mRNA - do nothing
                    pass
                elif feature == 'CDS':
                    # Explicitly skip CDS - do nothing  
                    pass
                elif feature == 'exon':
                    # ONLY accumulate exon features
                    exons.append(fields)

        # Write the last gene if it exists
        if current_gene and exons:
            write_gene_exons(outfile, current_gene, exons, gene_id, source)


# Write each gene entry with all the features
def write_gene_exons(outfile, gene, exons, gene_id, source):
    """
    Write all exon features of a gene as cDNA_match entries with shared gene ID.
    """
    # Get gene data
    seqname, _, feature, start, end, score, strand, frame, attributes = gene
    
    # Extract gene name from attributes (Name=xxx)
    gene_match = re.search(r'Name=([^;\s]+)', attributes)
    gene_name = gene_match.group(1) if gene_match else f'gene_{gene_id}'

    # Write ONLY exons as cDNA_match
    for exon in exons:
        seqname_t, _, feature_t, t_start, t_end, t_score, t_strand, t_frame, t_attributes = exon
        
        # Extract target name from exon attributes, fallback to gene name
        target_match = re.search(r'Name=([^;\s]+)', t_attributes)
        target_name = target_match.group(1) if target_match else gene_name
        
        # Use exon score (or . if missing)
        score_out = t_score if t_score and t_score != '.' else '.'
        
        outfile.write(f"{seqname_t}\t{source}\tcDNA_match\t{t_start}\t{t_end}\t{score_out}\t{t_strand}\t.\tID=match.{source}.{gene_id};Target={target_name}\n")


def convert_gmap_exon_to_EVM_all_files(input_dir, output_dir, source):
    '''
    Loop the convert_gmap_to_EVM function over all the filtered GMAP result annotation files for each genome
        - input_dir: Path to filtered-GFF3 GMAP result annotation files
        - output_dir: Path to filtered-EVM-GFF3 GMAP result annotation files
        - source: Name of the tool used to generate the annotation output.
            - 'gmapExon': run GMAP with exons as query but keep only exon features
    '''
    # Loop the main function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename
        output_name = output_name.rsplit('.', 1)[0] + '_EVM.gff3'
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files
        convert_gmap_to_EVM(Filepath_input, Filepath_output, source)
        print("Saved_" + output_name)
    