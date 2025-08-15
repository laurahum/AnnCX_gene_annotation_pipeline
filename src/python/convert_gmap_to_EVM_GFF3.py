#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 22:23:31 2024

@author: lahumada
"""


# This script converts filtered GFF3 output from GMAP (cDNA or exons) to a format compatible with EVM
#
# 1. Convert features
# The output from GMAP contains for features:
#     gene, mRNA, exon, CDS
# from all those features EVM (transcript alignment) accepts:
#     cDNA_match
# all the different features in GMAP filtered GFF3 file must be converted to cDNA_match features
#
# 2. Modify attributes column
# The attributes column of the GMAP filtered GFF3 files must be modified according to the example provided by EVM
# Give the entries associated to the same gene model, in order of entry, a counter in the attributes section:
# ex. ID=match.gmapcDNA.1 for the features that are associated to the first gene model entry.
#
# The following example file provided by EVM documentation was used as a guide to make the formatting changes
# https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/transcript_alignment.gff3
#
# Input: GMAP filtered GFF3 output (with features: gene, mRNA, exon, CDS) 
# Output: GMAP filtered EVM-compatible GFF3 file (with features: cDNA_match)

# Main function: convert_gmap_to_EVM_all_files(input_dir, output_dir, source)
# Note: source argument must be either:
#   'gmapcDNA'
#   'gmapExon'


# import sys
import os
import re


# Main function to convert the format
def convert_gmap_to_EVM(input_file, output_file, source):
    """
    Convert GMAP filtered-GFF3 output to GFF3 format compatible with EVM.
    gene, mRNA, exon, CDS -> cDNA_match
        - input_file (str): Path to the input GMAP output file
        - output_file (str): Path to the output GFF3 file compatible with EVM
        - source (str): Name of the tool used to generate the annotation output. Possible values:
            - gmapcDNA: run GMAP with cDNA as query
            - gmapExon: run GMAP with exons as query
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Counter for generating unique gene IDs
        gene_id = 1  
        # Store the current gene being processed
        current_gene = None  
        # List to accumulate features (mRNA, exon, CDS)
        mRNA = []
        exons = []
        cdss = []

        for line in infile:
            line = line.strip()
            
            # Check for the end of the current gene model
            if line == "###":  
                # If there is a current gene, write it to the output
                if current_gene:  
                    write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id, source)
                    gene_id += 1  # Increment the gene ID for the next gene
                
                # Reset current gene and features
                current_gene = None  
                mRNA = []
                exons = []
                cdss = []
            else:
                fields = line.split('\t')
                # Check if it is a valid annotation entry
                if len(fields) == 9:  
                    seqname, _, feature, start, end, score, strand, frame, attributes = fields
                    if feature == 'gene':  
                         # If there is a previous gene, write it
                        if current_gene: 
                            write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id, source)
                            gene_id += 1
                        # Start a new gene model
                        current_gene = fields  
                    # Accumulate features (mRNA, exon, CDS)
                    elif feature == 'mRNA':
                        mRNA.append(fields)  
                    elif feature == 'CDS':
                        cdss.append(fields)
                    elif feature == 'exon':
                        exons.append(fields)

        # Write the last gene if it exists (in case the file does not end with '###')
        if current_gene:
            write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id, source)


# Write each gene entry with all the features
def write_gene(outfile, gene, mRNA, exons, cdss, gene_id, source):
    """
    Write all features of a gene model as cDNA_match entries.
        - outfile (file): Output file object
        - gene (list): Gene feature fields
        - mRNA (list): List of mRNA
        - exons (list): list of exons
        - cdss (list): list of cdss
        - gene_id (int): Unique identifier for the gene
        - source (str): Name of the tool used to generate the annotation output. Possible values:
            - gmapcDNA: run GMAP with cDNA as query
            - gmapExon: run GMAP with exons as query
    """
    # Get data for each field
    seqname, _, feature, start, end, score, strand, frame, attributes = gene
    
    # Extract the gene name
    gene_match = re.search(r'Name=([^;]+)', attributes)
    gene_name = gene_match.group(1)

    # Write gene feature                                                           
    outfile.write(f"{seqname}\t{source}\tcDNA_match\t{start}\t{end}\t{score}\t{strand}\t.\tID=match.{source}.{gene_id};Target={gene_name}\n")
    
    # Write mRNA feature
    for i, mRNA in enumerate(mRNA, 1):
        seqname, _, feature, start, end, score, strand, frame, attributes = mRNA  
        outfile.write(f"{seqname}\t{source}\tcDNA_match\t{start}\t{end}\t{score}\t{strand}\t.\tID=match.{source}.{gene_id};Target={gene_name}\n")
    
    # Write exon features
    for i, exon in enumerate(exons, 1):
        seqname, _, feature, start, end, score, strand, frame, attributes = exon  
        
        # Get the coordinates annotations in the sequence used as query
        target_match = re.search(r'Target=\S+\s+(\S+)\s+(\S+)', attributes)
        target_start, target_end = target_match.groups()
        
        outfile.write(f"{seqname}\t{source}\tcDNA_match\t{start}\t{end}\t{score}\t{strand}\t.\tID=match.{source}.{gene_id};Target={gene_name} {target_start} {target_end}\n")
    
    # Write CDS features
    for i, cds in enumerate(cdss, 1):
        seqname, _, feature, start, end, score, strand, frame, attributes = cds
        
        # Get the coordinates annotations in the sequence used as query
        target_match = re.search(r'Target=\S+\s+(\S+)\s+(\S+)', attributes)
        target_start, target_end = target_match.groups()
        
        outfile.write(f"{seqname}\t{source}\tcDNA_match\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID=match.{source}.{gene_id};Target={gene_name} {target_start} {target_end}\n")


def convert_gmap_to_EVM_all_files(input_dir, output_dir, source):
    '''
    Loop the convert_gmap_to_EVM function over all the filtered GMAP result annotation files for each genome
        - input_dir: Path to filtered-GFF3 GMAP result annotation files
        - output_dir: Path to filtered-EVM-GFF3 GMAP result annotation files
        - source: Name of the tool used to generate the annotation output. Possible values:
            - gmapcDNA: run GMAP with cDNA as query
            - gmapExon: run GMAP with exons as query
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
        print("Saved_" + output_name + '\n')
    
