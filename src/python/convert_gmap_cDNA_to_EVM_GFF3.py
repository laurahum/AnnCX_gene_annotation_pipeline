#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 22:23:31 2024

@author: lahumada
"""

'''
This script converts filtered GFF3 output from GMAP (transcript) to a format compatible with EVM

1. Convert features
The output from GMAP contains for features:
    gene, mRNA, exon, CDS
from all those features EVM (gene predictions) accepts:
    gene, mRNA, exon, CDS
all the different features in GMAP filtered GFF3 file must be converted to the gene model
format compatible with EVM.

2. Modify attributes column
The attributes column of the GMAP filtered GFF3 files must be modified according to the example provided by EVM
Give the entries associated to the same gene model, in order of entry, a counter in the attributes section:
ex. ID=match.gmaptranscript.1 for the features that are associated to the first gene model entry.

The following example file provided by EVM documentation was used as a guide to make the formatting changes
https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/gene_predictions.gff3

Input: GMAP filtered GFF3 output (with features: gene, mRNA, exon, CDS) 
Output: GMAP filtered EVM-compatible GFF3 file (with features: gene, mRNA, exon, CDS) 

Main function: convert_gmap_cDNA_to_EVM_all_files(input_dir, output_dir, source)
Note: source argument must be:'gmaptranscript'
'''

# import sys
import os
import re

# Main function to convert the format
def convert_gmap_to_EVM(input_file, output_file, source):
    """
    Convert GMAP raw GFF3 output format to GFF3 format compatible with EVM.
    gene, mRNA, exon, CDS -> gene, mRNA, exon, CDS (EVM-standardized format)
        - input_file: Path to the raw GFF3 GMAP result annotation file
        - output_file: Path to the EVM GFF3 GMAP result annotation file
        - source (str): Name of the tool used to generate the annotation output. 
                        Must be 'gmaptranscript'
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        gene_id = 1
        current_gene = None 
        mRNAs = []
        exons = []
        cdss = []

        for line in infile:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            if len(fields) != 9:
                continue
                
            seqname, source_field, feature, start, end, score, strand, frame, attributes = fields
            
            if line.strip() == "###":
                # End of gene model
                if current_gene:
                    write_gene(outfile, current_gene, mRNAs, exons, cdss, gene_id, source)
                    gene_id += 1
                current_gene = None
                mRNAs, exons, cdss = [], [], []
                continue
            
            if feature == 'gene':
                # Write previous gene
                if current_gene:
                    write_gene(outfile, current_gene, mRNAs, exons, cdss, gene_id, source)
                    gene_id += 1
                current_gene = fields
                mRNAs, exons, cdss = [], [], []
            elif feature == 'mRNA':
                mRNAs.append(fields)
            elif feature == 'exon':
                exons.append(fields)
            elif feature == 'CDS':
                cdss.append(fields)

        # Write final gene
        if current_gene:
            write_gene(outfile, current_gene, mRNAs, exons, cdss, gene_id, source)


# Function to write a gene model
def write_gene(outfile, gene, mRNAs, exons, cdss, gene_id, source):
    seqname, _, feature, start, end, score, strand, frame, attributes = gene
    
    # Write gene feature
    outfile.write(f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t.\tID={seqname}.tPRED{gene_id:06d};Name={source} model {seqname}.m{gene_id:06d}\n")
    
    # Write mRNA feature (first one)
    if mRNAs:
        seqname_m, source_m, feature_m, m_start, m_end, m_score, m_strand, m_frame, m_attributes = mRNAs[0]
        outfile.write(f"{seqname}\t{source}\tmRNA\t{m_start}\t{m_end}\t{m_score}\t{m_strand}\t.\tID={seqname}.m000001;Parent={seqname}.tPRED{gene_id:06d}\n")
    
    # Write exon features
    for i, exon in enumerate(exons, 1):
        seqname_e, source_e, feature_e, e_start, e_end, e_score, e_strand, e_frame, e_attributes = exon
        outfile.write(f"{seqname_e}\t{source}\texon\t{e_start}\t{e_end}\t{e_score}\t{e_strand}\t.\tID={seqname}.e{i:06d};Parent={seqname}.m000001\n")
    
    # Write CDS features  
    for i, cds in enumerate(cdss, 1):
        seqname_c, source_c, feature_c, c_start, c_end, c_score, c_strand, c_frame, c_attributes = cds
        outfile.write(f"{seqname_c}\t{source}\tCDS\t{c_start}\t{c_end}\t{c_score}\t{c_strand}\t{c_frame}\tID=cds_of_{seqname}.m000001;Parent={seqname}.m000001\n")

# Main
def convert_gmap_cDNA_to_EVM_all_files(input_dir, output_dir, source):
    '''
    Loop the convert_gmap_to_EVM function over all the raw GMAP result annotation files for each genome
        - input_dir: Path to raw-GFF3 GMAP result annotation files
        - output_dir: Path to EVM-GFF3 GMAP result annotation files
        - source: 'gmaptranscript'
    '''
    # Loop the main function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename.rsplit('.', 1)[0] + '_EVM.gff3'
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files
        convert_gmap_to_EVM(Filepath_input, Filepath_output, source)
        print("Saved_" + output_name)
