#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 22:20:28 2024

@author: lahumada
"""

'''
This script converts raw GFF3 output from AUGUSTUS to a GFF3 format compatible with EVM

1. Convert features
The output from AUGUSTUS contains annotations for features:
    gene, transcript, stop_codon, CDS, start_codon
from all those features EVM (gene predictions) accepts only:
    gene, mRNA, exon and cds

- transcript in the augustus output will be now mRNA,
- If the AUGUSTUS output file lacks exon features, the CDS features will be used to create
exon features with the same coordinates. This is the approach that the EVM
example script to convert AUGUSTUS output takes when applied to AUGUSTUS output that lacks
exon features (augustus_GTF_to_EVM_GFF3.pl)
- Omit features stop_codon and start_codon in the EVM-compatible file

2. Modify attributes column
The attributes column of the AUGUSTUS raw GFF3 files must be modified according to the example provided by EVM.
The format of this column will still reflect the gene models.

The following example file provided by EVM documentation was used as a guide to make the formatting changes
https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/gene_predictions.gff3

Input: AUGUSTUS GFF3 (or GTF-like) output (with features: gene, transcript, stop_codon, CDS, start_codon)
Output: AUGUSTUS EVM-compatible GFF3 file (with features: gene, mRNA, exon, cds)

Main function: convert_augustus_to_EVM_all_files(input_dir, output_dir)
Note: source is omited as an argument and taken from the second column of Augustus GFF files (AUGUSTUS)
'''


import os

# Main function to convert the format
def convert_augustus_to_EVM(input_file, output_file):
    """
    Convert AUGUSTUS raw GFF3 output format to GFF3 format compatible with EVM.
    gene, [transcript/exon/CDS] -> gene, mRNA, exon and cds
    If no exon features present, creates exons from CDS coordinates
        - input_file: Path to the raw GFF3 AUGUSTUS result annotation file
        - output_file: Path to the EVM GFF3 AUGUSTUS result annotation file
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        gene_id = 1
        current_gene = None 
        mRNA = []
        exons = []
        cdss = []
        has_exons = False  # Track if we've seen real exons

        for line in infile:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            if len(fields) != 9:
                continue
                
            seqname, source, feature, start, end, score, strand, frame, attributes = fields
            
            if line.strip() == "###":
                # End of gene model
                if current_gene:
                    write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id, has_exons)
                    gene_id += 1
                current_gene = None
                mRNA, exons, cdss = [], [], []
                has_exons = False
                continue
            
            if feature == 'gene':
                # Write previous gene
                if current_gene:
                    write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id, has_exons)
                    gene_id += 1
                current_gene = fields
                mRNA, exons, cdss = [], [], []
                has_exons = False
            elif feature == 'transcript':
                mRNA.append(fields)
            elif feature == 'exon':
                exons.append(fields)
                has_exons = True
                cdss.append(fields) if feature == 'CDS' else None
            elif feature == 'CDS':
                cdss.append(fields)
                if not has_exons:  # Only create exons from CDS if no real exons
                    exons.append(fields)

        # Write final gene
        if current_gene:
            write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id, has_exons)



# Function to write a gene model
def write_gene(outfile, gene, mRNA, exons, cdss, gene_id, has_exons):
    seqname, source, feature, start, end, score, strand, frame, attributes = gene
    
    # Write gene feature
    outfile.write(f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t.\tID={seqname}.tPRED{gene_id:06d};Name={source} model {seqname}.m{gene_id:06d}\n")
    
    # Write mRNA feature (first one)
    if mRNA:
        seqname_m, source_m, feature_m, m_start, m_end, m_score, m_strand, m_frame, m_attributes = mRNA[0]
        outfile.write(f"{seqname}\t{source}\tmRNA\t{m_start}\t{m_end}\t{m_score}\t{m_strand}\t.\tID={seqname}.m000001;Parent={seqname}.tPRED{gene_id:06d}\n")
    
    # Write exon features
    for i, exon in enumerate(exons, 1):
        seqname_e, source_e, feature_e, e_start, e_end, e_score, e_strand, e_frame, e_attributes = exon
        # Use '.' phase for CDS-derived exons
        exon_phase = '.' if feature_e == 'CDS' else e_frame
        outfile.write(f"{seqname_e}\t{source}\texon\t{e_start}\t{e_end}\t{e_score}\t{e_strand}\t{exon_phase}\tID={seqname}.e{i:06d};Parent={seqname}.m000001\n")
    
    # Write CDS features  
    for i, cds in enumerate(cdss, 1):
        seqname_c, source_c, feature_c, c_start, c_end, c_score, c_strand, c_frame, c_attributes = cds
        outfile.write(f"{seqname_c}\t{source}\tCDS\t{c_start}\t{c_end}\t{c_score}\t{c_strand}\t{c_frame}\tID=cds_of_{seqname}.m000001;Parent={seqname}.m000001\n")

# Main
def convert_augustus_to_EVM_all_files(input_dir, output_dir):
    '''
    Loop the convert_augustus_to_EVM function over all the raw AUGUSTUS result annotation files for each genome
        - input_dir: Path to raw-GFF3 AUGUSTUS result annotation files
        - output_dir: Path to EVM-GFF3 AUGUSTUS result annotation files
    '''
    # Loop the main function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename
        output_name = output_name + '_EVM.gff3'
    
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files
        convert_augustus_to_EVM(Filepath_input, Filepath_output)
        print("Saved_" + output_name)
