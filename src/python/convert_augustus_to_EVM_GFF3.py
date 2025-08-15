#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 22:20:28 2024

@author: lahumada
"""


# This script converts raw GFF3 output from AUGUSTUS to a GFF3 format compatible with EVM
#
# 1. Convert features
# The output from AUGUSTUS contains annotations for features:
#     gene, transcript, stop_codon, CDS, start_codon
# from all those features EVM (gene predictions) accepts only:
#     gene, mRNA, exon and cds
#
# - transcript in the augustus output will be now mRNA,
# - Since the AUGUSTUS output file lacks exon features, the CDS features will be used to create
# exon features with the same coordinates. This is the approach that the EVM
# example script to convert AUGUSTUS output takes when applied to AUGUSTUS output that lacks
# exon features (augustus_GTF_to_EVM_GFF3.pl)
# - Omit features stop_codon and start_codon in the EVM-compatible file
#
# 2. Modify attributes column
# The attributes column of the AUGUSTUS raw GFF3 files must be modified according to the example provided by EVM.
# The format of this column will still reflect the gene models.
#
# The following example file provided by EVM documentation was used as a guide to make the formatting changes
# https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/gene_predictions.gff3
#
# Input: AUGUSTUS GFF3 (or GTF-like) output (with features: gene, transcript, stop_codon, CDS, start_codon)
# Output: AUGUSTUS EVM-compatible GFF3 file (with features: gene, mRNA, exon, cds)
#
# Note: source is omited as an argument and taken from the second column of Exonerate GFF files (AUGUSTUS)

# Main function: convert_augustus_to_EVM_all_files(input_dir, output_dir)



import os

# Main function to convert the format
def convert_augustus_to_EVM(input_file, output_file):
    """
    Convert AUGUSTUS raw GFF3 output format to GFF3 format compatible with EVM.
    gene, transcript, CDS -> gene, mRNA, exon and cds
        - input_file: Path to the raw GFF3 AUGUSTUS result annotation file
        - output_file: Path to the EVM GFF3 AUGUSTUS result annotation file
    """
     
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        
        # Counter for generating unique gene IDs       
        gene_id = 1
        # None = last gene model finished writing
        current_gene = None 
        
        # Lists to accumulate the exon and cdss for each gene model
        mRNA = []
        exons = []
        cdss = []

        for line in infile:
            if line.strip() == "###":
                # End of a gene model: write it to output file
                if current_gene:
                    write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id)
                    gene_id += 1
                    # outfile.write("\n") # Separate gene models with an empty line
                # Reset for next gene model
                current_gene = None 
                mRNA = []
                exons = []
                cdss = []
            else:
                # Process entry lines in the input annotation file
                fields = line.strip().split('\t')
                if len(fields) == 9:  # Annotation entry
                    seqname, source, feature, start, end, score, strand, frame, attributes = fields
                    if feature == 'gene':
                        # Write previous model in case it is still not written
                        if current_gene:
                            write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id)
                            gene_id += 1
                            # outfile.write("\n")
                        # Start a new gene model
                        current_gene = fields
                    # Accumulate mRNA feature 
                    elif feature == 'transcript':
                        mRNA.append(fields)
                    # Accumulate exon and cds features in the lists
                    elif feature == 'CDS':
                        cdss.append(fields)
                        exons.append(fields)  # The exon feature have the same coordinates as cds

        # Write the last gene if exists, in case the file does not end with '###'
        if current_gene:
            write_gene(outfile, current_gene, mRNA, exons, cdss, gene_id)


# Function to write a gene model
def write_gene(outfile, gene, mRNA, exons, cdss, gene_id):
    """
    Write a single gene model in GFF3 format compatible with EVM.
        - outfile (file): Output file object
        - gene (list): Gene feature fields
        - exons (list): List of exon feature fields
        - cdss (list): List of CDS feature fields
        - gene_id (int): Unique identifier for the gene
    """
    
    seqname, source, feature, start, end, score, strand, frame, attributes = gene
    # gene_name = f"g{gene_id}"
    
    # Write gene feature
    outfile.write(f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t.\tID={seqname}.tPRED{gene_id:06d};Name={source} model {seqname}.m{gene_id:06d}\n")
    
    # Write mRNA feature
    for i, mRNA in enumerate(mRNA, 1):
        seqname, source, feature, start, end, score, strand, frame, attributes = mRNA
        outfile.write(f"{seqname}\t{source}\tmRNA\t{start}\t{end}\t{score}\t{strand}\t.\tID={seqname}.m000001;Parent={seqname}.tPRED{gene_id:06d}\n")
    
    # Write exon features
    for i, exon in enumerate(exons, 1):
        seqname, source, feature, start, end, score, strand, frame, attributes = exon
        outfile.write(f"{seqname}\t{source}\texon\t{start}\t{end}\t{score}\t{strand}\t.\tID={seqname}.e{i:06d};Parent={seqname}.m000001\n")
    
    # Write CDS features
    for i, cds in enumerate(cdss, 1):
        seqname, source, feature, start, end, score, strand, frame, attributes = cds
        outfile.write(f"{seqname}\t{source}\tCDS\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID=cds_of_{seqname}.m000001;Parent={seqname}.m000001\n")

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
        print("Saved_" + output_name + '\n')
