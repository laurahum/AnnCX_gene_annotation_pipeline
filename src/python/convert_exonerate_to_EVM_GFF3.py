#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 02:18:09 2024

@author: lahumada
"""

# This script converts formatted filtered GFF output from Exonerate (cDNA) to a format compatible with EVM
#
# 1. Convert features
# The output from Exonerate contains for features:
#     gene, utr5, utr3, exon, splice5, splice3, intron, cds, similarity
# from all those features EVM (transcript alignment) accepts only:
#     cDNA_match
# - all the features (gene, exon, cds and similarity) must be converted to cDNA_match 
# - Omit features utr3, utr5, splice5, splice3, in the EVM-compatible file
#
# 2. Modify attributes column
# The attributes column of the Exonerate formatted filtered GFF3 files must be modified according to the example provided by EVM
# Give the entries associated to the same gene model, in order of entry, a counter in the attributes section:
# ex. ID=match.exonerate:cdna2genome.1 for the features that are associated to the first gene model entry.

# The following example file provided by EVM documentation was used as a guide to make the formatting changes
# https://github.com/EVidenceModeler/EVidenceModeler/blob/master/testing/transcript_alignment.gff3

# Input: Exonerate formatted filtered GFF output (with features: gene, utr5, utr3, exon, splice5, splice3, intron, cds, similarity)
# Output: Exonerate formatted filtered EVM-compatible GFF3 file (with features: DNA_match)
#
# Note: source is omited as an argument and taken from the second column of Exonerate GFF files (exonerate:cdna2genome)

# Main function: convert_exonerate_to_EVM_all_files(input_dir, output_dir)


import os
import re


# Function to convert the format
def convert_exonerate_to_EVM(input_file, output_file):
    """
    Convert Exonerate formatted-filtered-GFF output to GFF3 format compatible with EVM.
    gene, exon, cds, similarity -> cDNA_match
        - input_file (str): Path to the input Exonerate output file
        - output_file (str): Path to the output GFF3 file compatible with EVM
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:

        # Counter for generating unique gene IDs
        gene_id = 1 
        # None = last gene model finished writing
        current_gene = None 
        
        # Lists to accumulate the exon and cdss for each gene model
        similarity = []
        exons = []
        cdss = []

        for line in infile:
            if line.strip() == "###":
                # End of a gene model: write it to output file
                if current_gene:
                    write_gene(outfile, current_gene, similarity, exons, cdss, gene_id)
                    gene_id += 1
                    # outfile.write("\n") # Separate gene models with an empty line
                # Reset for next gene model
                current_gene = None 
                similarity = []
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
                            write_gene(outfile, current_gene, similarity, exons, cdss, gene_id)
                            gene_id += 1
                            # outfile.write("\n")
                        # Start a new gene model
                        current_gene = fields
                    # Accumulate similarity features
                    elif feature == 'similarity':
                        similarity.append(fields)
                    # Accumulate exon and cds features in the lists
                    elif feature == 'exon':
                        exons.append(fields)
                    elif feature == 'cds':
                        cdss.append(fields)

        # Write the last gene if exists, in case the file does not end with '###'
        if current_gene:
            write_gene(outfile, current_gene, similarity, exons, cdss, gene_id)


# Function to write a gene model
def write_gene(outfile, gene, similarity, exons, cdss, gene_id):
    """
    Write all features of a gene model as cDNA_match entries.
        - outfile (file): Output file object
        - gene (list): Gene feature fields
        - mRNA (list): List of mRNA
        - exons (list): list of exons
        - cdss (list): list of cdss
        - gene_id (int): Unique identifier for the gene
    """    
    # Get data for each field
    seqname, source, feature, start, end, score, strand, frame, attributes = gene
    
    # Extract the gene name
    gene_match = re.search(r'sequence\s*[=\s]\s*([^;\s]+)', attributes)
    gene_name = gene_match.group(1)
    
    # Write gene feature                                                                      
    outfile.write(f"{seqname}\t{source}\tcDNA_match\t{start}\t{end}\t{score}\t{strand}\t.\tID=match.{source}.{gene_id};Target={gene_name}\n")
    
    # Write similarity feature
    for i, similarity in enumerate(similarity, 1):
        seqname, source, feature, start, end, score, strand, frame, attributes = similarity  
        outfile.write(f"{seqname}\t{source}\tcDNA_match\t{start}\t{end}\t{score}\t{strand}\t.\tID=match.{source}.{gene_id};Target={gene_name}\n")
    
    # Write exon features
    for i, exon in enumerate(exons, 1):
        seqname, source, feature, start, end, score, strand, frame, attributes = exon  
        outfile.write(f"{seqname}\t{source}\tcDNA_match\t{start}\t{end}\t{score}\t{strand}\t.\tID=match.{source}.{gene_id};Target={gene_name}\n")
    
    # Write CDS features
    for i, cds in enumerate(cdss, 1):
        seqname, source, feature, start, end, score, strand, frame, attributes = cds
        outfile.write(f"{seqname}\t{source}\tcDNA_match\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID=match.{source}.{gene_id};Target={gene_name}\n")


# Main
def convert_exonerate_to_EVM_all_files(input_dir, output_dir):
    '''
    Loop the convert_exonerate_to_EVM function over all the formatted filtered Exonerate result annotation files for each genome
        - input_dir: Path to formatted-filtered-GFF Exonerate result annotation files
        - output_dir: Path to formatted-filtered-EVM-GFF Exonerate result annotation files
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
        convert_exonerate_to_EVM(Filepath_input, Filepath_output)
        print("Saved_" + output_name + '\n')
        
        
