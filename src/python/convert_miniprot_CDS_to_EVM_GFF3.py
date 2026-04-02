#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 05:42:23 2026

@author: lahumada
"""


'''
This script converts miniprot GFF3 output to a format compatible with EVM

Input: miniprot GFF3 output (with features: mRNA, CDS, stop_codon) 
Output: miniprot EVM-compatible GFF3 file (with features: nucleotide_to_protein_match)

Features to convert:
    CDS -> nucleotide_to_protein_match (mRNA provides grouping for unique IDs) for protein alignment

Main function: convert_miniprot_CDS_to_EVM_all_files(input_dir, output_dir, source)
Note: source must be 'miniprotCDS'
'''

import os
import re


def convert_miniprot_to_EVM(input_file, output_file, source):
    """
    Convert miniprot GFF3 output to GFF3 format compatible with EVM.
    Only CDS -> nucleotide_to_protein_match, using mRNA structure for ID grouping
        - input_file (str): Path to the input miniprot output file
        - output_file (str): Path to the output GFF3 file compatible with EVM
        - source: 'miniprotCDS'
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Counter for generating unique gene IDs
        gene_id = 1  
        # Store the current mRNA being processed
        current_mrna = None  
        # List to accumulate CDS features only
        cds_features = []

        for line in infile:
            line = line.strip()
            
            # Skip empty lines or comments
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            # Check if it is a valid annotation entry
            if len(fields) == 9:  
                seqname, source_field, feature, start, end, score, strand, frame, attributes = fields
                
                if feature == 'mRNA':  
                    # If there is a previous mRNA, write its CDS
                    if current_mrna and cds_features: 
                        write_mrna_cds(outfile, current_mrna, cds_features, gene_id, source)
                        gene_id += 1
                    # Start a new mRNA model
                    current_mrna = fields  
                    # Reset CDS list for new mRNA
                    cds_features = []
                    
                elif feature == 'CDS':
                    cds_features.append(fields)

        # Write the last mRNA if it exists
        if current_mrna and cds_features:
            write_mrna_cds(outfile, current_mrna, cds_features, gene_id, source)


def write_mrna_cds(outfile, mrna, cds_features, gene_id, source):
    """
    Write all CDS features of an mRNA as nucleotide_to_protein_match entries with shared gene ID.
    """
    # Get mRNA data
    seqname, _, feature, start, end, score, strand, frame, attributes = mrna
    
    # Extract target name from mRNA attributes (Target=rhKLRC1(NKG2A))
    target_match = re.search(r'Target=([^;]+)', attributes)
    target_name = target_match.group(1) if target_match else f'prot_{gene_id}'

    # Write ONLY CDS as nucleotide_to_protein_match
    for cds in cds_features:
        seqname_t, source_t, feature_t, t_start, t_end, t_score, t_strand, t_frame, t_attributes = cds
        
        # Extract Identity from CDS attributes
        identity_match = re.search(r'Identity=([\d.]+)', t_attributes)
        identity = identity_match.group(1) if identity_match else '0.0'
        
        # Use CDS score from column 6
        score_out = t_score if t_score and t_score != '.' else '.'
        
        outfile.write(f"{seqname_t}\t{source}\tnucleotide_to_protein_match\t{t_start}\t{t_end}\t{score_out}\t{t_strand}\t{t_frame}\tID=match.{source}.{gene_id};Target={target_name};Identity={identity}\n")


def convert_miniprot_CDS_to_EVM_all_files(input_dir, output_dir, source):
    """
    Process all miniprot GFF3 files in input directory.
    - source: 'miniprotCDS'
    """    
    # Loop the main function over the input files
    for filename in os.listdir(input_dir):
        # Get name for each input file and modify it for the output file names:
        output_name = filename
        output_name = output_name.rsplit('.', 1)[0] + '_EVM.gff3'
        
        # Directory to each input and output file:
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
    
        # Call function to process input files
        convert_miniprot_to_EVM(Filepath_input, Filepath_output, source)
        print("Saved_" + output_name)
      