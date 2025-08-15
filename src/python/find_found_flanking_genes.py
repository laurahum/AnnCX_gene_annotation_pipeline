#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 00:17:06 2024

@author: lahumada
"""

## This script takes as input a fasta file with two entries (the flanking genes)
## and the gff3 file output from GMAP and finds out if both entries in the fasta
## file were found in the gff3 file

## Main function: find_found_flanking_genes(fasta_input, gff3_input, txt_output)

## To-do: the output list of 'txt_output_yes.txt' in alphabetical order


import os
import re
from Bio import SeqIO



# Function to get fasta entries as a list
def read_fasta_entries(fasta_input):
    """Read entries from a FASTA file and return them as a list."""
    entries = []
    with open(fasta_input, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            entries.append(record.id)
    return entries

# Function to check if both entries are found in the GFF3 file
def check_entries_in_gff3(fasta_entries, gff3_file):
    """Check if both FASTA entries are present in the GFF3 file."""
    found_entries = set()
    with open(gff3_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                # Skip comment lines
                continue  
            fields = line.strip().split("\t")
            if len(fields) < 9:
                # Skip malformed lines
                continue  
            attributes = fields[8]
            for entry in fasta_entries:
                if entry in attributes:
                    found_entries.add(entry)
    return found_entries


# Function to write to file whether the genes were found
def write_to_txt(fasta_input, gff3_file, genome_name, txt_output_yes, txt_output_no):
    '''
    Write to a TXT file whether both flanking genes are found.
        - txt_output_yes = both genes were found 
        - txt_output_no = only one or no genes were found
    '''
    fasta_entries = read_fasta_entries(fasta_input)
    if len(fasta_entries) != 2:
        print("Error: The FASTA file with flanking genes should contain exactly two entries.")
        return

    found_entries = check_entries_in_gff3(fasta_entries, gff3_file)
    
    # List of the genome that do have both genes in the GFF3 file
    if len(found_entries) == 2:
        with open(txt_output_yes, 'a+') as file:
            file.write(genome_name)
            file.write("\n")
    
    # List of the genome that do not have both genes in the GFF3 file 
    elif len(found_entries) == 1:
        with open(txt_output_no, 'a+') as file:
            file.write(genome_name)
            # Reason:
            file.write(f"\nOnly one flanking gene ({found_entries.pop()}) was found in the GFF3 file.\n")
    else:
        with open(txt_output_no, 'a+') as file:
            file.write(genome_name)
            # Reason:
            file.write("\nNeither flanking gene was found in the GFF3 file.\n")


def find_found_flanking_genes(fasta_input, gff3_input, txt_output, name_genes):
    '''
    Find if both flanking genes are found by GMAP
        - fasta_input = FASTA file with two entries corresponding to the flanking genes
        - gff3_input = GFF3 file produced by GMAP with the annotations corresponding to the flanking genes 
        - txt_output = directory to write the results of this search
        	- txt_output_yes: List of the genome that do have both genes in the GFF3 file
        	- txt_output_no: List of the genome that do not have both genes in the GFF3 file (documented genome with only one flanking gene found or with none)

	- name_genes = name of the gene type to be annotated       	
    Return (txt_output_yes, txt_output_no): Paths to the files with the results of the search
    '''
    # Output file paths 
    txt_output_yes = os.path.join(txt_output, 'txt_output_yes.txt')
    txt_output_no = os.path.join(txt_output, 'txt_output_no.txt')
    
    # Create empty output files
    with open(txt_output_yes, 'w') as f_yes:
    	pass

    with open(txt_output_no, 'w') as f_no:
    	pass 
 
    # Loop the function over the GFF3 files created by gmap (genome vs ROI flanking regions)
    for file in os.listdir(gff3_input):
        file_name = file 
 	
 	# Get genome name
        pattern = rf'^(.*?)_{name_genes}_flanking'
        match = re.search(pattern, file_name)
        if match:
            genome_name = match.group(1)
    
        # GFF3 from each genome in the directory
        gff3_file = os.path.join(gff3_input, file)
    
        # Write to txt
        write_to_txt(fasta_input, gff3_file, genome_name, txt_output_yes, txt_output_no)
        
    return txt_output_yes, txt_output_no
