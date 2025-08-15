#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 03:13:35 2025

@author: lahumada
"""

import os
import re


def create_error_genome_name_report (errors, output_dir):
    ''' 
    Writes to a file all the errors found with the input genomes 
    (no FASTA files found, names have special characters, too long names)
    - errors (list): Error messages created in function get_name_fasta_files()
    - output_dir (str): Directory to save the error report
    '''
    output_file = os.path.join(output_dir, "genome_input_error.txt")
    
    with open(output_file, 'w') as file: 
        for err in errors:
            file.write(f"{err}\n")


def get_name_fasta_files(input_dir, output_dir, genomes_error_dir):
    '''
    Processes FASTA files in a given directory and creates a text file with their names (without extensions).
    Checks for special characters or dots and length (>60) in names. If any problems are found,
    prints errors and asks user to fix them.
    - input_dir (str): Directory containing FASTA files
    - output_dir (str): Directory to save the output file
    - genomes_error_dir (str): Directory to save errors found withe the input genomes
    '''
    
    # Number of characters allowed:
    num_char = 50
    
    # List to store file names without extensions
    file_names = []
    errors = [] # for any issues with the input genomes
    
    # Regex pattern to match fasta and adjacent endings
    fasta_pattern = re.compile(r'\.fasta$|\.fa$|\.fna$|\.ffn$|\.faa$|\.frn$', re.IGNORECASE)
    allowed_pattern = re.compile(r'^[A-Za-z0-9_-]+$')
    
    # Iterate through files in the folder
    for filename in os.listdir(input_dir):
        if fasta_pattern.search(filename):
            # Remove the extension
            name_without_extension = os.path.splitext(filename)[0]

            # Check for special characters or dot (.)
            if not allowed_pattern.match(name_without_extension):
                errors.append(f"[ERROR] Invalid characters in file name: '{filename}'. Remove special characters.")

            # Check for length > 60
            if len(name_without_extension) > num_char:
                errors.append(f"[ERROR] File name too long: '{filename}' ({len(name_without_extension)} characters). Maximum allowed is {num_char}.")

            file_names.append(name_without_extension)
        else:
            errors.append("No FASTA files found in input genome directory")

    if errors:
        create_error_genome_name_report(errors, genomes_error_dir)

    
    # Write file names to a TXT file
    output_file = os.path.join(output_dir, 'genome_file_names.txt')
    with open(output_file, 'w') as f:
        for name in file_names:
            f.write(f"{name}\n")
    
    print(f"Processed {len(file_names)} FASTA files. Results written to {output_file}")
    
    return output_file
