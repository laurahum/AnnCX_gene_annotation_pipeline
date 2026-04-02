#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 19:49:51 2026

@author: lahumada
"""


'''
This function takes as input a whole genome FASTA file with multiple entries
and processes it to extract each individual FASTA entry and saves them to individual
FASTA files. This is done for each genome name in the genome_file argument, which
is a TXT file that has one genome name per line.

main function: get_entries_fasta(input_file, output_dir, genome_file)
'''

import os
from pathlib import Path
from Bio import SeqIO

def get_entries_fasta(input_dir, output_dir, genome_file):
    '''
    Extracts entries (>) from a FASTA file and creates individual FASTA file per entry
    for each genome in genome_file.
    Args:
        - input_dir (str): Path to directory with multi-entry FASTA files to be processed
        - output_dir (str): Path to directory to store the extracted FASTA entries 
        - genome_file (str): Path to TXT file with the names of the genomes to be processed
    '''
    print("Processing ...")
    with open(genome_file, 'r') as file:
        genome_list = [line.strip() for line in file if line.strip()]

    for genome in genome_list:
        input_dir = Path(input_dir)
        input_file = [p for p in input_dir.iterdir() if f"{genome}" in p.name]
        
        output_genome_dir = os.path.join(output_dir, genome)
        output_genome_dir = Path(output_genome_dir)
        output_genome_dir.mkdir(parents=True, exist_ok=True)
        
        for i, record in enumerate(SeqIO.parse(str(input_file[0]), "fasta")):
            # Use record.id for meaningful names, write to output_dir
            output_file = os.path.join(output_genome_dir, f"{record.id}.fasta")
            with open(output_file, "w") as out_handle:  
                SeqIO.write(record, out_handle, "fasta")
    print(f"Wrote scaffolds to {output_dir}")
            