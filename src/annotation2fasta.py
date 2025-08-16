#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 03:23:08 2025

@author: lahumada
"""

# Standard library imports
import argparse
import sys
import os
from pathlib import Path

# Adjust Python path (preliminary)
sys.path.append(str(Path(__file__).parent.parent))

# Local application imports
from src.python import get_fasta_sequences_annotation

# Specific utility imports
from src.python.utils import (check_dir_path,
                              create_out_dir)


def main():
    parser = argparse.ArgumentParser(description='AnnCX feature annotate2fasta to convert annotations (GFF3) to genomic sequences (FASTA)')
    
    parser.add_argument('--annotation', type=check_dir_path, required=True, help='Directory where the annotation file (GFF3) to be converted to FASTA is located. Example: --annotation /path/to/annotation') # Directory
    parser.add_argument('--genome', type=check_dir_path, required=True, help='Directory where the genome or genomic region of interest file (FASTA) that was annotated is located. Example: --genome /path/to/genome') # Directory
    parser.add_argument('--namegenome', type=str, required=True, help='Name of the genome or genomic region of interest annotated. This name (ex: Macaca_mulatta) must be contained in the name of the genome (ex: Macaca_mulatta_genome.fasta) and annotation (ex. Macaca_mulatta_annotation.gff3) files. Example: --namegenome Macaca_mulatta') # String
    parser.add_argument('--nameproject', type=str, required=True, help='Name of the project. Example: --nameproject NKG2') # String
    parser.add_argument('--outdir', type=check_dir_path, required=True, help='Directory to save the output files. Example: --outdir /path/to/output') # Directory
    
    args = parser.parse_args()    

    # Make base output directory 
    output_dir_arg = os.path.abspath(args.outdir) # get absolute path
    output_dir = create_out_dir(output_dir_arg, f'annotation2fasta_{args.nameproject}', path=True)
    
    # Create a TXT file with the name of the genome (necessary due to the code of annotation2fasta.py)
    genome_txt_dir = create_out_dir(output_dir, 'genome_txt')
    genome_txt_file = os.path.join(genome_txt_dir, 'genome_list.txt')

    with open (genome_txt_file, 'w') as file:
        file.write(args.namegenome + '\n')
    print (f"Created txt file {genome_txt_file}")

    print ("Convert annotation to fasta")
    genome_txt_dir = create_out_dir(output_dir, 'fasta_sequences')
    get_fasta_sequences_annotation(str(genome_txt_file), 
                                   args.nameproject,
                                   os.path.abspath(str(args.genome)),
                                   os.path.abspath(str(args.annotation)),
                                   str(genome_txt_dir))


if __name__ == '__main__':
    main()
