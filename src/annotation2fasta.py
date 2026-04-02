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
			      check_file_path,
                              create_out_dir)


def main():
    parser = argparse.ArgumentParser(description='AnnCX feature annotate2fasta to convert annotations (GFF3) to genomic sequences (FASTA)')
    
    parser.add_argument('--annotation', type=check_dir_path, required=True, help='Directory where the annotation file (GFF3) to be converted to FASTA is located. Example: --annotation /path/to/annotation') # Directory
    parser.add_argument('--genome', type=check_dir_path, required=True, help='Directory where the genome or genomic region of interest file (FASTA) that was annotated is located. Example: --genome /path/to/genome') # Directory
    parser.add_argument('--txtgenome', type=check_file_path, required=True, help='TXT file containing the genomes to be processed with each name in one line. Example: --txtgenome /path/to/file.txt') # File
    parser.add_argument('--nameproject', type=str, required=True, help='Name of the project. Example: --nameproject NKG2') # String
    parser.add_argument('--outdir', type=check_dir_path, required=True, help='Directory to save the output files. Example: --outdir /path/to/output') # Directory
    
    args = parser.parse_args()    

    # Make base output directory 
    output_dir_arg = os.path.abspath(args.outdir) # get absolute path
    output_dir = create_out_dir(output_dir_arg, f'annotation2fasta_{args.nameproject}', path=True)

    print ("Convert annotation to fasta")
    genome_txt_dir = create_out_dir(output_dir, 'fasta_sequences')
    get_fasta_sequences_annotation(os.path.abspath(str(args.txtgenome)), 
                                   args.nameproject,
                                   os.path.abspath(str(args.genome)),
                                   os.path.abspath(str(args.annotation)),
                                   str(genome_txt_dir))


if __name__ == '__main__':
    main()
