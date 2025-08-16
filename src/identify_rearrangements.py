#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 06:21:55 2025

@author: lahumada
"""

# Standard library imports
import argparse
import os
import sys
from pathlib import Path

# Adjust Python path (preliminary)
sys.path.append(str(Path(__file__).parent.parent))

# Local application imports
from src.python import (run_bash_script, run_R_script)

# Specific utility imports
from src.python.utils import (check_dir_path,
                              check_file_path,
                              create_out_dir)


def main():
    parser = argparse.ArgumentParser(description='AnnCX feature identify_rearrangements to identify exon-level rearrangements')
    
    parser.add_argument('--subject', type=check_file_path, required=True, help='File (FASTA) with reference exon sequences (e.g header: > NKG2_exon_1). Example: --subject /path/to/reference/reference.fasta') # File
    parser.add_argument('--query', type=check_file_path, required=True, help='File (FASTA) with predicted exon sequences (e.g header: > NKG2_exon_1). Example: --query /path/to/query/query.fasta') # File
    parser.add_argument('--namegenes', type=check_file_path, required=True, help='File (TXT) with the names of the predicted genes. These names must be contained in the headers of the FASTA file for predicted exon sequences. Example: --namegenes /path/to/genes/genes.txt') # String
    parser.add_argument('--namegenome', type=str, required=True, help='Name of the genome in which the genes were predicted. Example: --namegenome Homo_sapiens') # String
    parser.add_argument('--outdir', type=check_dir_path, required=True, help='Directory to save the output files. Example: --outdir /path/to/output') # Directory
    parser.add_argument('--threads', default=-1, type=int, help='(OPTIONAL) Number of threads that can be used to run this feature. Example: --threads 3') # Int
    
    args = parser.parse_args()

    # Make base output directory 
    output_dir_arg = os.path.abspath(args.outdir) # get absolute path
    output_dir = create_out_dir(output_dir_arg, 'Identify_rearrangements', path=True)
    exonerate_dir = create_out_dir(output_dir, 'Exonerate_output')

    print ("1. Run Exonerate")
    exonerate_output = run_bash_script('exonerate_identify_artificial_rearrangements.sh', 'exonerate_run_identify_rearrangements', 
                                       os.path.abspath(args.subject),
                                       os.path.abspath(args.query),
                                       exonerate_dir,
                                       args.namegenome,
                                       args.threads)
    
    print ("2. Make heatmaps")
    exonerate_heatmaps_dir = create_out_dir(output_dir, 'Heatmaps')
    exonerate_file_path = exonerate_output.strip()

    run_R_script('Exonerate_heatmap_identify_artificial_rearrangements.R', 'Exonerate_heatmaps',
                 str(exonerate_file_path),
                 os.path.abspath(args.namegenes),
                 str(exonerate_heatmaps_dir),
                 args.namegenome)


if __name__ == '__main__':
    main()
