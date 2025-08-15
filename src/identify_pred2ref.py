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
from src.python import (run_bash_script, run_R_script)

# Specific utility imports
from src.python.utils import (check_dir_path,
                              check_file_path,
                              create_out_dir)


def main():
    parser = argparse.ArgumentParser(description='AnnCX feature identify_pred2ref to identify predicted genes')
    
    parser.add_argument('--subject', type=check_file_path, required=True, help='File (FASTA) with reference gene sequences (e.g. gene, exon_all, CDS_all, cDNA). Example: --subject /path/to/reference/reference.fasta') # File
    parser.add_argument('--query', type=check_file_path, required=True, help='File (FASTA) with predicted gene sequences (e.g. gene, exon_all, CDS_all, cDNA). Example: --query /path/to/query/query.fasta') # File
    parser.add_argument('--typeseq_query', type=str, required=True, help='Name of the type of query fasta sequences. Example: --typeseq_query cDNA') # String
    parser.add_argument('--typeseq_subject', type=str, required=True, help='Name of the type of subject fasta sequences. Example: --typeseq_subject cDNA') # String
    parser.add_argument('--namegenome', type=str, required=True, help='Name of the genome in which the genes were predicted. Example: --namegenome Homo_sapiens') # String
    parser.add_argument('--outdir', type=check_dir_path, required=True, help='Directory to save the output files. Example: --outdir /path/to/output') # Directory
    parser.add_argument('--threads', default=-1, type=int, help='(OPTIONAL) Number of threads that can be used to run this feature. Example: --threads 3') # Int
    
    args = parser.parse_args()

    # Make base output directory 
    output_dir = create_out_dir(args.outdir, 'Identify_pred2ref', path=True)
    blastn_dir = create_out_dir(output_dir, 'BLASTN_output')

    print ("1. Run BLASTN")
    blastn_output = run_bash_script('blastn_identify_predicted_genes.sh', 'blastn_run_identify_gene', 
                                    args.subject,
                                    args.query,
                                    blastn_dir,
                                    args.typeseq_query,
                                    args.typeseq_subject,
                                    args.namegenome,
                                    args.threads)
    
    print ("2. Make heatmaps")
    blastn_heatmaps_dir = create_out_dir(output_dir, 'Heatmaps')
    blast_file_path = blastn_output.strip()    
    
    run_R_script('BLASTN_heatmap_identify_predicted_annotations.R', 'BLASTN_heatmaps',
                 str(blast_file_path),
                 str(blastn_heatmaps_dir),
                 args.namegenome,
                 args.typeseq_query,
                 args.typeseq_subject)


if __name__ == '__main__':
    main()
