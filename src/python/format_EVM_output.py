#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 04:45:45 2026

@author: lahumada
"""

"""
Convert EVM to AnnCX format GFF3
"""

import re
import os

def convert_EVM_to_AnnCX_gff3(input_file, output_file):
    """Convert EVM format to AnnCX gene annotation format."""
    
    gene_counter = 0
    cds_counter = 0
    current_gene_id = None
    current_mrna_id = None
    first_header_written = False
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.rstrip()
            
            # Preserve the first ## header line exactly
            if line.startswith('##') and not first_header_written:
                outfile.write(line + '\n')
                first_header_written = True
                continue
            
            # Skip empty lines and other comments
            if not line or (line.startswith('#') and not line.startswith('##')):
                continue
            
            parts = line.split('\t')
            if len(parts) < 9:
                continue
                
            seqid, source, feature, start, end, score, strand, phase, attributes = parts
            source = "AnnCX"
                        
            # New gene - assign counter
            if feature == 'gene':
                gene_counter += 1
                cds_counter = 0
                current_gene_id = f"gene_{gene_counter}"
                new_attrs = f"ID={current_gene_id};Name={current_gene_id}"
                parts[1] = source
                parts[8] = new_attrs
                outfile.write('\t'.join(parts) + '\n')
            
            # mRNA
            elif feature == 'mRNA':
                current_mrna_id = f"{current_gene_id}.mrna1"
                new_attrs = f"ID={current_mrna_id};Parent={current_gene_id};Name={current_gene_id}"
                parts[1] = source
                parts[8] = new_attrs
                outfile.write('\t'.join(parts) + '\n')
            
            # exon
            elif feature == 'exon':
                exon_num_match = re.search(r'exon(\d+)', attributes)
                exon_num = int(exon_num_match.group(1)) if exon_num_match else 1
                exon_id = f"{current_mrna_id}.exon{exon_num}"
                new_attrs = f"ID={exon_id};Name={current_gene_id};Parent={current_mrna_id}"
                parts[1] = source
                parts[8] = new_attrs
                outfile.write('\t'.join(parts) + '\n')
            
            # CDS
            elif feature == 'CDS':
                cds_counter += 1
                cds_id = f"{current_mrna_id}.cds{cds_counter}"
                new_attrs = f"ID={cds_id};Name={current_gene_id};Parent={current_mrna_id}"
                parts[1] = source
                parts[8] = new_attrs
                outfile.write('\t'.join(parts) + '\n')
    

def convert_EVM_to_AnnCX_gff3_all_files(input_dir, output_dir):
    '''Loop over EVM gff3 files'''    
    for filename in os.listdir(input_dir):
        output_name = filename.rsplit('.', 1)[0] + '_AnnCX.gff3'
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
        
        convert_EVM_to_AnnCX_gff3(Filepath_input, Filepath_output)
        print(f"Formatted: {output_name}")
        