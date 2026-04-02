#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 20:56:45 2026
@author: lahumada
"""
"""
Convert minimap2 PAF to GFF3 gene models with gene/exon structure
- Parses CIGAR string to identify exons (M/D/I blocks separated by N introns) 
- Creates COMPLETE gene → mRNA → exon hierarchy PER PAF ENTRY
- Score in column 6 of GFF3 in mRNA entry = cm
- Stores all original PAF info in attributes of the corresponding mRNA entry

## EXTRACTED PAF FIELDS (mRNA attributes):
CORE FIELDS:
- qname, qlen, qstart, qend, tstart, tend, nmatch, alen, mapq

MINIMAP2 TAGS:
- NM:i - Edit distance (mismatches + gaps)
- ms:i - Matching score (chaining phase)  
- AS:i - Alignment score (DP score)
- nn:i - Number of ambiguous bases
- ts:A - Transcript strand (+/-)
- tp:A - Transcript type (P=protein, S=spliced)
- cm:i - Conserved minimizer count (Also in GFF3 score column 6)
- s1:i - Query minimizer count
- s2:i - Target minimizer count
- de:f - Deletion error rate
- dv:f - Sequence divergence
- rl:i - Repetitive seed length
- cg:Z - CIGAR string (full alignment)
- cs:Z - Compression sequence tag (mismatches/INDELs)

Main function: format_minimap2_output_to_gff_all_files(input_dir,output_dir)
"""


import os
import sys
import re

def parse_paf_line(line):
    '''Parse PAF line → dict with core fields + ALL tags'''
    fields = line.strip().split('\t')
    core = fields[:12]
    tags = fields[12:]
    
    record = {
        'qname': core[0], 'qlen': int(core[1]), 'qstart': int(core[2]), 'qend': int(core[3]),
        'strand': core[4], 'tname': core[5], 'tlen': int(core[6]), 
        'tstart': int(core[7]), 'tend': int(core[8]), 'nmatch': int(core[9]), 
        'alen': int(core[10]), 'mapq': int(core[11])
    }
    
    # **EXTRACT ALL TAGS** (not just selective ones)
    for tag in tags:
        if ':' in tag:
            parts = tag.split(':', 2)
            if len(parts) == 3:
                key, typ, value = parts
                record[key] = value
    
    cigar = record.get('cg', '')
    exons = parse_cigar_exons(cigar, record['tstart'], record['strand'])
    record['exons'] = exons
    
    return record

def parse_cigar_exons(cigar, tstart, strand):
    '''Parse CIGAR → list of exon coordinates (start, end)'''
    exons = []
    pos = int(tstart)
    
    for match in re.finditer(r'(\d+)([MDINSX=])', cigar):
        length = int(match.group(1))
        op = match.group(2)
        
        if op in 'MDIX=':  # Exonic regions
            ex_start = pos + 1  # GFF3 1-based
            ex_end = pos + length
            exons.append((ex_start, ex_end))
            pos += length
        elif op == 'N':  # Intron
            pos += length
    
    return exons

def format_minimap2_output_to_gff(Filepath_input, Filepath_output):
    '''Convert PAF → COMPLETE GFF3 gene models with ALL PAF info'''
    with open(Filepath_input, 'r') as file:
        all_gene_models = []
        
        for line_num, line in enumerate(file, 1):
            line = line.rstrip('\n')
            if line.startswith('#') or not line.strip():
                continue
                
            try:
                paf = parse_paf_line(line)
            except:
                print(f'# Warning: skipped malformed line {line_num}', file=sys.stderr)
                continue
            
            # Build ONE complete gene model from this PAF line
            tname = paf['tname']
            gene_name = paf['qname'].rsplit('_', 1)[0] if '_' in paf['qname'] else paf['qname']
            start = min(paf['tstart'], paf['tend']) + 1
            end = max(paf['tstart'], paf['tend'])
            strand = '+' if paf['strand'] == '+' else '-'
            score = int(paf.get('cm', 0))  # Use cm as GFF score
            
            gene_model = []
            
            # 1. GENE feature
            gene_id = f'{gene_name}.gene.{line_num}'
            gene_model.append(
                f'{tname}\tminimap2\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id};Name={gene_name}'
            )
            
            # 2. mRNA feature - **ALL PAF TAGS AS ATTRIBUTES (BASE-1 CONVERSION)**
            mrna_id = f'{gene_name}.mRNA.{line_num}'

            # Core PAF fields - **CONVERT base-0 → base-1 for GFF3**
            qstart_1based = paf['qstart'] + 1
            qend_1based = paf['qend']
            tstart_1based = paf['tstart'] + 1  
            tend_1based = paf['tend']

            core_attrs = [
                f'ID={mrna_id}', f'Parent={gene_id}',
                f'qname={paf["qname"]}', f'qlen={paf["qlen"]}',
                f'qstart={qstart_1based}', f'qend={qend_1based}',      # Base-1 labeled
                f'tstart={tstart_1based}', f'tend={tend_1based}',       # Base-1 labeled
                f'nmatch={paf["nmatch"]}', f'alen={paf["alen"]}',
                f'mapq={paf["mapq"]}'
            ]
            
            # **ALL TAG ATTRIBUTES** (dynamic extraction)
            tag_attrs = []
            for key, value in paf.items():
                if key not in ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 
                              'tstart', 'tend', 'nmatch', 'alen', 'mapq', 'exons']:
                    tag_attrs.append(f'{key}={value}')
            
            # Combine all attributes
            all_attrs = core_attrs + tag_attrs
            attributes = ';'.join(all_attrs)
            
            gene_model.append(
                f'{tname}\tminimap2\tmRNA\t{start}\t{end}\t{score}\t{strand}\t.\t{attributes}'
            )
            
            # 3. EXON features (sorted by genomic position)
            exons_sorted = sorted(paf['exons'], key=lambda x: x[0])
            for rank, (ex_start, ex_end) in enumerate(exons_sorted, 1):
                exon_id = f'{mrna_id}.exon.{rank}'
                exon_attrs = f'ID={exon_id};Parent={mrna_id};Rank={rank}'
                gene_model.append(
                    f'{tname}\tminimap2\texon\t{ex_start}\t{ex_end}\t.\t{strand}\t.\t{exon_attrs}'
                )
            
            all_gene_models.append(gene_model)
    
    # Write gene models in genomic order
    all_gene_models.sort(key=lambda x: int(x[0].split('\t')[3]))
    
    with open(Filepath_output, 'w') as output_file:
        output_file.write('##gff-version 3\n')
        for gene_model in all_gene_models:
            for line in gene_model:
                output_file.write(line + '\n')

def format_minimap2_output_to_gff_all_files(input_dir, output_dir):
    '''Loop over minimap2 PAF files → properly separated GFF3 gene models'''    
    for filename in os.listdir(input_dir):
        output_name = filename.rsplit('.', 1)[0] + '_FORMATTED.gff3'
        Filepath_input = os.path.join(input_dir, filename)
        Filepath_output = os.path.join(output_dir, output_name)
        
        format_minimap2_output_to_gff(Filepath_input, Filepath_output)
        print(f"Saved: {output_name}")
