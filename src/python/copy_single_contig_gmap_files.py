#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 22:27:55 2022

@author: lahumada
"""

## Copy to a new directory the GMAP data files (from the alignment of flanking 
# regions to ROI) that correspond only to the genome where the ROI is in the 
# same contig, as per the text file 'Single_contig_genome.txt'

# Main function: copy_single_contig_genomes(input_dir, output_dir, single_contig_list)


import shutil
import os



def copy_single_contig_genomes(input_dir, output_dir, single_contig_list):
    '''
    Copy to a new directory the GMAP GFF3 files, generated from the alignment of
    flanking genes to the FASTA genome files, that correspond to the genome where
    the ROI is in the same contig/scaffold/chromosome.
    	- input_dir: Directory to GMAP GFF3 annotation files with the flanking genes
    	- output_dir: Directory to copy GFF3 annotation files of single contig genome
    	- single_contig_list: TXT file list of genome with flanking genes in a single contig
    '''
    
    # Names of the files within the source directory
    files = os.listdir(input_dir)

    # List of genome with ROI in a single contig
    genome_file = open(single_contig_list,'r')
    genome_list = genome_file.readlines()

    ### Copy single contig genome files into the destination directory

    # Iterate over each element in the list of genome with ROI in single contig
    for genome_name in genome_list:
        # Remove the last character ('\n') on each of the genome name
        genome_name = genome_name.strip()
        print(genome_name)
    
        # Iterate over each element in the list of file names from source directory
        for file_name in files:
            # Only if the name of the file corresponds to a single contig genome
            if (genome_name in file_name):
                # Copy file to the destination directory
                shutil.copy2(os.path.join(input_dir,file_name), output_dir)
                print(file_name + ' copied')
                # Exit this for loop once the single contig genome file is found
                break
 
    
