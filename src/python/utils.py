#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 23:52:34 2025

@author: lahumada
"""

## Other functions needed along the pipeline execution

import os
import shutil
from pathlib import Path
import argparse


# Function to check whether the overlaps value is an int [0:7]
def check_overlaps_range(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Number of overlaps: {value} is not an integer")

    if ivalue < 0 or ivalue > 9:
        raise argparse.ArgumentTypeError(f"Number of overlaps: {value} is out of allowed range [0, 9]")
    return ivalue

# Function to check if a string is a directory
def check_dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

# Function to check if the argument passed is a directory
def check_file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)

## Function to create output directories 
def create_out_dir (base_dir, subdirs, path=False):
    '''Create an output directory converting the base_dir into a Path object if necessary
    If subdirs=='7_consensus_EVM' and exists, delete it first then recreate.'''
    if path:
        output_dir = Path(base_dir) / subdirs
    else:
        output_dir = base_dir / subdirs
    
    # SPECIAL CASE: Delete '2_run' (which is within '7_consensus_EVM') if it exists
    if str(subdirs) == '2_run' and output_dir.exists():
        shutil.rmtree(output_dir)
        print(f"Deleted existing: {output_dir}")
    
    # Create directory (fresh)
    output_dir.mkdir(parents=True, exist_ok=True)
        
    return output_dir

# Function to print the content of a TXT file
def print_results(path_to_file):
    '''Print the content of a TXT file'''
    with open(path_to_file, 'r') as file:
        print(file.read())
        
# Skippable prompt
def user_continue_prompt():
    while True:
        response = input("\nDo you want to continue with the annotation of the genomes with both flanking genes? (yes/no): ").lower().strip()
        if response in ['yes', 'y']:
            return True
        elif response in ['no', 'n']:
            return False
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")
