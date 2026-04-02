#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 22:56:25 2025

@author: lahumada
"""

## Custom read of GFF3 files and transform into a pandas dataframe

## Main function: read_gff3(file_path)


import pandas as pd

# Function to read the gff3 file and return a pandas dataframe
def read_gff3(file_path):
    '''
    Function to read a gff3 file and return a pandas dataframe
    Arg: - file_path: PATH to GFF3 file
    Return: pandas dataframe
    '''
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Keep lines that contain '###'
            if '###' in line:
                data.append(line)
            elif not line.startswith('#'): 
                # Skip comment lines
                fields = line.split('\t')
                # Only the fields that correspond to an annotation (have 9 fields)
                if len(fields) == 9:  
                    data.append(fields)
    
    # Create DataFrame only if there are data rows
    if not data:
        return None
    
    return pd.DataFrame(data)
