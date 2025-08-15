#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 04:14:01 2025

@author: lahumada
"""

## General script to run within Python the R scripts in the pipeline.

## Main function: run_R_script(script_name, *args)


import subprocess
import os
from pathlib import Path

def run_R_script(script_name, function_name, *args):
    '''
    Run R script from a Python script using subprocess
     - script_name : str, Name of the R script to be run
     - *args : str, Name of the different arguments required by the R script
    '''
    
    # Get the path to the root directory
    root_dir = Path(__file__).parent.parent.parent
    
    # Construct the full path to the bash script
    script_path = root_dir / 'src' / 'R' / script_name
    
    # Ensures every arg is quoted
    quoted_args = [f'"{arg}"' for arg in args]  
    
    # Run script
    command = ['Rscript', '-e', f'source("{script_path}"); {function_name}({", ".join(quoted_args)})']
    
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return result.stdout
        
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        print(e.stderr)
