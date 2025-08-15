#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 22:41:05 2025

@author: lahumada
"""

## General script to run within Python the Bash scripts in the pipeline.

## Main function: run_bash_script(script_name, *args)


import subprocess
import os
from pathlib import Path

def run_bash_script(script_name, function_name, *args):
    '''
    Run Bash script from a Python script using subprocess
     - script_name : str, Name of the bash script to be run
     - *args : str, Name of the different arguments required by the bash script
    '''
    
    # Get the path to the root directory
    root_dir = Path(__file__).parent.parent.parent
    
    # Construct the full path to the bash script
    script_path = root_dir / 'src' / 'bash' / script_name
    
    # Run script
    bash_command = ['/bin/bash', '-c', f'source {script_path} && {function_name} {" ".join(map(str, args))}']
    
    # Run with live streaming + capture
    process = subprocess.Popen(
    	['/bin/bash', '-c', f'source {script_path} && {function_name} {" ".join(map(str, args))}'],
    	stdout=subprocess.PIPE,
    	stderr=subprocess.STDOUT,
    	text=True,
    	bufsize=1,
    )
    
    # Print and capture output line by line
    output_lines = []
    for line in process.stdout:
        print(line, end='')  # Print live
        output_lines.append(line)
    
    # Wait for completion and check return code
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, bash_command)
    
    return ''.join(output_lines)  # Return captured output
    
    
#    try:
#        result = subprocess.run(command, check=True, capture_output=True, text=True)
#        return result.stdout
        
#    except subprocess.CalledProcessError as e:
#        print(f"Error running {script_name}: {e}")
#        print(e.stderr)
