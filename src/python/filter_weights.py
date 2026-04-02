#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 01:26:30 2026

@author: lahumada
"""

'''
Generates a filtered EVM weights file containing only the weight lines corresponding 
to user-selected gene annotation tools. This ensures EVM only considers evidence 
from tools that were actually run in the pipeline.

Main function: generate_tool_weights(input_weights_file, output_weights_file, selected_tools)
'''



def generate_tool_weights(input_weights_file, output_weights_file, selected_tools):
    """
    Creates a filtered weights file containing only lines where weight ID contains selected tool names.
    
    Args:
        input_weights_file (str): Path to the full weights TXT file
        output_weights_file (str): Path to generate filtered weights TXT file
        selected_tools (set): Set of selected tool names (e.g., {'gmap', 'minimap2'})
    
    Returns:
        bool: True if successful and lines found, False otherwise
    """
    matching_lines = []
    
    # Read input file line by line
    try:
        with open(input_weights_file, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty/comment lines
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    weight_id = parts[1]
                    
                    # Check if ANY selected tool name appears in weight_id (case-insensitive)
                    for tool in selected_tools:
                        if tool.lower() in weight_id.lower():
                            matching_lines.append(line)
                            break  # Found match, no need to check other tools
        
        # Write filtered weights
        with open(output_weights_file, 'w') as outfile:
            for line in matching_lines:
                outfile.write(line + '\n')
        
        if matching_lines:
            print(f"Generated {output_weights_file} with {len(matching_lines)} weight lines")
            print(f"  Selected tools: {selected_tools}")
        else:
            print(f"Warning: No matching weights found for tools: {selected_tools}")
            # Still create empty file so EVM doesn't fail
            with open(output_weights_file, 'w') as outfile:
                pass
            
    except FileNotFoundError:
        print(f"Error: Input weights file not found: {input_weights_file}")
    except Exception as e:
        print(f"Error processing weights file: {e}")
