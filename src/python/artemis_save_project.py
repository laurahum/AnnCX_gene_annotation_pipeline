#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 20:13:02 2022

Modified: Mon Feb 08 2022

@author: lahumada
"""

### Open automatically in Artemis (Project File Manager) the annotation data 
### for each genome in the annotation pipeline, and save each genome into 
### a Project (title, ROI and annotations).

### Artemis saves the Project File Manager data into a hidden TXT file and BAK
### file in the home directory:
### /home/lahumada/.artemis.project.properties
### /home/lahumada/.artemis.project.properties.bak (back up file)
### If this file is not created yet, this script creates the file.

### Each project looks like this:
### They are separated by a '#' character, which is also included at the start
### and end of the file:
### #
### project.project_name.title=project_name
### project.project_name.sequence=/path/to/genomic/sequence/annotated/genomic_sequence.fasta
### project.project_name.annotation=/path/to/annotation/file/annotation.gff3
### #

### The BAK file is created automatically by Artemis once the Artemis program
### is opened and closed by the user.

### Main function: artemis_save_project(single_contig_list, gene_to_annotate, ROI_dir, annotation_dirs)


import os
from pathlib import Path


# Function to get and sort files within a directory
def process_dir (directory_name):
    """
    Process a directory and return an alphabetically sorted list of its file contents.
        - directory_name (str): Path to the directory to process.

    Returns (list): Sorted list of filenames in the directory.
    """
    directory_files = os.listdir(directory_name)
    directory_files.sort() 
    return (directory_files)


## Make function to return each path to file given a directory and genome name
def give_path_to_feature (genome_name, directory_name):
    """
    Find and return the path to a file for a specific genome in a given directory.
        - genome_name (str): Name of the genome to search for.
        - directory_name (str): Path to the directory to search in.

    Returns (str): Full path to the file containing the genome name, or None if not found.
    """
    dir_files = process_dir (directory_name)
    for file_name in dir_files:
        # Only if the name of the file corresponds to the genome_name
        if (genome_name in file_name):
        # Return the filepath as a string
            filepath = os.path.join(directory_name, file_name)
            return (filepath)
            break


# Function to create a missed annotation report
def create_missed_annotation_report(errors, output_dir):
    """
    Write to a file for each genome any annotation missing files which can happen when
    one of the annotation tools does not work on a genomic region for any matter.
        - errors(list): List with the errors found
        - output_dir(str): Path to the directory where the report will be saved
    """
    
    # File name
    output_file = os.path.join(output_dir, "missed_annotation_report.txt")
    
    with open(output_file, "w") as file:
        for err in errors:
            file.write(f"{err}\n")


# Function to create a dictionary with all the data necessary to make an artemis project
def create_content_artemis_project(single_contig_list, gene_to_annotate, ROI_dir, annotation_dirs, output_dir_missed):
    
    """
    Create a list of dictionaries containing Artemis project information for a list of genome.
    Reads a list of genome from a file and create a dictionary for each genome
    containing information needed for an Artemis project. This includes the genome name,
    project title, sequence and annotation file paths. 
        - single_contig_list (str): Path to a file containing a list of genome names.
        - gene_to_annotate (str): Name of the gene to annotate. (Example: "NKG2")
        - ROI_dir: Path to the directory containing the annotated genomic FASTA files
        - annotation_dirs (list): List of paths to directories containing annotation files.
        - output_dir_missed (str): Path to save report created with create_missed_annotation_report()

    Returns (list): List of dictionaries, each containing project information for a single genome.
    """
    
    # Read genome file
    genome_file = open(single_contig_list,'r')
    genome_list = genome_file.readlines() 

    ## Make an empty dictionary for each genome project
    # Save the dictionaries within a list
    List_of_projects = []
    errors = [] # for any missing annotation files

    for genome_name in genome_list:
        genome_name = genome_name.strip()
    
        project_name = f"{genome_name}_{gene_to_annotate}"
    
        entry_sequence = f"project.{project_name}.sequence="
        entry_annotation = f"project.{project_name}.annotation="
        entry_title = f"project.{project_name}.title={project_name}"
        
        ROI_path = give_path_to_feature(genome_name, ROI_dir)
        annotation_path = []
        for each_dir in annotation_dirs:
            path = give_path_to_feature(genome_name, each_dir)
            if path != None:
                annotation_path.append(path)
            else:
                # Get parent dir to the annotation_dir 
                parent_dir = each_dir.rsplit('/', 1)[0]
                errors.append(f"{genome_name} is missing annotation for {parent_dir}")
        
        each_dict = {
            'genome_name':[genome_name],
            'title.name':[entry_title],
            'sequence.name':[entry_sequence],
            'annotation.name':[entry_annotation],
            'sequence.path':[ROI_path],
            'annotation.path': annotation_path
            }
    
        List_of_projects.append(each_dict)
    
        # Control
        print(genome_name + ' create_dictionary')
    
    if errors:
        create_missed_annotation_report(errors, output_dir_missed)        

    return (List_of_projects)
    
    
## Write data from the list of dictionaries (List_of_projects) to the Artemis TXT file
def artemis_save_project(single_contig_list, gene_to_annotate, ROI_dir, annotation_dirs, output_dir_missed):

    """
    Create and save Artemis project information for multiple genome to a the file '.artemis.project.properties' in the 
    user's home directory. This file is created by Artemis once the tool is opened and used. The function checks if the
    file does not exist yet and creates a new one. The function writes project information for each genome in the format
    required by Artemis, including the project title, sequence file path, and annotation file paths. Each project is 
    separated by a '#' character in the file.
        - single_contig_list (str): Path to a file containing a list of genome names.
        - gene_to_annotate (str): Name of the gene to annotate.
        - ROI_dir: Path to the directory containing the annotated genomic FASTA files
        - annotation_dirs (list): List of paths to directories containing annotation files.
        - output_dir_missed (str): Path to save report created with create_missed_annotation_report()
    """
    
    # 1. Find or create artemis TXT file where the projects are written
    # Get home directory
    home_dir = str(Path.home())
    artemis_file_path = os.path.join(home_dir, '.artemis.project.properties')

    # Check if the file exists and if not create it
    if not os.path.exists(artemis_file_path):
        Path(artemis_file_path).touch()

    # 2. Create list of dictionaries:
    List_of_projects = create_content_artemis_project(single_contig_list, gene_to_annotate, ROI_dir, annotation_dirs, output_dir_missed)

    # 3. Write each dictionary as an Artemis project: 
    with open(artemis_file_path, "a+") as new_file:
        # Move read cursor to the start of file.
        new_file.seek(0)
        # If file is not empty then append '\n'
        data = new_file.read()
        
        if len(data) > 0 :
            if data[-1] != '#\n':
                new_file.write('#' + '\n')
            else:
                new_file.write("\n")
        else:
            new_file.write('#')
        
        # Iterate over the list of dictionaries:            
        for dictionary in List_of_projects:
        
            # Make a list for each genome where I will be saving the data in the
            # format that the projects are written to the Artemis text file
            List_each_genome_data = []
       
            # Check that the loop goes over each genome
            each_genome = dictionary['genome_name'][0]
        
            # Fill the list with Project title, sequence and annotations
            ## Project title name
            project_title_name = dictionary['title.name'][0]
            List_each_genome_data.append(project_title_name)
            ## Project sequence name + file path
            project_sequence_name = dictionary['sequence.name'][0]
            project_sequence_path = dictionary['sequence.path'][0]
            project_sequence = (project_sequence_name + project_sequence_path)
            List_each_genome_data.append(project_sequence)
            # Project annotation name + annotation paths
            project_annotation_name = dictionary['annotation.name'][0]
            project_annotation_list = dictionary['annotation.path']
            project_annotation_path = ' '.join(p for p in project_annotation_list if p is not None) # Skip any annotation paths that are None
            project_annotation = (project_annotation_name + project_annotation_path)
            List_each_genome_data.append(project_annotation)
        
            # Write the data in List_each_genome_data to the Artemis text file:
            for element in List_each_genome_data:
                new_file.write(element + '\n')
        
            # Write the separator of projects: '#'
            new_file.write('#' + '\n')
            
            # Control
            print(each_genome + ' write_project')
        
