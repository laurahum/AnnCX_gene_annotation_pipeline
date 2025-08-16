#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 20:37:06 2025

@author: lahumada
"""

import sys
import os
import glob
import re
import shutil
import argparse
import subprocess
from contextlib import contextmanager # for the report

# Adjust Python path (preliminary)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Local application imports
from src.python import (
    get_name_fasta_files,
    find_found_flanking_genes,
    single_contig_genomes,
    copy_single_contig_genomes,
    extract_seqkit_data_all_file,
    create_N_annotation_gff3_all_files,
    format_genewise_output_to_gff_all_files,
    format_exonerate_output_to_gff_all_files,
    filter_blast_one_map_per_region,
    filter_gmap_one_map_per_region,
    filter_exonerate_one_map_per_region,
    filter_genewise_one_map_per_region,
    convert_augustus_to_EVM_all_files,
    convert_gmap_to_EVM_all_files,
    convert_blastn_to_EVM_all_files,
    convert_tblastn_to_EVM_all_files,
    convert_genewise_to_EVM_all_files,
    convert_exonerate_to_EVM_all_files,
    filter_EVM_results_get_only_genes,
    filter_EVM_results_overlap_number_tools,
    artemis_save_project,
    run_bash_script)

# Specific utility imports
from src.python import (check_overlaps_range,
                        check_dir_path,
                        check_file_path,
                        create_out_dir, 
                        print_results, 
                        user_continue_prompt)

# Post-processing script
from src.python import get_fasta_sequences_annotation


# Version
version = "1.0.0"


# Set WISECONFIGDIR necessary for running GeneWise 
wise_config_dir = os.path.join(os.environ.get("CONDA_PREFIX", ""), "share", "wise2", "wisecfg")
if os.path.isdir(wise_config_dir):
    os.environ["WISECONFIGDIR"] = wise_config_dir
else:
    sys.stderr.write(
        f"WARNING: WISECONFIGDIR directory '{wise_config_dir}' does not exist. "
        "GeneWise may not function properly.\n"
    )


# Context manager function to make a report
@contextmanager
def tee_output(report_file_path):
    class Tee:
        def __init__(self, file, terminal):
            self.file = file
            self.terminal = terminal
        def write(self, message):
            self.terminal.write(message)
            self.file.write(message)
        def flush(self):
            self.terminal.flush()
            self.file.flush()

    original_stdout = sys.stdout
    original_stderr = sys.stderr

    with open(report_file_path, 'w') as f:
        sys.stdout = Tee(f, original_stdout)
        sys.stderr = Tee(f, original_stderr)
        try:
            yield
        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr


def main():
    parser = argparse.ArgumentParser(description=f'AnnCX gene annotation pipeline {version}')
    
    # All arguments from user
    parser.add_argument('--genome', type=check_dir_path, required=True, help='Directory where the genomic FASTA files are located. Example: --genome /path/to/genomes') # Directory
    parser.add_argument('--namegenes', type=str, required=True, help='Name of the genes to be annotated. Example: --namegenes NKG2') # String
    parser.add_argument('--querycDNA', type=check_file_path, required=True, help='FASTA file with entries for cDNA sequences to use as query. Example: --querycDNA /path/to/query/querycDNA.fasta') # File
    parser.add_argument('--queryprot', type=check_file_path, required=True, help='FASTA file with entries for protein sequences to use as query. Example: --queryprot /path/to/query/queryprot.fasta') # File
    parser.add_argument('--queryexon', type=check_file_path, required=True, help='FASTA file with entries for exon sequences to use as query. Example: --queryexon /path/to/query/queryexon.fasta') # File
    parser.add_argument('--spsrepeatmasker', type=str, required=True, help='Name of the phylogenetic group to run RepeatMasker. Example: --spsrepeatmasker primates') # String
    parser.add_argument('--spsaugustus', type=str, required=True, help='Name of the species to run AUGUSTUS. Example: --spsaugustus human') # String
    parser.add_argument('--outdir', type=check_dir_path, required=True, help='Directory to save the output files. Example: --outdir /path/to/output') # Directory
    parser.add_argument('--flanking', type=check_file_path, help='(OPTIONAL) FASTA file with entries for gene sequences to use as flanking genes and extract a genomic region of interest from a larger genome file to annotate. Example: --flanking /path/to/flanking/geneflanking.fasta') # File
    parser.add_argument('--maxintron', default=7000, type=int, help='(OPTIONAL) Sets maximum intron size in basepairs  (Default = 7000). Example: --maxintron 6000') # Int
    parser.add_argument('--threads', default=-1, type=int, help='(OPTIONAL) Number of threads that can be used to run some tools. Example: --threads 3') # Int
    parser.add_argument('--overlapsEVM', type=check_overlaps_range, default=4, help='(OPTIONAL) Number of overlapping gene annotation tools to run EVM (0-7, Default = 4). Example: --overlapsEVM 2') # Int
    parser.add_argument('--skip-prompt', action='store_true', help='(OPTIONAL) If the argument is used, skips the prompt after finding the flanking genes that asks the user whether to continue with the annotation process. Example: --skip-prompt') # Boolean
    parser.add_argument('--version', action='version', version=f'AnnCX version {version}')

    args = parser.parse_args()
    

    # Get the directory of the main script
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Make base output directory
    output_dir_arg = os.path.abspath(args.outdir) # get absolute path
    output_dir = create_out_dir(output_dir_arg, f'AnnCX_{args.namegenes}', path=True)

    # Make a report file 
    report_dir = create_out_dir(output_dir, 'report')
    report_file = os.path.join(report_dir, "pipeline_report.txt")
    
    # Use tee_output context manager function 
    with tee_output(report_file):
        run_pipeline(current_dir, output_dir, args)



def run_pipeline(current_dir, output_dir, args):
    """ 
    Executes the AnnCX gene annotation pipeline using the provided configuration.
    - current_dir (str): Directory path to where the main script is located.
    - output_dir (str): Base directory path for the pipeline results
    - args (argparse.Namespace): Parsed command-line arguments containing all user-specified parameters.
    """
    
    print("=========================================================================================")
    print("=========================================================================================")
    print("Started AnnCX running")
    
    # Make directory for errors
    error_dir = create_out_dir(output_dir, 'errors')
    # Remove all contents in error_dir at pipeline start
    for item in error_dir.iterdir():
        if item.is_file() or item.is_dir():
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()

    # Path to file with the list of genomes given as input
    genomes_input_dir = create_out_dir(output_dir, 'genomes_input')
    genomes_error_dir = create_out_dir(error_dir, 'format_genomes_input')
    genomes_abs_path = os.path.abspath(args.genome) # get absolute path
    genomes_list_input_file = get_name_fasta_files(genomes_abs_path, str(genomes_input_dir), genomes_error_dir)
    
    # Look for any error TXT file in the error directory
    txt_files = glob.glob(os.path.join(genomes_error_dir, "*.txt"))

    # If there is at least one error TXT file, print its contents
    if txt_files:
        for file in txt_files:
            with open(file, "r") as f:
                print(f.read())
                
        print("\nPlease fix the above errors and re-run the pipeline.")
        sys.exit(1)


    # Directory with the genomic region to annotate (either directly the one input by the user or the extracted with flanking genes)
    roi_raw_dir = None
    # Path to file with the list of genomes to be annotated
    genomes_list_process_file = None
    
    # List of directories with the output annotations from the pipeline to visualize in Artemis
    artemis_annotations_dir = []
    

    ## Step 1. Extract genomic region of interest from a larger genomic sequence using two flanking genes
    if (os.path.abspath(args.flanking)):
        # Run gmap to annotate the flanking genes
        print("=========================================================================================")
        print("*** STEP 1. Extract ROI")
        print("=========================================================================================")

        roi_subdir = create_out_dir(output_dir, '1_extract_ROI')

        print("-----------------------------------------------------------------------------------------")       
        print("1.1. Find flanking genes")
        print("-----------------------------------------------------------------------------------------")
        print("1.1.1. Make GMAP database")
        genome_gmap_build_dir = create_out_dir(roi_subdir, '1_gmap_build')
        run_bash_script('gmap_build.sh', 'gmap_build_run', str(genomes_list_input_file), str(genome_gmap_build_dir), genomes_abs_path, "genome", args.threads)
        
        print("1.1.2. Run GMAP")
        genome_gmap_dir = create_out_dir(roi_subdir, '2_gmap')
        run_bash_script('gmap.sh', 'gmap_run', str(genomes_list_input_file), str(genome_gmap_build_dir), os.path.abspath(args.flanking), str(genome_gmap_dir), args.namegenes, "flanking", "genome", args.threads)
        
        # Find whether the flanking genes were found 
        print("-----------------------------------------------------------------------------------------")
        print("1.2. Found flanking genes:")
        print("-----------------------------------------------------------------------------------------")
        found_flanking_dir = create_out_dir(roi_subdir, '3_found_flanking_genes')
        yes_flanking, no_flanking = find_found_flanking_genes(os.path.abspath(args.flanking), str(genome_gmap_dir), str(found_flanking_dir), args.namegenes)
        
        print(f"Genomes with both flanking genes found: {yes_flanking}")
        print_results(yes_flanking)
        print(f"Other genomes: {no_flanking}")
        print_results(no_flanking) 
        
        # Ask user if they want to continue with the annotation
        if not args.skip_prompt:
            if not user_continue_prompt():
                print("Exiting AnnCX")
                sys.exit(0)

        print("Continuing running AnnCX...")
            
        # Find genomes in which the flanking genes are found in a single contig 
        print("-----------------------------------------------------------------------------------------")
        print("1.3. Find genomes in which the flanking genes are located in the same contig/scaffold/chr")
        print("-----------------------------------------------------------------------------------------")
        single_contig_dir = create_out_dir(roi_subdir, '4_find_single_contig_genomes')
        single_contig_sp = single_contig_genomes(str(genome_gmap_dir), str(single_contig_dir), str(yes_flanking), args.namegenes)
        genomes_list_process_file = single_contig_sp
        
        # Move GFF3 files from single contig genomes to another directory
        print("-----------------------------------------------------------------------------------------")
        print("1.4. Get genomes with flanking genes in a single contig/scaffold/chr")
        print("-----------------------------------------------------------------------------------------")
        copy_single_contig_dir = create_out_dir(roi_subdir, '5_found_single_contig_genomes')
        copy_single_contig_genomes(str(genome_gmap_dir), str(copy_single_contig_dir), str(genomes_list_process_file))
        
        # Get from the GFF3 files the data to be used by seqkit to extract ROI and save as BED
        print("-----------------------------------------------------------------------------------------")
        print("1.5. Get data to extract region of interest")
        print("-----------------------------------------------------------------------------------------")
        bed_for_seqkit_dir = create_out_dir(roi_subdir, '6_bed_files_for_seqkit')
        extract_seqkit_data_all_file(str(copy_single_contig_dir), str(bed_for_seqkit_dir))
    
        # Extract ROI with seqkit
        print("-----------------------------------------------------------------------------------------")
        print("1.6. Extract genomic region of interest with SeqKit")
        print("-----------------------------------------------------------------------------------------")
        roi_extracted_dir = create_out_dir(roi_subdir, '7_extracted_roi_raw')
        run_bash_script('seqkit_extract_ROI.sh', 'seqkit_run', str(bed_for_seqkit_dir), genomes_abs_path, str(roi_extracted_dir), args.namegenes, args.threads)
        
        roi_raw_dir = str(roi_extracted_dir)
    
    
    # Use the genomic region directly introduced by the user
    else:
        roi_raw_dir = genomes_abs_path
        genomes_list_process_file = genomes_list_input_file
    

    ## Step 2. Annotate potential gaps in the sequence
    print("=========================================================================================")
    print("*** STEP 2. Annotate gaps in the genomic region of interest")
    print("=========================================================================================")
    annotate_Ns_dir = create_out_dir(output_dir, '2_annotate_N_gaps')
    artemis_annotations_dir.append(str(annotate_Ns_dir))
    create_N_annotation_gff3_all_files(str(genomes_list_process_file), str(roi_raw_dir), str(annotate_Ns_dir))
    
    
    ## Step 3. Hard-masking of ROI
    print("=========================================================================================")
    print("*** STEP 3. Hard-masking genomic region of interest")
    print("=========================================================================================")
    repeatmasker_dir = create_out_dir(output_dir, '3_repeatmasker')
    output_repeatmasker = run_bash_script('RepeatMasker_hardmask_ROI.sh', 'repeatmasker_run', str(roi_raw_dir), str(repeatmasker_dir), args.spsrepeatmasker, args.threads)
    
    if output_repeatmasker:
        roi_hardmasked_folder = re.search(r'MASKED_FOLDER=(.*)', output_repeatmasker)
        repeat_annotations_folder = re.search(r'ANNOTATIONS_FOLDER=(.*)', output_repeatmasker)

        if roi_hardmasked_folder and repeat_annotations_folder:
            roi_hardmasked_dir = roi_hardmasked_folder.group(1)
            repeat_annotations_dir = repeat_annotations_folder.group(1)
            artemis_annotations_dir.append(str(repeat_annotations_dir))
            
    ## Step 4. Run gene annotation tools on genomic region of interest
    print("=========================================================================================") 
    print("*** STEP 4. Run gene annotation tools")
    print("=========================================================================================")
    
    # Create directories
    annotations_dir = create_out_dir(output_dir, '4_gene_annotation_tools')

    databases_dir = create_out_dir(annotations_dir, 'databases') 
    
    makeblastdb_dir = create_out_dir(databases_dir, 'makeblastdb')
    blast_base_dir = create_out_dir(annotations_dir, 'BLAST')
    blastn_base_dir = create_out_dir(blast_base_dir, 'BLASTN')
    blastn_raw_dir = create_out_dir(blastn_base_dir, 'raw')
    tblastn_base_dir = create_out_dir(blast_base_dir, 'TBLASTN')
    tblastn_raw_dir = create_out_dir(tblastn_base_dir, 'raw')   
   
    gmap_build_dir = create_out_dir(databases_dir, 'gmap_build')
    gmap_base_dir = create_out_dir(annotations_dir, 'GMAP')
    gmap_cDNA_base_dir = create_out_dir(gmap_base_dir, 'GMAP_cDNA')
    gmap_cDNA_raw_dir = create_out_dir(gmap_cDNA_base_dir, 'raw')
    gmap_exon_base_dir = create_out_dir(gmap_base_dir, 'GMAP_exon')
    gmap_exon_raw_dir = create_out_dir(gmap_exon_base_dir, 'raw')
    
    genewise_base_dir = create_out_dir(annotations_dir, 'GeneWise')
    genewise_raw_dir = create_out_dir(genewise_base_dir, 'raw')
    
    exonerate_base_dir = create_out_dir(annotations_dir, 'Exonerate')
    exonerate_raw_dir = create_out_dir(exonerate_base_dir, 'raw')

    protprof_dir = create_out_dir(output_dir, 'protein_profile')
    augustus_base_dir = create_out_dir(annotations_dir, 'AUGUSTUS')
    augustus_raw_dir = create_out_dir(augustus_base_dir, 'raw')

    # Run tools
    print("-----------------------------------------------------------------------------------------")
    print("4.1. Run BLAST")
    print("-----------------------------------------------------------------------------------------")
    # Make blast database with makeblastdb
    print("4.1.1. Make BLAST database")
    run_bash_script('blast_makeblastdb_ROI_hardmasked.sh', 
                    'makeblastdb_run', 
                    str(roi_hardmasked_dir), 
                    str(makeblastdb_dir), 
                    str(genomes_list_process_file))

    print("4.1.2. Run BLASTN")
    run_bash_script('blastn_query_ROI_hardmasked_out6.sh', 
                    'blastn_run', 
                    str(makeblastdb_dir), 
                    os.path.abspath(args.querycDNA), 
                    str(blastn_raw_dir), 
                    str(genomes_list_process_file), 
                    args.namegenes, 
                    args.threads)

    print("4.1.3. Run TBLASTN")
    run_bash_script('tblastn_query_ROI_hardmasked_out6.sh', 
                    'tblastn_run', 
                    str(makeblastdb_dir), 
                    os.path.abspath(args.queryprot), 
                    str(tblastn_raw_dir), 
                    str(genomes_list_process_file), 
                    args.namegenes, 
                    args.threads)
 
    print("-----------------------------------------------------------------------------------------")
    print("4.2. Run GMAP")
    print("-----------------------------------------------------------------------------------------")
    # Make GMAP database with gmap_build
    print("4.2.1. Make GMAP database")
    run_bash_script('gmap_build.sh', 
                    'gmap_build_run', 
                    str(genomes_list_process_file), 
                    str(gmap_build_dir), 
                    str(roi_hardmasked_dir), 
                    "ROI", 
                    args.threads)

    print("4.2.2. Run GMAP cDNA query")
    run_bash_script('gmap.sh', 
                    'gmap_run', 
                    str(genomes_list_process_file), 
                    str(gmap_build_dir), 
                    os.path.abspath(args.querycDNA), 
                    str(gmap_cDNA_raw_dir), 
                    args.namegenes, 
                    "cDNA", 
                    "ROI", 
                    args.threads)

    print("4.2.3. Run GMAP exon query")
    run_bash_script('gmap.sh', 
                    'gmap_run', 
                    str(genomes_list_process_file), 
                    str(gmap_build_dir), 
                    os.path.abspath(args.queryexon), 
                    str(gmap_exon_raw_dir), 
                    args.namegenes, 
                    "exon", 
                    "ROI", 
                    args.threads)

    print("-----------------------------------------------------------------------------------------")
    print("4.3. Run GeneWise")
    print("-----------------------------------------------------------------------------------------")
    run_bash_script('genewise_query_ROI_hardmasked.sh', 
                    'genewise_run', 
                    str(roi_hardmasked_dir), 
                    os.path.abspath(args.queryprot), 
                    str(genewise_raw_dir), 
                    str(genomes_list_process_file), 
                    args.namegenes)

    print("-----------------------------------------------------------------------------------------")
    print("4.4. Run Exonerate")
    print("-----------------------------------------------------------------------------------------")
    run_bash_script('exonerate_query_ROI_hardmasked.sh', 
                    'exonerate_run', 
                    str(roi_hardmasked_dir), 
                    os.path.abspath(args.querycDNA), 
                    str(exonerate_raw_dir), 
                    str(genomes_list_process_file), 
                    args.namegenes, 
                    args.maxintron, 
                    args.threads)

    print("-----------------------------------------------------------------------------------------")
    print("4.5. Run AUGUSTUS")
    print("-----------------------------------------------------------------------------------------")

    # Make protein profile
    print("4.5.1. Make protein profile")
    output_protprof = run_bash_script('make_augustus_protprof.sh', 
                                      'make_augustus_protprof', 
                                      os.path.abspath(args.queryprot), 
                                      str(protprof_dir), 
                                      args.namegenes)
    
    if output_protprof:
        protprof_match = re.search(r'PROTPROF_FOLDER=(.*)', output_protprof)
        if protprof_match:
            protprof_file = protprof_match.group(1)
    
    # Run Augustus on raw ROI
    print("4.5.2. Run AUGUSTUS")
    run_bash_script('augustus_query_ROI_unmasked_protprof.sh', 'augustus_run', str(roi_raw_dir), str(protprof_file), str(augustus_raw_dir), str(genomes_list_process_file), args.namegenes, args.spsaugustus)
    artemis_annotations_dir.append(augustus_raw_dir)
    
    
    ## Step 5. Formatting gene annotation output
    print("=========================================================================================")
    print("*** STEP 5. Formatting gene annotation output")
    print("=========================================================================================")
    
    # Make directories
    blastn_format_base_dir = create_out_dir(blastn_base_dir, 'formatted')   
    blastn_format_1_dir = create_out_dir(blastn_format_base_dir, 'gff')   
    blastn_format_2_dir = create_out_dir(blastn_format_base_dir, 'gff_formatted')   

    tblastn_format_base_dir = create_out_dir(tblastn_base_dir, 'formatted')
    tblastn_format_1_dir = create_out_dir(tblastn_format_base_dir, 'gff')   
    tblastn_format_2_dir = create_out_dir(tblastn_format_base_dir, 'gff_formatted')   
    
    genewise_format_dir = create_out_dir(genewise_base_dir, 'formatted')   
    exonerate_format_dir = create_out_dir(exonerate_base_dir, 'formatted')   

    print("-----------------------------------------------------------------------------------------")
    print("5.1. Formatting BLAST")
    print("-----------------------------------------------------------------------------------------")
    run_bash_script('format_blast_outfmt6_to_gff.sh', 'format_blastn_out6_to_gff', str(blastn_raw_dir), str(blastn_format_1_dir), str(blastn_format_2_dir), str(genomes_list_process_file), args.namegenes)
    run_bash_script('format_blast_outfmt6_to_gff.sh', 'format_tblastn_out6_to_gff', str(tblastn_raw_dir), str(tblastn_format_1_dir), str(tblastn_format_2_dir), str(genomes_list_process_file), args.namegenes)

    print("-----------------------------------------------------------------------------------------")
    print("5.2. Formatting GeneWise")
    print("-----------------------------------------------------------------------------------------")
    format_genewise_output_to_gff_all_files(str(genewise_raw_dir), str(genewise_format_dir))

    print("-----------------------------------------------------------------------------------------")
    print("5.3. Formatting Exonerate")
    print("-----------------------------------------------------------------------------------------")
    format_exonerate_output_to_gff_all_files(str(exonerate_raw_dir), str(exonerate_format_dir))
    
    
    ## Step 6. Filtering gene annotation output
    print("=========================================================================================")
    print("*** STEP 6. Filtering gene annotation output")
    print("=========================================================================================")
    
    # Make directories 
    blastn_filter_dir = create_out_dir(blastn_base_dir, 'filtered')
    tblastn_filter_dir = create_out_dir(tblastn_base_dir, 'filtered')
    gmap_cDNA_filter_dir = create_out_dir(gmap_cDNA_base_dir, 'filtered')
    gmap_exon_filter_dir = create_out_dir(gmap_exon_base_dir, 'filtered')
    genewise_filter_dir = create_out_dir(genewise_base_dir, 'filtered')
    exonerate_filter_dir = create_out_dir(exonerate_base_dir, 'filtered')

    print("-----------------------------------------------------------------------------------------")
    print("6.1. Filtering BLAST")
    print("-----------------------------------------------------------------------------------------")
    filter_blast_one_map_per_region(str(blastn_format_2_dir), str(blastn_filter_dir)) # BLASTN
    filter_blast_one_map_per_region(str(tblastn_format_2_dir), str(tblastn_filter_dir)) # TBLASTN

    print("-----------------------------------------------------------------------------------------")
    print("6.2. Filtering GMAP")
    print("-----------------------------------------------------------------------------------------")
    filter_gmap_one_map_per_region(str(gmap_cDNA_raw_dir), str(gmap_cDNA_filter_dir)) # GMAP cDNA
    filter_gmap_one_map_per_region(str(gmap_exon_raw_dir), str(gmap_exon_filter_dir)) # GMAP exon

    print("-----------------------------------------------------------------------------------------")
    print("6.3. Filtering GeneWise")
    print("-----------------------------------------------------------------------------------------")
    filter_genewise_one_map_per_region(str(genewise_format_dir), str(genewise_filter_dir))

    print("-----------------------------------------------------------------------------------------")
    print("6.4. Filtering Exonerate")
    print("-----------------------------------------------------------------------------------------")
    filter_exonerate_one_map_per_region(str(exonerate_format_dir), str(exonerate_filter_dir))
    
    # Add to artemis directories 
    artemis_annotations_dir.append(blastn_filter_dir)
    artemis_annotations_dir.append(tblastn_filter_dir)
    artemis_annotations_dir.append(gmap_cDNA_filter_dir)
    artemis_annotations_dir.append(gmap_exon_filter_dir)
    artemis_annotations_dir.append(genewise_filter_dir)
    artemis_annotations_dir.append(exonerate_filter_dir)


    ## Step 7. Combine gene annotation tool output into a consensus annotation with EVM
    print("=========================================================================================")
    print("*** STEP 7. Combine gene annotation tool output into a consensus annotation with EVM")
    print("=========================================================================================")
    
    # Make directories
    blastn_evm_dir = create_out_dir(blastn_base_dir, 'evm')
    tblastn_evm_dir = create_out_dir(tblastn_base_dir, 'evm')
    gmap_cDNA_evm_dir = create_out_dir(gmap_cDNA_base_dir, 'evm')
    gmap_exon_evm_dir = create_out_dir(gmap_exon_base_dir, 'evm')
    genewise_evm_dir = create_out_dir(genewise_base_dir, 'evm')
    exonerate_evm_dir = create_out_dir(exonerate_base_dir, 'evm')   
    augustus_evm_dir = create_out_dir(augustus_base_dir, 'evm')   
    
    evm_base_dir = create_out_dir(output_dir, '7_consensus_EVM')
    evm_concatenate_dir = create_out_dir(evm_base_dir, '1_concatenate')
    evm_run_base_dir = create_out_dir(evm_base_dir, '2_run')
    evm_run_raw_dir = create_out_dir(evm_run_base_dir, 'raw')
    evm_run_gff3_dir = create_out_dir(evm_run_base_dir, 'gff3')
    evm_filter_dir = create_out_dir(evm_base_dir, '3_filter')
    evm_filter_1_dir = create_out_dir(evm_filter_dir, 'filter_1')
    evm_filter_2_dir = create_out_dir(evm_filter_dir, 'filter_2')
    evm_filter_3_dir = create_out_dir(evm_filter_dir, 'filter_3')

    print("-----------------------------------------------------------------------------------------")    
    print("7.1. Convert gene annotation output to EVM")
    print("-----------------------------------------------------------------------------------------")
    print("7.1.1. Convert BLAST")
    convert_blastn_to_EVM_all_files(str(blastn_filter_dir), str(blastn_evm_dir), 'blastn') # BLASTN
    convert_tblastn_to_EVM_all_files(str(tblastn_filter_dir), str(tblastn_evm_dir), 'tblastn') #TBLASTN

    print("7.1.2. Convert GMAP")
    convert_gmap_to_EVM_all_files(str(gmap_cDNA_filter_dir), str(gmap_cDNA_evm_dir), 'gmapcDNA') # GMAP cDNA
    convert_gmap_to_EVM_all_files(str(gmap_exon_filter_dir), str(gmap_exon_evm_dir), 'gmapExon') # GMAP exon

    print("7.1.3. Convert GeneWise")
    convert_genewise_to_EVM_all_files(str(genewise_filter_dir), str(genewise_evm_dir), 'genewise')
    
    print("7.1.4. Convert Exonerate")
    convert_exonerate_to_EVM_all_files(str(exonerate_filter_dir), str(exonerate_evm_dir))
    
    print("7.1.5. Convert AUGUSTUS")
    convert_augustus_to_EVM_all_files(str(augustus_raw_dir), str(augustus_evm_dir))
    
    print("-----------------------------------------------------------------------------------------")    
    print("7.2. Concatenate converted files")
    print("-----------------------------------------------------------------------------------------")
    run_bash_script('EVM_concatenate_GFF3.sh', 'EVM_concatenate_GFF3', 
                    str(genewise_evm_dir),
                    str(tblastn_evm_dir),
                    str(blastn_evm_dir),
                    str(gmap_cDNA_evm_dir),
                    str(exonerate_evm_dir),
                    str(gmap_exon_evm_dir),
                    str(evm_concatenate_dir),
                    str(genomes_list_process_file))

    print("-----------------------------------------------------------------------------------------")
    print("7.3. Run EVM")
    print("-----------------------------------------------------------------------------------------")
    # Get weights path
    weights_path = os.path.join(current_dir, 'weights', 'weights_recommended.txt')
    # Run EVM
    run_bash_script('EVM_all_tools_ROI.sh', 'EVM_run', 
                    str(roi_hardmasked_dir),
                    str(augustus_evm_dir),
                    str(evm_concatenate_dir),
                    str(evm_run_raw_dir),
                    str(evm_run_gff3_dir),
                    str(weights_path),
                    str(genomes_list_process_file),
                    args.threads)

    print("-----------------------------------------------------------------------------------------")    
    print("7.4. Filtering EVM")
    print("-----------------------------------------------------------------------------------------")
    # Filter 1/3
    print("7.4.1. Filter 1/3")
    filter_EVM_results_get_only_genes(str(evm_run_gff3_dir), str(evm_filter_1_dir))
    # Filter 2/3
    print("7.4.2. Filter 2/3")
    run_bash_script('filter_EVM_results_find_overlaps.sh', 'filter_EVM_results_find_overlaps',
                    str(evm_filter_1_dir),
                    str(augustus_raw_dir),
                    str(blastn_filter_dir),
                    str(tblastn_filter_dir),
                    str(gmap_cDNA_filter_dir),
                    str(gmap_exon_filter_dir),
                    str(genewise_filter_dir),
                    str(exonerate_filter_dir),
                    str(evm_filter_2_dir),
                    str(genomes_list_process_file))
    # Filter 3/3
    print("7.4.3. Filter 3/3")
    filter_EVM_results_overlap_number_tools(str(evm_run_gff3_dir),
                                            str(evm_filter_2_dir),
                                            str(evm_filter_3_dir),
                                            str(genomes_list_process_file),
                                            args.overlapsEVM,
                                            version)

    # Append results of EVM to Artemis 
    artemis_annotations_dir.append(evm_filter_3_dir)

    print("=========================================================================================")    
    print("*** STEP Post-processing: Convert annotations to FASTA")
    print("=========================================================================================")
    fasta_out_dir = create_out_dir(output_dir, 'output_raw_FASTA_annotations')
    get_fasta_sequences_annotation(str(genomes_list_process_file), 
                                   args.namegenes,
                                   str(roi_raw_dir),
                                   str(evm_filter_3_dir),
                                   str(fasta_out_dir))

    print("=========================================================================================")    
    print("*** STEP 8. Prepare visualization in Artemis")
    print("=========================================================================================")
    missed_annotation_base_dir = create_out_dir(error_dir, 'missed_annotations')
    
    artemis_save_project(str(genomes_list_process_file),
                         args.namegenes,
                         roi_raw_dir,
                         artemis_annotations_dir, 
                         missed_annotation_base_dir)
    
    # Look for any error TXT file in the error directory
    txt_files = glob.glob(os.path.join(missed_annotation_base_dir, "*.txt"))

    # If there is at least one .txt file, print its contents
    if txt_files:
        print("\n Missing annotations were detected:\n")
        for file in txt_files:
            print(f"{os.path.basename(file)}:")
            with open(file, "r") as f:
                print(f.read())
    else:
        print("No missing annotation files found.") 
    
    print("=========================================================================================")
    print("=========================================================================================")
    print("Finished AnnCX running")
    print("Open Artemis to start manual gene curation:")
    
    # Open Artemis automatically
    subprocess.run(['art'])
    
    
if __name__ == "__main__":
    main()
