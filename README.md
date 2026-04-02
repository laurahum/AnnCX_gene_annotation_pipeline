# AnnCX_gene_annotation_pipeline (v1.0.0)

AnnCX is gene annotation pipeline designed for the analysis of gene-rich complex genomic regions. AnnCX automates all intermediate steps to provide comprehensive results with minimal user input. Key novel features include:
- Iterative extraction of target regions from whole-genome sequences to focus on complex genomic areas (OPTIONAL).
- Incorporate a 7+ diverse array of individual annotation tools with an emphasis on the accurate identification of exon-intron boundaries.
- Genome-wide identification of query genes to detect possible genes located outside the region of interest (OPTIONAL).
- Detection and annotation of ambiguous nucleotide regions (assembly gaps) to provide a clear representation of problematic genome assemblies that commonly affect complex regions regions.
- Support for manual curation by implementing an automatic visualization of the raw annotation results once the pipeline has completed its execution.
<br>

Custom supplementary tools to assess the quality of genome assembly: 
- **identify_pred2ref**: to help aid on the gene identification process
- **identify_rearrangements**: to help identify exon-level rearrangements.

Other custom supplementarty tools:
- **annotation2fasta**: convert annotation GFF3 files into FASTA sequences
<br>

![Alt text](https://github.com/laurahum/AnnCX_gene_annotation_pipeline/blob/main/figures/overview_AnnCX.png?raw=true)
&emsp;*Figure 1: Overview of AnnCX's structure*
###### *User-required input files are highlighted in red.*
<br>

![Alt text](https://github.com/laurahum/AnnCX_gene_annotation_pipeline/blob/main/figures/artemis_example_output.png?raw=true)
&emsp;*Figure 2: AnnCX annotation output (Artemis)*
###### *Annotation outputs produced by AnnCX are depicted within grey lines arranged in a mirrored fashion above and below the genomic sequence, according to the orientation of each annotated feature. The annotation entry corresponding to the consensus annotation produced by AnnCX is highlighted in green.*
<br>

## Installation

1. Clone this repository:
   
   `git clone https://github.com/laurahum/AnnCX_gene_annotation_pipeline.git`     # Link to repository
   
   `cd path/to/pipeline/folder`     # Directory of the repository folder

2. Installation

   `chmod +x install_AnnCX.sh`    # Give rights to the installation script
   
   `bash install_AnnCX.sh`        # Run installation script
<br>

## Usage for:
### - Genome annotation

1. Activate the conda environment where the pipeline is installed:
   
   `conda activate AnnCX`

2. Run the pipeline:

    *NOTE: For best compatibility, name FASTA files and sequence headers using simple names without dots (.), special characters, or spaces*
    
    *NOTE: Please ensure that Artemis is closed before running the pipeline. If Artemis is open, the genome files will be annotated, but the new annotations may not appear in the Artemis interface*

   Example using the example data provided in the repository:
   
   ```
   (AnnCX) user@computer:~/path/to/pipeline/folder$ ./AnnCX.py \
   --genome examples/genomic_sequences/genome \
   --namegenes NKG2 \
   --querytranscript examples/Genome_annotation/transcript_sequences.fasta \
   --queryprot examples/Genome_annotation/protein_sequences.fasta \
   --queryexon examples/Genome_annotation/exon_sequences.fasta \
   --spsrepeatmasker primates \
   --spsaugustus human \
   --outdir /path/to/output/folder/ \
   --flanking examples/Genome_annotation/flanking_regions.fasta \  # Step 1 - Optional extract region of interest
   --WGannotation \  # Step 0 - Optional to search for query genes genome-wide
   --threads your_threads
   ```
   <br>
   
   Step 1: extract ROI = OPTIONAL - If the user provides flanking regions, AnnCX checks how many genome files contain both flanking genes and, by default, prompts whether to continue extract ROI and annotate on those genomes. Use `--skip-prompt` to proceed automatically with the extraction without asking for confirmation.
   
   Example skip step 1:
      ```
   (AnnCX) user@computer:~/path/to/pipeline/folder$ ./AnnCX.py \
   --genome examples/genomic_sequences/genomic_ROI \
   --namegenes NKG2 \
   --querytranscript examples/Genome_annotation/transcript_sequences.fasta \
   --queryprot examples/Genome_annotation/protein_sequences.fasta \
   --queryexon examples/Genome_annotation/exon_sequences.fasta \
   --spsrepeatmasker primates \
   --spsaugustus human \
   --outdir /path/to/output/folder/ \
   --threads your_threads
   ```
   <br>
   
   
   Step 3: Hardmasking ROI = OPTIONAL - If the user provides an already hardmasked genomic sequence and a repeat annotation file.
   
   Example skip step 3:
   ```
   (AnnCX) user@computer:~/path/to/pipeline/folder$ ./AnnCX.py \
   --genome examples/genomic_sequences/genomic_ROI \
   --genome examples/genomic_sequences/genomic_ROI_masked \
   --namegenes NKG2 \
   --querytranscript examples/Genome_annotation/transcript_sequences.fasta \
   --queryprot examples/Genome_annotation/protein_sequences.fasta \
   --queryexon examples/Genome_annotation/exon_sequences.fasta \
   --skipRepeatmasker \
   --repeatannotations examples\Genome_annotation\repeat_annotations \ # Optional
   --spsaugustus human \
   --outdir /path/to/output/folder/ \
   --threads your_threads
   ```
   <br>
   
3. Open results in Artemis: 

   Artemis is opened automatically for the visualization of the annotation output produced by AnnCX. Use --skipCreateArtemis to skip the creation of an Artemis project and --skipOpenArtemis to skip openning Artemis automatically at the end.
   
   After the pipeline has finished running and you have closed Artemis, you can reopen Artemis within the conda environment to visualize the annotation results again:
  
   `(AnnCX) user@computer:~/path/to/pipeline/folder$ art`


<br>

### - Identify predicted genes:   identify_pred2ref feature

This feature produces a heatmap plot (SVG) and report for the set of genes predicted in an annotated genome.

1. Activate the conda environment:
   
   `conda activate AnnCX`

2. Run the tool:
   
   Example using the example data provided in the repository:
   
   ```
   (AnnCX) user@computer:~/path/to/pipeline/folder$ ./src/identify_pred2ref.py \
   --subject examples/Genome_annotation/cDNA_sequences.fasta \
   --query examples/Genome_annotation/cDNA_sequences.fasta \
   --typeseq_query cDNA \
   --typeseq_subject cDNA \
   --namegenome Macaca_mulatta \
   --outdir /path/to/output/folder \
   ```
<br>

![Alt text](https://github.com/laurahum/AnnCX_gene_annotation_pipeline/blob/main/figures/Identify_pred2ref.png?raw=true)
&emsp;*Figure 3: identify_pred2ref feature output heatmap plot*
###### *Predicted genes (x-axis) against reference genes (y-axis). Columns show alignment results for each predicted gene. Individual cells represent the match quality for each pair-wise alignment and indicate BLAST’s percentage identity. Each column is colored based on a composite match score using percentage of identity, coverage, bit-score and penalized by gap openings, in a gradient red→white→blue (from lowerst to highest). Cells highlighted in black show which reference gene is most similar.*
<br>

### - Identify rearrangements:   identify_rearrangements feature

This feature produces one heatmap plot (SVG) per gene predicted in an annotated genome.

1. Activate the conda environment:
   
   `conda activate AnnCX`

2. Run the tool:
   
   Example using the example data provided in the repository:
   
   ```
   (AnnCX) user@computer:~/path/to/pipeline/folder$ ./src/identify_rearrangements.py \
   --subject examples/Genome_annotation/exon_sequences.fasta \
   --query examples/Genome_annotation/exon_sequences.fasta \
   --namegenes examples/Identify_rearrangements/gene_names.txt
   --namegenome Macaca_mulatta \
   --outdir /path/to/output/folder \
   --threads your_threads
   ```
<br>

![Alt text](https://github.com/laurahum/AnnCX_gene_annotation_pipeline/blob/main/figures/Identify_rearrangements.png?raw=true)
&emsp;*Figure 4: identify_rearrangements feature output heatmap plot*
###### *Exons of one of the predicted genes (x-axis) against exons of reference genes (y-axis). Columns represent the alignment results of each predicted exon and are colored using Exonerate’s percentage of identity to generate a red→white→blue gradient (from lowest to highest). Exons from each reference genes are separated by horizontal lines. Each cell displays pair-wise Exonerate identity score and those with the highest identity per exon are highlighted with a black square.*
<br>

### - Convert annotation to fasta:   annotation2fasta feature   

1. Activate the conda environment:
   
   `conda activate AnnCX`

2. Run the tool:

    *NOTE: Input annotation files (--annotation) must be in GFF3 format with annotation features: gene, mRNA, exon, CDS*

    *NOTE: Input FASTA file with the genomic region (--genome) that was annotated. If the annotations are produced by AnnCX using flanking genes, give the extracted ROI, not the whole genome FASTA (/path/to/output/AnnCX_genenames/1_extract_ROI/7_extracted_roi_raw)*
   
   Example using the example data provided in the repository:
   
   ```
   (AnnCX) user@computer:~/path/to/pipeline/folder$ ./src/annotation2fasta.py \
   --annotation examples/annotate2fasta/GFF3 \
   --genome examples/genomic_sequences/genome \
   --txtgenome examples/annotate2fasta/TXT/txt_genome.txt \
   --nameproject extract_annotation_NKG2 \
   --outdir /path/to/output/folder
   ```
<br>

## Requirements

- Conda (Miniconda or Anaconda)
- Linux operating system
- AnnCX has many dependencies, so please allow sufficient time for installation

<br>

________________________________________________________________________________
Full documentation including output files explanation: 
[AnnCX_documentation](https://laurahum.github.io/AnnCX_gene_annotation_pipeline/)
