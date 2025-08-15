####
# Created on: 03.08.2025
# Author: lahumada

# Coloring the cells using the ranking of the values of each column
# https://stackoverflow.com/questions/15962291/using-rank-across-columns-to-create-new-variable
#####


# Dependencies
library(stringr)
library(data.table)
library(ggplot2)
library(dplyr)
library(magrittr)


#' Generate heatmaps from identify_rearrangements feature
#'
#' Generates heatmaps to visualize which of the predicted annotated 
#' exons per gene corresponds to which exon of the reference genes. 
#' pident values are displayed in each corresponding cell of the plot. 
#' Cell colors are based on the rank of values within each column, where:
#'   - Blue indicates the highest rank (best pident) in that column
#'   - Red indicates the lowest rank (worst pident)
#' 
#' @param File_input_alignment (str): File path exonerate output file.
#' @param File_input_genes (str): File path TXT file with names of the genes predicted
#' @param Dir_output (str): Directory path to save SVG heatmap file
#' @param genome_name (str) Name of the genome being processed (Example: 'Macaca_mulatta')
#'
#' @details
#' 1. Load exonerate output as a dataframe
#' 2. Subset relevant columns
#'   - V1 = Annotation of prediction (query)
#'   - V2 = Annotation of reference (subject)
#'   - V3 = Percentage of identity from running exonerate: query V1 vs subject V2
#' 3. Split dataframe per gene
#' 4. Generate and save heatmap as SVG
#' @return
#' Invisibly returns NULL. The main output is the saved SVG heatmap file.
#'
#' @import stringr
#' @import data.table
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#'
#' @examples
#' \dontrun{
#' Exonerate_heatmaps("exonerate_output.txt", "output_directory/", "Homo_sapiens", 
#'                 "predicted_gene_exons", "reference_gene_exons")
#' }
Exonerate_heatmaps <- function(File_input_alignment, File_input_genes, Dir_output, genome_name){
  #### 1. Load the dataframe output from exonerate
  df = read.delim(File_input_alignment, header = FALSE, stringsAsFactors = FALSE)
  
  #### 2. Make a subset
  # V1 = Annotation of prediction
  # V2 = Annotation of reference
  # V3 = Percentage of identity from running exonerate: query V1 vs subject V2
  df_subset_raw = df[,c(1:3)]
  
  # For the rows with the same query and subject, only keep the alignment with max percent identity
  df_subset <- df_subset_raw %>%
    arrange(desc(V3)) %>%
    distinct(V1, V2, .keep_all = TRUE)
  
  # Get only the second decimal place, otherwise it doesn't fit into the squares in the heatmap
  df_subset$V3 = round(df_subset$V3, 2)
  
  #### 3. Load file with the names of genes
  genes = read.delim(File_input_genes, header = FALSE, stringsAsFactors = FALSE)
  
  #### 4. Split the dataframe per gene
  # list of dataframes
  list_of_df = list()
  
  for (i in (1:length(rownames(genes)))) {
    # Boolean for each gene name
    vector = grepl(genes[i,], df_subset[,1])
    # Split the dataframe using the boolean vectors and store dataframes in list
    list_of_df[[i]] = subset(df_subset, vector)
  }
  
  # Function to make a plot for each gene and save it in File_output directory
  heatmap_plotting = function(species_and_gene, dataframe, intercepts, score_index){
    
    # Make the plot
    heatmap <- ggplot(dataframe, aes(x = V1, y = V2)) +
      geom_tile(aes(fill = rank), colour = "black") +  # fill the cells with colours based on the rank values
      scale_fill_gradient2(low = "steelblue2", 
                           high = "palevioletred2", 
                           mid = "white", 
                           midpoint = (max(dataframe$rank)+min(dataframe$rank))/2,
                           name = "Exonerate\nPercentage\nof Identity \n",
                           breaks = c(min(dataframe$rank), 0.5, max(dataframe$rank)),
                           labels = c("Max", "0.5", "Min"),
                           limits = c(min(dataframe$rank), max(dataframe$rank))) +
      geom_text(aes(label= V3))  # display in the cells the corresponding percentage number

    # Add the theme formatting
    base_size <- 9
    heatmap + theme_bw() + 
      theme(panel.grid = element_blank()) + # draw a rectangle around the whole graph
      geom_hline(yintercept = intercepts, color = "black") +  # draws horizontal lines to separate the genes in the y axis
      geom_tile(data=dataframe[score_index,], fill="transparent", colour = "black", size = 1) +  # makes a black square around the values with the highest score per column in the x axis
      labs(title = paste(species_and_gene, "identify_rearrangements", sep = "_"), x = "Query:    Predicted", y = "Subject:    Reference") + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0)) + 
      theme(plot.title = element_text(size = 15, hjust = 0.5), 
            axis.ticks = element_blank(), 
            legend.text = element_text(size= 11),
            axis.title.x = element_text(size = 12,
                                        margin = margin(t=12)),
            axis.title.y = element_text(size = 12),
            axis.text.y = element_text(size = 10, 
                                       angle = 0, hjust = 0.5, 
                                       colour = "grey50",
                                       margin = margin(t=12)),
            axis.text.x = element_text(size = 10, 
                                       angle = 45,  # Change this value to adjust the tilt angle
                                       hjust = 1,   # Adjust horizontal justification
                                       vjust = 1,   # Adjust vertical justification
                                       colour = "grey50"))
    
  }
  
  # Function to get the intercepts to draw horizontal lines between the different genes in y axis
  get_intercepts = function(dataframe){
    column = dataframe$V2
    exon_names = levels(as.factor(column))
    exon_names = gsub("_e.*", "", exon_names, ignore.case = TRUE)
    exon_unique = unique(exon_names)
    
    intercepts = c()
    for (i in (1:length(exon_unique))){
      position = grep(exon_unique[i], exon_names)[1]
      intercepts[i] = position - 0.5
    }
    return(intercepts)
  }

  # Function to get the row index of the max raw scores per exon
  get_max_scores = function(dataframe){
    max_score = dataframe[, list(max.score=max(V3)), by=V1]
    index = list()
    for (i in (1:length(max_score$V1))){
      position = which(dataframe$V1==max_score$V1[i] & dataframe$V3==max_score$max.score[i])
      index = append(index, list(position))
    }
    index = as.numeric(unlist(index))
    index = unique(index)
    return(index)
  }
  
  #### 7. Run heatmap function
  for (i in (1:length(list_of_df))){
    # Remove duplicated values
    list_of_df[[i]] = unique(list_of_df[[i]])
    
    # Transform the dataframes in the list_of_genes into data.table format
    list_of_df[[i]] = data.table(V1 = list_of_df[[i]][,1],
                                 V2 = list_of_df[[i]][,2], 
                                 V3 = list_of_df[[i]][,3])
    # Add a column at the end of each data.table with the ranking of the values per column
    list_of_df[[i]][, rank := rank(-V3), by = V1]
    
    # Create the file directory for each gene
    gene_name = genes$V1[sapply(genes$V1, function(g) any(grepl(g, list_of_df[[i]]$V1, fixed = TRUE)))]
    species_and_gene = paste(genome_name, gene_name, sep = "_")
    File_name = paste(species_and_gene, "identify_rearrangements.svg", sep = "_")
    File_output = paste(Dir_output, File_name, sep = "/")
    
    intercepts = get_intercepts(dataframe = list_of_df[[i]])
    score_index = get_max_scores(dataframe = list_of_df[[i]])
    
    # Open svg to save the plot
    svg(filename = File_output)
    
    # Run function to make the heatmap
    plot = heatmap_plotting(species_and_gene = species_and_gene, dataframe=list_of_df[[i]], intercepts = intercepts, score_index = score_index)
    print(plot)
    
    # Close the plotting window to finish saving the plot
    dev.off() 
  }
}