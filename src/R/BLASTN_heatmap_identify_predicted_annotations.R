####
# Created on: 03.08.2025
# Author: lahumada
# Edited by: Merve Nida Bastuerk

# Coloring the cells using the ranking of the values of each column
# https://stackoverflow.com/questions/15962291/using-rank-across-columns-to-create-new-variable
#####


# Dependencies
library(stringr)
library(data.table)
library(ggplot2)
library(dplyr)
library(magrittr)


#' Generate heatmaps from identify_pred2ref feature
#'
#' Generates heatmaps to visualize which of the predicted annotated 
#' genes corresponds to which kind of reference gene. 
#' pident values are displayed in each corresponding cell of the plot. 
#' Cell colors are based on the rank of values within each column, where:
#'   - Blue indicates the highest rank (best match) in that column
#'   - Red indicates the lowest rank (worst match)
#' 
#' @param File_input_blast (str): File path blastn output file.
#' @param Dir_output_heatmaps (str): Directory path to save SVG heatmap file
#' @param Dir_output_report (str): Directory path to save TXT report file
#' @param genome_name (str) Name of the genome being processed (Example: 'Macaca_mulatta')
#' @param gene_level_query (str) Gene level descriptor for the query (predicted genes) (Example: cDNA)
#' @param gene_level_subject (str) Gene level descriptor for the subject (reference genes) (Example: cDNA)
#'
#' @details
#' 1. Load blastn output as a dataframe
#' 2. Subset relevant columns
#'   - V1 = Annotation of prediction (query)
#'   - V2 = Annotation of reference (subject)
#'   - V3 = Percentage of identity from running blastn: query V1 vs subject V2
#'   - V6 = gapopen
#'   - V12 = Bitscore
#'   - V13 = qcovs
# 3. Compute rankings and scores to identify best matches.
# 4. Generate and save heatmap as SVG
# 5. Generate and save report as TXT
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
#' BLASTN_heatmaps("blast_output.txt", "output_directory_heatmaps/", "output_directory_report/", "Homo_sapiens", 
#'                 "cDNA", "cDNA")
#' }
BLASTN_heatmaps <- function(File_input_blast, Dir_output_heatmaps, Dir_output_report, genome_name, gene_level_query, gene_level_subject){
  
  #### Load the dataframe output from blastn
  df = read.delim(File_input_blast, header = FALSE, stringsAsFactors = FALSE)
  
  #### Make a subset
  # V1 = Annotation of prediction
  # V2 = Annotation of reference
  # V3 = Percentage of identity from running blastn: query V1 vs subject V2
  # V6 = gapopen
  # V12 = Bitscore
  # V13 = qcovs
  df_subset_raw = df[,c(1:3, 6, 12, 13)]
  df_subset_raw$original_order <- seq_len(nrow(df_subset_raw))
  
  # For the rows with the same query and subject, only keep the alignment with max bitscore
  df_subset <- df_subset_raw %>%
    arrange(desc(V3), V6, desc(V12), original_order) %>%
    distinct(V1, V2, .keep_all = TRUE)
  
  # Arrange rows in order
  df_subset <- df_subset %>%
    arrange(original_order)
  
  # Get only the second decimal place, otherwise it doesn't fit into the squares in the heatmap
  df_subset$V3 = round(df_subset$V3, 2)
  
  #### Function to make heatmap
  # Make a plot for each genome and save it in File_output directory
  heatmap_plotting = function(species, dataframe, best_matches) {
    # Calculate the number of rows and columns
    n_rows <- length(unique(dataframe$V2))
    n_cols <- length(unique(dataframe$V1))
    
    # Make the plot
    heatmap_plot <- ggplot(dataframe, aes(x = V1, y = V2)) +
      geom_tile(aes(fill = rank), colour = "black") +
      scale_fill_gradient2(low = "steelblue2", 
                           high = "palevioletred2", 
                           mid = "white", 
                           midpoint = mean(dataframe$rank),
                           name = "BLASTN\nPercentage\nof Identity\n(Cell value)\n\nOverall Quality\n(Cell color)\n",
                           breaks = c(min(dataframe$rank), max(dataframe$rank)),
                           labels = c("Best", "Worst"),
                           limits = c(min(dataframe$rank), max(dataframe$rank))) +
      geom_text(aes(label = round(V3,2))) +
      geom_tile(data = best_matches, aes(x = V1, y = V2), 
                fill = "transparent", colour = "black", size = 1)
    
    # Add the theme formatting
    base_size <- 9
    heatmap_plot <- heatmap_plot + theme_bw() + 
      theme(aspect.ratio = 1) +  # This adjusts the overall plot aspect ratio
      theme(panel.border = element_blank(), panel.grid = element_blank()) +
      labs(title = species, x = paste("Query:    Predicted", paste("(",gene_level_query,")", sep=""), sep=" "), y = paste("Subject:    Reference", paste("(",gene_level_subject,")", sep=""), sep=" ")) + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_discrete(expand = c(0, 0)) + 
      theme(plot.title = element_text(size = 15, hjust = 0.5), 
            axis.ticks = element_blank(), 
            legend.text = element_text(size= 11),
            axis.title.x = element_text(size = 12,
                                        margin = margin(t=12)),
            axis.title.y = element_text(size = 12),
            axis.text.y = element_text(size = 10, 
                                       angle = 0, hjust = 0.5, colour = "grey50",
                                       margin = margin(t=12)),
            axis.text.x = element_text(size = 10, 
                                       angle = 45,
                                       hjust = 1,
                                       vjust = 1,
                                       colour = "grey50"))
    
    # Return the plot
    return(heatmap_plot)
  }
  
  # Function to get the row index of the max raw scores per match
  get_max_scores = function(dataframe) {
    dataframe[, score := V3 + V13/100 - V6/max(V6) + V12/max(V12)]
    best_matches = dataframe[, .SD[score == max(score)], by = V1]
    return(best_matches)
  }
  
  # Transform the dataframe into data.table format
  df_data_table = data.table(V1 = df_subset[,1], V2 = df_subset[,2], 
                             V3 = df_subset[,3], V6 = df_subset[,4], 
                             V12 = df_subset[,5], V13 =df_subset[,6])
  
  # Calculate rank based on all three criteria
  df_data_table[, rank := rank(-V3 -V13/100 + V6/max(V6) - V12/max(V12)), by = V1]
  
  # Get the best matches
  best_matches = get_max_scores(df_data_table)
  
  # Create file directory
  File_name = paste(genome_name, "_identify_pred2ref.svg", sep = "")
  File_output = paste(Dir_output_heatmaps, File_name, sep = "/")
  
  # Open svg to save the plot
  svg(filename = File_output)
  
  #### Run heatmap function
  plot = heatmap_plotting(species = genome_name, dataframe = df_data_table, best_matches = best_matches)
  print(plot)
  
  # Close the plotting window to finish saving the plot
  dev.off()
  
  ##########################
  # Generate genotype report
  best_out_file <- file.path(
    Dir_output_report,
    paste0(genome_name, "_identify_pred2ref_best_matches.txt")
  )
  
  best_matches_to_write <- best_matches[, .(query = V1, match = V2, pident = V3, gapscore = V6, bitscore = V12, qcovs = V13, score = score)]
  
  fwrite(
    best_matches_to_write,
    file = best_out_file,
    sep = "\t",
    col.names = FALSE
  )
}