# ******************************************************************************
# Creation of pie charts for glycosylation data analysis
# Author: Pedro Fortes González (refactored)
# Date: 2024-09-16
# ******************************************************************************

# ______________________________________________________________________________
# 0. Environment setup ----
# ______________________________________________________________________________

# Clean environment and set working directory
rm(list = ls())
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


# Load required libraries
library(ggplot2)      # Data visualization
library(ggrepel)      # Text label positioning
library(utils)        # Utilities

# Set seed for reproducibility
set.seed(0156576321)



# ______________________________________________________________________________
# 1. Define color palette and category order ----
# ______________________________________________________________________________

# Standard color palette for glycosylation types
glyco_colors <- c("Fucosylated" = "red", "Fucosialylated" = "cyan2", 
                  "Sialylated" = "purple", "Oligomannose" = "green4", 
                  "Other" = "orange")

# Standard order of categories
glyco_categories <- c("Fucosylated", "Fucosialylated", "Sialylated", "Oligomannose", "Other")


# ______________________________________________________________________________
# 2. Helper functions ----
# ______________________________________________________________________________

#' Load and prepare data for pie charts
#'
#' @param file_path Path to the CSV file
#' @return Processed dataframe with percentages and ordered factors
#' 
prepare_pie_data <- function(file_path) {
  # Read data
  df <- read.csv(file_path)
  
  # Calculate percentages
  df <- df %>%
    mutate(Value = factor(Value, levels = glyco_categories), 
           Percentage = (Count / sum(Count)) * 100)
  
  return(df)
  }

# ______________________________________________________________________________

#' Create pie chart with or without labels
#'
#' @param data Processed dataframe
#' @param technique Name of the isolation technique
#' @param entity_type Either "Peptides" or "Proteins"
#' @param show_labels Whether to show percentage labels
#' @param size Font size for labels (default: 10)
#' @return ggplot object with the pie chart
create_pie_chart <- function(data, technique, entity_type, show_labels = TRUE, size = 10) {
  
  # Base plot
  p <- ggplot(
    data, aes(x = "", y = Count, fill = Value)) + 
    geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y") +
    scale_fill_manual(values = glyco_colors) + theme_void() +
    labs(title = paste0("Glycosylated ", entity_type, " detected by ", 
                        technique, " (in Vesiclepedia)"), 
         subtitle = "PTM type and percentage", fill = "PTM Type") + 
    theme(legend.text = element_blank(), legend.title = element_blank(), legend.position = "none")
  
  # Add labels if requested
  if (show_labels) {p <- p + 
    geom_label_repel(aes(label = paste0(round(Percentage, 1), "% (n=", Count, ")")), 
                     position = position_stack(vjust = 0.5), size = size, 
                     direction = "both", box.padding = 2, label.padding = 1, 
                     force = 3, show.legend = FALSE)
  }
  return(p)
}

# ______________________________________________________________________________

#' Generate and save pie charts for a given dataset
#'
#' @param technique Isolation technique name
#' @param entity_type Either "Peptides" or "Proteins"
#' @param data_path Base path for data files
#' @param output_path Base path for output files
#' @param chart_index Chart number for file naming
generate_charts <- function(technique, entity_type, data_path, output_path, chart_index) {
  
  # Construct file paths
  if (entity_type == "Peptides") {
    file_suffix <- "PTM cluster.csv"
  } else {
    file_suffix <- "PTM_types_by_protein.csv"
  }
  
  input_file <- paste0(data_path, technique, "_POOLS_123_filtered_vcp_filtered_glyc_", file_suffix)
  
  # Load and prepare data
  chart_data <- prepare_pie_data(input_file)
  
  # Create charts with and without labels
  chart_with_labels <- create_pie_chart(data = chart_data, technique = technique, 
                                        entity_type = entity_type, show_labels = TRUE)
  
  chart_without_labels <- create_pie_chart(data = chart_data, technique = technique, 
                                           entity_type = entity_type, show_labels = FALSE)
  
  # Save charts
  output_file_labels <- paste0(output_path, "textbox/", 
                               sprintf("%02d%s_sectors_vcp_%s_glyc_%s_textbox.png", 
                                       chart_index, 
                                       ifelse(entity_type == "Proteins", "a", "b"), 
                                       tolower(substr(entity_type, 1, 4)), technique)
                               )
  
  output_file_no_labels <- paste0(output_path, "no_text/", 
                                  sprintf("%02d%s_sectors_vcp_%s_glyc_%s_notext.png", 
                                          chart_index, 
                                          ifelse(entity_type == "Proteins", "a", "b"), 
                                          tolower(substr(entity_type, 1, 4)), technique)
                                  )
  
  ggsave(filename = output_file_labels, plot = chart_with_labels, width = 16, height = 9)
  ggsave(filename = output_file_no_labels, plot = chart_without_labels, width = 16, height = 9)
  
  return(list(with_labels = chart_with_labels, without_labels = chart_without_labels))
}



# ______________________________________________________________________________
# 3. Define paths and techniques ----
# ______________________________________________________________________________

# Base paths
peptide_data_path <- "./data/input_processed/2_value_counts/vesiclepedia_glycosylated/peptides/"
protein_data_path <- "./data/input_processed/2_value_counts/vesiclepedia_glycosylated/proteins/"
output_base_path <- "./data/input_processed/3_figures/sector_diagrams/"

# Isolation techniques
techniques <- c("ExoGAG", "IP_CD9", "SEC", "UC")



# ______________________________________________________________________________
# 4. Generate all charts ----
# ______________________________________________________________________________
cat("\n===== GENERATING CHARTS FOR GLYCOSYLATED PEPTIDES =====\n")

# Process peptides charts
peptide_charts <- list()
for (i in seq_along(techniques)) {
  technique <- techniques[i]
  cat(paste0("\nProcessing charts for ", technique, " peptides...\n"))
  
  peptide_charts[[technique]] <- generate_charts(
    technique = technique, entity_type = "Peptides", data_path = peptide_data_path,
    output_path = output_base_path, chart_index = i
    )
  
  # Print chart to console
  print(peptide_charts[[technique]]$with_labels)
}

# ______________________________________________________________________________

cat("\n\n===== GENERATING CHARTS FOR GLYCOSYLATED PROTEINS =====\n")

# Process proteins charts
protein_charts <- list()
for (i in seq_along(techniques)) {
  technique <- techniques[i]
  cat(paste0("\nProcessing charts for ", technique, " proteins...\n"))
  
  protein_charts[[technique]] <- generate_charts(
    technique = technique, entity_type = "Proteins", data_path = protein_data_path,
    output_path = output_base_path, chart_index = i
    )
  
  # Print chart to console
  print(protein_charts[[technique]]$with_labels)
}

# ______________________________________________________________________________
cat("\n\nAll charts saved in the folders 'output/3_figures/sector_diagrams/no_text/ & 'output/3_figures/sector_diagrams/textbox/'\n")
# ______________________________________________________________________________
