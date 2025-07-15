# ******************************************************************************
# Statistical analysis and visualization for mass spectrometry data
# Version: 1.0.0
# Last updated: 2025-07-15
#
# This script plots sector diagrams on the main glycosylation groups that
# populate the samples in our analysis.
#
# Author: Pedro Fortes González
# ******************************************************************************
# ******************************************************************************

# ______________________________________________________________________________
# 0. Environment setup ----
# ______________________________________________________________________________

# Clean WEnv
rm(list = ls())

# Prevent automatic creation of Rplots.pdf
pdf(NULL)

# Helper function for user interaction w/ the console
read_input <- function(prompt) {
  cat(prompt)
  return(readLines("stdin", n=1))
}

# Set WD
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  # RStudio
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
  script_dir <- dirname(getActiveDocumentContext()$path)
  repo_root <- dirname(script_dir)
  
} else {
  # CLI
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
  script_dir <- dirname(script_path)
  repo_root <- dirname(script_dir)
  setwd(script_dir)
  
}

cat("\nScript's dir:", script_dir, "\n")
cat("Repo's dir:", repo_root, "\n")


# ______________________________________________________________________________
# Check and install dependencies if needed ----
# ______________________________________________________________________________
check_and_install_dependencies <- function(required_packages) {
  cat("\nChecking R dependencies...\n")
  
  missing_packages <- c()
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      missing_packages <- c(missing_packages, package)
    }
  }
  
  if (length(missing_packages) > 0) {
    cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
    
    # Ask user if they want to install missing packages
    ## opt 1
    #install_choice <- readline(prompt = "Do you want to install these packages? (y/n): ")
    ## opt 2
    install_choice <- read_input("Do you want to install these packages? (y/n): ")
    print(paste("You selected:", install_choice))
    
    if (tolower(install_choice) == "y") {
      cat("Installing packages...\n")
      for (pkg in missing_packages) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
      cat("Installation completed. Loading packages...\n")
    } else {
      cat("Required packages were not installed. The script may fail.\n")
      return(FALSE)
    }
  } else {
    cat("All dependencies are installed.\n")
  }
  
  return(TRUE)
}

# Define required packages for this script
required_packages <- c("dplyr", "ggplot2", "ggrepel", "utils")

# Check and install dependencies
dependencies_ok <- check_and_install_dependencies(required_packages)
if (!dependencies_ok) {
  cat("ERROR: Could not load all dependencies.\n")
  quit(status = 1)
}


# Load required libraries
library(ggplot2)      # Data visualization
library(ggrepel)      # Text label positioning
library(utils)        # Utilities
library(dplyr)        # Data manipulation

# Set seed for reproducibility
set.seed(0156576321)

# ______________________________________________________________________________
# 1. Analysis Type Selection ----
# ______________________________________________________________________________

# Function to get analysis choice
get_analysis_choice <- function() {
  # Check if config file exists
  if (file.exists("analysis_config.txt")) {
    # Read choice from file
    saved_choice <- readLines("analysis_config.txt", n=1)
    cat("Found saved analysis choice:", saved_choice, "\n")
    return(saved_choice)
  }
  
  # If no file exists
  cat("\n==========================================")
  cat("\nSelect analysis type:")
  cat("\n1. Technique analysis (default)")
  cat("\n2. Pooled analysis")
  cat("\n==========================================\n")
  return(read_input("Enter your choice (1/2): "))
}

# Get the analysis choice
analysis_choice <- get_analysis_choice()

# Set analysis parameters based on user choice
if (analysis_choice == "2") {
  ANALYSIS_TYPE <- "Pooled"
  GROUPING_VAR <- "Pool"
  cat("\nSelected POOLED analysis - Grouping by Pool\n")
} else {
  ANALYSIS_TYPE <- "Technique"
  GROUPING_VAR <- "Technique"
  cat("\nSelected TECHNIQUE analysis - Grouping by Technique\n")
}


# ______________________________________________________________________________
# 2. Define color palette and category order ----
# ______________________________________________________________________________

# Standard color palette for glycosylation types
glyco_colors <- c("Fucosylated" = "red", "Fucosialylated" = "cyan2", 
                  "Sialylated" = "purple", "Oligomannose" = "green4", 
                  "Other" = "orange")

# Standard order of categories
glyco_categories <- c("Fucosylated", "Fucosialylated", "Sialylated", "Oligomannose", "Other")


# ______________________________________________________________________________
# 3. Helper functions ----
# ______________________________________________________________________________

#' Load and prepare data for pie charts
#'
#' @param file_path Path to the CSV file
#' @return Processed dataframe with percentages and ordered factors
#' 
prepare_pie_data <- function(file_path) {
  # Check if file exists
  if (!file.exists(file_path)) {
    cat("\nWARNING: File not found:", file_path, "\n")
    # Try to find a similar file
    base_dir <- dirname(file_path)
    base_name <- basename(file_path)
    pattern <- gsub("_POOLS_\\d+_", ".*", basename(file_path))
    matching_files <- list.files(base_dir, pattern = pattern, full.names = TRUE)
    
    if (length(matching_files) > 0) {
      file_path <- matching_files[1]
      cat("Using alternative file:", file_path, "\n")
    } else {
      cat("No alternative files found. Cannot proceed.\n")
      return(NULL)
    }
  }
  
  # Read data
  df <- tryCatch({
    read.csv(file_path, sep=";")
  }, error = function(e) {
    cat("\nERROR reading file:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
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
#' @param grouping_var Name of the grouping_var derived from analysis_config.txt
#' @param entity_type Either "Peptides" or "Proteins"
#' @param show_labels Whether to show percentage labels
#' @param size Font size for labels (default: 10)
#' @return ggplot object with the pie chart
create_pie_chart <- function(data, grouping_var, entity_type, show_labels = TRUE, size = 10) {
  # Check if data is valid
  if (is.null(data) || nrow(data) == 0) {
    cat("\nNo data available for", grouping_var, entity_type, "\n")
    return(NULL)
  }
  
  # Base plot
  p <- ggplot(
    data, aes(x = "", y = Count, fill = Value)) + 
    geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y") +
    scale_fill_manual(values = glyco_colors) + theme_void() +
    labs(title = paste0("Glycosylated ", entity_type, " detected by ", 
                        grouping_var, " (in Vesiclepedia)"), 
         subtitle = "PTM type and percentage", fill = "PTM Type") + 
    theme(legend.text = element_blank(), legend.title = element_blank(), legend.position = "none")
  
  # Add labels if requested
  if (show_labels) {
    p <- p + 
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
#' @param grouping_var grouping_var name
#' @param entity_type Either "Peptides" or "Proteins"
#' @param data_path Base path for data files
#' @param output_path Base path for output files
#' @param chart_index Chart number for file naming
generate_charts <- function(grouping_var, entity_type, data_path, output_path, chart_index) {
  
  # Construct file paths
  if (entity_type == "Peptides") {
    file_suffix <- "PTM_Cluster.csv"
  } else {
    file_suffix <- "PTM_types_by_protein.csv"
  }
  
  # UPDATED: Use a more flexible approach to find the right files
  # Try the most specific pattern first
  pattern <- paste0(".*", grouping_var, "_filtered_vcp_filtered_glyc_", file_suffix)
  matching_files <- list.files(data_path, pattern = pattern, full.names = TRUE)
  
  if (length(matching_files) > 0) {
    input_file <- matching_files[1]
    cat("\nFound matching file:", basename(input_file), "\n")
  } else {
    # Fall back to original pattern (which will likely fail but preserve original behavior)
    input_file <- paste0(data_path, grouping_var, "_POOLS_123_filtered_vcp_filtered_glyc_", file_suffix)
    cat("\nTrying default file pattern:", basename(input_file), "\n")
  }
  
  # Rest of function remains unchanged
  chart_data <- prepare_pie_data(input_file)
  
  # Skip if data loading failed
  if (is.null(chart_data)) {
    cat("\nSkipping chart generation for", grouping_var, entity_type, "- data loading failed\n")
    return(NULL)
  }
  
  # Create charts with and without labels
  chart_with_labels <- create_pie_chart(data = chart_data, grouping_var = grouping_var, 
                                        entity_type = entity_type, show_labels = TRUE)
  
  chart_without_labels <- create_pie_chart(data = chart_data, grouping_var = grouping_var, 
                                           entity_type = entity_type, show_labels = FALSE)
  
  # Skip if chart creation failed
  if (is.null(chart_with_labels) || is.null(chart_without_labels)) {
    cat("\nSkipping chart saving for", grouping_var, entity_type, "- chart creation failed\n")
    return(NULL)
  }
  
  # Create output directories if needed
  dir.create(file.path(output_path, "textbox"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_path, "no_text"), recursive = TRUE, showWarnings = FALSE)
  
  # Save charts
  output_file_labels <- paste0(output_path, "textbox/", 
                               sprintf("%02d%s_sectors_vcp_%s_glyc_%s_textbox.png", 
                                       chart_index, 
                                       ifelse(entity_type == "Proteins", "a", "b"), 
                                       tolower(substr(entity_type, 1, 4)), grouping_var)
  )
  
  output_file_no_labels <- paste0(output_path, "no_text/", 
                                  sprintf("%02d%s_sectors_vcp_%s_glyc_%s_notext.png", 
                                          chart_index, 
                                          ifelse(entity_type == "Proteins", "a", "b"), 
                                          tolower(substr(entity_type, 1, 4)), grouping_var)
  )
  
  tryCatch({
    ggsave(filename = output_file_labels, plot = chart_with_labels, width = 16, height = 9)
    cat("\nSaved chart with labels to:", output_file_labels)
    
    ggsave(filename = output_file_no_labels, plot = chart_without_labels, width = 16, height = 9)
    cat("\nSaved chart without labels to:", output_file_no_labels)
  }, error = function(e) {
    cat("\nERROR saving chart:", e$message, "\n")
  })
  
  return(list(with_labels = chart_with_labels, without_labels = chart_without_labels))
}



# ______________________________________________________________________________
# 4. Define paths and grouping_levels (values of grouping_var) ----
# ______________________________________________________________________________

# Base paths
peptide_data_path <- "../output/2_value_counts/vesiclepedia_glycosylated/peptides/"
protein_data_path <- "../output/2_value_counts/vesiclepedia_glycosylated/proteins/"
output_base_path <- "../output/3_figures/sector_diagrams/"

# grouping_levels (from analysis results)
grouping_file_path <- "../output/1_filtered_dfs/vesiclepedia/summary_all.csv"
grouping_df <- read.csv(grouping_file_path, sep = ";", fileEncoding = "UTF-8")

# Convert columns to factor
#grouping_df[, GROUPING_VAR] <- lapply(grouping_df[, GROUPING_VAR], factor)
#grouping_levels <- unique(grouping_df[, GROUPING_VAR])
# Convert to factor directly
grouping_df[[GROUPING_VAR]] <- factor(grouping_df[[GROUPING_VAR]])
grouping_levels <- levels(grouping_df[[GROUPING_VAR]])



# ______________________________________________________________________________
# 5. Generate all charts ----
# ______________________________________________________________________________
cat("\n===== GENERATING CHARTS FOR GLYCOSYLATED PEPTIDES =====\n")

# Process peptides charts
peptide_charts <- list()
for (i in seq_along(grouping_levels)) {
  grouping_level <- grouping_levels[i]
  cat(paste0("\nProcessing charts for ", grouping_level, " peptides...\n"))
  
  peptide_charts[[grouping_level]] <- generate_charts(
    grouping_var = grouping_level, 
    entity_type = "Peptides", 
    data_path = peptide_data_path,
    output_path = output_base_path, 
    chart_index = i
  )
  
  # Print chart to console
  if (!is.null(peptide_charts[[grouping_level]])) {
    print(peptide_charts[[grouping_level]]$with_labels)
  }
}

# ______________________________________________________________________________

cat("\n\n===== GENERATING CHARTS FOR GLYCOSYLATED PROTEINS =====\n")

# Process proteins charts
protein_charts <- list()
for (i in seq_along(grouping_levels)) {
  grouping_level <- grouping_levels[i]
  cat(paste0("\nProcessing charts for ", grouping_level, " proteins...\n"))
  
  protein_charts[[grouping_level]] <- generate_charts(
    grouping_var = grouping_level, 
    entity_type = "Proteins", 
    data_path = protein_data_path,
    output_path = output_base_path, 
    chart_index = i
  )
  
  # Print chart to console
  if (!is.null(protein_charts[[grouping_level]])) {
    print(protein_charts[[grouping_level]]$with_labels)
  }
}

# ______________________________________________________________________________
cat("\n\nAll charts saved in the folders 'output/3_figures/sector_diagrams/no_text/ & 'output/3_figures/sector_diagrams/textbox/'\n")
# ______________________________________________________________________________