# ******************************************************************************
# Statistical analysis and visualization for mass spectrometry data
# Version: 1.0.0
# Last updated: 2025-07-15
#
# This script performs non-parametric statistical analysis (Kruskal-Wallis test
# followed by Dunn's post-hoc tests with FDR correction) on processed mass 
# spectrometry data to compare glycomic profiles across isolation techniques.
# Non-parametric tests were selected due to the typical non-normal distribution
# of peptide/protein count data.
#
# Author: Pedro Fortes González
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
required_packages <- c("dplyr", "ggplot2", "ggpubr", "rstatix", "envalysis", "gridExtra")

# Check and install dependencies
dependencies_ok <- check_and_install_dependencies(required_packages)
if (!dependencies_ok) {
  cat("ERROR: Could not load all dependencies.\n")
  quit(status = 1)
}

# Load required libraries
library(dplyr)      # Data manipulation
library(ggplot2)    # Visualization
library(ggpubr)     # Extensions for ggplot2
library(rstatix)    # Statistical tests
library(envalysis)  # Publication theme
library(gridExtra)  # Graphics organization

# Set seed for reproducibility
set.seed(0156576321)

# ______________________________________________________________________________
# 1. Analysis Setting ----
# ______________________________________________________________________________

# Analysis & Plot Type Selection ----

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

# Set analysis parameters and plot type based on user choice
if (analysis_choice == "2") {
  ANALYSIS_TYPE <- "Pooled"
  GROUPING_VAR <- "Pool"
  PLOT_TYPE <- "both"
  cat("\nSelected POOLED analysis - Grouping by Pool")
  cat("\nVisualization: BOXPLOTS, BARPLOTS and SECTOR DIAGRAMS")
} else {
  ANALYSIS_TYPE <- "Technique"
  GROUPING_VAR <- "Technique"
  PLOT_TYPE <- "boxplots"
  cat("\nSelected TECHNIQUE analysis - Grouping by Technique\n")
  cat("\nVisualization: BOXPLOTS and SECTOR DIAGRAMS")
}

# Create output directory for barplots if needed
if (PLOT_TYPE != "boxplots") {
  dir.create("../output/3_figures/barplots", recursive = TRUE, showWarnings = FALSE)
}

# ______________________________________________________________________________
# Function to reorder factor levels moving the last to first position
# ______________________________________________________________________________
reorder_technique_factor <- function(df) {
  # Only apply if we're doing technique analysis and analysis_choice is "1"
  if (ANALYSIS_TYPE == "Technique" && analysis_choice == "1") {
    # Get the levels of the Technique factor
    tech_levels <- levels(df$Technique)
    
    if (length(tech_levels) > 1) {
      # Debug output
      cat("\nOriginal Technique levels:", paste(tech_levels, collapse=", "), "\n")
      
      # Move last level to first position
      new_order <- c(tech_levels[length(tech_levels)], tech_levels[-length(tech_levels)])
      
      # Apply the new order
      df$Technique <- factor(df$Technique, levels = new_order)
      
      # Debug output
      cat("Reordered Technique levels:", paste(levels(df$Technique), collapse=", "), "\n")
    }
  }
  return(df)
}

# ______________________________________________________________________________
# 2. Helper functions ----
# ______________________________________________________________________________

#' Perform Kruskal-Wallis test and visualize result
#'
#' @param data DataFrame with data
#' @param formula Formula for the test (ex: n.Peptides.Total ~ Technique)
#' @return Test result as a table
perform_kruskal_test <- function(data, formula) {
  test_result <- kruskal_test(data, formula)
  table_result <- tableGrob(test_result)
  grid.arrange(table_result)
  return(test_result)
}

# ______________________________________________________________________________

#' Create boxplot with pairwise comparisons
#'
#' @param data DataFrame with data
#' @param y_var Name of the y-axis variable (string)
#' @param grouping_var Name of the grouping variable (string) - "Technique" or "Pool"
#' @param title Plot title
#' @param y_label Y-axis label
#' @param y_breaks Vector with y-axis breaks (optional)
#' @param y_limits Vector with y-axis limits (optional)
#' @param output_path Path to save the plot
#' @param export_tests Whether to export test results (default: TRUE)
#' @param export_name Custom name for the exported test results dataframe (optional)
#' @return ggplot object
create_boxplot_with_stats <- function(data, y_var, grouping_var = GROUPING_VAR, title, y_label, 
                                      y_breaks = NULL, y_limits = NULL, output_path,
                                      export_tests = TRUE, export_name = NULL) {
  
  # Extract y variable data
  y_data <- data[[y_var]]
  data_max <- max(y_data, na.rm = TRUE)
  
  # Calculate y limits if not provided
  if (is.null(y_limits)) {
    y_min <- 0  # Starting from 0 for count data
    
    # Scale padding based on the magnitude of the data
    if (data_max < 70) {
      y_max <- data_max * 1.5
    } else if (data_max < 100) {
      y_max <- data_max * 1.4
    } else if (data_max < 500) {
      y_max <- data_max * 1.3
    } else if (data_max < 1000) {
      y_max <- data_max + 200
    } else if (data_max < 2000) {
      # Reduce padding for values > 1000
      y_max <- data_max + 300  # Was 500, now 300 (reduced by 200 units)
    } else {
      # Reduce padding for values > 2000
      y_max <- data_max + 400  # Was data_max * 1.2, now using fixed addition
    }
    
    # Round y_max to a nice number
    y_max <- ceiling(y_max / 10) * 10  # Round to nearest 10
    if (y_max > 100) y_max <- ceiling(y_max / 50) * 50  
    if (y_max > 500) y_max <- ceiling(y_max / 100) * 100  # Round to nearest 100 for larger values
    if (y_max > 1000) y_max <- ceiling(y_max / 500) * 500
    #if (y_max > 3000) y_max <- ceiling(y_max / 1000) * 1000  # Round to nearest 1000 for even larger values
    
    y_limits <- c(y_min, y_max)
  }
  
  # Calculate y breaks if not provided
  if (is.null(y_breaks)) {
    y_range <- y_limits[2] - y_limits[1]
    
    # Determine number of breaks based on range
    if (y_range <= 20) {
      break_step <- 5
    } else if (y_range <= 50) {
      break_step <- 5
    } else if (y_range <= 100) {
      break_step <- 10
    } else if (y_range <= 500) {
      break_step <- 50
    } else if (y_range <= 1000) {
      break_step <- 100
    } else if (y_range <= 5000) {
      break_step <- 500
    } else {
      break_step <- 1000
    }
    
    # Generate sequence of breaks
    y_breaks <- seq(y_limits[1], y_limits[2], by = break_step)
  }
  
  # Create formula for statistical tests
  formula_str <- paste(y_var, "~", grouping_var)
  formula_obj <- as.formula(formula_str)
  
  # Try to run Kruskal-Wallis test first to check if it's valid
  kw_test <- kruskal_test(data, formula_obj)
  
  # Auto-adjust plot dimensions based on number of groups and label lengths
  num_groups <- length(unique(data[[grouping_var]]))
  max_label_length <- max(nchar(as.character(data[[grouping_var]])))
  
  # Calculate dynamic plot dimensions
  plot_height <- max(9, 6 + (num_groups * 0.5))  # Base height + additional height per group
  plot_width <- max(12, 10 + (max_label_length * 0.2))  # Base width + additional width for long labels
  
  # Check if Kruskal-Wallis test has valid results
  if (is.na(kw_test$p[1]) || kw_test$p[1] > 0.05) {
    # If not valid or not significant, create basic boxplot without p-values
    p <- ggplot(data, aes_string(x = grouping_var, y = y_var)) + 
      geom_boxplot() + geom_jitter(width = 0, height = 0, alpha = 0.5) + 
      scale_x_discrete(limits = rev(levels(data[[grouping_var]]))) + 
      labs(title = title, x = grouping_var, y = y_label, 
           caption = "No significant differences found between groups (Kruskal-Wallis test)") +  
      theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
      theme(plot.caption = element_text(size = 20),
            axis.text.y = element_text(size = max(20, 36 - max_label_length * 0.5))) + # Auto-adjust axis text size
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
    
    # Save plot with dynamic dimensions
    ggsave(filename = output_path, plot = p, width = plot_width, height = plot_height)
    return(p)
  }
  
  # If Kruskal-Wallis test is valid, proceed with Dunn test
  tryCatch({
    # Dunn test for pairwise comparisons
    posthoc_tests <- data %>% 
      dunn_test(formula_obj, p.adjust.method = "fdr") %>%
      mutate(p_value_rounded = sprintf("p = %.3f", p.adj),
             p.adj.signif = as.character(p.adj.signif)) %>%
      add_xy_position(x = grouping_var)
    
    # Export the posthoc_tests dataframe to the global environment with user-defined name
    if (export_tests && nrow(posthoc_tests) > 0) {
      # Use custom export name if provided, or default to y_var + "_posthoc_tests"
      df_name <- ifelse(!is.null(export_name), 
                        export_name, 
                        paste0(y_var, "_posthoc_tests"))
      
      assign(df_name, posthoc_tests, envir = .GlobalEnv)
      cat(paste0("Exported test results to dataframe: ", df_name, "\n"))
    }
    
    # Check if posthoc_tests has valid results
    if (nrow(posthoc_tests) > 0) {
      
      # Filter to show only significant p-values
      posthoc_tests_filtered <- posthoc_tests %>% 
        #filter(p.adj.signif != "ns") %>%
        # Sort the comparisons to group similar ones
        arrange(group1, group2)
      
      # Create text for p-values
      p_values_text <- paste0(
        "Adjusted p-values for pairwise comparisons (Dunn's Test):\n",
        paste(posthoc_tests_filtered$group1, "vs", posthoc_tests_filtered$group2, 
              ": p=", round(posthoc_tests_filtered$p.adj, 3), collapse = "\n")
      )
      
      # Calculate y position and parameters based on data scale
      if (data_max < 70) {
        y_position <- data_max * 1.3  # Start bar at 15% above max
        step_increase <- 0.6           # Small step between bars
        bracket_nudge <- -0.05 * data_max  # Small nudge relative to data
      } else if (data_max < 100) {
        y_position <- data_max * 1.25  # Start bar at 13% above max
        step_increase <- 0.5
        bracket_nudge <- -0.04 * data_max
      } else if (data_max < 500) {
        y_position <- data_max * 1.2  # Start bar at 10% above max
        step_increase <- 0.4
        bracket_nudge <- -0.03 * data_max
      } else if (data_max < 1000) {
        y_position <- data_max * 1.15  # Start bar at 8% above max
        step_increase <- 0.3
        bracket_nudge <- -0.02 * data_max
      } else if (data_max < 2000) {
        y_position <- data_max * 1.1  # Start bar at 7% above max
        step_increase <- 0.2
        bracket_nudge <- -0.01 * data_max
      } else {
        y_position <- data_max * 1.08  # Start bar at 6% above max
        step_increase <- 0.1
        bracket_nudge <- -0.008 * data_max
      }
      
      # Create custom positions for multiple comparisons
      if (nrow(posthoc_tests_filtered) > 1) {
        # Calculate positions for multiple bars
        posthoc_tests_filtered <- posthoc_tests_filtered %>%
          mutate(
            y.position = y_position + ((row_number() - 1) * step_increase * data_max)
          )
      } else {
        # For a single comparison, just use the calculated y_position
        posthoc_tests_filtered$y.position <- y_position
      }
      
      # Create boxplot with p-values
      p <- ggplot(data, aes_string(x = grouping_var, y = y_var)) + 
        geom_boxplot() + 
        geom_jitter(width = 0, height = 0, alpha = 0.5) + 
        scale_x_discrete(limits = rev(levels(data[[grouping_var]]))) + 
        stat_pvalue_manual(
          data = posthoc_tests_filtered, 
          label = "p.adj.signif",
          coord.flip = TRUE, 
          angle = 0, 
          hjust = 0, 
          xmin = "group1", 
          xmax = "group2",
          y.position = "y.position",  # Use our manually calculated positions
          size = 8, 
          bracket.size = 0.6,
          tip.length = 0.01,
          bracket.nudge.y = bracket_nudge
        ) +
        labs(title = title, x = grouping_var, y = y_label, caption = p_values_text) +  
        theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
        theme(plot.caption = element_text(size = 20),
              axis.text.y = element_text(size = max(20, 36 - max_label_length * 0.5))) +
        scale_y_continuous(breaks = y_breaks, limits = y_limits) + 
        coord_flip()
    } else {
      # No valid posthoc tests, create basic boxplot
      p <- ggplot(data, aes_string(x = grouping_var, y = y_var)) + 
        geom_boxplot() + geom_jitter(width = 0, height = 0, alpha = 0.5) + 
        scale_x_discrete(limits = rev(levels(data[[grouping_var]]))) + 
        labs(title = title, x = grouping_var, y = y_label, 
             caption = "Kruskal-Wallis test significant but no valid pairwise comparisons") +  
        theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
        theme(plot.caption = element_text(size = 20),
              axis.text.y = element_text(size = max(20, 36 - max_label_length * 0.5))) + 
        scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
    }
  }, error = function(e) {
    # Handle any error in the Dunn test by creating a basic boxplot
    message("Error in Dunn test: ", e$message)
    p <- ggplot(data, aes_string(x = grouping_var, y = y_var)) + 
      geom_boxplot() + geom_jitter(width = 0, height = 0, alpha = 0.5) + 
      scale_x_discrete(limits = rev(levels(data[[grouping_var]]))) + 
      labs(title = title, x = grouping_var, y = y_label, 
           caption = "Could not perform statistical tests") +  
      theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
      theme(plot.caption = element_text(size = 20),
            axis.text.y = element_text(size = max(20, 36 - max_label_length * 0.5))) + 
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
  })
  
  # Save plot with dynamic dimensions
  ggsave(filename = output_path, plot = p, width = plot_width, height = plot_height)
  
  return(p)
}

# ______________________________________________________________________________

#' Create barplot with pairwise comparisons
#'
#' @param data DataFrame with data
#' @param y_var Name of the y-axis variable (string)
#' @param grouping_var Name of the grouping variable (string) - "Technique" or "Pool"
#' @param title Plot title
#' @param y_label Y-axis label
#' @param error_type Type of error bars: "se" (standard error), "sd" (standard deviation), or "ci" (95% confidence interval)
#' @param y_breaks Vector with y-axis breaks (optional)
#' @param y_limits Vector with y-axis limits (optional)
#' @param output_path Path to save the plot
#' @param export_tests Whether to export test results (default: TRUE)
#' @param export_name Custom name for the exported test results dataframe (optional)
#' @return ggplot object
create_barplot_with_stats <- function(data, y_var, grouping_var = GROUPING_VAR, title, y_label,
                                      y_breaks = NULL, y_limits = NULL, output_path,
                                      export_tests = TRUE, export_name = NULL) {
  
  # Calculate summary statistics (simplified)
  summary_stats <- data %>%
    group_by(!!sym(grouping_var)) %>%
    summarise(
      n = n(),
      mean = mean(!!sym(y_var), na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate y limits if not provided
  if (is.null(y_limits)) {
    y_min <- 0
    y_max <- max(summary_stats$mean, na.rm = TRUE)
    
    # Add some padding (20% of range)
    y_padding <- 0.2 * y_max
    y_max <- y_max + y_padding
    
    # Round y_max to a nice number
    y_max <- ceiling(y_max / 10) * 10
    if (y_max > 100) y_max <- ceiling(y_max / 100) * 100
    if (y_max > 1000) y_max <- ceiling(y_max / 1000) * 1000
    
    y_limits <- c(y_min, y_max)
  }
  
  # Calculate y breaks if not provided
  if (is.null(y_breaks)) {
    y_range <- y_limits[2] - y_limits[1]
    
    # Determine number of breaks based on range
    if (y_range <= 20) {
      break_step <- 5
    } else if (y_range <= 50) {
      break_step <- 5
    } else if (y_range <= 100) {
      break_step <- 10
    } else if (y_range <= 500) {
      break_step <- 50
    } else if (y_range <= 1000) {
      break_step <- 100
    } else if (y_range <= 5000) {
      break_step <- 500
    } else {
      break_step <- 1000
    }
    
    # Generate sequence of breaks
    y_breaks <- seq(y_limits[1], y_limits[2], by = break_step)
  }
  
  # Create formula for statistical tests
  formula_str <- paste(y_var, "~", grouping_var)
  formula_obj <- as.formula(formula_str)
  
  # Try to run Kruskal-Wallis test first to check if it's valid
  kw_test <- kruskal_test(data, formula_obj)
  
  # Auto-adjust plot dimensions based on number of groups and label lengths
  num_groups <- length(unique(data[[grouping_var]]))
  max_label_length <- max(nchar(as.character(data[[grouping_var]])))
  
  # Calculate dynamic plot dimensions
  plot_height <- max(9, 6 + (num_groups * 0.5))  # Base height + additional height per group
  plot_width <- max(12, 10 + (max_label_length * 0.2))  # Base width + additional width for long labels
  
  # Check if Kruskal-Wallis test has valid results
  if (is.na(kw_test$p[1]) || kw_test$p[1] > 0.05) {
    # If not valid or not significant, create basic barplot without p-values
    p <- ggplot(summary_stats, aes_string(x = grouping_var, y = "mean")) + 
      geom_bar(stat = "identity", fill = "steelblue", width = 0.7, color = "black") +
      scale_x_discrete(limits = rev(levels(data[[grouping_var]]))) + 
      labs(title = title, x = grouping_var, y = y_label, 
           caption = paste0("No significant differences found between groups (Kruskal-Wallis test)")) +
      theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
      theme(plot.caption = element_text(size = 20),
            axis.text.y = element_text(size = max(20, 36 - max_label_length * 0.5))) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
    
    # Save plot with dynamic dimensions
    ggsave(filename = output_path, plot = p, width = plot_width, height = plot_height)
    return(p)
  }
  
  # If Kruskal-Wallis test is valid, proceed with Dunn test
  tryCatch({
    # Dunn test for pairwise comparisons
    posthoc_tests <- data %>% 
      dunn_test(formula_obj, p.adjust.method = "fdr") %>%
      mutate(p_value_rounded = sprintf("p = %.3f", p.adj),
             p.adj.signif = as.character(p.adj.signif)) %>%
      add_xy_position(x = grouping_var)
    
    # Export the posthoc_tests dataframe to the global environment with user-defined name
    if (export_tests && nrow(posthoc_tests) > 0) {
      # Use custom export name if provided, or default to y_var + "_posthoc_tests"
      df_name <- ifelse(!is.null(export_name), 
                        export_name, 
                        paste0(y_var, "_barplot_tests"))
      
      assign(df_name, posthoc_tests, envir = .GlobalEnv)
      cat(paste0("Exported test results to dataframe: ", df_name, "\n"))
    }
    
    # Check if posthoc_tests has valid results
    if (nrow(posthoc_tests) > 0) {
      
      # Filter to show only significant p-values
      posthoc_tests_filtered <- posthoc_tests %>%
        filter(p.adj.signif != "ns") %>% 
        # Sort the comparisons to group similar ones
        arrange(group1, group2)
      
      # Calculate maximum allowed y position (95% of y-axis maximum)
      max_allowed_y <- y_limits[2] * 0.95
      
      # Add manual positioning for the significance bars
      if (nrow(posthoc_tests_filtered) > 0) {
        # Create positions for each comparison that stay within the plot boundaries
        posthoc_tests_filtered <- posthoc_tests_filtered %>%
          mutate(
            # Set each bar position to be within the y-limits
            y.position = pmin(
              # Start at 70% of data max and increase by 7% for each bar
              data_max * (0.7 + (row_number() - 1) * 0.07),
              # But never exceed 95% of the y-axis maximum
              max_allowed_y
            )
          )
      }
      
      # Create text for p-values
      p_values_text <- paste0(
        "Adjusted p-values for pairwise comparisons (Dunn's Test):\n",
        paste(posthoc_tests_filtered$group1, "vs", posthoc_tests_filtered$group2, 
              ": p=", round(posthoc_tests_filtered$p.adj, 3), collapse = "\n")
      )
      
      # Increase upper limit by 50% to make room for significance bars
      y_limits[2] <- y_limits[2] * 1.5
      
      # Create barplot with p-values
      p <- ggplot(summary_stats, aes_string(x = grouping_var, y = "mean")) + 
        geom_bar(stat = "identity", fill = "steelblue", width = 0.7, color = "black") +
        scale_x_discrete(limits = rev(levels(data[[grouping_var]]))) + 
        stat_pvalue_manual(data = posthoc_tests_filtered, label = "p.adj.signif",
                           coord.flip = TRUE, angle = 0, hjust = 0, xmin = "group1", xmax = "group2",
                           size = 8, bracket.size = 0.6, tip.length = 0.01) +
        labs(title = title, x = grouping_var, y = paste0(y_label), 
             caption = p_values_text) +  
        theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
        theme(plot.caption = element_text(size = 20),
              axis.text.y = element_text(size = max(20, 36 - max_label_length * 0.5))) +
        scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
    } else {
      # No valid posthoc tests, create basic barplot
      p <- ggplot(summary_stats, aes_string(x = grouping_var, y = "mean")) + 
        geom_bar(stat = "identity", fill = "steelblue", width = 0.7, color = "black") +
        scale_x_discrete(limits = rev(levels(data[[grouping_var]]))) + 
        labs(title = title, x = grouping_var, y = paste0(y_label), 
             caption = "Kruskal-Wallis test significant but no valid pairwise comparisons") +  
        theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
        theme(plot.caption = element_text(size = 20),
              axis.text.y = element_text(size = max(20, 36 - max_label_length * 0.5))) +
        scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
    }
  }, error = function(e) {
    # Handle any error in the Dunn test by creating a basic barplot
    message("Error in Dunn test: ", e$message)
    p <- ggplot(summary_stats, aes_string(x = grouping_var, y = "mean")) + 
      geom_bar(stat = "identity", fill = "steelblue", width = 0.7, color = "black") +
      scale_x_discrete(limits = rev(levels(data[[grouping_var]]))) + 
      labs(title = title, x = grouping_var, y = paste0(y_label), 
           caption = "Could not perform statistical tests") +  
      theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
      theme(plot.caption = element_text(size = 20),
            axis.text.y = element_text(size = max(20, 36 - max_label_length * 0.5))) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
  })
  
  # Save plot with dynamic dimensions
  ggsave(filename = output_path, plot = p, width = plot_width, height = plot_height)
  
  return(p)
}

# ______________________________________________________________________________

#' Load and prepare data
#'
#' @param file_path Path to CSV file
#' @param analysis_type Type of analysis ('Technique' or 'Pooled')
#' @return Processed DataFrame
#' Load and prepare data
#'
#' @param file_path Path to CSV file
#' @param analysis_type Type of analysis ('Technique' or 'Pooled')
#' @return Processed DataFrame
load_and_prepare_data <- function(file_path, analysis_type = ANALYSIS_TYPE) {
  df <- read.csv(file_path, sep = ";", fileEncoding = "UTF-8")
  
  # Convert categorical columns to factor
  as_factor <- c("Sample", "Technique", "Pool")
  df[, as_factor] <- lapply(df[, as_factor], factor)
  
  # Assign unique ID and filter data based on analysis type
  df <- df %>% mutate(ID = row_number()) 
  
  # Apply different filtering strategies based on analysis type
  if (analysis_type == "Technique") {
    # For technique analysis, focus on specific pools
    if (any(c("POOL_1", "POOL_2", "POOL_3") %in% df$Pool)) {
      df <- df %>% filter(Pool %in% c("POOL_1", "POOL_2", "POOL_3"))
      cat("\nFiltered data for technique analysis using pools 1-3\n")
    } else if ("ALL_POOLS" %in% df$Pool) {
      # Alternative if individual pools aren't available
      df <- df %>% filter(Pool == "ALL_POOLS")
      cat("\nFiltered data for technique analysis using ALL_POOLS\n")
    }
    # Recode technique if necessary
    if("IP_CD9" %in% levels(df$Technique)) {
      df <- df %>% mutate(Technique = recode(Technique, "IP_CD9" = "IP CD9"))
    } 
    if ("Raw_milk" %in% levels(df$Technique)){
      df <- df %>% mutate(Technique = recode(Technique, "Raw_milk" = "Raw milk"))
    }
    
    # Apply the reordering function to move the last technique (Raw_milk/Raw milk) to first position
    df <- reorder_technique_factor(df)
  }
  return(df)
}


# ______________________________________________________________________________
# 3. TOTAL DATA ANALYSIS ----
# ______________________________________________________________________________

cat("\n\n===== TOTAL DATA ANALYSIS =====\n\n")

# Load data
data_path <- "../output/1_filtered_dfs/vesiclepedia/summary_all.csv"
df_total <- load_and_prepare_data(data_path, ANALYSIS_TYPE)

# reorder techniques
levels(df_total$Technique)
df_total$Technique <- factor(df_total$Technique, 
                             levels = c(
                               "Raw milk", "ExoGAG", "IP CD9", "SEC", "UC"))

# Verify data structure
cat("Total data summary:\n")
print(summary(df_total))

# ______________________________________________________________________________

## 3.1 Protein Analysis ----
cat("\n--- Total Protein Analysis ---\n")

# Kruskal-Wallis test
kw_formula <- as.formula(paste("n.Proteins.Total ~", GROUPING_VAR))
test_total_prot <- perform_kruskal_test(df_total, kw_formula)
print(test_total_prot)

# Define common variables for plots
suffix <- ifelse(ANALYSIS_TYPE == "Technique", "per_technique", "per_pool")
protein_title <- paste("Proteins per", ANALYSIS_TYPE)

# Create plots based on user selection
if (PLOT_TYPE == "boxplots" || PLOT_TYPE == "both") {
  boxplot_total_pept <- create_boxplot_with_stats(
    data = df_total,
    y_var = "n.Proteins.Total",
    grouping_var = GROUPING_VAR,
    title = protein_title,
    y_label = "Proteins (n)",
    export_name = paste0("total_prot_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/boxplots/01a1_boxplot_total_prot_", suffix, ".png")
  )
}

if (PLOT_TYPE == "barplots" || PLOT_TYPE == "both") {
  barplot_total_pept <- create_barplot_with_stats(
    data = df_total,
    y_var = "n.Proteins.Total",
    grouping_var = GROUPING_VAR,
    title = protein_title,
    y_label = "Proteins (n)",
    export_name = paste0("total_prot_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/barplots/01a2_barplot_total_prot_", suffix, ".png")
  )
}

# ______________________________________________________________________________

## 3.2 Peptide Analysis ----
cat("\n--- Total Peptide Analysis ---\n")

# Kruskal-Wallis test formula
kw_formula <- as.formula(paste("n.Peptides.Total ~", GROUPING_VAR))
test_total_pept <- perform_kruskal_test(df_total, kw_formula)
print(test_total_pept)

# Define common variables for plots
suffix <- ifelse(ANALYSIS_TYPE == "Technique", "per_technique", "per_pool")
peptide_title <- paste("Peptides per", ANALYSIS_TYPE)

# Create plots based on user selection
if (PLOT_TYPE == "boxplots" || PLOT_TYPE == "both") {
  boxplot_total_pept <- create_boxplot_with_stats(
    data = df_total,
    y_var = "n.Peptides.Total",
    grouping_var = GROUPING_VAR,
    title = peptide_title,
    y_label = "Peptides (n)",
    export_name = paste0("total_pept_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/boxplots/01b1_boxplot_total_pept_", suffix, ".png")
  )
}

if (PLOT_TYPE == "barplots" || PLOT_TYPE == "both") {
  barplot_total_pept <- create_barplot_with_stats(
    data = df_total,
    y_var = "n.Peptides.Total",
    grouping_var = GROUPING_VAR,
    title = peptide_title,
    y_label = "Peptides (n)",
    export_name = paste0("total_pept_barplot_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/barplots/01b2_barplot_total_pept_", suffix, ".png")
  )
}


# ______________________________________________________________________________
# 4. VESICLEPEDIA ANALYSIS ----
# ______________________________________________________________________________

cat("\n\n===== VESICLEPEDIA ANALYSIS =====\n\n")

# We already have the data loaded in df_total, no need to load it again

# ______________________________________________________________________________

## 4.1 Filtered Protein Analysis ----
cat("\n--- Vcp Protein Analysis ---\n")

# Kruskal-Wallis test
kw_formula <- as.formula(paste("n.Proteins.Filtered ~", GROUPING_VAR))
test_vcp_prot <- perform_kruskal_test(df_total, kw_formula)
print(test_vcp_prot)

# Define common variables for plots
suffix <- ifelse(ANALYSIS_TYPE == "Technique", "per_technique", "per_pool")
vcp_prot_title <- paste("Vesiclepedia Proteins")

# Create plots based on user selection
if (PLOT_TYPE == "boxplots" || PLOT_TYPE == "both") {
  boxplot_vcp_prot <- create_boxplot_with_stats(
    data = df_total,
    y_var = "n.Proteins.Filtered",
    grouping_var = GROUPING_VAR,
    title = vcp_prot_title,
    y_label = "Proteins (n)",
    export_name = paste0("vcp_prot_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/boxplots/06a1_boxplot_vcp_prot_", suffix, ".png")
  )
}

if (PLOT_TYPE == "barplots" || PLOT_TYPE == "both") {
  barplot_vcp_prot <- create_barplot_with_stats(
    data = df_total,
    y_var = "n.Proteins.Filtered",
    grouping_var = GROUPING_VAR,
    title = vcp_prot_title,
    y_label = "Proteins (n)",
    export_name = paste0("vcp_prot_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/barplots/06a2_barplot_vcp_prot_", suffix, ".png")
  )
}

# ______________________________________________________________________________

## 4.2 Filtered Peptide Analysis ----
cat("\n--- Vcp Peptide Analysis ---\n")

# Kruskal-Wallis test formula
kw_formula <- as.formula(paste("n.Peptides.Filtered ~", GROUPING_VAR))
test_vcp_pept <- perform_kruskal_test(df_total, kw_formula)
print(test_vcp_pept)

# Define common variables for plots
suffix <- ifelse(ANALYSIS_TYPE == "Technique", "per_technique", "per_pool")
vcp_pept_title <- paste("Vesiclepedia Peptides")


# Create plots based on user selection
if (PLOT_TYPE == "boxplots" || PLOT_TYPE == "both") {
  boxplot_vcp_pept <- create_boxplot_with_stats(
    data = df_total,
    y_var = "n.Peptides.Filtered",
    grouping_var = GROUPING_VAR,
    title = vcp_pept_title,
    y_label = "Peptides (n)",
    export_name = paste0("vcp_pept_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/boxplots/06b1_boxplot_vcp_pept_", suffix, ".png")
  )
}

if (PLOT_TYPE == "barplots" || PLOT_TYPE == "both") {
  barplot_vcp_pept <- create_barplot_with_stats(
    data = df_total,
    y_var = "n.Peptides.Filtered",
    grouping_var = GROUPING_VAR,
    title = vcp_pept_title,
    y_label = "Peptides (n)",
    export_name = paste0("vcp_pept_barplot_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/barplots/06b2_barplot_vcp_pept_", suffix, ".png")
  )
}


# ______________________________________________________________________________
# 5. GLYCOSYLATION ANALYSIS ----
# ______________________________________________________________________________

cat("\n\n===== GLYCOSYLATION ANALYSIS =====\n\n")

# Load glycosylation data
glyc_path <- "../output/1_filtered_dfs/vesiclepedia_glycosylated/summary_all.csv"
df_glyc <- load_and_prepare_data(glyc_path, ANALYSIS_TYPE)

# reorder techniques
levels(df_glyc$Technique)
df_glyc$Technique <- factor(df_glyc$Technique, 
                             levels = c(
                               "Raw milk", "ExoGAG", "IP CD9", "SEC", "UC"))

# Verify data structure
cat("Glycosylated data summary:\n")
print(summary(df_glyc))



# ______________________________________________________________________________

## 5.1 Glycosylated Protein Analysis ----
cat("\n--- Glycosylated Protein Analysis ---\n")

# Kruskal-Wallis test
kw_formula <- as.formula(paste("n.Proteins.Filtered ~", GROUPING_VAR))
test_glyc_prot <- perform_kruskal_test(df_glyc, kw_formula)
print(test_glyc_prot)

# Define common variables for plots
suffix <- ifelse(ANALYSIS_TYPE == "Technique", "per_technique", "per_pool")
glyc_prot_title <- paste("Vesiclepedia Glycoproteins")

# Create plots based on user selection
if (PLOT_TYPE == "boxplots" || PLOT_TYPE == "both") {
  boxplot_vcp_prot_glyc <- create_boxplot_with_stats(
    data = df_glyc,
    y_var = "n.Proteins.Filtered",
    grouping_var = GROUPING_VAR,
    title = glyc_prot_title,
    y_label = "Proteins (n)",
    export_name = paste0("vcp_prot_glyc_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/boxplots/06c1_boxplot_vcp_prot_glyc_", suffix, ".png")
  )
}

if (PLOT_TYPE == "barplots" || PLOT_TYPE == "both") {
  barplot_vcp_prot_glyc <- create_barplot_with_stats(
    data = df_glyc,
    y_var = "n.Proteins.Filtered",
    grouping_var = GROUPING_VAR,
    title = glyc_prot_title,
    y_label = "Proteins (n)",
    export_name = paste0("vcp_prot_glyc_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/barplots/06c2_barplot_vcp_prot_glyc_", suffix, ".png")
  )
}

# ______________________________________________________________________________

## 5.2 Glycosylated Peptide Analysis ----
cat("\n--- Glycosylated Peptide Analysis ---\n")

# Kruskal-Wallis test formula
kw_formula <- as.formula(paste("n.Peptides.Filtered ~", GROUPING_VAR))
test_glyc_pept <- perform_kruskal_test(df_glyc, kw_formula)
print(test_glyc_pept)

# Define common variables for plots
suffix <- ifelse(ANALYSIS_TYPE == "Technique", "per_technique", "per_pool")
glyc_pept_title <- paste("Vesiclepedia Glycopeptides")

# Create plots based on user selection
if (PLOT_TYPE == "boxplots" || PLOT_TYPE == "both") {
  boxplot_vcp_pept_glyc <- create_boxplot_with_stats(
    data = df_glyc,
    y_var = "n.Peptides.Filtered",
    grouping_var = GROUPING_VAR,
    title = glyc_pept_title,
    y_label = "Peptides (n)",
    export_name = paste0("vcp_pept_glyc_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/boxplots/06d1_boxplot_vcp_pept_glyc_", suffix, ".png")
  )
}

if (PLOT_TYPE == "barplots" || PLOT_TYPE == "both") {
  barplot_vcp_pept_glyc <- create_barplot_with_stats(
    data = df_glyc,
    y_var = "n.Peptides.Filtered",
    grouping_var = GROUPING_VAR,
    title = glyc_pept_title,
    y_label = "Peptides (n)",
    export_name = paste0("vcp_pept_glyc_barplot_pwcomps_", suffix),
    output_path = paste0("../output/3_figures/barplots/06d2_barplot_vcp_pept_glyc_", suffix, ".png")
  )
}

# ______________________________________________________________________________
if (PLOT_TYPE == "both") {
  cat(paste0("\n\nAnalysis completed. Figures have been saved to:\n",
             "- Boxplots: 'output/3_figures/boxplots/'\n", 
             "- Barplots: 'output/3_figures/barplots/'\n", 
             "All with suffix '_", suffix, "'\n"))
} else if (PLOT_TYPE == "barplots") {
  cat(paste0("\n\nAnalysis completed. All figures have been saved to 'output/3_figures/barplots/' with suffix '_", suffix, "'\n"))
} else {
  cat(paste0("\n\nAnalysis completed. All figures have been saved to 'output/3_figures/boxplots/' with suffix '_", suffix, "'\n"))
}
# ______________________________________________________________________________