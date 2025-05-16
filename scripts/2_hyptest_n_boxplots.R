# ******************************************************************************
# Statistical analysis and visualization for mass spectrometry data
# Author: Pedro Fortes González (refactored)
# Date: 2024-09-16
# ******************************************************************************

# ______________________________________________________________________________
# 0. Environment setup ----
# ______________________________________________________________________________

# Clean WEnv
rm(list = ls())

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
    install_choice <- readline(prompt = "Do you want to install these packages? (y/n): ")
    
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
# 1. Helper functions ----
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
#' @param title Plot title
#' @param y_label Y-axis label
#' @param y_breaks Vector with y-axis breaks (optional)
#' @param y_limits Vector with y-axis limits (optional)
#' @param output_path Path to save the plot
#' @param export_tests Whether to export test results (default: TRUE)
#' @param export_name Custom name for the exported test results dataframe (optional)
#' @return ggplot object
create_boxplot_with_stats <- function(data, y_var, title, y_label, 
                                      y_breaks = NULL, y_limits = NULL, output_path,
                                      export_tests = TRUE, export_name = NULL) {
  
  # Extract y variable data
  y_data <- data[[y_var]]
  
  # Calculate y limits if not provided
  if (is.null(y_limits)) {
    y_min <- 0  # Starting from 0 for count data
    y_max <- max(y_data, na.rm = TRUE)
    
    # Add some padding (10% of range)
    y_padding <- 0.1 * y_max
    y_max <- y_max + y_padding
    
    # Round y_max to a nice number
    y_max <- ceiling(y_max / 10) * 10  # Round to nearest 10
    if (y_max > 100) y_max <- ceiling(y_max / 100) * 100  # Round to nearest 100 for larger values
    if (y_max > 1000) y_max <- ceiling(y_max / 1000) * 1000  # Round to nearest 1000 for even larger values
    
    y_limits <- c(y_min, y_max)
  }
  
  # Calculate y breaks if not provided
  if (is.null(y_breaks)) {
    y_range <- y_limits[2] - y_limits[1]
    
    # Determine number of breaks based on range
    if (y_range <= 20) {
      # For small ranges, use more granular breaks
      break_step <- 2
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
  
  # Create formula for dunn test
  formula_str <- paste(y_var, "~ Technique")
  formula_obj <- as.formula(formula_str)
  
  # Try to run Kruskal-Wallis test first to check if it's valid
  kw_test <- kruskal_test(data, formula_obj)
  
  # Check if Kruskal-Wallis test has valid results
  if (is.na(kw_test$p[1]) || kw_test$p[1] > 0.05) {
    # If not valid or not significant, create basic boxplot without p-values
    p <- ggplot(data, aes_string(x = "Technique", y = y_var)) + 
      geom_boxplot() + geom_jitter(width = 0, height = 0, alpha = 0.5) + 
      scale_x_discrete(limits = rev(levels(data$Technique))) + 
      labs(title = title, x = "Techniques", y = y_label, 
           caption = "No significant differences found between groups (Kruskal-Wallis test)") +  
      theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
      theme(plot.caption = element_text(size = 20)) + 
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
    
    # Save plot
    ggsave(filename = output_path, plot = p, width = 12, height = 9)
    return(p)
  }
  
  # If Kruskal-Wallis test is valid, proceed with Dunn test
  tryCatch({
    # Dunn test for pairwise comparisons
    posthoc_tests <- data %>% 
      dunn_test(formula_obj, p.adjust.method = "fdr") %>%
      mutate(p_value_rounded = sprintf("p = %.3f", p.adj),
             p.adj.signif = as.character(p.adj.signif)) %>%
      add_xy_position(x = "Technique")
    
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
      # Create text for p-values
      p_values_text <- paste0(
        "Adjusted p-values for pairwise comparisons (Dunn's Test):\n",
        paste(posthoc_tests$group1, "vs", posthoc_tests$group2, 
              ": p=", round(posthoc_tests$p.adj, 3), collapse = "; ")
      )
      
      # Filter to show only significant p-values
      posthoc_tests_filtered <- posthoc_tests %>% filter(p.adj.signif != "ns")
      
      # Create boxplot with p-values
      p <- ggplot(data, aes_string(x = "Technique", y = y_var)) + 
        geom_boxplot() + geom_jitter(width = 0, height = 0, alpha = 0.5) + 
        scale_x_discrete(limits = rev(levels(data$Technique))) + 
        stat_pvalue_manual(data = posthoc_tests_filtered, label = "p.adj.signif",
                           coord.flip = TRUE, angle = 0, hjust = 0, xmin = "group1", xmax = "group2", 
                           size = 8) + 
        labs(title = title, x = "Techniques", y = y_label, caption = p_values_text) +  
        theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
        theme(plot.caption = element_text(size = 20)) + 
        scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
    } else {
      # No valid posthoc tests, create basic boxplot
      p <- ggplot(data, aes_string(x = "Technique", y = y_var)) + 
        geom_boxplot() + geom_jitter(width = 0, height = 0, alpha = 0.5) + 
        scale_x_discrete(limits = rev(levels(data$Technique))) + 
        labs(title = title, x = "Techniques", y = y_label, 
             caption = "Kruskal-Wallis test significant but no valid pairwise comparisons") +  
        theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
        theme(plot.caption = element_text(size = 20)) + 
        scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
    }
  }, error = function(e) {
    # Handle any error in the Dunn test by creating a basic boxplot
    message("Error in Dunn test: ", e$message)
    p <- ggplot(data, aes_string(x = "Technique", y = y_var)) + 
      geom_boxplot() + geom_jitter(width = 0, height = 0, alpha = 0.5) + 
      scale_x_discrete(limits = rev(levels(data$Technique))) + 
      labs(title = title, x = "Techniques", y = y_label, 
           caption = "Could not perform statistical tests") +  
      theme_publish(base_size = 36, base_family = "Times", base_linewidth = 0.7) + 
      theme(plot.caption = element_text(size = 20)) + 
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + coord_flip()
  })
  
  # Save plot
  ggsave(filename = output_path, plot = p, width = 12, height = 9)
  
  return(p)
}

# ______________________________________________________________________________

#' Load and prepare data
#'
#' @param file_path Path to CSV file
#' @return Processed DataFrame
load_and_prepare_data <- function(file_path) {
  df <- read.csv(file_path, sep = ";", fileEncoding = "UTF-8")
  
  # Convert categorical columns to factor
  as_factor <- c("Sample", "Technique", "Pool")
  df[, as_factor] <- lapply(df[, as_factor], factor)
  
  # Assign unique ID and filter data
  df <- df %>%  
    mutate(ID = row_number()) %>%  
    filter(Pool %in% c("POOL_1", "POOL_2", "POOL_3"))
  
  # Recode technique if necessary
  if("IP_CD9" %in% levels(df$Technique)) {
    df <- df %>% mutate(Technique = recode(Technique, "IP_CD9" = "IP CD9"))
  }
  
  return(df)
}



# ______________________________________________________________________________
# 2. TOTAL DATA ANALYSIS ----
# ______________________________________________________________________________
cat("\n\n===== TOTAL DATA ANALYSIS =====\n\n")

# Load data
data_path <- "../output/1_filtered_dfs/vesiclepedia/summary_all.csv"
df_total <- load_and_prepare_data(data_path)

# Verify data structure
cat("Total data summary:\n")
print(summary(df_total))



# ______________________________________________________________________________

## 2.1 Peptide Analysis ----

cat("\n--- Total Peptide Analysis ---\n")

# Kruskal-Wallis test
test_total_pept <- perform_kruskal_test(df_total, n.Peptides.Total ~ Technique)
print(test_total_pept)

# Create boxplot
boxplot_total_pept <- create_boxplot_with_stats(
  data = df_total,
  y_var = "n.Peptides.Total",
  title = "Peptides per Technique",
  y_label = "Peptides (n)",
  export_name = "total_pept_pwcomps",
  output_path = "../output/3_figures/boxplots/01a1_boxplot_total_pept.png"
)



# ______________________________________________________________________________

## 2.2 Protein Analysis ----
cat("\n--- Total Protein Analysis ---\n")

# Kruskal-Wallis test
test_total_prot <- perform_kruskal_test(df_total, n.Proteins.Total ~ Technique)
print(test_total_prot)

# Create boxplot
boxplot_total_prot <- create_boxplot_with_stats(
  data = df_total,
  y_var = "n.Proteins.Total",
  title = "Proteins per Technique",
  y_label = "Proteins (n)",
  export_name = "total_prot_pwcomps",
  output_path = "../output/3_figures/boxplots/01b1_boxplot_total_prot.png"
)




# ______________________________________________________________________________
# 3. VESICLEPEDIA ANALYSIS ----
# ______________________________________________________________________________

cat("\n\n===== VESICLEPEDIA ANALYSIS =====\n\n")

# We already have the data loaded in df_total, no need to load it again

# ______________________________________________________________________________

## 3.1 Filtered Protein Analysis ----
cat("\n--- Filtered Protein Analysis (Vesiclepedia) ---\n")

# Kruskal-Wallis test
test_vcp_prot <- perform_kruskal_test(df_total, n.Proteins.Filtered ~ Technique)
print(test_vcp_prot)

# Create boxplot
boxplot_vcp_prot <- create_boxplot_with_stats(
  data = df_total,
  y_var = "n.Proteins.Filtered",
  title = "Vesiclepedia Proteins",
  y_label = "Proteins (n)",
  export_name = "vcp_prot_pwcomps",
  output_path = "../output/3_figures/boxplots/06a1_boxplot_vcp_prot.png"
)



# ______________________________________________________________________________

## 3.2 Filtered Peptide Analysis ----
cat("\n--- Filtered Peptide Analysis (Vesiclepedia) ---\n")

# Kruskal-Wallis test
test_vcp_pept <- perform_kruskal_test(df_total, n.Peptides.Filtered ~ Technique)
print(test_vcp_pept)

# Create boxplot
boxplot_vcp_pept <- create_boxplot_with_stats(
  data = df_total,
  y_var = "n.Peptides.Filtered",
  title = "Vesiclepedia Peptides",
  y_label = "Peptides (n)",
  export_name = "vcp_pept_pwcomps",
  output_path = "../output/3_figures/boxplots/06b1_boxplot_vcp_pept.png"
)


# ______________________________________________________________________________
# 4. GLYCOSYLATION ANALYSIS ----
# ______________________________________________________________________________

cat("\n\n===== GLYCOSYLATION ANALYSIS =====\n\n")

# Load glycosylation data
glyc_path <- "../output/1_filtered_dfs/vesiclepedia_glycosylated/summary_all.csv"
df_glyc <- load_and_prepare_data(glyc_path)

# Verify data structure
cat("Glycosylated data summary:\n")
print(summary(df_glyc))



# ______________________________________________________________________________

## 4.1 Glycosylated Protein Analysis ----
cat("\n--- Glycosylated Protein Analysis ---\n")

# Kruskal-Wallis test
test_glyc_prot <- perform_kruskal_test(df_glyc, n.Proteins.Filtered ~ Technique)
print(test_glyc_prot)

# Create boxplot
boxplot_glyc_prot <- create_boxplot_with_stats(
  data = df_glyc,
  y_var = "n.Proteins.Filtered",
  title = "Vesiclepedia Glycoproteins",
  y_label = "Glycoproteins (n)",
  export_name = "glyc_prot_pwcomps",
  output_path = "../output/3_figures/boxplots/06a2_boxplot_vcp_prot_glyc_pwcomp.png"
)


# ______________________________________________________________________________

## 4.2 Glycosylated Peptide Analysis ----
cat("\n--- Glycosylated Peptide Analysis ---\n")

# Kruskal-Wallis test
test_glyc_pept <- perform_kruskal_test(df_glyc, n.Peptides.Filtered ~ Technique)
print(test_glyc_pept)

# Create boxplot
boxplot_glyc_pept <- create_boxplot_with_stats(
  data = df_glyc,
  y_var = "n.Peptides.Filtered",
  title = "Vesiclepedia Glycopeptides",
  y_label = "Glycopeptides (n)",
  export_name = "glyc_pept_pwcomps",
  output_path = "../output/3_figures/boxplots/06b2_boxplot_vcp_pept_glyc_pwcomp.png"
)

# ______________________________________________________________________________
cat("\n\nAnalysis completed. All figures have been saved to 'output/3_figures/boxplots/'\n")
# ______________________________________________________________________________
