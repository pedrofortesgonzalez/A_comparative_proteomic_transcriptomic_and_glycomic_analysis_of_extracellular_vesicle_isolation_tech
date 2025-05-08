# ******************************************************************************
# Statistical analysis and visualization for mass spectrometry data
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
#' @param y_breaks Vector with y-axis breaks
#' @param y_limits Vector with y-axis limits
#' @param output_path Path to save the plot
#' @return ggplot object
create_boxplot_with_stats <- function(data, y_var, title, y_label, 
                                      y_breaks, y_limits, output_path) {
  # Create formula for dunn test
  formula_str <- paste(y_var, "~ Technique")
  formula_obj <- as.formula(formula_str)
  
  # Dunn test for pairwise comparisons
  posthoc_tests <- data %>% 
    dunn_test(formula_obj, p.adjust.method = "fdr") %>%
    mutate(p_value_rounded = sprintf("p = %.3f", p.adj),
           p.adj.signif = as.character(p.adj.signif)) %>%
    add_xy_position(x = "Technique")
  
  # Create text for p-values
  p_values_text <- paste0(
    "Adjusted p-values for pairwise comparisons (Dunn's Test):\n",
    paste(posthoc_tests$group1[1:3], "vs", posthoc_tests$group2[1:3], 
          ": p=", round(posthoc_tests$p.adj[1:3], 3), collapse = "; "), "\n",
    paste(posthoc_tests$group1[4:length(posthoc_tests$group1)], "vs", 
          posthoc_tests$group2[4:length(posthoc_tests$group2)], 
          ": p=", round(posthoc_tests$p.adj[4:length(posthoc_tests$p.adj)], 3), 
          collapse = "; ")
  )
  
  # Filter to show only significant p-values
  posthoc_tests_filtered <- posthoc_tests %>% filter(p.adj.signif != "ns")
  
  # Create boxplot with aes_string instead of quasiquotation
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
  df <- read.csv(file_path, sep = ",", fileEncoding = "UTF-8")
  
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
  y_breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000),
  y_limits = c(0, 3000),
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
  y_breaks = c(0, 100, 200, 300, 400, 500, 600),
  y_limits = c(0, 600),
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
  y_breaks = c(0, 50, 100, 150, 200, 250, 300),
  y_limits = c(0, 320),
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
  y_breaks = c(0, 500, 1000, 1500, 2000),
  y_limits = c(0, 2100),
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
  y_breaks = c(0, 5, 10, 15),
  y_limits = c(0, 17),
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
  y_breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45),
  y_limits = c(0, 45),
  output_path = "../output/3_figures/boxplots/06b2_boxplot_vcp_pept_glyc_pwcomp.png"
)

# ______________________________________________________________________________
cat("\n\nAnalysis completed. All figures have been saved to 'output/3_figures/boxplots/'\n")
# ______________________________________________________________________________
