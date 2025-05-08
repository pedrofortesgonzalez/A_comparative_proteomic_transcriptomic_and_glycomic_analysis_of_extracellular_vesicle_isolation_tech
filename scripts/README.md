# Scripts
This folder contains all scripts needed to perform the analysis of glycosylations in EV proteins and peptides. The scripts are organized in numerical order according to the analysis workflow.

## Folder structure

scripts/

├── __pycache__/                       # Python compiled files (automatically generated)
├── 1_ptm_detection.py                 # Initial processing of mass spectrometry data
├── 2_hyptest_&_boxplots.R             # Statistical analysis and boxplot generation
├── 3_sector_diagrams.R                # Creation of pie charts
├── user_defined_funcs.py              # Helper functions for data processing with Python
├── cleaning_1_python.py               # Cleanup of unused Python libraries
└── cleaning_2_R.R                     # Cleanup of unused R libraries

## Script Descriptions
### 1_ptm_detection.py
This script performs the initial processing of mass spectrometry data. Its main functions are:

Creating the folder structure in output/
Processing input CSV files
Filtering proteins and peptides present in Vesiclepedia
Identifying proteins and peptides with glycosylations of interest
Generating count summaries for subsequent analyses

#### Main dependencies: pandas, numpy, os
#### Input files:
CSV files in input_data/input/
input_data/glycosylation_list.csv
input_data/vesiclepedia_proteins_240712.csv

#### Output files: Multiple CSV files in output/ subfolders

### user_defined_funcs.py
Module containing custom functions used by 1_ptm_detection.py. It separates the main logic from the script for better code organization and maintenance.
cleaning_1_python.py
Utility script that uses the vulture module to analyze 1_ptm_detection.py and identify imported libraries that are not used, helping to keep the code clean and efficient.
Dependencies: vulture
2_hyptest_&_boxplots.R
This script performs statistical analyses on the processed data and generates boxplot visualizations. Its main functions are:

Performing Kruskal-Wallis tests to compare isolation techniques
Running Dunn's post-hoc tests for pairwise comparisons
Generating boxplots that show:

Total peptide and protein counts
Counts of peptides and proteins present in Vesiclepedia
Counts of glycosylated peptides and proteins present in Vesiclepedia



Main dependencies: dplyr, ggplot2, ggpubr, rstatix, envalysis, gridExtra
Input files: Processed CSV files in output/processed_data/
Output files: Boxplots in output/figures/boxplots/
3_sector_diagrams.R
This script creates pie charts to visualize the distribution of glycosylation types in proteins and peptides. The main functions are:

Processing count data for different glycosylation types
Generating pie charts by isolation technique
Creating versions with and without text labels for each chart

Main dependencies: ggplot2, ggrepel
Input files: Count CSV files in output/processed_data/value_counts/
Output files: Pie charts in output/figures/sector_diagrams/
cleaning_2_R.R
Utility script that uses the nolock package to identify imported but unused libraries in R scripts, helping to optimize the code.
Dependencies: nolock, rstudioapi
Workflow
The complete analysis follows this workflow:

Data Processing: Run 1_ptm_detection.py to process raw data
Statistical Analysis: Run 2_hyptest_&_boxplots.R to generate statistical analyses and boxplots
Distribution Visualization: Run 3_sector_diagrams.R to create pie charts

## Requirements and Dependencies
Python

pandas
numpy
os
vulture (for code cleanup)

R

dplyr
ggplot2
ggpubr
rstatix
envalysis
gridExtra
ggrepel
nolock (for code cleanup)
rstudioapi

## Usage Notes

Execute the scripts in numerical order to ensure correct workflow.
The cleaning scripts (cleaning_1_python.py and cleaning_2_R.R) are optional and can be used to optimize code.
Ensure all dependencies are installed before running the scripts.
