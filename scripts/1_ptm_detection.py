#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PTM Detection: Analysis of glycosylated proteins from mass spectrometry data.

This script processes mass spectrometry data to identify post-translational
modifications, with a focus on protein glycosylation.

Author: Pedro Fortes Gonzalez (refactored)
"""

###############################################################################
# 1. Environment Setup
######################
# Import system libraries and set working directory
import os
import sys
import pandas as pd
import importlib
from tqdm.notebook import tqdm
import session_info


# %%

# Detecting repo's root dir and set working dir
try:
    # Try to obtain the actual script's directory (works for executing with bash)
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.dirname(SCRIPT_DIR)
    print(f"\nRepository authomatically detected in:\n{REPO_ROOT}\n")
    
except Exception as e:
    # If we work from an IDE/notebook where __file__ is not defined
    print(f"\nPath could not be authomatically detected:\n{str(e)}")
    REPO_ROOT = input("\nPlease, re-introduce the absolute path to the repository in your computer: ")
    REPO_ROOT = REPO_ROOT.strip()
    SCRIPT_DIR = REPO_ROOT + "/scripts"
    print(f"\nUsing the given path:\n{REPO_ROOT}\n")
    
    # Validate that the path exists
    if not os.path.exists(REPO_ROOT):
        print(f"\nERROR: The given path was not found:\n{REPO_ROOT}\n")
        sys.exit(1)

# Verify if we're in a bash environment
# If executed from bash, determine the relative path to current directory
if os.getcwd() != REPO_ROOT and os.path.basename(os.getcwd()) == "scripts":
    REPO_ROOT = os.path.dirname(os.getcwd())
    SCRIPT_DIR = REPO_ROOT + "/scripts"
    print(f"\nExecution from bash detected. Repo's path:\n{REPO_ROOT}\n")

# Print depuration info
print(f"\nCurrent script's directory:\n{SCRIPT_DIR}\n")
print(f"\nRepository root dir:\n{REPO_ROOT}\n")

# %%

# Show package list
session_info.show(dependencies=False, html=False, std_lib=False)

# Configure progress bar for pandas operations
tqdm.pandas(desc="Processing")

# Import custom functions module
sys.path.append(SCRIPT_DIR)
import user_defined_funcs as udf

# Expand dataframe display for better visualization
udf.expand_dataframe_view()

# %%

###############################################################################
# 2. Data Import
################

# Define file's paths relative to REPO_ROOT
GLYCOSYLATION_LIST_PATH = os.path.join(REPO_ROOT, 'input_data', 'glycosylation_list.csv')
VESICLEPEDIA_PATH = os.path.join(REPO_ROOT, 'input_data', 'Vesiclepedia_proteins_240712.csv')

# Verify that the files exist
for file_path in [GLYCOSYLATION_LIST_PATH, VESICLEPEDIA_PATH]:
    if not os.path.exists(file_path):
        print(f"\nERROR: Required file not found:\n{file_path}\n")
        print(f"\nSearch path:\n{REPO_ROOT}\n")
        print("Please verify that the files are in the correct path.\n")
        sys.exit(1)

# Import list of PTMs (glycosylations) of interest
ptm_df = pd.read_csv(GLYCOSYLATION_LIST_PATH, usecols=[0], sep=";")  # Select only first column
ptms_of_interest = list(set(ptm_df.iloc[:, 0].dropna().unique()))  # Remove NaNs and duplicates

# Verify PTM list import
print(f"\nNumber of PTMs imported: {len(ptms_of_interest)}")

# Import Vesiclepedia database (downloaded on 2024-07-12)
vesiclepedia = pd.read_csv(VESICLEPEDIA_PATH, sep=";")

# Verify Vesiclepedia database import
print(f"\nNumber of entries in Vesiclepedia: {vesiclepedia.shape[0]}")

# %%

###############################################################################
# 3. File Preprocessing
#######################
# Reload the updated functions module to ensure we have the latest version
importlib.reload(udf)

# Define input directory and technique identifiers
RAW_DATA_DIR = os.path.join(REPO_ROOT, str(input("\nEnter the name of your input folder: ")))
INPUT_DIR = os.path.join(REPO_ROOT, 'input_data/input')
TECHNIQUES = ["ExoGAG", "SEC", "IP_CD9", "UC"]

# Step 3.1: Create copy of the input directory for processing while preserving the original data
udf.create_input_directory(RAW_DATA_DIR, INPUT_DIR)
    
# Setp 3.2.: Ensure all files are properly named
udf.rename_pool_files(INPUT_DIR, TECHNIQUES)  # Erase files without pool identifier and correct the incorrect ones
udf.simplify_filenames(INPUT_DIR, TECHNIQUES) # Clean up filenames by removing prefixes and standardizing format
    
# Step 3.3: Detect CSV delimiters and standardize all CSV files to use semicolon (;) as delimiter
udf.standardize_csv_delimiters(INPUT_DIR)

# Step 3.4: Create combined pool files from individual pools
udf.combine_pool_files(INPUT_DIR, TECHNIQUES)  # Combine pool1 + pool2 + pool3 into POOLS_123 files

# %%

###############################################################################
# 4. Column Cleanup and Extraction
##################################

# Step 4.1: Extract protein accession numbers from the Accession column
importlib.reload(udf)

# Define inputs for protein name extraction
ACCESSION_COL = 'Accession'
PROTEIN_NAME_COL = "Prot_Name"

# Regular expression pattern to extract standardized protein accession numbers
# Matches patterns like UniProt accessions (e.g., P12345, Q9XYZ0) and other common formats
PROTEIN_PATTERN = r"([A-Z]\d[A-Z0-9]{3}[0-9]-?\d*|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z]?[A-Z0-9]+[0-9])"

# Extract protein names from accession information and create a new column
udf.extract_protein_names(INPUT_DIR, ACCESSION_COL, PROTEIN_NAME_COL, PROTEIN_PATTERN)

# Step 4.2: Clean peptide sequences by removing modifications in parentheses
importlib.reload(udf)
PEPTIDE_COL = 'Peptide'
PEPTIDE_SEQ_COL = "Peptide_Sequence"
PEPTIDE_PATTERN = r'\([^)]*\)'  # Pattern to remove content inside parentheses

# Extract clean peptide sequences without modification annotations
udf.extract_peptide_sequences(INPUT_DIR, PEPTIDE_COL, PEPTIDE_SEQ_COL, PEPTIDE_PATTERN)

# Step 4.3: Classify PTMs into standardized categories
importlib.reload(udf)
PTM_COL = 'PTM'
PTM_CLUSTER_COL = "PTM_Cluster"

# Create a new column with standardized PTM classifications
udf.classify_ptm_types(INPUT_DIR, PTM_COL, PTM_CLUSTER_COL)

# %%

###############################################################################
# 5. Create Output Directory Structure
###############################################################################

# Reload the functions module to ensure we're using the latest version
importlib.reload(udf)

# Define input directory
OUTPUT_DIR = os.path.join(REPO_ROOT, 'output')

# Create standardized output directory structure with subdirectories
# This creates folders for filtered data, value counts, and figures
output_dir = udf.create_output_directories(OUTPUT_DIR)

# Display confirmation
print(f"\nCreated output directory structure at: {output_dir}")

# %%


###############################################################################
# 6. Data Filtering
###################

# Define common variables for filtering operations
TECHNIQUES = ["ExoGAG", "SEC", "IP_CD9", "UC"]
POOLS = ["POOL_1", "POOL_2", "POOL_3", "NO_POOL", "POOLS_123"]
INTEREST_COLS = ['Peptide_Sequence', "Prot_Name", 'PTM', "Accession", "Peptide", "PTM_Cluster"]

# Base directories
### INPUT_DIR
### OUTPUT_DIR

#----------
# Step 6.1: Filter by Vesiclepedia
#-------------------------------
importlib.reload(udf)

# Define output directory for Vesiclepedia filtered data
VCP_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "1_filtered_dfs/vesiclepedia")

# Filter proteins by presence in Vesiclepedia database
vcp_summary = udf.filter_vesiclepedia_proteins(
    INPUT_DIR, VCP_OUTPUT_DIR, TECHNIQUES, POOLS, ptms_of_interest, 
    vesiclepedia, INTEREST_COLS
    )

# Display summary of Vesiclepedia filtering
print("\nVesiclepedia filtering summary:")
print(vcp_summary)

#----------
# Step 6.2: Filter by glycosylation
#-------------------------------
importlib.reload(udf)

# Define input/output directories for glycosylation filtering
# Input is the Vesiclepedia-filtered data from previous step
GLYC_INPUT_DIR = VCP_OUTPUT_DIR
GLYC_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "1_filtered_dfs/vesiclepedia_glycosylated")

# Filter Vesiclepedia proteins by presence of glycosylation PTMs
#glyc_summary = udf.filter_trial(
glyc_summary = udf.filter_glycosylated_proteins(
    GLYC_INPUT_DIR, GLYC_OUTPUT_DIR, TECHNIQUES, POOLS, ptms_of_interest, INTEREST_COLS)

    
    
# Display summary of glycosylation filtering
print("\nGlycosylation filtering summary:")
print(glyc_summary)

# Show only individual pools (excluding combined pools and non-pooled samples)
individual_pools = glyc_summary.loc[(glyc_summary["Pool"] != "POOLS_123") & 
                                    (glyc_summary["Pool"] != "NO_POOL"), :]
print("\nSummary for individual pools only:")
print(individual_pools)

# %%

###############################################################################
# 7. Peptide Value Counts Analysis
##################################

# Define common variables for count operations
TECHNIQUES = ["ExoGAG", "SEC", "IP_CD9", "UC"]
POOLS = ["POOL_1", "POOL_2", "POOL_3", "NO_POOL", "POOLS_123"]
INTEREST_COLS = ["Prot_Name", 'PTM', 'Peptide_Sequence', "PTM_Cluster"]

# Base output directory
### OUTPUT_DIR

#----------
# Step 7.1: Count peptides in total dataset (unfiltered)
#-------------------------------
importlib.reload(udf)

# Define input/output directories
TOTAL_INPUT_DIR = INPUT_DIR
TOTAL_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "2_value_counts/total/peptides")

print("\nAnalyzing peptide distributions in total dataset...")
udf.count_peptides_by_category(TOTAL_INPUT_DIR, TOTAL_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COLS)

#----------
# Step 7.2: Count peptides in Vesiclepedia-filtered dataset
#-------------------------------
importlib.reload(udf)

# Define input/output directories
VCP_INPUT_DIR = os.path.join(OUTPUT_DIR, "1_filtered_dfs/vesiclepedia/")
VCP_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "2_value_counts/vesiclepedia/peptides")

print("\nAnalyzing peptide distributions in Vesiclepedia-filtered dataset...")
udf.count_peptides_by_category(VCP_INPUT_DIR, VCP_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COLS)

#----------
# Step 7.3: Count peptides in glycosylated Vesiclepedia proteins
#-------------------------------
importlib.reload(udf)

# Define input/output directories
GLYC_INPUT_DIR = os.path.join(OUTPUT_DIR, "1_filtered_dfs/vesiclepedia_glycosylated/")
GLYC_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "2_value_counts/vesiclepedia_glycosylated/peptides")

print("\nAnalyzing peptide distributions in glycosylated Vesiclepedia proteins...")
udf.count_peptides_by_category(GLYC_INPUT_DIR, GLYC_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COLS)

# %%

###############################################################################
# 8. Protein PTM Analysis
#########################

# Define common variables for protein analysis
TECHNIQUES = ["ExoGAG", "SEC", "IP_CD9", "UC"]
POOLS = ["POOL_1", "POOL_2", "POOL_3", "NO_POOL", "POOLS_123"]
INTEREST_COL = "PTM_Cluster"
GROUP_BY_COL = "Prot_Name"

# Base output directory
### OUTPUT_DIR

#----------
# Step 8.1: Analyze PTM distribution across proteins in total dataset
#-------------------------------
importlib.reload(udf)

# Define input/output directories
TOTAL_INPUT_DIR = INPUT_DIR
TOTAL_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "2_value_counts/total/proteins/")

print("\nAnalyzing PTM distribution across all proteins...")
udf.count_proteins_by_ptm(
    TOTAL_INPUT_DIR, TOTAL_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COL, 
    group_by_col=GROUP_BY_COL
    )

#----------
# Step 8.2: Analyze PTM distribution across Vesiclepedia proteins
#-------------------------------
importlib.reload(udf)

# Define input/output directories
VCP_INPUT_DIR = os.path.join(OUTPUT_DIR, "1_filtered_dfs/vesiclepedia/")
VCP_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "2_value_counts/vesiclepedia/proteins")

print("\nAnalyzing PTM distribution across Vesiclepedia proteins...")
udf.count_proteins_by_ptm(
    VCP_INPUT_DIR, VCP_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COL, 
    group_by_col=GROUP_BY_COL
    )

#----------
# Step 8.3: Analyze PTM distribution across glycosylated Vesiclepedia proteins
#-------------------------------
importlib.reload(udf)

# Define input/output directories
GLYC_INPUT_DIR = os.path.join(OUTPUT_DIR, "1_filtered_dfs/vesiclepedia_glycosylated/")
GLYC_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "2_value_counts/vesiclepedia_glycosylated/proteins")

print("\nAnalyzing PTM distribution across glycosylated Vesiclepedia proteins...")
udf.count_proteins_by_ptm(
    GLYC_INPUT_DIR, GLYC_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COL, 
    group_by_col=GROUP_BY_COL
    )

# %%

###############################################################################
# Erase temporary processing dir
################################

# Erase temporary input dir this directory
udf.delete_directory(INPUT_DIR)

# %%

print("\n\nAll tables have been saved to 'output/1_filtered_dfs/' & 'output/2_value_counts/'\n\nPython pre-processing finished\nTime to move on to R!")
###############################################################################
