#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PTM Detection: Analysis of glycosylated proteins from mass spectrometry data.

This script processes mass spectrometry data to identify post-translational
modifications, with a focus on protein glycosylation.

Author: Pedro Fortes Gonzalez
Date:
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

# Set working directory
WORKING_DIR = r"/Users/pedrofortesgonzalez/Desktop/NefroCHUS/projects/240614_detectar_prots_glucosiladas_bravo-susana_real/script_for_repo/scripts/"
os.chdir(WORKING_DIR)

# Show package list
session_info.show(dependencies=False, html=False, std_lib=False)

# Configure progress bar for pandas operations
tqdm.pandas(desc="Processing")

# Import custom functions module
sys.path.append(WORKING_DIR)
import user_defined_funcs as udf

# Expand dataframe display for better visualization
udf.expand_dataframe_view()

# %%

###############################################################################
# 2. Data Import
################
# Import list of PTMs (glycosylations) of interest
GLYCOSYLATION_LIST_PATH = '../input_data/glycosylation_list.csv'
ptm_df = pd.read_csv(GLYCOSYLATION_LIST_PATH, usecols=[0], sep=";")  # Select only first column
ptms_of_interest = list(set(ptm_df.iloc[:, 0].dropna().unique()))  # Remove NaNs and duplicates

# Verify PTM list import
print(f"\nNumber of PTMs imported: {len(ptms_of_interest)}")

# Import Vesiclepedia database (downloaded on 2024-07-12)
VESICLEPEDIA_PATH = "../input_data/Vesiclepedia_proteins_240712.csv"
vesiclepedia = pd.read_csv(VESICLEPEDIA_PATH)

# Verify Vesiclepedia database import
print(f"\nNumber of entries in Vesiclepedia: {vesiclepedia.shape[0]}")

# %%

###############################################################################
# 3. File Preprocessing
#######################
# Reload the updated functions module to ensure we have the latest version
importlib.reload(udf)

# Define input directory and technique identifiers
DATA_DIR = '../input_data/input'
TECHNIQUES = ["ExoGAG", "SEC", "IP_CD9", "UC"]

# Step 3.1: Ensure all files are properly named with pool indicators
udf.rename_pool_files(DATA_DIR, TECHNIQUES)  # Add "_NO_POOL" to files without pool identifier

# Step 3.2: Create combined pool files from individual pools
udf.combine_pool_files(DATA_DIR, TECHNIQUES)  # Combine pool1 + pool2 + pool3 into POOLS_123 files

# Step 3.3: Clean up filenames by removing prefixes and standardizing format
udf.simplify_filenames(DATA_DIR, TECHNIQUES)

# %%

###############################################################################
# 4. Column Cleanup and Extraction
##################################

# Step 4.1: Extract protein accession numbers from the Accession column
importlib.reload(udf)

# Define inputs for protein name extraction
DATA_DIR = '../input_data/input'
ACCESSION_COL = 'Accession'
PROTEIN_NAME_COL = "Prot Name"

# Regular expression pattern to extract standardized protein accession numbers
# Matches patterns like UniProt accessions (e.g., P12345, Q9XYZ0) and other common formats
PROTEIN_PATTERN = r"([A-Z]\d[A-Z0-9]{3}[0-9]-?\d*|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z]?[A-Z0-9]+[0-9])"

# Extract protein names from accession information and create a new column
udf.extract_protein_names(DATA_DIR, ACCESSION_COL, PROTEIN_NAME_COL, PROTEIN_PATTERN)

# Step 4.2: Clean peptide sequences by removing modifications in parentheses
importlib.reload(udf)
PEPTIDE_COL = 'Peptide'
PEPTIDE_SEQ_COL = "Peptide Sequence"
PEPTIDE_PATTERN = r'\([^)]*\)'  # Pattern to remove content inside parentheses

# Extract clean peptide sequences without modification annotations
udf.extract_peptide_sequences(DATA_DIR, PEPTIDE_COL, PEPTIDE_SEQ_COL, PEPTIDE_PATTERN)

# Step 4.3: Classify PTMs into standardized categories
importlib.reload(udf)
PTM_COL = 'PTM'
PTM_CLUSTER_COL = "PTM cluster"

# Create a new column with standardized PTM classifications
udf.classify_ptm_types(DATA_DIR, PTM_COL, PTM_CLUSTER_COL)

# %%

###############################################################################
# 5. Create Output Directory Structure
###############################################################################

# Reload the functions module to ensure we're using the latest version
importlib.reload(udf)

# Define input directory
DATA_DIR = '../output'

# Create standardized output directory structure with subdirectories
# This creates folders for filtered data, value counts, and figures
output_dir = udf.create_output_directories(DATA_DIR)

# Display confirmation
print(f"\nCreated output directory structure at: {output_dir}")

# %%


###############################################################################
# 6. Data Filtering
###################

# Define common variables for filtering operations
TECHNIQUES = ["ExoGAG", "SEC", "IP_CD9", "UC"]
POOLS = ["POOL_1", "POOL_2", "POOL_3", "NO_POOL", "POOLS_123"]
INTEREST_COLS = ['Peptide Sequence', "Prot Name", 'PTM', "Accession", 
                 "Peptide", "PTM cluster"]

# Base directories
DATA_DIR = '../input_data/input'
OUTPUT_BASE_DIR = '../output'

#----------
# Step 6.1: Filter by Vesiclepedia
#-------------------------------
importlib.reload(udf)

# Define output directory for Vesiclepedia filtered data
VCP_OUTPUT_DIR = f"{OUTPUT_BASE_DIR}/1_filtered_dfs/vesiclepedia"

# Filter proteins by presence in Vesiclepedia database
vcp_summary = udf.filter_vesiclepedia_proteins(
    DATA_DIR, VCP_OUTPUT_DIR, TECHNIQUES, POOLS, ptms_of_interest, 
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
GLYC_OUTPUT_DIR = f"{OUTPUT_BASE_DIR}/1_filtered_dfs/vesiclepedia_glycosylated"

# Filter Vesiclepedia proteins by presence of glycosylation PTMs
glyc_summary = udf.filter_glycosylated_proteins(
    GLYC_INPUT_DIR, GLYC_OUTPUT_DIR, TECHNIQUES, POOLS, ptms_of_interest, 
    INTEREST_COLS
    )

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
INTEREST_COLS = ["Prot Name", 'PTM', 'Peptide Sequence', "PTM cluster"]

# Base output directory
OUTPUT_BASE_DIR = '../output'

#----------
# Step 7.1: Count peptides in total dataset (unfiltered)
#-------------------------------
importlib.reload(udf)

# Define input/output directories
TOTAL_INPUT_DIR = '../input_data/input/'
TOTAL_OUTPUT_DIR = f"{OUTPUT_BASE_DIR}/2_value_counts/total/peptides"

print("\nAnalyzing peptide distributions in total dataset...")
udf.count_peptides_by_category(TOTAL_INPUT_DIR, TOTAL_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COLS)

#----------
# Step 7.2: Count peptides in Vesiclepedia-filtered dataset
#-------------------------------
importlib.reload(udf)

# Define input/output directories
VCP_INPUT_DIR = f"{OUTPUT_BASE_DIR}/1_filtered_dfs/vesiclepedia/"
VCP_OUTPUT_DIR = f"{OUTPUT_BASE_DIR}/2_value_counts/vesiclepedia/peptides"

print("\nAnalyzing peptide distributions in Vesiclepedia-filtered dataset...")
udf.count_peptides_by_category(VCP_INPUT_DIR, VCP_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COLS)

#----------
# Step 7.3: Count peptides in glycosylated Vesiclepedia proteins
#-------------------------------
importlib.reload(udf)

# Define input/output directories
GLYC_INPUT_DIR = f"{OUTPUT_BASE_DIR}/1_filtered_dfs/vesiclepedia_glycosylated/"
GLYC_OUTPUT_DIR = f"{OUTPUT_BASE_DIR}/2_value_counts/vesiclepedia_glycosylated/peptides"

print("\nAnalyzing peptide distributions in glycosylated Vesiclepedia proteins...")
udf.count_peptides_by_category(GLYC_INPUT_DIR, GLYC_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COLS)

# %%

###############################################################################
# 8. Protein PTM Analysis
#########################

# Define common variables for protein analysis
TECHNIQUES = ["ExoGAG", "SEC", "IP_CD9", "UC"]
POOLS = ["POOL_1", "POOL_2", "POOL_3", "NO_POOL", "POOLS_123"]
INTEREST_COL = "PTM cluster"
GROUP_BY_COL = "Prot Name"

# Base output directory
OUTPUT_BASE_DIR = '../output'

#----------
# Step 8.1: Analyze PTM distribution across proteins in total dataset
#-------------------------------
importlib.reload(udf)

# Define input/output directories
TOTAL_INPUT_DIR = '../input_data/input/'
TOTAL_OUTPUT_DIR = f"{OUTPUT_BASE_DIR}/2_value_counts/total/proteins/"

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
VCP_INPUT_DIR = f"{OUTPUT_BASE_DIR}/1_filtered_dfs/vesiclepedia/"
VCP_OUTPUT_DIR = f"{OUTPUT_BASE_DIR}/2_value_counts/vesiclepedia/proteins"

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
GLYC_INPUT_DIR = f"{OUTPUT_BASE_DIR}/1_filtered_dfs/vesiclepedia_glycosylated/"
GLYC_OUTPUT_DIR = f"{OUTPUT_BASE_DIR}/2_value_counts/vesiclepedia_glycosylated/proteins"

print("\nAnalyzing PTM distribution across glycosylated Vesiclepedia proteins...")
udf.count_proteins_by_ptm(
    GLYC_INPUT_DIR, GLYC_OUTPUT_DIR, TECHNIQUES, POOLS, INTEREST_COL, 
    group_by_col=GROUP_BY_COL
    )

# %%

print("\n\nAll tables have been saved to 'output/1_filtered_dfs/' & 'output/2_value_counts/'\nPython pre-processing finished\nTime to move on to R!")
###############################################################################
