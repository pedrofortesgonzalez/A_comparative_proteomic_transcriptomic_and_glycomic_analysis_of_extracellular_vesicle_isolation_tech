#!/usr/bin/env python
# coding: utf-8
"""
PTM Detection from Mass Spectrometry Data.

This script analyzes Excel files from mass spectrometry to identify
proteins with post-translational modifications (PTMs).

Author: Pedro Fortes González
"""
# import libraries
import os
import pandas as pd
import numpy as np
import re
from pathlib import Path
import shutil
from tqdm.notebook import tqdm
import session_info
import csv

# Show package list
session_info.show(dependencies=False, html=False, std_lib=False)



##############################################################################
# DataFrame visualization functions


def expand_dataframe_view():
    """
    Display DataFrames with maximum visibility.
    """
    pd.set_option('display.max_rows', None)       # Show all rows
    pd.set_option('display.max_columns', None)    # Show all columns
    pd.set_option('display.width', None)          # No width restriction
    pd.set_option('display.max_colwidth', None)   # Show full cell content


def reset_dataframe_view():
    """
    Reset all display settings modified by expand_dataframe_view()
    back to pandas defaults.
    """
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')
    pd.reset_option('display.max_colwidth')


    
##############################################################################
# File preprocessing functions

def create_input_directory(RAW_DATA_DIR, PROCESSED_DATA_DIR):
    """
    Create a copy of the raw data directory in the processed data location.
    
    Args:
        RAW_DATA_DIR: Path to the source directory containing raw data
        PROCESSED_DATA_DIR: Path where the copy should be created
    
    Returns:
        Path: Path to the created processed data directory
    
    This function copies the entire directory structure from RAW_DATA_DIR
    to PROCESSED_DATA_DIR, including all subdirectories and files.
    """
    # Convert to Path objects if they are strings
    raw_dir = Path(RAW_DATA_DIR)
    processed_dir = Path(PROCESSED_DATA_DIR)
    
    # Check if source directory exists
    if not raw_dir.exists():
        raise FileNotFoundError(f"\nSource directory does not exist: {raw_dir}")
    
    # Create the target directory if it doesn't exist
    processed_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy the entire directory structure
    print(f"\nCopying data from {raw_dir} to {processed_dir}...")
    
    # Count files for progress reporting
    total_files = sum(1 for _ in raw_dir.rglob('*') if _.is_file())
    copied_files = 0
    
    # Iterate through all items in the source directory
    for item in raw_dir.rglob('*'):
        # Calculate relative path from RAW_DATA_DIR
        relative_path = item.relative_to(raw_dir)
        # Create target path
        target_path = processed_dir / relative_path
        
        if item.is_dir():
            # Create directory if it doesn't exist
            target_path.mkdir(parents=True, exist_ok=True)
        else:
            # Copy the file
            shutil.copy2(item, target_path)
            copied_files += 1
            if copied_files % 10 == 0 or copied_files == total_files:
                print(f"\nCopied {copied_files}/{total_files} files...")
    
    print(f"\nSuccessfully copied {copied_files} files to {processed_dir}")
    return processed_dir


def rename_pool_files(directory, techniques):
    """
    Rename mass spectrometry files to indicate pool status.
    
    Args:
        directory: Path to directory containing files to rename
        techniques: List of technique identifiers used in filenames
        
    This function adds "_NO_POOL" suffix to files without "POOL" in the name
    and corrects specific misnamed files.
    """
    files = [file for file in Path(directory).rglob('*')]
    for file in files:
        if "POOL" not in file.name: 
            file.unlink()
            print(f"\nFile {file.name} deleted")
            
        # Correct specific misnamed file
        if "DB search psm_ID_CD9_NO_POOL.csv" in file.name:
            new_name = file.with_name("DB search psm_IP_CD9_NO_POOL.csv")
            file.rename(new_name)
            file = new_name
            print(f"\nFile renamed: {file.name}")
            

def simplify_filenames(directory, techniques):
    """
    Remove extraneous prefixes from mass spectrometry filenames.
    
    Args:
        directory: Path to directory containing files to rename
        techniques: List of technique identifiers used in filenames
        
    Removes "DB search psm_" prefix and replaces spaces with underscores
    in filenames to create cleaner, more consistent file naming.
    """
    files = [file for file in Path(directory).rglob('*')]
    counter = 0
    
    for file in files:
        if "DB search psm" in file.name:  # Remove prefix from filename
            counter += 1
            new_name = file.with_name(
                file.name.replace(' ', '_')
                         .replace('__', '_')
                         .lstrip('DB_search_psm_')
            )
            file.rename(new_name)
            file = new_name
    
    print(f"Renamed {counter} files")
            
            
            
def standardize_csv_delimiters(input_directory, file_pattern="*.csv", backup=False):
    """
    Detect CSV delimiters and standardize all CSV files to use semicolon (;) as delimiter.
    
    Args:
        input_directory: Path to directory containing CSV files
        file_pattern: Pattern to match CSV files (default: "*.csv")
        backup: Whether to create backup of original files (default: True)
        
    Returns:
        List of files that were modified
        
    This function scans all CSV files in the specified directory (and subdirectories),
    detects their delimiter, and if not already a semicolon, converts them to use semicolon
    as the delimiter. Original files can be backed up with a .bak extension.
    """
    
    input_dir = Path(input_directory)
    modified_files = []
    
    # Find all CSV files in the directory and subdirectories
    csv_files = list(input_dir.rglob(file_pattern))
    print(f"\nFound {len(csv_files)} CSV files to check")
    
    for file_path in csv_files:
        try:
            # Detect the delimiter using csv.Sniffer
            with open(file_path, 'r', newline='') as f:
                sample = f.read(4096)
                if not sample:  # Skip empty files
                    print(f"\nSkipping empty file: {file_path}")
                    continue
                
                sniffer = csv.Sniffer()
                dialect = sniffer.sniff(sample)
                detected_delimiter = dialect.delimiter
            
            # If delimiter is already semicolon, skip this file
            if detected_delimiter == ';':
                print(f"\nFile already uses semicolon delimiter: {file_path}")
                continue
                
            # Read the CSV with the detected delimiter
            df = pd.read_csv(file_path, sep=detected_delimiter)
            
            # Verify that the DataFrame has more than one column (sanity check)
            if len(df.columns) <= 1 and detected_delimiter != ';':
                # Try with semicolon as a fallback
                try:
                    test_df = pd.read_csv(file_path, sep=';')
                    if len(test_df.columns) > 1:
                        print(f"\nSniffer failed for {file_path}, but semicolon works better")
                        continue  # File already uses semicolon effectively
                except:
                    pass  # Stick with the original detection
            
            # Create backup if requested
            if backup:
                backup_path = str(file_path) + '.bak'
                shutil.copy2(file_path, backup_path)
                print(f"\nCreated backup: {backup_path}")
            
            # Write the file back with semicolon delimiter
            df.to_csv(file_path, sep=';', index=False)
            print(f"\nConverted {file_path} from '{detected_delimiter}' to ';' delimiter")
            modified_files.append(file_path)
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    print(f"\nConverted {len(modified_files)} files to use semicolon delimiter")
    return modified_files            


def combine_pool_files(directory, techniques):
    """
    Combine individual pool files (1, 2, 3) into a single combined file.
    
    Args:
        directory: Path to directory containing pool files
        techniques: List of technique identifiers to process
        
    For each technique, finds all pool files and combines them into a single CSV
    with the naming pattern "{technique}_POOLS_123.csv"
    """
    for technique in techniques:
        print(f"\n{technique}")
        # Create search pattern with regex to find pool files
        pattern = re.compile(f"{technique}.*POOL(?!S).*\d\.csv")
        # Find files matching the regex pattern
        files = [
            file for file in Path(directory).rglob('*') 
            if pattern.search(file.name)
        ]
        
        # Process each file for the current technique
        file_count = 0
        for file in files:
            file_count += 1
            df = pd.read_csv(file, sep=";")
            print(f"\nDimensions {technique} pool {file_count} df: {df.shape}")
        
        # Combine all pool files into one DataFrame
        combined_pool_df = pd.concat(
            [pd.read_csv(file, sep=";") for file in files], 
            ignore_index=True
        )
        
        # Save the resulting DataFrame
        pool_output_path = directory + f"/{technique}_POOLS_123.csv"
        combined_pool_df.to_csv(pool_output_path, index=False, sep=";")
        print(f"Dimensions {technique} combined pools: {combined_pool_df.shape}")






##############################################################################
# Data extraction and processing functions


def extract_protein_names(input_dir, target_col, new_col, pattern):
    """
    Extract protein accession numbers from a target column using regex pattern.
    
    Args:
        input_dir: Directory containing CSV files to process
        target_col: Column containing protein information to extract from
        new_col: Name of the new column to store extracted protein names
        pattern: Regular expression pattern to identify protein accessions
        
    This function processes all CSV files in the input directory, extracts protein
    names using the provided regex pattern, and adds them as a new column.
    """
    dataframes = list(Path(input_dir).glob('*.csv'))
    print("Processing:")
    regex = re.compile(pattern)
    file_counter = 0
    
    for csv_file in tqdm(dataframes, desc="Processing CSV files"):
        file_counter += 1
        df = pd.read_csv(csv_file, encoding='latin1', sep=";", low_memory=False)
        
        if target_col in df.columns:
            # Extract and clean desired patterns
            df[new_col] = df[target_col].apply(
                lambda x: extract_and_clean_accessions(x, regex) if pd.notna(x) else []
            )
            # Expand rows with multiple accessions
            df = df.explode(new_col)
            # Save the modified DataFrame and overwrite original file
            df.to_csv(csv_file, index=False, sep=";")
            print(f"File #{file_counter}: {csv_file.name} overwritten with new column: '{new_col}'.")
        else:
            print(f"Column {target_col} not found in {csv_file}")


def extract_and_clean_accessions(accession_str, regex):
    """
    Extract protein accessions from a string using the provided regex pattern.
    
    Args:
        accession_str: String containing protein accession information
        regex: Compiled regular expression pattern to use for extraction
        
    Returns:
        List of matched accessions
    """
    matches = regex.findall(accession_str)  # Find all matches in the string
    return matches


def extract_peptide_sequences(input_dir, target_col, new_col, pattern):
    """
    Extract peptide sequences from a target column using regex pattern.
    
    Args:
        input_dir: Directory containing CSV files to process
        target_col: Column containing peptide information to extract from
        new_col: Name of the new column to store extracted peptide sequences
        pattern: Regular expression pattern to remove from peptide sequences
        
    This function processes all CSV files in the input directory, cleans peptide
    sequences by removing the specified pattern, and adds them as a new column.
    """
    dataframes = list(Path(input_dir).glob('*.csv'))
    print("Processing:")
    regex = re.compile(pattern)
    
    for csv_file in tqdm(dataframes, desc="Processing CSV files"):
        df = pd.read_csv(csv_file, sep=";")
        
        if target_col in df.columns:
            # Extract desired pattern and assign to new column
            df[new_col] = df[target_col].str.replace(regex, '', regex=True)
            # Save the modified DataFrame and overwrite original file
            df.to_csv(csv_file, sep=";", index=False)
            print(f"{csv_file.name} overwritten with '{new_col}'")
        else:
            print(f"Column {target_col} not found in {csv_file}")
            
def classify_ptm_types(input_dir, target_col, new_col):
    """
    Classify PTMs into categories based on pattern matching in glycosylation markers.
    
    Args:
        input_dir: Directory containing CSV files to process
        target_col: Column containing PTM information to classify
        new_col: Name of the new column to store the PTM classification
        
    This function categorizes PTMs into: Fucosialylated, Fucosylated, Sialylated,
    Oligomannose, No PTM, or Other based on specific glycosylation markers.
    """
    dataframes = list(Path(input_dir).glob('*.csv'))
    print("Processing:")
    
    for csv_file in tqdm(dataframes, desc="Processing CSV files"):
        df = pd.read_csv(csv_file, sep=";")
        
        if target_col in df.columns:
            print(f"\n{csv_file}: creating {new_col}: ")
            df[target_col] = df[target_col].astype(str)  # Convert column to string for pattern matching
            
            # Define conditions for PTM classification
            conditions = [
                # Fucosialylated: contains fucose/dHex markers AND sialic acid markers
                df[target_col].str.contains(r'dHex|Fucos|Biantennary', regex=True) & 
                df[target_col].str.contains(r'NeuGc|NeuAc|Kdn|Neuraminic', regex=True),
                
                # Fucosylated: contains fucose/dHex markers BUT NO sialic acid markers
                df[target_col].str.contains(r'dHex|Fucos|Biantennary', regex=True) & 
                ~df[target_col].str.contains(r'NeuGc|NeuAc|Kdn|Neuraminic', regex=True),
                
                # Sialylated: contains sialic acid markers BUT NO fucose/dHex markers
                ~df[target_col].str.contains(r'dHex|Fucos|Biantennary', regex=True) & 
                df[target_col].str.contains(r'NeuGc|NeuAc|Kdn|Neuraminic', regex=True),
                
                # Oligomannose: contains HexNAc markers but no fucose or sialic acid markers
                df[target_col].str.contains(r'HexNAc|N\-linked\sglycan\score', regex=True) & 
                ~df[target_col].str.contains(r'dHex|Fucos|NeuGc|NeuAc|Kdn|Neuraminic|Biantennary', regex=True),
                
                # No PTM: explicitly "nan" values
                df[target_col] == "nan"
            ]
            
            # Define values corresponding to conditions (in order)
            choices = ['Fucosialylated', 'Fucosylated', 'Sialylated', 'Oligomannose', "No PTM"]
            
            # Apply conditions to create the new column
            df[new_col] = np.select(conditions, choices, default='Other')
            
            # Handle empty cells in target_col, marking them as "No PTM"
            empty_row_catcher = (
                df[target_col].isna() | 
                (df[target_col] == '') | 
                (df[target_col] == " ") | 
                (df[target_col] == "nan")
            )
            df.loc[empty_row_catcher, new_col] = 'No PTM'
        else:
            df[new_col] = "No PTM"
            
        # Check result
        print(df[new_col].value_counts(dropna=False))
        
        # Save the modified DataFrame and overwrite original file
        df.to_csv(csv_file, index=False, sep=";")
        print(f"{csv_file.name} overwritten with {new_col}")



######################################################################################################################
# Function for creation of output directories

def create_output_directories(input_directory):
    """
    Create standardized directory structure for analysis outputs.
    
    Args:
        input_directory: Path to the input directory containing analysis files
        
    Creates a structured output directory with subdirectories for:
    1. Filtered DataFrames (vesiclepedia, glycosylated)
    2. Value counts (total, vesiclepedia, vesiclepedia_glycosylated)
    3. Figures for visualizations
    
    The output directory is created at the same level as the input and script
    directories.
    """
    # Create main output directory
    input_dir = Path(input_directory)
    parent_dir = input_dir.parent
    output_dir = parent_dir / ('output')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories for filtered DataFrames
    filtered_dir = output_dir / "1_filtered_dfs"
    filtered_dir.mkdir(parents=True, exist_ok=True)
    (filtered_dir / "vesiclepedia").mkdir(parents=True, exist_ok=True)
    (filtered_dir / "vesiclepedia_glycosylated").mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories for value counts
    counts_dir = output_dir / "2_value_counts"
    counts_dir.mkdir(parents=True, exist_ok=True)
    
    # Total counts directories
    total_dir = counts_dir / "total"
    total_dir.mkdir(parents=True, exist_ok=True)
    (total_dir / "peptides").mkdir(parents=True, exist_ok=True)
    (total_dir / "proteins").mkdir(parents=True, exist_ok=True)
    
    # Vesiclepedia counts directories
    vcp_dir = counts_dir / "vesiclepedia"
    vcp_dir.mkdir(parents=True, exist_ok=True)
    (vcp_dir / "peptides").mkdir(parents=True, exist_ok=True)
    (vcp_dir / "proteins").mkdir(parents=True, exist_ok=True)
    
    # Vesiclepedia glycosylated counts directories
    vcp_glyc_dir = counts_dir / "vesiclepedia_glycosylated"
    vcp_glyc_dir.mkdir(parents=True, exist_ok=True)
    (vcp_glyc_dir / "peptides").mkdir(parents=True, exist_ok=True)
    (vcp_glyc_dir / "proteins").mkdir(parents=True, exist_ok=True)
    
    # Create directory for figures
    figures_dir = output_dir / "3_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    (figures_dir/ "boxplots").mkdir(parents=True, exist_ok=True)
    (figures_dir/ "sector_diagrams/no_text").mkdir(parents=True, exist_ok=True)
    (figures_dir/ "sector_diagrams/textbox").mkdir(parents=True, exist_ok=True)
    
    # Return the main output directory path for convenience
    return output_dir


    
##############################################################################
# Functions for data filtering


def filter_vesiclepedia_proteins(input_dir, output_dir, techniques, pools, ptms_of_interest, vesiclepedia, interest_cols):
    """
    Filter protein data to keep only those present in Vesiclepedia database.
    
    Args:
        input_dir: Directory containing CSV files to process
        output_dir: Directory where filtered files will be saved
        techniques: List of isolation techniques in the filenames
        pools: List of pool identifiers in the filenames
        ptms_of_interest: List of PTMs to analyze
        vesiclepedia: DataFrame containing Vesiclepedia protein database
        interest_cols: Columns to include in filtered output
        
    Returns:
        DataFrame summarizing the results of filtering across all files
    """
    output_dir = Path(output_dir)
    input_dir = Path(input_dir)

    # Create summary DataFrame
    summary_df = pd.DataFrame()

    # Process all CSV files in the provided directory
    csv_files = sorted(list(Path(input_dir).glob('*.csv')), key=lambda x: x.name)

    for csv_file in tqdm(csv_files, desc="Processing CSV files"):
        df = pd.read_csv(csv_file, sep=";")
        
        # Extract sample information from filename
        sample_name = csv_file.stem
        sample_name = sample_name.replace(' ', '_').replace('__', '_').rstrip('.csv').lstrip('DB_search_psm_')
        print(sample_name)

        # Identify technique and pool
        technique = next((technique for technique in techniques if re.search(technique, sample_name)), "")
        pool = next((pool for pool in pools if re.search(pool, sample_name)), "")
    
        # Filter proteins by presence in Vesiclepedia
        peptides_count = df["Peptide Sequence"].value_counts()
        print(f"Total peptides: {len(peptides_count)}")
        
        condition = df["Prot Name"].isin(vesiclepedia["Accession"])
        df_filtered = df.loc[condition, interest_cols]
        
        # Save filtered data
        output_file = output_dir / f"{sample_name}_filtered_vcp.csv"
        df_filtered.to_csv(output_file, index=False, sep=";")
        
        filtered_peptides_count = df_filtered["Peptide Sequence"].value_counts()
        print(f"Peptides in Vesiclepedia: {len(filtered_peptides_count)}")
    
        # Calculate statistics
        n_prots = len(df["Prot Name"].value_counts())
        n_ptms = len(df["PTM"].value_counts())
        n_peptides = len(df["Peptide Sequence"].value_counts())
       
        n_prots_filt = len(df_filtered["Prot Name"].value_counts())
        n_ptms_filt = len(df_filtered["PTM"].value_counts())
        n_peptides_filt = len(df_filtered["Peptide Sequence"].value_counts())
        
        ratio_ptms_prots = (100 * n_prots_filt) / n_prots if n_prots > 0 else 0
        ratio_ptms_peptides = (100 * n_peptides_filt) / n_peptides if n_prots > 0 else 0
        
        # Add results to summary DataFrame
        temp_df = pd.DataFrame({
            "Sample": [sample_name], 
            "Technique": [technique],
            "Pool": [pool], 
            "n Proteins Total": [n_prots],
            "n PTM Total": [n_ptms], 
            "n Peptides Total": [n_peptides],
            "n Proteins Filtered": [n_prots_filt], 
            "n PTM Filtered": [n_ptms_filt],
            "n Peptides Filtered": [n_peptides_filt], 
            "% Proteins with PTM / Total": [ratio_ptms_prots], 
            "% Peptides with PTM / Total": [ratio_ptms_peptides]
        })

        summary_df = pd.concat([summary_df, temp_df], ignore_index=True)
        
    # Sort summary by technique and pool
    summary_df = summary_df.sort_values(by=['Technique', "Pool"])
    
    # Export summary results in different configurations
    # All samples
    summary_file_path = output_dir / 'summary_all.csv'
    summary_df.to_csv(summary_file_path, index=False, sep=";")

    # POOLS_123 only
    sf = summary_df.loc[summary_df["Pool"]=="POOLS_123"]
    filepath = output_dir / 'summary_pools_123.csv'
    sf.to_csv(filepath, index=False, sep=";")

    # POOLS_123 & NO_POOL
    sf = summary_df.loc[(summary_df["Pool"]=="POOLS_123") | (summary_df["Pool"]=="NO_POOL")]
    filepath = output_dir / 'summary_123nopool.csv'
    sf.to_csv(filepath, index=False, sep=";")

    # Individual pools (POOL 1, POOL 2, POOL 3)
    sf = summary_df.loc[(summary_df["Pool"]!="POOLS_123") & (summary_df["Pool"]!="NO_POOL")]
    filepath = output_dir / 'summary_individual_pools.csv'
    sf.to_csv(filepath, index=False, sep=";")
    
    return summary_df

            

def filter_glycosylated_proteins(input_dir, output_dir, techniques, pools, ptms_of_interest, interest_cols):
    """
    Filter protein data to keep only those with glycosylation PTMs of interest.
    
    Args:
        input_dir: Directory containing CSV files to process
        output_dir: Directory where filtered files will be saved
        techniques: List of isolation techniques in the filenames
        pools: List of pool identifiers in the filenames
        ptms_of_interest: List of glycosylation PTMs to filter for
        interest_cols: Columns to include in filtered output
        
    Returns:
        DataFrame summarizing the results of filtering across all files
    """
    output_dir = Path(output_dir)
    input_dir = Path(input_dir)

    # Create summary DataFrame
    summary_df = pd.DataFrame()

    # Process all CSV files in the provided directory
    csv_files = sorted(list(Path(input_dir).glob('*.csv')), key=lambda x: x.name)

    for csv_file in tqdm(csv_files, desc="Processing CSV files"):
        if str("summary") not in csv_file.name:
            df = pd.read_csv(csv_file, sep=";")
            
            # Extract sample information from filename
            sample_name = csv_file.stem
            sample_name = sample_name.replace(' ', '_').replace('__', '_').rstrip('.csv').lstrip('DB_search_psm_')
            print(f"\n{sample_name}")

            # Identify technique and pool
            technique = next((technique for technique in techniques if re.search(technique, sample_name)), "")
            pool = next((pool for pool in pools if re.search(pool, sample_name)), "")
            # Filter proteins by glycosylation PTMs
            peptides_count = df["Peptide Sequence"].value_counts()
            print(f"\nTotal peptides: {len(peptides_count)}")
            
            condition = df["PTM"].isin(ptms_of_interest)
            df_filtered = df.loc[condition, interest_cols]
            
            # Save filtered data
            output_file = output_dir / f"{sample_name}_filtered_glyc.csv"
            df_filtered.to_csv(output_file, index=False, sep=";")
            
            filtered_peptides_count = df_filtered["Peptide Sequence"].value_counts()
            print(f"\nGlycosylated peptides: {len(filtered_peptides_count)}")

            # Filter proteins by glycosylation PTMs
            peptides_count = df["Peptide Sequence"].value_counts()
            print(f"\nTotal peptides: {len(peptides_count)}")
            
            condition = df["PTM"].isin(ptms_of_interest)
            df_filtered = df.loc[condition, interest_cols]
            
            # Save filtered data
            output_file = output_dir / f"{sample_name}_filtered_glyc.csv"
            df_filtered.to_csv(output_file, index=False, sep=";")
            
            filtered_peptides_count = df_filtered["Peptide Sequence"].value_counts()
            print(f"\nGlycosylated peptides: {len(filtered_peptides_count)}")
    
            # Calculate statistics
            n_prots = len(df["Prot Name"].value_counts())
            n_ptms = len(df["PTM"].value_counts())
            n_peptides = len(df["Peptide Sequence"].value_counts())
       
            n_prots_filt = len(df_filtered["Prot Name"].value_counts())
            n_ptms_filt = len(df_filtered["PTM"].value_counts())
            n_peptides_filt = len(df_filtered["Peptide Sequence"].value_counts())
        
            ratio_ptms_prots = (100 * n_prots_filt) / n_prots if n_prots > 0 else 0
            ratio_ptms_peptides = (100 * n_peptides_filt) / n_peptides if n_prots > 0 else 0
        
            # Add results to summary DataFrame
            temp_df = pd.DataFrame({
                "Sample": [sample_name], 
                "Technique": [technique],
                "Pool": [pool], 
                "n Proteins Total": [n_prots],
                "n PTM Total": [n_ptms], 
                "n Peptides Total": [n_peptides],
                "n Proteins Filtered": [n_prots_filt], 
                "n PTM Filtered": [n_ptms_filt],
                "n Peptides Filtered": [n_peptides_filt], 
                "% Proteins with PTM / Total": [ratio_ptms_prots], 
                "% Peptides with PTM / Total": [ratio_ptms_peptides]
            })

            summary_df = pd.concat([summary_df, temp_df], ignore_index=True)
        
    # Sort summary by technique and pool
    summary_df = summary_df.sort_values(by=['Technique', "Pool"])
    
    # Export summary results in different configurations
    # All samples
    summary_file_path = output_dir / 'summary_all.csv'
    summary_df.to_csv(summary_file_path, index=False, sep=";")

    # POOLS_123 only
    sf = summary_df.loc[summary_df["Pool"]=="POOLS_123"]
    filepath = output_dir / 'summary_pools_123.csv'
    sf.to_csv(filepath, index=False, sep=";")

    # POOLS_123 & NO_POOL
    sf = summary_df.loc[(summary_df["Pool"]=="POOLS_123") | (summary_df["Pool"]=="NO_POOL")]
    filepath = output_dir / 'summary_123nopool.csv'
    sf.to_csv(filepath, index=False, sep=";")

    # Individual pools (POOL 1, POOL 2, POOL 3)
    sf = summary_df.loc[(summary_df["Pool"]!="POOLS_123") & (summary_df["Pool"]!="NO_POOL")]
    filepath = output_dir / 'summary_individual_pools.csv'
    sf.to_csv(filepath, index=False, sep=";")
    
    return summary_df



######################################################################################################################
# Functions for counting


def count_peptides_by_category(input_dir, output_dir, techniques, pools, interest_cols):
    """
    Count and analyze peptides by different categories of interest.
    
    Args:
        input_dir: Directory containing CSV files to process
        output_dir: Directory where count files will be saved
        techniques: List of isolation techniques in the filenames
        pools: List of pool identifiers in the filenames
        interest_cols: Columns to perform value counts on
        
    For each column of interest, creates a CSV file with value counts,
    showing the distribution of values and their percentages.
    """
    # Convert path arguments to Path objects
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    
    # Process all CSV files in the input directory
    csv_files = sorted(list(Path(input_dir).glob('*.csv')), key=lambda x: x.name)
    print("Processing:")
    
    for csv_file in tqdm(csv_files, desc="Processing CSV files"):
        if "summary" not in csv_file.name:
            # Read dataframe
            df = pd.read_csv(csv_file, sep=";")
            
            # Extract sample information from filename
            sample_name = csv_file.stem
            print(sample_name)
            
            # Identify technique and pool (not used currently but kept for future use)
            technique = [t for t in techniques if re.search(t, sample_name)]
            pool = [p for p in pools if re.search(p, sample_name)]
            
            # Count values for each column of interest
            for col in interest_cols:
                if col in df.columns:
                    # Get value counts
                    value_counts = df[col].value_counts()
                    
                    # Create summary DataFrame with counts and percentages
                    summary = pd.DataFrame({
                        'Value': value_counts.index, 
                        'Count': value_counts.values
                    })
                    
                    # Calculate percentage of each value
                    total_count = summary['Count'].sum()
                    summary['Percentage'] = (summary['Count'] / total_count) * 100
                    
                    # Export the counts to CSV
                    output_file = output_dir / f"{sample_name}_{col}.csv"
                    summary.to_csv(output_file, index=False, sep=";")
                else:
                    print(f"Warning: Column '{col}' not found in file {sample_name}.")


def count_proteins_by_ptm(input_dir, output_dir, techniques, pools, interest_col, group_by_col="Prot Name"):
    """
    Count proteins by PTM types, grouping at the protein level.
    
    Args:
        input_dir: Directory containing CSV files to process
        output_dir: Directory where count files will be saved
        techniques: List of isolation techniques in the filenames
        pools: List of pool identifiers in the filenames
        interest_col: Column to count values from (typically "PTM cluster")
        group_by_col: Column to group by (default "Prot Name")
        
    Groups proteins and counts the different types of PTMs for each protein,
    then summarizes the distribution of PTM types across all proteins.
    """
    # Convert path arguments to Path objects
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    
    # Process all CSV files in the input directory
    csv_files = sorted(list(Path(input_dir).glob('*.csv')), key=lambda x: x.name)
    print("Processing:")
    
    for csv_file in tqdm(csv_files, desc="Processing CSV files"):
        if "summary" not in csv_file.name:
            # Read dataframe
            df = pd.read_csv(csv_file, sep=";")
            
            # Extract sample information from filename
            sample_name = csv_file.stem
            print(sample_name)
            
            # Identify technique and pool (not used currently but kept for future use)
            technique = [t for t in techniques if re.search(t, sample_name)]
            pool = [p for p in pools if re.search(p, sample_name)]
            
            # Group by protein name and count PTM types
            if group_by_col in df.columns and group_by_col == "Prot Name":
                # Group by protein and count PTM types
                df_grouped = df.groupby(group_by_col)[interest_col].value_counts()
                df_grouped = df_grouped.reset_index(name="Count")
                
                # Count the occurrence of each PTM type across all proteins
                ptm_counts = df_grouped[interest_col].value_counts()
                
                # Create summary DataFrame with counts and percentages
                summary = pd.DataFrame({
                    'Value': ptm_counts.index, 
                    'Count': ptm_counts.values
                })
                
                # Clean up value names
                summary["Value"] = summary["Value"].astype(str)
                summary['Value'] = summary['Value'].str.replace(r"\W", "", regex=True)
                
                # Calculate percentage of each PTM type
                total_count = summary['Count'].sum()
                summary['Percentage'] = (summary['Count'] / total_count) * 100
                
                # Export the counts to CSV
                output_file = output_dir / f"{sample_name}_PTM_types_by_protein.csv"
                summary.to_csv(output_file, index=False, sep=";")
            else:
                print(f"Warning: Required grouping column '{group_by_col}' not found in file {sample_name}.")
                continue
            


######################################################################################################################
# Functions deleting the temporary dir. created for processing mass spectrometry results


def delete_directory(directory_path):
    """
    Delete a directory and all its contents without confirmation.
    
    Args:
        directory_path (str): Path to the directory to delete
    
    Returns:
        bool: True if directory was successfully deleted, False otherwise
    """
    try:
        # Check if directory exists
        if os.path.exists(directory_path):
            # Delete the directory and all its contents
            shutil.rmtree(directory_path)
            print(f"\nDeleted directory: {directory_path}")
            return True
        else:
            print(f"\nDirectory does not exist: {directory_path}")
            return False
            
    except Exception as e:
        print(f"\nError deleting directory: {str(e)}")
        return False            