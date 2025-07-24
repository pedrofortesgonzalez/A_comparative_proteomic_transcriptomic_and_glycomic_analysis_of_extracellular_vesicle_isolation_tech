#!/usr/bin/env python
# coding: utf-8
"""
PTM Detection from Mass Spectrometry Data.

This script contains user-defined functions (UDFs) that enable
script no. 1 to run.

Author: Pedro Fortes González (refactored)
"""
# import libraries
import os
import sys
import pandas as pd
import numpy as np
import re
from pathlib import Path
import shutil
from tqdm.notebook import tqdm
import csv


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
# Basic variables functions

def set_pools():
    # Ask if this is a pooled analysis
    pool_analysis_confirmation = input("\nIs this a pooled analysis? (y/n): ").strip().lower()
    # Check if answer is affirmative
    affirmative_answers = ['y', 'yes', 'sí', 'si', 's']
    is_pooled = pool_analysis_confirmation in affirmative_answers
    # Return
    return is_pooled


def set_techs():   
    # Ask user to specify techniques in the analysis
    print(f"\n-----------------------------------------------", f"\n-----------------------------------------------")
    print("Choose a set of analysed techniques: ", f"\n  0 --> None", f"\n  1 --> Default ('Raw_milk', 'ExoGAG', 'SEC', 'IP_CD9', 'UC')")
    print(f"  2 --> Custom", f"\n-----------------------------------------------")
    choose_techniques = int(input("\nType an option (0/1/2): "))

    if choose_techniques == 0:
        TECHNIQUES = []
        print(f"\nYour chose to run the analysis without specified techniques")
    elif choose_techniques == 1:
        TECHNIQUES = ['Raw_milk', 'ExoGAG', 'SEC', 'IP_CD9', 'UC']
        print(f"\nYour chosen set of techniques is the Default set:", TECHNIQUES)
    else:
        TECHNIQUES = re.split(r",\s*", input("\nEnter technique names sepparated by commas: "))
        print(f"\nYour chosen set of techniques is the Custom set: ", TECHNIQUES)

    # has_techniques is defined as a bool
    has_techniques = len(TECHNIQUES) > 0
    # other ways to define it
    #has_techniques = True if len(TECHNIQUES) > 0 else False
    #has_techniques = bool(TECHNIQUES)

    # Return
    return has_techniques, TECHNIQUES


def set_analysis(var1, var2, TECHNIQUES):
    print("\nYour analysis is: ")
    if not var1 and not var2:
        print("Not pooled and has no techniques\n\n!!!!!!!!!!!!!!!!!!!!\nUNSUPPORTED ANALYSIS\nHalting execution...")
        sys.exit(2)
    elif not var1:
        print(f"Unpooled technique analysis: comparing {len(TECHNIQUES)} different techniques")
    elif not var2:
        print(f"Pool-only analysis: comparing only pools")
    else:
        print(f"Pooled technique analysis: comparing pools of {len(TECHNIQUES)} different techniques")
    
    
    
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


def rename_pool_files(directory, techniques, is_pooled, has_techniques):
    """
    Rename mass spectrometry files to indicate pool status.
    
    Args:
        directory: Path to directory containing files to rename
        techniques: List of technique identifiers used in filenames

    Returns:
        has_techniques: Boolean indicating if this analisys is for comparing techniques
        
    This function adds "_NO_POOL" suffix to files without "POOL" in the name
    and corrects specific misnamed files.
    """

    # Iterate over files
    files = [file for file in Path(directory).rglob('*')]

    for file in files:
        sep = "*"
        starnum = len(file.name)
        
        print(f"\n{sep*starnum*2}")
        print(f"{sep*(int(starnum/2) - 2)} {file.name} {sep*(int(starnum/2) - 2)}")
        print(f"{sep*starnum*2}")   
        
        # Always correct specific misnamed file (common to both branches)
        if "DB search psm_ID_CD9_NO_POOL.csv" in file.name:
            new_name = file.with_name("DB search psm_IP_CD9_NO_POOL.csv")
            file.rename(new_name)
            file = new_name
            print(f"\nFile renamed: {file.name}")
        
        # If analysis IS pooled, handle files without "POOL"
        if is_pooled and "POOL" not in file.name:
            print(f"\n *********** WARNING ***********")
            print(f"Files without specified 'POOL' in their names will be removed from pooled analyses")
            print(f"\nFile without specified pool: {file.name}")
            
            answer = input("\nDo you wish to rename it? (y/n): ").strip().lower()
            affirmative_answers = ['y', 'yes', 'sí', 'si', 's']
            if answer in affirmative_answers:
                addendum = input("\nPlease specify pool identifier (e.g., 'pool 1' / 'all pools'): ").upper()
                addendum = re.sub(r"\s+", "_", addendum)
                
                # Create new name by inserting addendum before file extension
                stem = file.stem  # filename without extension
                suffix = file.suffix  # file extension
                new_name = file.with_name(f"{stem}_{addendum}{suffix}")
                file.rename(new_name)
                file = new_name  # Update file reference
                print(f"\nFile renamed to: {file.name}")
            else:
                file.unlink()
                print(f"\nFile {file.name} deleted")
            

def simplify_filenames(directory, garbage):
    """
    Remove extraneous prefixes from mass spectrometry filenames.
    
    Args:
        directory: Path to directory containing files to rename
        techniques: List of technique identifiers used in filenames
        
    Removes elements in garbage list and replaces spaces with underscores
    in filenames to create cleaner, more consistent file naming.
    """
    files = [file for file in Path(directory).rglob('*')]
    counter = 0
    
    for file in files:
        if any(prefix in file.name for prefix in garbage): 
            counter += 1
            new_name_str = file.name

            # Now we erase ALL prefixes (any() doesn't save prefix)
            for prefix in garbage:
                new_name_str = re.sub(prefix, "", new_name_str) # Remove prefixes from filename
                
            new_name_str = re.sub(r"\s+", "_", new_name_str) # spaces -> underscores
            new_name_str = re.sub(r"_+", "_", new_name_str) # multiple underscores -> 1 underscore
            new_name_str = new_name_str.strip("_") # erase underscores at end/beginning
            
            new_name = file.with_name(new_name_str)
            file.rename(new_name)
            file = new_name
            print(f"\nRenamed: {file.name}")
    
    print(f"\nRenamed {counter} files")
            
            
            
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


def extract_pool_number(file_path):
    """Extract pool number from file name"""
    match = re.search(r'POOL_(\d+)', file_path.name)
    return int(match.group(1)) if match else 0

    
def combine_pool_files(directory, techniques, has_techniques):
    """
    Combine individual pool files (1, 2, 3) into a single combined file.
    
    Args:
        directory: Path to directory containing pool files
        techniques: List of technique identifiers to process
        has_techniques: Boolean indicating if this analisys is for comparing techniques
        
    For each technique, finds all pool files and combines them into a single CSV
    with the naming pattern "{technique}_POOLS_{pool_len_elements_concatenated}.csv"
    """
    
    # has_techniques == True
    ########################
    if has_techniques:
        for technique in techniques:
            print(f"\n{technique}")
            # Create search pattern with regex to find pool files
            pattern = re.compile(f"{technique}.*POOL(?!S).*\d\.csv")
            # Find files matching the regex pattern
            files = [
                file for file in Path(directory).rglob('*') 
                if pattern.search(file.name)
            ]
        
            # Sort files per pool number
            files = sorted(files, key=extract_pool_number)
        
            # Process each file for the current technique
            for i, file in enumerate(files, 1):
                df = pd.read_csv(file, sep=";")
                # Extract pool number from file name
                pool_number = extract_pool_number(file)
                print(f"\nDimensions {technique} pool {pool_number} df: {df.shape}")
        
            # Combine all pool files into one DataFrame
            combined_pool_df = pd.concat(
                [pd.read_csv(file, sep=";") for file in files], 
                ignore_index=True
            )
        
            # Save the resulting DataFrame
            pool_output_path = directory + f"/{technique}_POOLS_123.csv"
            combined_pool_df.to_csv(pool_output_path, index=False, sep=";")
            print(f"Dimensions {technique} combined pools: {combined_pool_df.shape}")


    # has_techniques == False
    #########################
    else:
        # Create search pattern with regex to find pool files
        pattern = re.compile(r".*POOL(?!S).*\d\.csv")  # Corregido: añadido 'r' para raw string
        # Find files matching the regex pattern
        files = [
            file for file in Path(directory).rglob('*') 
            if pattern.search(file.name)
        ]

        # Sort files per pool number
        files = sorted(files, key=extract_pool_number)
        
        # Process each file
        for i, file in enumerate(files, 1):
            df = pd.read_csv(file, sep=";")
            pool_number = extract_pool_number(file)
            print(f"\nDimensions pool {pool_number} df: {df.shape}")
            
        # Combine all pool files into one DataFrame
        combined_pool_df = pd.concat(
        [pd.read_csv(file, sep=";") for file in files], 
        ignore_index=True
        )
        
        # Save the resulting DataFrame
        file_count = len(files)
        file_count_name = ''.join(str(x) for x in range(1, file_count+1))
        pool_output_path = directory + f"/ALL_POOLS.csv"
        combined_pool_df.to_csv(pool_output_path, index=False, sep=";")
        print(f"Dimensions combined pools: {combined_pool_df.shape}")
        

def separate_pool_files(directory):
    """
    Separate combined pool files into individual pool files.
    
    Args:
        directory: Path to directory containing combined pool files
        
    For each file that contains pool data, separates them by pool
    based on the TARGET column values.
    """
    # Look for files that might be combined (either explicitly marked as ALL_POOLS or potential combined files)
    all_files = list(Path(directory).rglob('*.csv'))
    
    # First priority: files explicitly marked as ALL_POOLS
    all_pools_files = [f for f in all_files if 'ALL_POOLS' in f.name]
    
    # Second priority: files that don't have individual pool indicators
    potential_combined_files = [
        f for f in all_files 
        if not re.search(r'POOL_\d', f.name) and 'ALL_POOLS' not in f.name
    ]
    
    # Combine both lists, prioritizing ALL_POOLS files
    files_to_process = all_pools_files + potential_combined_files
    
    if not files_to_process:
        print("No files found to separate.")
        return
        
    for file in files_to_process:
        print(f"\nProcessing: {file.name}")
        
        try:
            df = pd.read_csv(file, sep=";")
                
            # Check if TARGET column exists
            TARGET = 'Source File Mod'
            if TARGET not in df.columns:
                print(f"Warning: {TARGET} column not found in {file.name}")
                continue
                
            # Extract pool numbers from TARGET column
            df['Pool_Number'] = df['Source File Mod'].str.extract(r'Pool_(\d+)\.mgf')
                
            # Check if we found any pool numbers
            pool_numbers = df['Pool_Number'].dropna().unique().tolist()
            if len(pool_numbers) == 0:
                print(f"No pool information found in {file.name}")
                continue
                
            # Ensure consistent type for pool numbers (convert to strings)
            pool_numbers = sorted(pool_numbers, key=lambda x: int(x) if x.isdigit() else x)
            print(f"Found pools: {pool_numbers}")
                
            # Separate by pools
            for pool_num in pool_numbers:
                pool_data = df[df['Pool_Number'] == pool_num].copy()
                pool_data = pool_data.drop('Pool_Number', axis=1)  # Remove helper column
                    
                # Generate output filename
                file_stem = file.stem
                # Remove _ALL_POOLS if it exists to get clean base name
                clean_stem = file_stem.replace('_ALL_POOLS', '')
                output_file = Path(directory) / f"{clean_stem}_POOL_{pool_num}.csv"
                    
                # Save individual pool file
                pool_data.to_csv(output_file, index=False, sep=";")
                print(f"Created: {output_file.name}")
                print(f"Dimensions pool {pool_num}: {pool_data.shape}")
                
            # Rename original to indicate it contains all pools (only if not already marked)
            if 'ALL_POOLS' not in file.name:
                all_pools_name = file.with_name(f"{file.stem}_ALL_POOLS.csv")
                # Check if destination exists before renaming
                if all_pools_name.exists():
                    print(f"Warning: {all_pools_name.name} already exists, skipping rename")
                else:
                    try:
                        file.rename(all_pools_name)
                        print(f"Original file renamed to: {all_pools_name.name}")
                    except Exception as e:
                        print(f"Error renaming file: {e}")
            else:
                print(f"File already marked as ALL_POOLS: {file.name}")
                
        except Exception as e:
            print(f"Error processing {file.name}: {e}")

def handle_pool_files(directory, techniques, is_pooled, has_techniques):
    """
    Handle pool files: either combine individual pools or separate combined pools.
    
    Args:
        directory: Path to directory containing files
        techniques: List of technique identifiers used in filenames
        is_pooled: Boolean indicating if this is a pooled analysis
        has_techniques: Boolean indicating if this analysis is for comparing techniques
    """
    if not is_pooled:
        print("\n=== NON-POOLED ANALYSIS ===")
        print("Pool handling skipped for non-pooled analysis.")
        return
        
    print("\n=== POOLED ANALYSIS DETECTED ===")
    
    # Ensure directory is a Path object for consistency
    directory_path = Path(directory) if not isinstance(directory, Path) else directory
    
    files = list(directory_path.glob('*.csv'))
    
    # Check for individual pool files (POOL_1, POOL_2, etc.)
    individual_pool_files = [f for f in files if re.search(r'POOL_\d', f.name)]
    
    # Check for combined files - either ALL_POOLS or files without pool indicators
    combined_files = [f for f in files if 'ALL_POOLS' in f.name or 
                     (not re.search(r'POOL_\d', f.name) and f.name.endswith('.csv'))]
    
    has_individual_pools = len(individual_pool_files) > 0
    has_combined_files = len(combined_files) > 0
    
    if has_individual_pools:
        print(f"Found individual pool files: {len(individual_pool_files)} files")
    if has_combined_files:
        print(f"Found potential combined files: {len(combined_files)} files")
    
    # Decide what action to take
    if has_individual_pools and not has_combined_files:
        print("\n→ Action: Combining individual pool files...")
        try:
            # Convert Path to string for backward compatibility
            combine_pool_files(str(directory_path), techniques, has_techniques)
        except Exception as e:
            print(f"Error combining pool files: {e}")
            
    elif has_combined_files and not has_individual_pools:
        print("\n→ Action: Separating combined files into individual pools...")
        try:
            # separate_pool_files already accepts Path objects
            separate_pool_files(directory_path)
        except Exception as e:
            print(f"Error separating pool files: {e}")
            
    elif has_individual_pools and has_combined_files:
        print("\n⚠️  Warning: Found both individual and combined files!")
        action = input("What would you like to do? (combine/separate/auto): ").strip().lower()
        
        try:
            if action.startswith('c'):
                combine_pool_files(str(directory_path), techniques, has_techniques)
            elif action.startswith('s'):
                separate_pool_files(directory_path)
            elif action.startswith('a'):
                # Auto mode: choose based on counts (go with the operation that would process more files)
                if len(individual_pool_files) >= len(combined_files):
                    print("Auto-selected: combining files (more individual files found)")
                    combine_pool_files(str(directory_path), techniques, has_techniques)
                else:
                    print("Auto-selected: separating files (more combined files found)")
                    separate_pool_files(directory_path)
            else:
                print("No action taken.")
        except Exception as e:
            print(f"Error during pool file handling: {e}")
    else:
        print("\n⚠️  No pool files found.")


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
    Extract Peptide_Sequences from a target column using regex pattern.
    
    Args:
        input_dir: Directory containing CSV files to process
        target_col: Column containing peptide information to extract from
        new_col: Name of the new column to store extracted Peptide_Sequences
        pattern: Regular expression pattern to remove from Peptide_Sequences
        
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
            
            # Create clean_target_col copy to remove garbage from it
            clean_target_col = "Clean_PTM"
            df[clean_target_col] = df[target_col].str.replace(r'[A-Z]\d{1,2}:', '', regex=True)
            df[clean_target_col] = df[clean_target_col].str.replace(r'(\([\D]{2,}\))?:(true|false)', '', regex=True)

            # Explode clean_target_col since it may contain multiple values/row
            #df[clean_target_col] = df[clean_target_col].str.split(';') # split turns row into a list with specified sep
            #df[clean_target_col] = df[clean_target_col].explode(clean_target_col) # explode turns list-like rows into multiple rows
            #df[clean_target_col] = df[clean_target_col].str.replace('^\s', "", regex=True)# clean initial spaces after removing separator (";")
            
            # Define conditions for PTM classification
            conditions = [
                # Fucosialylated: contains fucose/dHex markers AND sialic acid markers
                df[clean_target_col].str.contains(r'dHex|Fucos|Biantennary', regex=True) & 
                df[clean_target_col].str.contains(r'NeuGc|NeuAc|Kdn|Neuraminic', regex=True),
                
                # Fucosylated: contains fucose/dHex markers BUT NO sialic acid markers
                df[clean_target_col].str.contains(r'dHex|Fucos|Biantennary', regex=True) & 
                ~df[clean_target_col].str.contains(r'NeuGc|NeuAc|Kdn|Neuraminic', regex=True),
                
                # Sialylated: contains sialic acid markers BUT NO fucose/dHex markers
                ~df[clean_target_col].str.contains(r'dHex|Fucos|Biantennary', regex=True) & 
                df[clean_target_col].str.contains(r'NeuGc|NeuAc|Kdn|Neuraminic', regex=True),
                
                # Oligomannose: contains HexNAc markers but no fucose or sialic acid markers
                df[clean_target_col].str.contains(r'HexNAc|N\-linked\sglycan\score', regex=True) & 
                ~df[clean_target_col].str.contains(r'dHex|Fucos|NeuGc|NeuAc|Kdn|Neuraminic|Biantennary', regex=True),
                
                # No PTM: explicitly "nan" values
                df[clean_target_col] == "nan"
            ]
            
            # Define values corresponding to conditions (in order)
            choices = ['Fucosialylated', 'Fucosylated', 'Sialylated', 'Oligomannose', "No PTM"]
            
            # Apply conditions to create the new column
            df[new_col] = np.select(conditions, choices, default='Other')
            
            # Handle empty cells in clean_target_col, marking them as "No PTM"
            empty_row_catcher = (
                df[clean_target_col].isna() | 
                (df[clean_target_col] == '') | 
                (df[clean_target_col] == " ") | 
                (df[clean_target_col] == "nan")
            )
            df.loc[empty_row_catcher, new_col] = 'No PTM'
        else:
            df[new_col] = "No PTM"
        
        # Save the modified DataFrame and overwrite original file
        df.to_csv(csv_file, index=False, sep=";")
        print(f"{csv_file.name} overwritten with {new_col}")

        # Check results
        print(f"\n{df[new_col].value_counts(dropna=False).head(5)}")

        print(f"\n{df[clean_target_col].value_counts(dropna=False).head(5)}")



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


def filter_vesiclepedia_proteins(input_dir, output_dir, techniques, pools, ptms_of_interest, vesiclepedia, interest_cols, is_pooled, has_techniques, analysis_choice):
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

        # Identify technique and pool based on analysis type
        
        ## has_techniques = True
        if has_techniques:
            # Original logic: extract technique from filename
            technique = next((technique for technique in techniques if re.search(technique, sample_name)), "")

        ## has_techniques = False
        else:
            # Generic case: no technique-specific processing
            technique = "Generic"

        ## is_pooled = True
        if is_pooled:
            if has_techniques:
                # Original logic: extract pool from filename using pools list
                pool = next((pool for pool in pools if re.search(pool, sample_name)), "")
            else:
                # Generic pooled case: try to extract pool info from filename or use generic
                if re.search(r'POOL_(\d+)', sample_name):
                    pool_match = re.search(r'POOL_(\d+)', sample_name)
                    pool = f"POOL_{pool_match.group(1)}"
                elif 'ALL_POOLS' in sample_name:
                    pool = "ALL_POOLS"
                else:
                    pool = "Unknown_Pool"

        ## is_pooled = False
        else:
            pool = "NO_POOL"
    
        # Filter proteins by presence in Vesiclepedia
        peptides_count = df["Peptide_Sequence"].value_counts()
        print(f"Total peptides: {len(peptides_count)}")
        
        condition = df["Prot_Name"].isin(vesiclepedia["Accession"])
        df_filtered = df.loc[condition, interest_cols]
        
        # Save filtered data
        output_file = output_dir / f"{sample_name}_filtered_vcp.csv"
        df_filtered.to_csv(output_file, index=False, sep=";")
        
        filtered_peptides_count = df_filtered["Peptide_Sequence"].value_counts()
        print(f"Peptides in Vesiclepedia: {len(filtered_peptides_count)}")
    
        # Calculate statistics
        n_prots = len(df["Prot_Name"].value_counts())
        n_ptms = len(df["Clean_PTM"].value_counts())
        n_peptides = len(df["Peptide_Sequence"].value_counts())
       
        n_prots_filt = len(df_filtered["Prot_Name"].value_counts())
        n_ptms_filt = len(df_filtered["Clean_PTM"].value_counts())
        n_peptides_filt = len(df_filtered["Peptide_Sequence"].value_counts())
        
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
    # If working w/ the default set of techniques, move "Whole_milk" to the first position
    if analysis_choice.startswith("1"):
        factor_tech_vals = list(summary_df['Technique'].unique())  # Convert to a list first
        factor_tech_vals = [factor_tech_vals[-1]] + factor_tech_vals[:-1]  # Move last to first position
        summary_df['Technique'] = pd.Categorical(summary_df['Technique'], categories=factor_tech_vals, ordered=True)
        summary_df = summary_df.sort_values(by=['Technique', "Pool"])
    
    # Export summary results in different configurations
    # Always export all samples
    summary_file_path = output_dir / 'summary_all.csv'
    summary_df.to_csv(summary_file_path, index=False, sep=";")
    
    if is_pooled and has_techniques:
        # Original logic for pooled analysis with techniques
        
        # POOLS_123 only
        sf = summary_df.loc[summary_df["Pool"]=="POOLS_123"]
        if not sf.empty:
            filepath = output_dir / 'summary_pools_123.csv'
            sf.to_csv(filepath, index=False, sep=";")

        # POOLS_123 & NO_POOL
        sf = summary_df.loc[(summary_df["Pool"]=="POOLS_123") | (summary_df["Pool"]=="NO_POOL")]
        if not sf.empty:
            filepath = output_dir / 'summary_123nopool.csv'
            sf.to_csv(filepath, index=False, sep=";")

        # Individual pools (POOL 1, POOL 2, POOL 3)
        sf = summary_df.loc[(summary_df["Pool"]!="POOLS_123") & (summary_df["Pool"]!="NO_POOL")]
        if not sf.empty:
            filepath = output_dir / 'summary_individual_pools.csv'
            sf.to_csv(filepath, index=False, sep=";")
    
    elif is_pooled and not has_techniques:
        # Generic pooled analysis - group by pool types
        
        # All pools combined
        sf = summary_df.loc[summary_df["Pool"]=="ALL_POOLS"]
        if not sf.empty:
            filepath = output_dir / 'summary_all_pools.csv'
            sf.to_csv(filepath, index=False, sep=";")
        
        # Individual pools (POOL_1, POOL_2, etc.)
        sf = summary_df.loc[summary_df["Pool"].str.contains("POOL_", na=False)]
        if not sf.empty:
            filepath = output_dir / 'summary_individual_pools.csv'
            sf.to_csv(filepath, index=False, sep=";")
    
    else:
        # Non-pooled analysis - group by technique if available
        if has_techniques:
            # Export by technique
            for technique in summary_df["Technique"].unique():
                sf = summary_df.loc[summary_df["Technique"]==technique]
                if not sf.empty:
                    filepath = output_dir / f'summary_{technique}.csv'
                    sf.to_csv(filepath, index=False, sep=";")
        else:
            # Generic non-pooled - just the main summary
            print("Generic non-pooled analysis - only summary_all.csv exported")
    
    return summary_df

            

def filter_glycosylated_proteins(input_dir, output_dir, techniques, pools, ptms_of_interest, 
                                interest_cols, is_pooled, has_techniques, analysis_choice):
    """
    Filter protein data to keep only those with glycosylation PTMs of interest.
    
    Args:
        input_dir: Directory containing CSV files to process
        output_dir: Directory where filtered files will be saved
        techniques: List of isolation techniques in the filenames
        pools: List of pool identifiers in the filenames
        ptms_of_interest: List of glycosylation PTMs to filter for
        interest_cols: Columns to include in filtered output
        is_pooled: Boolean indicating if this is a pooled analysis
        has_techniques: Boolean indicating if files follow technique naming patterns
    
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

            # Identify technique and pool based on analysis type
            if has_techniques:
                # Original logic: extract technique from filename
                technique = next((technique for technique in techniques if re.search(technique, sample_name)), "")
            else:
                # Generic case: no technique-specific processing
                technique = "Generic"
            
            if is_pooled:
                if has_techniques:
                    # Original logic: extract pool from filename using pools list
                    pool = next((pool for pool in pools if re.search(pool, sample_name)), "")
                else:
                    # Generic pooled case: try to extract pool info from filename
                    if re.search(r'POOL_(\d+)', sample_name):
                        pool_match = re.search(r'POOL_(\d+)', sample_name)
                        pool = f"POOL_{pool_match.group(1)}"
                    elif 'ALL_POOLS' in sample_name:
                        pool = "ALL_POOLS"
                    else:
                        pool = "Unknown_Pool"
            else:
                # Non-pooled analysis
                pool = "NO_POOL"
            
            # Filter proteins by glycosylation PTMs
            peptides_count = df["Peptide_Sequence"].value_counts()
            print(f"Total peptides: {len(peptides_count)}")
            
            condition = df["Clean_PTM"].isin(ptms_of_interest)
            df_filtered = df.loc[condition, interest_cols]
            
            # Save filtered data
            output_file = output_dir / f"{sample_name}_filtered_glyc.csv"
            df_filtered.to_csv(output_file, index=False, sep=";")
            
            filtered_peptides_count = df_filtered["Peptide_Sequence"].value_counts()
            print(f"Glycosylated peptides: {len(filtered_peptides_count)}")
    
            # Calculate statistics
            n_prots = len(df["Prot_Name"].value_counts())
            n_ptms = len(df["Clean_PTM"].value_counts())
            n_peptides = len(df["Peptide_Sequence"].value_counts())
       
            n_prots_filt = len(df_filtered["Prot_Name"].value_counts())
            n_ptms_filt = len(df_filtered["Clean_PTM"].value_counts())
            n_peptides_filt = len(df_filtered["Peptide_Sequence"].value_counts())
        
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
    # If working w/ the default set of techniques, move "Whole_milk" to the first position
    if analysis_choice.startswith("1"):
        factor_tech_vals = list(summary_df['Technique'].unique())  # Convert to a list first
        factor_tech_vals = [factor_tech_vals[-1]] + factor_tech_vals[:-1]  # Move last to first position
        summary_df['Technique'] = pd.Categorical(summary_df['Technique'], categories=factor_tech_vals, ordered=True)
        summary_df = summary_df.sort_values(by=['Technique', "Pool"])
    
    # Export summary results - adapt based on analysis type
    # Always export all samples
    summary_file_path = output_dir / 'summary_all.csv'
    summary_df.to_csv(summary_file_path, index=False, sep=";")
    
    if is_pooled and has_techniques:
        # Original logic for pooled analysis with techniques
        
        # POOLS_123 only
        sf = summary_df.loc[summary_df["Pool"]=="POOLS_123"]
        if not sf.empty:
            filepath = output_dir / 'summary_pools_123.csv'
            sf.to_csv(filepath, index=False, sep=";")

        # POOLS_123 & NO_POOL
        sf = summary_df.loc[(summary_df["Pool"]=="POOLS_123") | (summary_df["Pool"]=="NO_POOL")]
        if not sf.empty:
            filepath = output_dir / 'summary_123nopool.csv'
            sf.to_csv(filepath, index=False, sep=";")

        # Individual pools (POOL 1, POOL 2, POOL 3)
        sf = summary_df.loc[(summary_df["Pool"]!="POOLS_123") & (summary_df["Pool"]!="NO_POOL")]
        if not sf.empty:
            filepath = output_dir / 'summary_individual_pools.csv'
            sf.to_csv(filepath, index=False, sep=";")
    
    elif is_pooled and not has_techniques:
        # Generic pooled analysis - group by pool types
        
        # All pools combined
        sf = summary_df.loc[summary_df["Pool"]=="ALL_POOLS"]
        if not sf.empty:
            filepath = output_dir / 'summary_all_pools.csv'
            sf.to_csv(filepath, index=False, sep=";")
        
        # Individual pools (POOL_1, POOL_2, etc.)
        sf = summary_df.loc[summary_df["Pool"].str.contains("POOL_", na=False)]
        if not sf.empty:
            filepath = output_dir / 'summary_individual_pools.csv'
            sf.to_csv(filepath, index=False, sep=";")
    
    else:
        # Non-pooled analysis - group by technique if available
        if has_techniques:
            # Export by technique
            for technique in summary_df["Technique"].unique():
                sf = summary_df.loc[summary_df["Technique"]==technique]
                if not sf.empty:
                    filepath = output_dir / f'summary_{technique}.csv'
                    sf.to_csv(filepath, index=False, sep=";")
        else:
            # Generic non-pooled - just the main summary
            print("\nGeneric non-pooled analysis - only summary_all.csv exported")
    
    return summary_df



######################################################################################################################
# Functions for counting


def count_peptides_by_category(input_dir, output_dir, techniques, pools, interest_cols, 
                              is_pooled, has_techniques):
    """
    Count and analyze peptides by different categories of interest.
    
    Args:
        input_dir: Directory containing CSV files to process
        output_dir: Directory where count files will be saved
        techniques: List of isolation techniques in the filenames
        pools: List of pool identifiers in the filenames
        interest_cols: Columns to perform value counts on
        is_pooled: Boolean indicating if this is a pooled analysis
        has_techniques: Boolean indicating if files follow technique naming patterns
        
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
            
            # Identify technique and pool based on analysis type
            if has_techniques:
                # Original logic: extract technique from filename
                technique = [t for t in techniques if re.search(t, sample_name)]
                if is_pooled:
                    # Extract pool from filename using pools list
                    pool = [p for p in pools if re.search(p, sample_name)]
                else:
                    pool = ["NO_POOL"]
            else:
                # Generic case: no technique-specific processing
                technique = ["Generic"]
                if is_pooled:
                    # Generic pooled case: try to extract pool info from filename
                    if re.search(r'POOL_(\d+)', sample_name):
                        pool_match = re.search(r'POOL_(\d+)', sample_name)
                        pool = [f"POOL_{pool_match.group(1)}"]
                    elif 'ALL_POOLS' in sample_name:
                        pool = ["ALL_POOLS"]
                    else:
                        pool = ["Unknown_Pool"]
                else:
                    pool = ["NO_POOL"]
            
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
                    
                    # Create output filename that reflects the analysis type
                    if has_techniques and technique:
                        tech_str = '_'.join(technique)
                        if is_pooled and pool:
                            pool_str = '_'.join(pool)
                            output_file = output_dir / f"{sample_name}_{col}.csv"
                        else:
                            output_file = output_dir / f"{sample_name}_{col}.csv"
                    else:
                        if is_pooled and pool:
                            pool_str = '_'.join(pool)
                            output_file = output_dir / f"{sample_name}_{col}.csv"
                        else:
                            output_file = output_dir / f"{sample_name}_{col}.csv"
                    
                    # Export the counts to CSV
                    summary.to_csv(output_file, index=False, sep=";")
                else:
                    print(f"Warning: Column '{col}' not found in file {sample_name}.")


def count_proteins_by_ptm(input_dir, output_dir, techniques, pools, interest_col, 
                         is_pooled, has_techniques, group_by_col="Prot_Name"):
    """
    Count proteins by PTM types, grouping at the protein level.
    
    Args:
        input_dir: Directory containing CSV files to process
        output_dir: Directory where count files will be saved
        techniques: List of isolation techniques in the filenames
        pools: List of pool identifiers in the filenames
        interest_col: Column to count values from (typically "PTM_Cluster")
        is_pooled: Boolean indicating if this is a pooled analysis
        has_techniques: Boolean indicating if files follow technique naming patterns
        group_by_col: Column to group by (default "Prot_Name")
        
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
            
            # Identify technique and pool based on analysis type
            if has_techniques:
                # Original logic: extract technique from filename
                technique = [t for t in techniques if re.search(t, sample_name)]
                if is_pooled:
                    # Extract pool from filename using pools list
                    pool = [p for p in pools if re.search(p, sample_name)]
                else:
                    pool = ["NO_POOL"]
            else:
                # Generic case: no technique-specific processing
                technique = ["Generic"]
                if is_pooled:
                    # Generic pooled case: try to extract pool info from filename
                    if re.search(r'POOL_(\d+)', sample_name):
                        pool_match = re.search(r'POOL_(\d+)', sample_name)
                        pool = [f"POOL_{pool_match.group(1)}"]
                    elif 'ALL_POOLS' in sample_name:
                        pool = ["ALL_POOLS"]
                    else:
                        pool = ["Unknown_Pool"]
                else:
                    pool = ["NO_POOL"]
            
            # Group by protein name and count PTM types
            if group_by_col in df.columns and group_by_col == "Prot_Name":
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
                
                # Create output filename that reflects the analysis type
                if has_techniques and technique:
                    tech_str = '_'.join(technique)
                    if is_pooled and pool:
                        pool_str = '_'.join(pool)
                        output_file = output_dir / f"{sample_name}_PTM_types_by_protein.csv"
                    else:
                        output_file = output_dir / f"{sample_name}_PTM_types_by_protein.csv"
                else:
                    if is_pooled and pool:
                        pool_str = '_'.join(pool)
                        output_file = output_dir / f"{sample_name}_PTM_types_by_protein.csv"
                    else:
                        output_file = output_dir / f"{sample_name}_PTM_types_by_protein.csv"
                
                # Export the counts to CSV
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