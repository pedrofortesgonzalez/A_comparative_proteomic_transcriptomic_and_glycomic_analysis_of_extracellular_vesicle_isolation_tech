# Input Data

This folder contains all input files required to run the analysis of glycosylations in extracellular vesicle proteins and peptides.

## Folder Structure

```
input_data/
├── glycosylation_list.csv     # List of glycosylations of interest
├── vesiclepedia_proteins_240712.csv  # Extracellular vesicle proteins (updated 07/12/2024)
└── input/                     # Mass spectrometry data for processing
```

## File Descriptions

### glycosylation_list.csv
This file contains all chemical nomenclatures of glycosylations of interest that we aim to detect in the mass spectrometry data. It includes different types such as:
- Fucosylations
- Fucosialylations
- Sialylations
- Oligomannoses
- Other glycosylations

### vesiclepedia_proteins_240712.csv
List of proteins associated with extracellular vesicles according to the Vesiclepedia database, updated as of July 12, 2024. This file is used to filter the proteins detected in our experiments and determine which ones are associated with extracellular vesicles according to the scientific literature.

### input/ Folder
Contains CSV files with raw mass spectrometry data to be processed by the scripts. Each file represents a different experiment with a specific extracellular vesicle isolation technique, including:
- ExoGAG
- IP_CD9
- SEC (Size Exclusion Chromatography)
- UC (Ultracentrifugation)

The data is organized in multiple sample pools (POOL_1, POOL_2, POOL_3) to allow for comparison and statistical analysis.

## CSV Input File Format

The mass spectrometry input files must have the following format:
- The first column must contain protein/peptide identifiers
- Must include columns for peptide sequence
- Must include information about post-translational modifications (PTMs)

## Important Notes

1. Do not modify the original files in this folder, as the scripts use these files as a baseline reference.
2. If you need to use an updated version of Vesiclepedia, make sure to update the filename in the corresponding scripts.
3. The glycosylation list can be expanded following the same format if additional types need to be searched.
