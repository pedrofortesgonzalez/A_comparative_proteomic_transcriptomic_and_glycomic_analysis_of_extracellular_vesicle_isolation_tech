# Data

This folder contains all input files required to run the analysis of glycosylations in EV proteins and peptides.

## Folder Structure

```
data/
├── glycosylation_list.csv            # List of glycosylations of interest
├── vesiclepedia_proteins_240712.csv  # EV proteins (update from July 12th, 2024)
└── input/                            # Temporary folder with mass spectrometry data for processing
```


## Important Notes

1. Do not modify the original files in this folder, as the scripts use these files as a baseline reference.
2. If you need to use an updated version of Vesiclepedia, make sure to update the filename in the corresponding scripts.
3. The glycosylation list can be expanded following the same format if additional types need to be searched.
   
# The **mass spectrometry input files must contain the following columns at least**:
   
    - A protein/peptide identifier in UniProtKB Accession format (Accession column)
    - A column for peptide sequence (Peptide column)
    - A column for post-translational modifications (PTM column). Examples of these may be located in `glycosylation_list.csv`



## File Descriptions

### `glycosylation_list.csv`
This file contains all chemical nomenclatures of glycosylations of interest that we aim to detect in the mass spectrometry data. It includes different types such as:

- Fucosylations
- Fucosialylations
- Sialylations
- Oligomannoses
- Other glycosylations



### `vesiclepedia_proteins_240712.csv`
List of proteins associated with EVs according to the Vesiclepedia database, updated as of July 12, 2024. This file is used to filter the proteins detected in our experiments and determine which ones are associated with EVs according to the scientific literature.


### `input/` Folder
Contains CSV files with mass spectrometry data to be processed. It is a copy of the input folder located in the main folder that allows the scripts to work without modifying the original data, and is both created and removed by the script `1_ptm_detection.py`.
