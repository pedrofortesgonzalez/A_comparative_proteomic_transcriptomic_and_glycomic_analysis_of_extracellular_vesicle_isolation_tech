# Input Data

This folder contains all input files required to run the analysis of glycosylations in EV proteins and peptides.

## Folder Structure

```
input_data/
├── glycosylation_list.csv            # List of glycosylations of interest
├── vesiclepedia_proteins_240712.csv  # EV proteins (updated 07/12/2024)
└── input/                            # Mass spectrometry data for processing
```


## Important Notes

1. Do not modify the original files in this folder, as the scripts use these files as a baseline reference.
2. If you need to use an updated version of Vesiclepedia, make sure to update the filename in the corresponding scripts.
3. The glycosylation list can be expanded following the same format if additional types need to be searched.
   
# The **mass spectrometry input files must have the following format**:
   
    - The first column must contain protein/peptide identifiers
    - Must include columns for peptide sequence
    - Must include information about post-translational modifications (PTMs)



## File Descriptions

### <ins>glycosylation_list.csv</ins>
This file contains all chemical nomenclatures of glycosylations of interest that we aim to detect in the mass spectrometry data. It includes different types such as:

- Fucosylations
- Fucosialylations
- Sialylations
- Oligomannoses
- Other glycosylations



### <ins>vesiclepedia_proteins_240712.csv</ins>
List of proteins associated with EVs according to the Vesiclepedia database, updated as of July 12, 2024. This file is used to filter the proteins detected in our experiments and determine which ones are associated with EVs according to the scientific literature.


### <ins>input/</ins> Folder
Contains CSV files with mass spectrometry data to be processed. Each file represents a different experiment with a specific EV isolation technique, including:

- ExoGAG
- IP_CD9
- SEC
- UC
