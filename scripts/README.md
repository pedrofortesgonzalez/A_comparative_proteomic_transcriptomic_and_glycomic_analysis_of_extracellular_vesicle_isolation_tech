# Scripts
This folder contains all scripts needed to perform the analysis of glycosylations in EV proteins and peptides. The scripts are organized in numerical order according to the analysis workflow.

## Folder structure
```{bash}
scripts/
├── 1_ptm_detection.py                 # Initial processing of mass spectrometry data
├── 2_hyptest_&_boxplots.R             # Statistical analysis and boxplot generation
├── 3_sector_diagrams.R                # Creation of pie charts
├── user_defined_funcs.py              # Helper functions for data processing with Python
├── cleaning_1_python.py               # Cleanup of unused Python libraries
├── cleaning_2_R.R                     # Cleanup of unused R libraries
└── __pycache__/                       # Python compiled files (automatically generated)

```
## Workflow
The complete analysis follows this workflow:

- Data Processing: Run <ins>1_ptm_detection.py</ins> to process raw data
  - The console will ask you to provide the name of the folder where your data is located. This folder should be located in the main folder
- Boxplots and Hypothesis Testing: Run <ins>2_hyptest_&_boxplots.R</ins> to generate statistical analyses and boxplots
- Glycosylation Proportions Visualization: Run <ins>3_sector_diagrams.R</ins> to create pie charts

### Usage Notes
Execute the scripts in numerical order to ensure correct workflow.
The cleaning scripts (<ins>cleaning_1_python.py</ins> and <ins>cleaning_2_R.R</ins>) are optional and can be used to optimize code.
Ensure all dependencies are installed before running the scripts.

### Requirements and Dependencies

| **Python** 3.10.16 | **R** 4.4.2 |
| :----: | :----: |
| pandas == 2.2.3 | dplyr_1.1.4 |
| session_info == 1.0.0 | ggplot2_3.5.1 |
| tqdm == 4.67.1 | ggpubr_0.6.0 |
| numpy == 1.26.4 | rstatix_0.7.2 |
|  | envalysis_0.7.0 |
|  | gridExtra_2.3 |
|  | ggrepel_0.9.6 |
|  | rstudioapi_0.17.1 |

### Code Cleanup Libraries

| **Python** | **R** |
| :----: | :----: |
| vulture == 1.0 | nolock_1.1.0 |


***
## Main Scripts
### <ins>1_ptm_detection.py</ins>
This script performs the initial processing of mass spectrometry data. Its main functions are:

- Creating the folder structure in output/
- Processing input CSV files
- Filtering proteins and peptides present in Vesiclepedia (EV-related)
- Identifying proteins and peptides with glycosylations of interest
- Generating count summaries for subsequent analyses

**Input files**:

- CSV files in data/input/
- data/glycosylation_list.csv
- data/vesiclepedia_proteins_240712.csv

**Output files**: Multiple CSV files in

- output/1_filtered_dfs/
- output/2_value_counts/

***
### <ins>2_hyptest_&_boxplots.R</ins>
This script performs statistical analyses on the processed data and generates boxplot visualizations. Its main functions are:

- Performing Kruskal-Wallis tests to compare isolation techniques
- Running Dunn's post-hoc tests for pairwise comparisons
- Generating boxplots that show:
  - Total peptide and protein counts
  - Counts of peptides and proteins present in Vesiclepedia
  - Counts of glycosylated peptides and proteins present in Vesiclepedia

**Input files**:
Processed CSV files in output/processed_data/

**Output files**:
Boxplots in output/figures/boxplots/

***
### <ins>3_sector_diagrams.R</ins>
This script creates pie charts to visualize the distribution of glycosylation types in proteins and peptides. The main functions are:

- Processing count data for different glycosylation types
- Generating pie charts by isolation technique
- Creating versions with and without text labels for each chart

**Input files**:
Count CSV files in output/processed_data/value_counts/

**Output files**:
Pie charts in output/figures/sector_diagrams/



***
# Additional Scripts
## <ins>user_defined_funcs.py</ins>
Module containing custom functions used by <ins>1_ptm_detection.py</ins>. It separates the main logic from the script for better code organization and maintenance.

## <ins>cleaning_1_python.py</ins>
Utility script that uses the vulture module to analyze <ins>1_ptm_detection.py</ins> and identify imported libraries that are not used, helping to keep the code clean and efficient.

## <ins>cleaning_2_R.R</ins>
Utility script that uses the nolock package to identify imported but unused libraries in R scripts, helping to optimize the code.
