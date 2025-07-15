# Scripts
This folder contains all scripts needed to perform the analysis of glycosylations in EV proteins and peptides. The scripts are organized in numerical order according to the analysis workflow.


***
## Workflow and Usage Notes
#### <ins>Step 1</ins>
Drag your mass spectrometry data to the main repository folder. Make sure:

   - your files are in .csv format, and contain the following columns:  **'Peptide', "Accession"** and **"PTM"**
   - the folder containing your data has no spaces in its name

#### <ins>Step 2</ins>
Execute the numbered scripts inside the `scripts/` folder. You can do it:

   - **Automatically (recommended)**: if you are familiar with bash command line, then execute
```{bash}
# Navigate to the repository directory
cd /path/to/repository

# Execute the pipeline script
scripts/run_analysis.sh
```

   - **Manually**: execute each script individually in the correct sequence in order to ensure correct workflow. You can do so either with an IDE (like Spyder, Visual Studio Code, Jupyter Notebook, etc.) or through the following bash command line:
```{bash}
# First, run the Python preprocessing script:
python3 scripts/1_ptm_detection.py

# Next, run the R statistical analysis script:
Rscript scripts/2_hyptest_n_boxplots.R

# Finally, run the visualization script:
Rscript scripts/3_sector_diagrams.R
```

#### Usage notes

    - Note 1: Script 1 (`1_ptm_detection.py`) will prompt you to type the name of your input folder. For the data analysed in the article, you would type mass_spectrometry_data`.
    - Note 2: Using an IDE will give you more control over the steps in each script and allow you to directly inspect code and intermediate results.
    
#### <ins>Step 3</ins>
The cleaning scripts (`cleaning_1_python.py` and `cleaning_2_R.R`) are optional and can be used to optimize code to remove unused libraries.

#### <ins>Step 4</ins>
All results will be stored in the `../output` directory, organized into subdirectories for filtered data, counts, and visualizations. 


***
## Main Scripts
### `run_analysis.sh`
This bash scripts runs scripts number 1, 2 and 3. Previous to that, it detects possible missing dependencies essential for the execution.


### `1_ptm_detection.py`
This script performs the initial processing of mass spectrometry data. Its main functions are:

- Creating the folder structure in `../output/`
- Processing input CSV files
- Filtering proteins and peptides present in Vesiclepedia (EV-related)
- Identifying proteins and peptides with glycosylations of interest
- Generating count summaries for subsequent analyses

**Input files**:
- CSV files in `./input_data/input/` --> a direct temporary copy of your input folder
- `../input_data/glycosylation_list.csv`
- `../input_data/vesiclepedia_proteins_240712.csv`

**Output files**: Multiple CSV files are generated in
- `../output/1_filtered_dfs/`
- `../output/2_value_counts/`

**Usage notes**
- The script will prompt for the folder name containing your data
- *Note: When running this script in an IDE environment (like PyCharm, VSCode, etc.), you may be prompted to provide the absolute path to the repository, as automatic directory detection may not work properly in IDEs*


### `2_hyptest_n_boxplots.R`
This script performs statistical analyses on the processed data and generates boxplot visualizations. Its main functions are:

- Performing Kruskal-Wallis tests to compare isolation techniques
- Running Dunn's post-hoc tests for pairwise comparisons
- Generating boxplots that show:
  - Total peptide and protein counts
  - Counts of peptides and proteins present in Vesiclepedia
  - Counts of glycosylated peptides and proteins present in Vesiclepedia

**Input files**:
Processed CSV files in `../3_figures/filtered_dfs/`

**Output files**:
- Boxplots in `../3_figures/boxplots/`
- Boxplots are created with automatically adjusted y-axes


### `3_sector_diagrams.R`
This script creates pie charts to visualize the distribution of glycosylation types in proteins and peptides. The main functions are:

- Processing count data for different glycosylation types
- Generating pie charts by isolation technique
- Creating versions with and without text labels for each chart

**Input files**:
Count CSV files in `../output/2_value_counts/`

**Output files**:
- Pie charts in `../output/3_figures/sector_diagrams/`
- Generates both plain and annotated versions



***
## Additional Scripts
### `user_defined_funcs.py`
Module containing custom functions used by `1_ptm_detection.py`. It separates the main logic from the script for better code organization and maintenance.

### `cleaning_1_python.py`
Utility script that uses the vulture module to analyze `1_ptm_detection.py` and identify imported libraries that are not used, helping to keep the code clean and efficient.

### `cleaning_2_R.R`
Utility script that uses the nolock package to identify imported but unused libraries in R scripts, helping to optimize the code.
