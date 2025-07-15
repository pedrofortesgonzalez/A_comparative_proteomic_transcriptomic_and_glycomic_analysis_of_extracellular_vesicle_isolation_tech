***Please go through this README file for a comprehensive description of the features of this repository.***

# Glycosylation Analysis in Extracellular Vesicles

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📋 Overview
This repository contains a computational pipeline for analyzing glycomic profiles of extracellular vesicles (EVs) isolated through four different techniques: ExoGAG, IP-CD9, SEC and UC. The workflow processes mass spectrometry data to identify post-translational modifications (PTMs), with a focus on protein glycosylation patterns. It serves as supplementary material for the publication:
```
Pereira-Hernández *et al*. `A comparative proteomic, transcriptomic and glycomic analysis of extracellular vesicle isolation techniques highlights ExoGAG efficiency for a more complete identification of breast milk molecular signaling pathways`
```

## 📁 Repository Structure
```
.
├── input_data/                     # Input data directory
│   ├── glycosylation_list.csv      # List of glycosylation PTMs of interest
│   └── Vesiclepedia_proteins_*.csv # Reference database of EV proteins
├── output/                         # Generated during analysis (not in repo)
│   ├── 1_filtered_dfs/             # Filtered datasets
│   ├── 2_value_counts/             # Count summaries by category
│   └── 3_figures/                  # Generated visualizations
└── scripts/
    ├── 1_ptm_detection.py          # Python script for data preprocessing
    ├── 2_hyptest_n_boxplots.R      # R script for statistical testing and boxplots
    ├── 3_sector_diagrams.R         # R script for pie chart generation
    ├── run_analysis.sh             # Main pipeline execution script
    ├── user_defined_funcs.py       # Custom Python functions
    └── analysis_config.txt         # Configuration file for analysis type
```


## How to use

### Getting started
- Python 3.8+
- R 4.0+
- Required Python packages:
  - pandas
  - tqdm
  - session_info
- Required R packages:
  - dplyr
  - ggplot2
  - ggpubr
  - rstatix
  - envalysis
  - gridExtra
  - ggrepel

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/pedrofortesgonzalez/Comparative_glycomic_analysis_of_EV_isolation_techniques.git
   cd Comparative_glycomic_analysis_of_EV_isolation_techniques
   ```

2. Install required dependencies:
   ```bash
   # Python dependencies
   pip install -r requirements.txt
   
   # R dependencies will be installed automatically when running the scripts
   ```

**Alternatively**, download it as a folder and place it in a location you can find on your hard drive (you may need to copy and paste its location in case auto-detection fails).

***
### Usage
1. Place your mass spectrometry data files in a directory inside the repository folder
2. Run the complete analysis pipeline:
   ```bash
   ./scripts/run_analysis.sh /path/to/your/data
   ```

3. Alternatively, run each step individually:
   ```bash
   # Step 1: Data preprocessing
   python scripts/1_ptm_detection.py
   
   # Step 2: Statistical analysis and boxplots
   Rscript scripts/2_hyptest_n_boxplots.R
   
   # Step 3: Generate sector diagrams
   Rscript scripts/3_sector_diagrams.R
   ```

***
### Analysis Workflow
1. **Data Preprocessing** (`1_ptm_detection.py`):
   - Standardize CSV file formats and delimiters
   - Extract protein accessions and clean peptide sequences
   - Classify PTMs into standardized categories
   - Filter proteins by presence in Vesiclepedia database
   - Filter by glycosylation status

2. **Statistical Analysis** (`2_hyptest_n_boxplots.R`):
   - Perform Kruskal-Wallis tests to identify significant differences
   - Apply Dunn's post-hoc tests with FDR correction
   - Generate boxplots and barplots with statistical annotations

3. **Visualization** (`3_sector_diagrams.R`):
   - Create pie charts showing distribution of glycosylation types
   - Generate sector diagrams with and without text labels


## Results
The analysis produces several types of visualizations:

1. **Boxplots** comparing protein/peptide counts across isolation techniques
2. **Sector diagrams** showing the distribution of glycosylation types
3. **Summary statistics** for hypothesis testing results

   
## License
This project is licensed under the MIT License - see the LICENSE file for details.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## Citation
If you use this code in your research, please cite:
```
Pereira-Hernández, M. *et al*. (2025). `A comparative proteomic, transcriptomic and glycomic analysis of extracellular vesicle isolation techniques highlights ExoGAG efficiency for a more complete identification of breast milk molecular signaling pathways`. [Journal information pending]
```


## Contact
- GitHub Profile: [Pedro Fortes González](https://github.com/pedrofortesgonzalez)
- LinkedIN: [Pedro Fortes González](www.linkedin.com/in/pedrofortesgonzalez)
- ORCID: [Pedro Fortes González](https://orcid.org/0009-0008-7016-0292)
