# Please go through the README files for a comprehensive description of the features of this repository.

## Glycosylation Analysis in Extracellular Vesicles
This repository contains the code and data used for the detection and analysis of glycosylations in extracellular vesiclen (EV) proteins and peptides, as supplementary material for the publication:

[Publication Title] (Authors, Journal, Year)
--doi--

# Workflow
*The complete analysis follows this workflow:*

- *Data Processing: Run <ins>1_ptm_detection.py</ins> to process raw data*
  - *The console will ask you to provide the name of the folder where your data is located. This folder should be located in the main folder*
- *Boxplots and Hypothesis Testing: Run <ins>2_hyptest_&_boxplots.R</ins> to generate statistical analyses and boxplots*
- *Glycosylation Proportions Visualization: Run <ins>3_sector_diagrams.R</ins> to create pie charts*

## Workflow and Usage Notes

1. Drag your mass spectrometry data to the main folder. Make sure:
    - your files are in .csv format, and contain columns for **rellenar**
    - the folder containing your data has no spaces in its name
2. Execute the scripts inside the <ins>scripts/</ins> folder. You can do it:
    - Manually: opening up each script with an IDE or **whatever**. Execute the scripts in numerical order to ensure correct workflow.
    - Automatically: if you are familiar with bash command line, then execute
3. The cleaning scripts (<ins>cleaning_1_python.py</ins> and <ins>cleaning_2_R.R</ins>) are optional and can be used to optimize code.
4. Ensure all dependencies are installed before running the scripts.

## Description
This project analyzes the performance of different EV isolation techniques (ExoGAG, IP CD9, SEC, UC) for detecting glycosylated proteins and peptides. The analysis focuses on identifying:

- The total number of proteins and peptides detected by each technique
- Proteins and peptides present in the Vesiclepedia database
- Proteins and peptides with specific glycosylations

The results include statistical analyses (Kruskal-Wallis and Dunn's tests) and visualizations (boxplots and pie charts) that allow comparing the performance of the different techniques.
