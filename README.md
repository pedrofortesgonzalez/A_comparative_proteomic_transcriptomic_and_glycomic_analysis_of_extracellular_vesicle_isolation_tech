## Please go through the README files for a comprehensive description of the features of this repository.

# Glycosylation Analysis in Extracellular Vesicles
This repository contains the code and data used for the detection and analysis of glycosylations in extracellular vesiclen (EV) proteins and peptides, as supplementary material for the publication:

[Publication Title] (Authors, Journal, Year)
--doi--

# Workflow
The complete workflow and required dependencies for this analysis are described in <ins>./scripts/README.md</ins>.


# Description
This project analyzes the performance of different EV isolation techniques (ExoGAG, IP CD9, SEC, UC) for detecting glycosylated proteins and peptides. The analysis focuses on identifying:

- The total number of proteins and peptides detected by each technique
- Proteins and peptides present in the Vesiclepedia database
- Proteins and peptides with specific glycosylations

The results include statistical analyses (Kruskal-Wallis and *post hoc* Dunn's tests if omnibus test is significant) and visualizations (boxplots and pie charts) that allow comparing the performance of the different techniques.
