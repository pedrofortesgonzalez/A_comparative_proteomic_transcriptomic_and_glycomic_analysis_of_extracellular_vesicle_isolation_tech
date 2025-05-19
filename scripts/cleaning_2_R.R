# ******************************************************************************
# Cleaning R scripts
# Author: Pedro Fortes González (refactored)
# ******************************************************************************

# Set working environment
rm(list=ls())                                    # Empty working environment
library(rstudioapi)                              # Set Working Dir. authomatically to repository root folder
setwd(dirname(getActiveDocumentContext()$path))
library(nolock)                                  # This library detects imported but unused packages


# Detect unused packages in script n. 2
script_path <- "./2_hyptest_n_boxplots.R"        # Set script path
libr_unused(script = script_path)                # Detect unused packages
libr_used(script = script_path)                  # Detect used packages


# Detect unused packages in script n. 3
script_path <- "./3_sector_diagrams.R"           # Set script path
libr_unused(script = script_path)                # Detect unused packages
libr_used(script = script_path)                  # Detect used packages
