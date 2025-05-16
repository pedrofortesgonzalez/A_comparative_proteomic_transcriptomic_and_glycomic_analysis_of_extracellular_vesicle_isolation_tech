# set WEnv
rm(list=ls())                                    # clean WEnv
library(rstudioapi)                              # set Working Dir. authomatically to repository root folder
setwd(dirname(getActiveDocumentContext()$path))
library(nolock)                                  # for managing dirs


# detect unused packages in script n. 2
script_path <- "./2_hyptest_n_boxplots.R"        # set script path
libr_unused(script = script_path)                # detect unused packages
libr_used(script = script_path)                  # detect used packages


# detect unused packages in script n. 3
script_path <- "./3_sector_diagrams.R"           # set script path
libr_unused(script = script_path)                # detect unused packages
libr_used(script = script_path)                  # detect used packages
