# set WEnv
rm(list=ls())                                    # clean WEnv
library(rstudioapi)                              # configurar wdirectory automáticamente a la carpeta raíz del proyecto
setwd(dirname(getActiveDocumentContext()$path))
library(nolock)                                  # para gestionar los directorios


# detect unused packages in script n. 2
script_path <- "./2_hyptest_&_boxplots.R"        # set script path
libr_unused(script = script_path)                # detect unused packages
libr_used(script = script_path)                  # detect used packages


# detect unused packages in script n. 3
script_path <- "./3_sector_diagrams.R"        # set script path
libr_unused(script = script_path)                # detect unused packages
libr_used(script = script_path)                  # detect used packages
