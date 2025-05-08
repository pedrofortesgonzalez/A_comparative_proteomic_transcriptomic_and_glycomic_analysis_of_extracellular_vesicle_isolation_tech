#!/usr/bin/env python
# coding: utf-8

# Import system libraries
import os
import vulture

# Set working directory
WORKING_DIR = r"/Users/pedrofortesgonzalez/Desktop/NefroCHUS/projects/240614_detectar_prots_glucosiladas_bravo-susana_real/script_for_repo/scripts"
os.chdir(WORKING_DIR)

# Show package list
session_info.show(dependencies=False, html=False, std_lib=False)


# Clean script of unused libraries
!vulture ./1_ptm_detection.py

