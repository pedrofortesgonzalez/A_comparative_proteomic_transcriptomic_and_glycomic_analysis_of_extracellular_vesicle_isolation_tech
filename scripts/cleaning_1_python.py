#!/usr/bin/env python
# coding: utf-8

# Import system libraries
import os
import vulture

# Set working 
# Detecting repo's root dir
try:
    # Try to obtain the actual script's directory (works for executing with bash)
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(SCRIPT_DIR)
    
except Exception as e:
    # If we work from an IDE/notebook where __file__ is not defined
    print(f"\nPath could not be authomatically detected:\n{str(e)}")
    REPO_ROOT = input("\nPlease, re-introduce the absolute path to the repository in your computer: ")
    REPO_ROOT = REPO_ROOT.strip()
    SCRIPT_DIR = REPO_ROOT + "/scripts"
    print(f"\nUsing the given path:\n{REPO_ROOT}\n")
    
    # Validate that the path exists
    if not os.path.exists(REPO_ROOT):
        print(f"\nERROR: The given path was not found:\n{REPO_ROOT}\n")
        sys.exit(1)



# Show package list
session_info.show(dependencies=False, html=False, std_lib=False)


# Clean script of unused libraries
!vulture ./1_ptm_detection.py

