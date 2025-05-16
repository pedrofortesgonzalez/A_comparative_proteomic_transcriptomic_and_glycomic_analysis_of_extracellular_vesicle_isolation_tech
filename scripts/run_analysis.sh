#!/bin/bash
# Script for executing sequentially the numbered scripts

####################
# Pre-execution code
####################

# Obtain the dir. where this bash script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo -e "\nThis script is running in the following directory:\n${SCRIPT_DIR}\n"

# Define colors for a better visualisation
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "\n${YELLOW}Starting analysis...${NC}"

# Check Python dependencies
echo -e "\n${YELLOW}Checking Python dependencies...${NC}"
python3 -c "
import importlib.util
import sys

required_packages = ['pandas', 'numpy']  # Añade aquí todas las dependencias
missing_packages = []

for package in required_packages:
    if importlib.util.find_spec(package) is None:
        missing_packages.append(package)

if missing_packages:
    print(f'ERROR: Missing Python packages: {\", \".join(missing_packages)}')
    print('Please install them using:')
    print(f'pip install {\" \".join(missing_packages)}')
    print('or')
    print(f'conda install {\" \".join(missing_packages)}')
    sys.exit(1)
else:
    print('All Python dependencies are installed!')
"

if [ $? -ne 0 ]; then
    echo -e "\n${RED}Please install the missing dependencies and try again.${NC}"
    exit 1
fi

###############################
# 1. Execute script n1 (Python)
###############################
echo -e "\n${YELLOW}Executing pre-processing script in Python...${NC}"
python3 "${SCRIPT_DIR}/1_ptm_detection.py"
if [ $? -ne 0 ]; then
    echo -e "\n${RED}Error in script number 1. Aborting the process.${NC}"
    exit 1
fi
echo -e "\n${GREEN}Pre-processing was completed successfully!.${NC}"

################################
# 2. Execute script number 2 (R)
################################
# Execute script n2 (R)
echo -e "\n${YELLOW}Executing first analysis in R...${NC}"
Rscript "${SCRIPT_DIR}/2_hyptest_n_boxplots.R"
if [ $? -ne 0 ]; then
    echo -e "\n${RED}Error in script number 2. Aborting the process.${NC}"
    exit 1
fi
echo -e "\n${GREEN}First analysis was completed successfully!.${NC}"

################################
# 3. Execute script number 3 (R)
################################
echo -e "\n${YELLOW}Executing second analysis en R...${NC}"
Rscript "${SCRIPT_DIR}/3_sector_diagrams.R"
if [ $? -ne 0 ]; then
    echo -e "\n${RED}Error in script number 3. Aborting the process.${NC}"
    exit 1
fi
echo -e "\n${GREEN}Second analysis was completed successfully!.${NC}"


#############################
# Pipeline completion message
#############################
echo "=========================================="
echo "====== Pipeline Execution Completed ======"
echo "=========================================="
echo -e "${GREEN}Every script was run successfully!\nGo to './output' folder to check your results!${NC}"
echo "=========================================="