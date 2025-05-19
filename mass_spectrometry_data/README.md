# Input Data
This folder contains the mass spectrometry data files for analysis. 

## Requirements
- All files must be in CSV format
- Each file must contain at least the following columns:
  - `Accession` - Protein/peptide identifier in UniProtKB Accession format
  - `Peptide` - Peptide sequence
  - `PTM` - Post-translational modifications

## Important Notes
- This folder is will have a temporary copy created during processing (`../input_data/input`)
- The script `1_ptm_detection.py` will ask for the name of your input folder (in this case, it is the name of this folder, `mass_spectrometry_data`)
- Files will be processed in their original format - no renaming is necessary
- Make sure your input folder name does not contain spaces

## Example
A properly formatted input CSV file should have headers and data similar to:

```csv
Accession,Peptide,PTM,OtherColumn1,OtherColumn2
P12345,PEPTIDESEQUENCE,HexNAc(1)Hex(1)Fuc(1),value1,value2
Q67890,ANOTHERPEPTIDE,HexNAc(2)Hex(2)NeuAc(1),value3,value4
