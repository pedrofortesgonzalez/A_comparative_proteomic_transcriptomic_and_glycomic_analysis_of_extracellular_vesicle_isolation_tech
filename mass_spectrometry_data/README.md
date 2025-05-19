# Input Data
This folder is where you should place your mass spectrometry data files for analysis.

## Requirements
- All files must be in CSV format
- Each file must contain at least the following columns:
  - `Accession` - Protein/peptide identifier in UniProtKB Accession format
  - `Peptide` - Peptide sequence
  - `PTM` - Post-translational modifications

## Important Notes
- This folder is automatically created as a temporary copy of your original data during processing
- The script `1_ptm_detection.py` will ask for the name of your input folder
- Files will be processed in their original format - no renaming is necessary
- Make sure your input folder name does not contain spaces

## Example
A properly formatted input CSV file should have headers and data similar to:

```csv
Accession,Peptide,PTM,OtherColumn1,OtherColumn2
P12345,PEPTIDESEQUENCE,HexNAc(1)Hex(1)Fuc(1),value1,value2
Q67890,ANOTHERPEPTIDE,HexNAc(2)Hex(2)NeuAc(1),value3,value4
