### example.py
### Zian Liu
### This file shows one way to import the xlsx files generated from DNAshapeR through my scripts. 

# Libraries
import numpy as np
import pandas as pd

# Designate folder; I like to put the DNAshapeR outputs in a separate folder
dir_shape = "../DNAshapeR_reference/"
filename_ref = "ref_7mers_structure_cpg.xlsx"

# Import file with pandas
F_strucref = pd.ExcelFile(dir_shape + filename_ref)
# These are the names for the shape features
Sheetnames = ['HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Buckle', 'Opening', 'ProT', 'Shear', 
              'Stagger', 'Stretch', 'MGW', 'EP']
# One by one import the sheets
DF_strucref = pd.DataFrame()
for item in Temp_sheetname:
    temp_df = F_strucref.parse(sheet_name=item, index_col=0)
    if np.shape(temp_df)[1] == 3:
        temp_df.rename({'V1': item+'_L', 'V2': item+'_C', 'V3': item+'_R'}, axis='columns', inplace=True)
    else:
        temp_df.rename({'V1': item+'_L', 'V2': item+'_CL', 'V3': item+'_CR', 'V4': item+'_R'}, axis='columns', inplace=True)
    DF_strucref = pd.concat([DF_strucref, temp_df], axis=1)

# This will then result in one table where the headers are designated "L (left)", "C (center)", "R (right)", or "CL (center left)" and "CR (center right)".
# The following is one way to parse the table based on my predictors:

Seq_list = np.array(['AAAAAAA', 'AAAAATT', 'CCAACTT', 'CCCCCCG'])
Pred = DF_strucref.loc[Seq_list, :]
print(Pred)

# Done!
