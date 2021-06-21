### generate_DNAhapes_5mer.R
### Zian Liu
### 4/6/2020

## This script generates the DNA shapes table 5-mers. For baseline analyses. 


# Library
library(tidyverse)
library(openxlsx)
library(DNAshapeR)   # Note: install from the author's github


# Build the dictionary and write it to a text file
dict_4 <- c('A','C','G','T')
a = c(rep('A', 256), rep('C', 256), rep('G', 256), rep('T', 256))
b = rep(c(rep('A', 64), rep('C', 64), rep('G', 64), rep('T', 64)), 4)
c = rep(c(rep('A', 16), rep('C', 16), rep('G', 16), rep('T', 16)), 16)
d = rep(c(rep('A', 4), rep('C', 4), rep('G', 4), rep('T', 4)), 64)
e = rep(c(rep('A', 1), rep('C', 1), rep('G', 1), rep('T', 1)), 256)

dict_fa_1 <- paste(a, b, c, d, e, sep = '')
dict_fa_title <- paste(rep('>seq', 1024), 1:1024, sep = '_')

# Interweave the previous two vectors into one vector & write it into text file
write.table(as.vector(rbind(dict_fa_title, dict_fa_1)), 
            file = 'fasta_dict_5mer.fa', row.names = FALSE, col.names = FALSE, quote = FALSE)


# The package doesn't have an 'All' option, define all shapes used
shape_type_list = c('HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Buckle', 'Opening', 'ProT', 'Shear', 'Stagger', 'Stretch', 'MGW', 'EP')

# Reload the text file & run it through DNAshapeR
pred_1 <- list()
for (shape_type in shape_type_list){
  pred_1[[shape_type]] <- getShape(
    paste(getwd(), "fasta_dict_5mer.fa", sep = '/'), shapeType = shape_type
  )
}

# Note that there are no edge predictions, so remove the NAs
# Also, add rownames to the values
for (z in 1:length(pred_1)){
  pred_1[[z]] <- as.data.frame(t(na.omit(t(pred_1[[z]]))))
  rownames(pred_1[[z]]) <- dict_fa_1
}

# Make an output
write.xlsx(pred_1, file = 'ref_5mers_structure.xlsx', rowNames = TRUE)

# Done. 
