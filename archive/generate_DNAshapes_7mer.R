### generate_DNAhapes_7mer.R
### Zian Liu
### 5/19/2020

## This script generates the DNA shapes table in 7-mer format. The output can be used for DNA shape deviation analyses. 


# Library
library(tidyverse)
library(openxlsx)
library(DNAshapeR)   # Note: install from the author's github for latest ver


# Build the dictionary and write it to a text file
dict_4 <- c('A','C','G','T')
a = c(rep('A', 4096), rep('C', 4096), rep('G', 4096), rep('T', 4096))
b = rep(c(rep('A', 1024), rep('C', 1024), rep('G', 1024), rep('T', 1024)), 4)
c = rep(c(rep('A', 256), rep('C', 256), rep('G', 256), rep('T', 256)), 16)
d = rep(c(rep('A', 64), rep('C', 64), rep('G', 64), rep('T', 64)), 64)
e = rep(c(rep('A', 16), rep('C', 16), rep('G', 16), rep('T', 16)), 256)
f = rep(c(rep('A', 4), rep('C', 4), rep('G', 4), rep('T', 4)), 1024)
g = rep(c('A', 'C', 'G', 'T'), 4096)

dict_fa_1 <- paste(a, b, c, d, e, f, g, sep = '')
dict_fa_title <- paste(rep('>seq', 16384), 1:16384, sep = '_')

# Interweave the previous two vectors into one vector & write it into text file
write.table(as.vector(rbind(dict_fa_title, dict_fa_1)), 
            file = 'DNAshapeR_reference/fasta_dict_7mer.fa', row.names = FALSE, col.names = FALSE, quote = FALSE)


# The package doesn't have an 'All' option, define all shapes used
shape_type_list = c('HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Buckle', 
                    'Opening', 'ProT', 'Shear', 'Stagger', 'Stretch', 'MGW', 'EP')

# Reload the text file & run it through DNAshapeR
pred_1 <- list()
for (shape_type in shape_type_list){
  pred_shape <- getShape(
    paste(getwd(), "DNAshapeR_reference/fasta_dict_7mer.fa", sep = '/'), shapeType = shape_type
  )
  pred_1[[shape_type]] <- pred_shape[[1]]
}

# Note that there are no edge predictions, so remove the NAs
# Also, add rownames to the values
for (z in 1:length(pred_1)){
  pred_1[[z]] <- as.data.frame(t(na.omit(t(pred_1[[z]]))))
  rownames(pred_1[[z]]) <- dict_fa_1
}

# Make an output
write.xlsx(pred_1, file = 'DNAshapeR_reference/ref_7mers_structure.xlsx', rowNames = TRUE)

# Done. 
