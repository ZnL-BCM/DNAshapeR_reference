### generate_DNAhapes_7mer_csv.R
### Zian Liu
### 4/22/2021

## This script generates the DNA shapes table in 7-mer format. 
#  The output can be used for DNA shape deviation analyses. 


# Library
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
ShapeTypeList = c('HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Buckle', 
                  'Opening', 'ProT', 'Shear', 'Stagger', 'Stretch', 'MGW', 'EP')

# Import the fastas and run shapes
ShapeTF = data.frame()

All_colnames = c()
for (ShType in ShapeTypeList){
  # Get DNA shape
  TempDF <- getShape(
    "fasta_dict_7mer.fa", shapeType = ShType
  )[[ShType]]
  # Remove all NA columns
  TempDF <- TempDF[, colSums(is.na(TempDF)) < 1]
  # Change column names
  if (length(TempDF) / length(dict_fa_1) < 2) {
    All_colnames <- c(All_colnames, paste(ShType, 1, sep="_"))
  } else {
    All_colnames <- c(All_colnames, paste(ShType, 1:ncol(TempDF), sep="_"))
  }
  # Add data
  if (ncol(ShapeTF) == 0) {
    ShapeTF <- TempDF
  } else {
    ShapeTF = cbind(ShapeTF, TempDF)
  }
}

rownames(ShapeTF) <- dict_fa_1
colnames(ShapeTF) <- All_colnames

# Make an output
write.csv(ShapeTF, file = "ref_7mers_structure.csv", row.names = TRUE)


# Done. 