### generate_DNAshapes_9mer_csv_cpg.R
### Zian Liu
### 5/11/2021

## This script generates the DNA shapes table 9-mers both with and without methylation 


# Library
library(DNAshapeR)   # Note: install from the author's github


# Build the dictionary and write it to a text file
dict_4 <- c('A','C','G','T')
a = c(rep('A', 65536), rep('C', 65536), rep('G', 65536), rep('T', 65536))
b = rep(c(rep('A', 16384), rep('C', 16384), rep('G', 16384), rep('T', 16384)), 4)
c = rep(c(rep('A', 4096), rep('C', 4096), rep('G', 4096), rep('T', 4096)), 16)
d = rep(c(rep('A', 1024), rep('C', 1024), rep('G', 1024), rep('T', 1024)), 64)
e = rep(c(rep('A', 256), rep('C', 256), rep('G', 256), rep('T', 256)), 256)
f = rep(c(rep('A', 64), rep('C', 64), rep('G', 64), rep('T', 64)), 1024)
g = rep(c(rep('A', 16), rep('C', 16), rep('G', 16), rep('T', 16)), 4096)
h = rep(c(rep('A', 4), rep('C', 4), rep('G', 4), rep('T', 4)), 16384)
i = rep(c(rep('A', 1), rep('C', 1), rep('G', 1), rep('T', 1)), 65536)

dict_fa_1 <- paste(a, b, c, d, e, f, g, h, i, sep = '')
dict_fa_title <- paste(rep('>seq', 262144), 1:262144, sep = '_')

# Create dict_fa_2 which replaces all occurrences of "CG" to "Mg"
dict_fa_2 <- gsub("CG", "Mg", dict_fa_1)

# Interweave the previous vectorss & write them into text file
write.table(as.vector(rbind(dict_fa_title, dict_fa_1)), 
            file = 'fasta_dict_9mer.fa', row.names = FALSE, col.names = FALSE, quote = FALSE)
# Interweave the previous two vectors into one vector & write it into text file
write.table(as.vector(rbind(dict_fa_title, dict_fa_2)), 
            file = 'fasta_dict_9mer_cpg.fa', 
            row.names = FALSE, col.names = FALSE, quote = FALSE)



# The package doesn't have an 'All' option, define all shapes used
ShapeTypeList = c('HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Buckle', 
                  'Opening', 'ProT', 'Shear', 'Stagger', 'Stretch', 'MGW', 'EP')
ShapeTypeListM = c("ProT", "HelT", "Roll", "MGW")


# Import the fastas and run shapes
ShapeTF = data.frame()

All_colnames = c()

for (ShType in ShapeTypeList){
  # Get DNA shape
  TempDF <- getShape(
    "fasta_dict_9mer.fa", shapeType = ShType
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

for (ShType in ShapeTypeListM){
  ShTypeM = paste(ShType, "m", sep="_")
  # Get DNA shape
  TempDF <- getShape(
    "fasta_dict_9mer_cpg.fa", shapeType = ShType, methylate=TRUE
  )[[ShType]]
  # Remove all NA columns
  TempDF <- TempDF[, colSums(is.na(TempDF)) < 1]
  # Change column names
  if (length(TempDF) / length(dict_fa_1) < 2) {
    All_colnames <- c(All_colnames, paste(ShType, "m", 1, sep="_"))
  } else {
    All_colnames <- c(All_colnames, paste(ShType, "m", 1:ncol(TempDF), sep="_"))
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
write.csv(ShapeTF, file = "ref_9mers_structure_cpg.csv", row.names = TRUE)


# Done. 