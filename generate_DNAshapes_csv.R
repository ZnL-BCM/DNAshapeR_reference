#!/usr/bin/env Rscript
### generate_DNAhapes_csv.R
### Zian Liu
### 6/21/2021

## This script generates the DNA shapes table for k-mers given any odd number k.


# Library
library(DNAshapeR)   # Note: install from the author's github


# Read inputs
args = commandArgs(trailingOnly=TRUE)
  #defaults=c(5, "fasta_dict_5mer.fa", "ref_5mers_structure.csv"))
k = as.integer(args[1])
FaName = args[2]
OutputFile = args[3]
FlagCPG = as.logical(args[4])

print("Length of sequences k: ")
print(k)
print("Name of output fasta file: ")
print(FaName)
print("Name of final output csv file: ")
print(OutputFile)
print("Whether to consider CpG status (logical): ")
print(FlagCPG)


# Build the basic letters dictionary
dict_4 <- c('A','C','G','T')
# For each value in k, make the iteration
TotalIt = 4^(k-1)
Seqs = matrix(nrow=4^k, ncol=k)
for (i in 1:k) {
  Div = 4^(i-1)
  TempIt = TotalIt/Div
  Seqs[, i] <- rep(
    c(rep('A', TempIt), rep('C', TempIt), rep('G', TempIt), rep('T', TempIt)), Div)
}
# Collapse the matrix to form individual sequences
dict_fa_1 <- apply(format(Seqs), 1, paste, collapse="")
# Add titles
dict_fa_title <- paste(rep('>seq', 4^k), 1:(4^k), sep = '_')
# Interweave the previous two vectors into one vector & write it into text file
write.table(as.vector(rbind(dict_fa_title, dict_fa_1)), 
            file = FaName, row.names = FALSE, col.names = FALSE, quote = FALSE)
# Notification
print("Fasta sequences generated and saved.")
# If CpG is toggled, save a CpG-specific fasta
if (FlagCPG) {
  # Create dict_fa_2 which replaces all occurrences of "CG" to "Mg"
  dict_fa_2 <- gsub("CG", "Mg", dict_fa_1)
  # Interweave the previous two vectors into one vector & write it into text file
  write.table(as.vector(rbind(dict_fa_title, dict_fa_2)), 
              file = paste(FaName, "cpg", sep="."), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# Define all shapes used
ShapeTypeList = c('HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Buckle', 
                  'Opening', 'ProT', 'Shear', 'Stagger', 'Stretch', 'MGW', 'EP')
# If CpG is toggled, define CpG specific shapes
if (FlagCPG) {
  ShapeTypeListM = c("ProT", "HelT", "Roll", "MGW")
}


# Import the fastas and run shapes
ShapeTF = data.frame()
All_colnames = c()
for (ShType in ShapeTypeList){
  # Get DNA shape
  TempDF <- getShape(
    FaName, shapeType = ShType
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
print("DNA shapes generated.")


# If CpG is toggled, run CpG-specific shapes:
if (FlagCPG) {
  for (ShType in ShapeTypeListM){
    ShTypeM = paste(ShType, "m", sep="_")
    # Get DNA shape
    TempDF <- getShape(
      paste(FaName, "cpg", sep="."), shapeType = ShType, methylate=TRUE
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
  print("CpG-specific DNA shapes generated.")
}

# Change names
rownames(ShapeTF) <- dict_fa_1
colnames(ShapeTF) <- All_colnames
# Make an output
write.csv(ShapeTF, file = OutputFile, row.names = TRUE)
# Notification
print("CSV file generated.")