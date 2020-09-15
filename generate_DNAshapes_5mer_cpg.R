### generate_DNAhapes_7mer_cpg.R
### Zian Liu
### 6/29/2020

## This script generates the DNA shapes table in 7-mer format. The output can be used for DNA shape deviation analyses. 
## This script specifically assumes middle "CG" are all CpG: Annotate them as "Mg" instead

# Library
library(tidyverse)
library(openxlsx)
library(DNAshapeR)   # Note: install from the author's github for latest ver


# Build the dictionary and write it to a text file
a = c(rep('A', 256), rep('C', 256), rep('G', 256), rep('T', 256))
b = rep(c(rep('A', 64), rep('C', 64), rep('G', 64), rep('T', 64)), 4)
c = rep(c(rep('A', 16), rep('C', 16), rep('G', 16), rep('T', 16)), 16)
d = rep(c(rep('A', 4), rep('C', 4), rep('G', 4), rep('T', 4)), 64)
e = rep(c('A', 'C', 'G', 'T'), 256)
dict_fa_1 <- paste(a, b, c, d, e, sep = '')

# Create dict_fa_2 which replaces 3 and 4 to "Mg"
dict_fa_2_tmp <- strsplit(dict_fa_1, "")
dict_fa_2 <- character(length(dict_fa_1))
for (i in 1:length(dict_fa_2_tmp)){
  if (dict_fa_2_tmp[[i]][3] == 'C' & dict_fa_2_tmp[[i]][4] == 'G'){
    dict_fa_2_tmp[[i]][3] = 'M'
    dict_fa_2_tmp[[i]][4] = 'g'
  }
  dict_fa_2[i] <- paste(dict_fa_2_tmp[[i]], collapse = "")
}


dict_fa_title <- paste(rep('>seq', 1024), 1:1024, sep = '_')

# Interweave the previous two vectors into one vector & write it into text file
write.table(as.vector(rbind(dict_fa_title, dict_fa_2)), 
            file = 'fasta_dict_5mer_cpg.fa', 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# The package doesn't have an 'All' option, define all shapes used
shape_type_list = c('HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Buckle', 
                    'Opening', 'ProT', 'Shear', 'Stagger', 'Stretch', 'MGW', 'EP')

# Reload the text file & run it through DNAshapeR
pred_2 <- list()
for (shape_type in shape_type_list){
  if (shape_type %in% c('HelT', 'Roll', 'ProT', 'MGW')){
    pred_shape <- getShape(
      paste(getwd(), "fasta_dict_5mer_cpg.fa", sep = '/'), shapeType = shape_type, 
      methylate=TRUE
    )
  } else {
    pred_shape <- getShape(
      paste(getwd(), "fasta_dict_5mer.fa", sep = '/'), shapeType = shape_type
    )
  }
  pred_2[[shape_type]] <- pred_shape[[1]]
}

# Note that there are no edge predictions, so remove the NAs
# Also, add rownames to the values
for (z in 1:length(pred_2)){
  pred_2[[z]] <- as.data.frame(t(na.omit(t(pred_2[[z]]))))
  rownames(pred_2[[z]]) <- dict_fa_1   # Note that we will use the regular 5-mer names
}

# Make an output
write.xlsx(pred_2, file = 'ref_5mers_structure_cpg.xlsx', rowNames = TRUE)

# Done. 
