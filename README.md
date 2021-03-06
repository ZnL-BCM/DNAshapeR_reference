# DNAshapeR_reference
### Zian Liu
#### Last updated: 6/21/2021

## Description

This is a collection of convenient R scripts for generating reference excel tables based on the DNAshape method (Zhou et al. 2013) and the DNAshapeR package (Chiu et al. 2016). 

The generated csv tables will have 14 types of shapes corresponding to the 14 shape features available in the package. The rows represent unique 7-mers, while the columns represent shape type and position. If CpG is considerd, there will be 18 types of shapes instead.

This can be easily imported to Python by something like ``pandas.read_csv(filename, header=0, index_col=0)`` or into other programming languages.

If you have any questions, please go to the DNAshapeR vignette curated by the original authors with *browseVignettes("DNAshapeR")* in R, or check out Tsu-Pei's GitHub page (http://tsupeichiu.github.io/DNAshapeR/).

## Usage

For usage, run the following:
```
Rscript generate_DNAshapes_csv.R k fasta_name.fa output_name.csv <TRUE/FALSE>
```
The four commands correspond to:
* k, the length of k-mers (odd integer)
* Output prefix for intermediate fasta files (recommend to put in a separate directory!)
* Output csv file name
* Whether to consider CpG methylation (logical variable)

Note that the **archive/** folder contains scripts which include hard-coded inputs and outputs, and these scripts are NOT GUARANTEED to work. Please use them with caution; using the *generate_DNAshapes_csv.R* script is recommended. 


## Citations

I do not claim any rights to these scripts, all rights belong to the original authors; feel free to reuse/redistribute. Please remember to cite the original authors if you are using these scripts:

#### All scripts:

Chiu T, Comoglio F, Zhou T, Yang L, Paro R, Rohs R (2016). “DNAshapeR: an R/Bioconductor package for DNA shape prediction and feature encoding.” *Bioinformatics*, 32, 1211-1213.

Chiu T, Rao S, Mann R, Honig B, Rohs R (2017). “Genome-wide prediction of minor-groove electrostatic potential enables biophysical modeling of protein-DNA binding.” *Nucleic Acids Res.*, 45(21), 12565-12576.

Li J, Sagendorf J, Chiu T, Pasi M, Perez A, Rohs R (2017). “Expanding the repertoire of DNA shape features for genome-scale studies of transcription factor binding.” *Nucleic Acids Res.*, 45(22), 12877-12887.

#### If methylation shapes are used, also cite this article:
Rao S, Chiu T, Kribelbauer J, Mann R, Bussemaker H, Rohs R (2018). “Systematic prediction of DNA shape changes due to CpG methylation explains epigenetic effects on protein-DNA binding.” *Epigenetics Chromatin*, 11:6.


## Thank you!
