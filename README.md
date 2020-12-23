# DNAshapeR_reference
### Zian Liu
#### Last updated: 12/23/2020

## Description

This is a collection of convenient R scripts for generating reference xlsx tables based on the DNAshape method (Zhou et al. 2013) and the DNAshapeR package (Chiu et al. 2016). 

The generated xlsx tables will have 14 tabs corresponding to the 14 shape features available in the package. In each tab, the rows represent unique 7-mers, while the columns represent position. 

This can then be easily imported to Python/other programming languages for downstream analyses. Example script is provided (*example.py*).

Note that DNAshapeR provides CpG-based shapes as CpG context can drastically affect DNA shape features. They have only provided four shapes, so the CpG shape features use CpG shapes for the four features with CpG alterations, while the other 10 shapes will use regular shape features.

If you have any questions, you can go to the DNAshapeR vignette curated by the original authors with *browseVignettes("DNAshapeR")* in R, or check out Tsu-Pei's awesome Github page (http://tsupeichiu.github.io/DNAshapeR/).


## Citations

I do not claim any rights to these scripts, all rights belong to the original authors; feel free to reuse/redistribute. Please remember to cite the original authors if you are using these scripts:

#### All scripts:

Chiu T, Comoglio F, Zhou T, Yang L, Paro R, Rohs R (2016). “DNAshapeR: an R/Bioconductor package for DNA shape prediction and feature encoding.” *Bioinformatics*, 32, 1211-1213.

Chiu T, Rao S, Mann R, Honig B, Rohs R (2017). “Genome-wide prediction of minor-groove electrostatic potential enables biophysical modeling of protein-DNA binding.” *Nucleic Acids Res.*, 45(21), 12565-12576.

Li J, Sagendorf J, Chiu T, Pasi M, Perez A, Rohs R (2017). “Expanding the repertoire of DNA shape features for genome-scale studies of transcription factor binding.” *Nucleic Acids Res.*, 45(22), 12877-12887.

#### If methylation shapes are used, also cite this article:
Rao S, Chiu T, Kribelbauer J, Mann R, Bussemaker H, Rohs R (2018). “Systematic prediction of DNA shape changes due to CpG methylation explains epigenetic effects on protein-DNA binding.” *Epigenetics Chromatin*, 11:6.


## Thank you!
