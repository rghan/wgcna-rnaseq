#################################################################################
## packages.R - script to load & update packages for analysis
#################################################################################
## update packages
update.packages()
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("BiocUpgrade")

## load packages
library(WGCNA)
library(cluster)
library(flashClust)

library(dplyr)
library(tidyr)

library(DESeq2)
library(BiocParallel)
