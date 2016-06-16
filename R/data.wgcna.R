## load data
## set multicore for DESeq2
register(MulticoreParam(22))
#===============================================================================##
## count data (filtered version removing rows with zero counts in all libraries)
count <- read.csv("./data/all.ww.bard.2012_gene_counts.csv",
                  header = TRUE,
                  row.names=1,
                  stringsAsFactors = FALSE)
# filter out rows (genes) containing only zero counts, reducing it to 27,926 from 29,971
count <- count[rowSums(count^2)>0,]
head(count)

## Create metadata to describe experimental conditions:
meta <- data.frame(
  row.names = colnames(count),
  brix = rep(c("20","22","24","26"), each=3),
  cultivar = rep(c("CS","ME","PN","CD","SM","CF","SB"), each=12),
  skin = rep(c(rep("Red",36),rep("White",24),rep("Red",12),rep("White",12))))

## Make design matrix
#brix <- relevel(meta$brix, ref="20")
#cultivar <- relevel(meta$cultivar, ref="CS")
modelMeta <- data.frame(model.matrix(~0 + brix + cultivar + skin, meta)) %>% 
  # add sample and skinRed columns
  mutate(skinRed = ifelse(skinWhite == 0, 1, 0),
         cultivarCD = ifelse(cultivarCF == 0 & cultivarCS == 0 & cultivarME == 0
                             & cultivarPN == 0 & cultivarSB == 0 
                             & cultivarSM == 0, 1 , 0)) %>% 
  # reorder columns
  select(brix20:cultivarPN, cultivarCD, cultivarSB:skinRed)
# assign rownames from meta
rownames(modelMeta) <- rownames(meta)
write.csv(modelMeta, "./data/modelMeta.csv", row.names = TRUE)


## build an DESeqDataSet from count and meta information
dds <- DESeqDataSetFromMatrix(
  countData = count,
  colData = meta,
  design = ~ brix * cultivar)

## Assign gene names to dds, not sure what to use
#rownames(dds) <- rownames(count)
mcols(dds) <- cbind(mcols(dds), count) # from https://support.bioconductor.org/p/62374/
all(rownames(dds) == rownames(count))

## Pre-filtering the dataset
nrow(dds)

## remove all genes with counts < 15 in more than 75% of samples (84*0.75=63)
## suggested by WGCNA on RNAseq FAQ
## 75% - maintains the TPS gene in dataset. balances low counts/false positives
dds75 <- dds[ rowSums(counts(dds) >= 15) >= 63, ]
nrow(dds75)     ## 16,606 genes in data set

#------------------------------------------------------------------------------
## some other options I played with for filtering counts in all samples
#------------------------------------------------------------------------------
## remove all genes with counts < 20 in more than 90% of samples (84*0.9=75.6)
## suggested by WGCNA on RNAseq FAQ
dds90 <- dds[ rowSums(counts(dds) >= 20) >= 75.6, ]
nrow(dds90)
## remove all genes with counts < 5in more than 80% of samples (84*0.8=67.2)
## suggested by another WGCNA pipeline
dds80 <- dds[ rowSums(counts(dds) >= 10) >= 67.2, ]
nrow(dds80)
#------------------------------------------------------------------------------


## Perform the DESeq2 normalization, required before WGCNA analysis
dds2 <- DESeq(dds75, betaPrior = FALSE, parallel = TRUE)
## Perform a variance-stabilizing transformation
vsd <- getVarianceStabilizedData(dds2)
## Many functions expect the matrix to be transposed
datExpr <- t(vsd) 
## check rows/cols
nrow(datExpr)
ncol(datExpr)

## save a copy 
write.csv(datExpr, "./data/datExpr.brix_08Apr16.csv", row.names = TRUE)
## proceed to wgcna_analysis.R

## create a list of V1 genes used in WGCNA that can be used for other purposes
geneNames <- rownames(vsd)
geneNames.array <- as.data.frame(geneNames)
colnames(geneNames.array) <- "V1"
write.csv(geneNames.array, "./data/filtered.genes4wgcna.csv", row.names = FALSE)