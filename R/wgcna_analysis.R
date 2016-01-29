#################################################################################
## wgcna_analysis.R - analyze the RNAseq data from data.R
#################################################################################
library(WGCNA)
# To allow multi-threading within WGCNA with all available cores, use 
allowWGCNAThreads(nThreads = 22) 
# or let WGCNA pick total number
enableWGCNAThreads()
# to disable threading if necessary.
disableWGCNAThreads() 
# Alternatively, set the following environment variable on your system:
ALLOW_WGCNA_THREADS=14

#################################################################################
# Load variance stabilized DESeq2 data (vsd2) into WGCNA and reformat
#################################################################################
#set your current working directory (where all your files are)
setwd(".")

#form a data frame analogous to expression data that will hold the clinical traits.
datTraits <- modelMeta
rownames(datTraits) <- modelMeta$sample
datTraits$sample <- NULL
table(rownames(datTraits)==rownames(vsd2)) # should return TRUE if datasets align correctly

# calculate the cluster tree using flashClust or hclust
wgcnaTree <- flashClust(dist(vsd2), method = "average")
traitColors <- numbers2colors(datTraits, signed = TRUE)
# Plot the sample tree
pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(wgcnaTree,
                    traitColors,
                    groupLabels= names(datTraits),
                    main="Outlier dendrogram and trait heatmap")
dev.off()

#output data from R
save(vsd2, datTraits, file="SamplesAndTraits.RData")
#################################################################################
# soft-thresholding
#################################################################################
## choosing a set of soft-thresholding powers
powers= c(c(1:10), seq(from = 12, to = 20, by = 1))
## call network topology analysis function
sft = pickSoftThreshold(vsd2, powerVector = powers, verbose = 5, 
                        blockSize = 20000, networkType="signed hybrid",
                        corFnc = "bicor", 
                        corOptions = list(use = 'p', maxPOutliers = 0.1))
## plot the results:
sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", 
     type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels=powers, cex=cex1, col="red")
## this line corresponds to using an R^2 cut-off of h
abline(h=0.90, col="red")
## Mean connectivity as a funciton of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", 
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

## from this plot, we would choose a power of 10 becuase it's the lowest power for 
## which the scale free topology index reaches 0.90

#################################################################################
## Block-wise network construction and module detection, attempting to look at
## network connectivity of all genes at once.
bwnet = blockwiseModules(vsd2, maxBlockSize = 4000, # half of the # of genes
                         power = 15, TOMType = "signed", minModuleSize = 30,
                         corType = "bicor",
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "berrybrixTOM-blockwiseSigned",
                         verbose = 3)

# Load the results of single-block analysis
#load(file = "FemaleLiver-02-networkConstruction-auto.RData");
# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels);
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

## To see how many modules were identified and what the module sizes, one can use
## table(bwLabels):
table(bwLabels)

# open a graphics window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Comparing the single block and block-wise network analysis
sizeGrWindow(12,9)
plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors,dynamicColors, mergedColors),
                    c("Single block", "2 blocks", "dynamic colors", "merged colors"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#################################################################################
# Construct a gene co-expression matrix and generate modules
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
enableWGCNAThreads()
#################################################################################
## Co-expression similarity and adjacency using assigned softpower
softPower=13
adjacency=adjacency(vsd2, power=softPower, type="signed hybrid")

#################################################################################
## Topological Overlap Matrix (TOM)
# Turn adjacency into topological overlap, i.e. translate the adjacency into 
# topological overlap matrix and calculate the corresponding dissimilarity:
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

## generate a clustered gene tree
geneTree <- flashClust(as.dist(dissTOM), method="average")

sizeGrWindow(8, 6)
plot(geneTree, xlab = "", sub = "", main = "Gene Clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# This sets the minimum number of genes to cluster into a module
minModuleSize <- 30 
# generating modules and assigning them colors
# Module identification using dynamic tree cut: 
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = minModuleSize)
table(dynamicMods)

## Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
## Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#################################################################################
## Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList <- moduleEigengenes(vsd2, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-bicor(MEs)
# Cluster module eigengenes
METree <- flashClust(as.dist(MEDiss), method= "average")
#plots tree showing how the eigengenes cluster together
sizeGrWindow(10,7)
plot(METree, main= "Clustering of module eigengenes", xlab = "", sub = "")

## Note: the tutorial chooses a height cut of 0.25, corresonding to correlation
## of 0.75, to merge

# set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres <- 0.3
# plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
# Call an automatic merging function
merge <- mergeCloseModules(vsd2, dynamicColors, cutHeight = MEDissThres, 
                           corFnc = "bicor", verbose = 3)
# The merged module colors
mergedColors <- merge$colors
length(unique(mergedColors))
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs


# To see the merging upon module colors, plot the gene dendrogram again, with
# the original and merged colors underneath
sizeGrWindow(12, 9)
pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

## In subsequent analysis, the merged module colors in mergedColors will be used
# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresonding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "brix.vsd.filtered.networkConstruction-stepByStep.RData")

#################################################################################
# Quantifying module-trait associations
#################################################################################
## Correlate eigengenes with external traits and look for the most significant associations
# Define numbers of genes and samples
nGenes <- ncol(vsd2)
nSamples <- nrow(vsd2)
# Recalculate MEs with color labels
MEsO <- moduleEigengenes(vsd2, moduleColors)$eigengenes
MEs <- orderMEs(MEsO)
moduleTraitCor <- bicor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

## Color code each association by correlation value
sizeGrWindow(12, 9)
# Will display correlations and their p-values
pdf(file = "Plots/module-trait.signed.relationships.pdf", width = 12, height = 9)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#################################################################################
## Gene relationship to trait & important modules: Gene Significance and Module Membership

# Note: see pg 16 of the ch12 tutorial for details on plotting multiple scatter
# plots.

# Define variable brix20 containing the brix20 column of datTrait
brix20 = as.data.frame(datTraits$brix20);
names(brix20) = "skinWhite"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(bicor(vsd2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(bicor(vsd2, brix20, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(brix20), sep="");
names(GSPvalue) = paste("p.GS.", names(brix20), sep="");

module = "lightcyan2"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#=====================================================================================
# Visualizing of networks within R
#=====================================================================================
## Visualizing the gene network
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(vsd2, power = 13);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^12; # color is faint at power=7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# remove the grey module
length(moduleColors=="grey")
restGenes <- (moduleColors != "grey")
# Call the plot function
sizeGrWindow(12,9)
TOMplot(plotTOM, geneTree, restGenes, main = "Network heatmap plot, all genes")
