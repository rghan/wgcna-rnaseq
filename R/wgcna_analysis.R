#===============================================================================
## wgcna_analysis.R - analyze the RNAseq data from data.R
#===============================================================================
library(WGCNA)
# To allow multi-threading within WGCNA with all available cores, use 
allowWGCNAThreads(nThreads = 22) 
# or let WGCNA pick total number
# enableWGCNAThreads() # not advised as it assigns 23 cores causing problems
# to disable threading if necessary.
disableWGCNAThreads() 
# Alternatively, set the following environment variable on your system:
ALLOW_WGCNA_THREADS=14

#===============================================================================
# Load variance stabilized DESeq2 data (datExpr) into WGCNA and reformat
#===============================================================================
datExpr <- read.csv("./data/datExpr.brix.csv", header = TRUE,
                    stringsAsFactors = FALSE, row.names = 1)
# assign vsd rownames to geneNames vector, used for TOM below
geneNames <- rownames(vsd)
# form a data frame analogous to expression data that will hold the clinical traits.
datTraits <- read.csv("./data/modelMeta.csv", header = TRUE, 
                      stringsAsFactors = FALSE, row.names = 1)

# Check if expression and trait data row names match
table(rownames(datTraits)==rownames(datExpr)) # should return TRUE if datasets align correctly

#===============================================================================
# Begin WGCNA
#===============================================================================
#set your current working directory (where all your files are)
setwd(".")

# calculate the cluster tree using flashClust or hclust
wgcnaTree <- flashClust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = TRUE)
# Plot the sample tree
pdf(file = "Plots/sampleClustering.20counts.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(wgcnaTree,
                    traitColors,
                    groupLabels= names(datTraits),
                    main="Outlier dendrogram and trait heatmap")
dev.off()

#output data from R
save(datExpr, datTraits, file="SamplesAndTraits.RData")
#===============================================================================
# soft-thresholding
#===============================================================================
## choosing a set of soft-thresholding powers
powers= c(c(1:10), seq(from = 12, to = 20, by = 1))
## call network topology analysis function
## may need to set this if an error is called
## allowWGCNAThreads(nThreads = 22)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, 
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


#===============================================================================
# Construct a gene co-expression matrix and generate modules
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
#enableWGCNAThreads() # This seems to cause problems w/ various functions
allowWGCNAThreads(nThreads = 22) 

#===============================================================================
## Co-expression similarity and adjacency using assigned softpower
softPower=11
adjacency=adjacency(datExpr, power=softPower, type="signed hybrid")

#===============================================================================
## Topological Overlap Matrix (TOM)
# Turn adjacency into topological overlap, i.e. translate the adjacency into 
# topological overlap matrix and calculate the corresponding dissimilarity:
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
colnames(TOM) = rownames(TOM) = geneNames # from PKLab Harvard 
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
#===============================================================================
## Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
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
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, 
                           corFnc = "bicor", verbose = 3)
# The merged module colors
mergedColors <- merge$colors
length(unique(mergedColors))
table(mergedColors)
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs


# To see the merging upon module colors, plot the gene dendrogram again, with
# the original and merged colors underneath
sizeGrWindow(12, 9)
pdf(file = "Plots/geneDendro.20counts.pdf", wi = 9, he = 6)
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

#===============================================================================
# Quantifying module-trait associations
#===============================================================================
## Correlate eigengenes with external traits and look for the most significant associations
# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
# Recalculate MEs with color labels
MEsO <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEsO)
moduleTraitCor <- bicor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

## Color code each association by correlation value
sizeGrWindow(12, 9)
# Will display correlations and their p-values
pdf(file = "Plots/mod-trait.relationships20.counts.pdf", width = 12, height = 9)
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

#===============================================================================
## Gene relationship to trait & important modules: Gene Significance and Module Membership

# Note: see pg 16 of the ch12 tutorial for details on plotting multiple scatter
# plots.

# Define variable brix20 containing the brix20 column of datTrait
brix20 = as.data.frame(datTraits$brix20);
names(brix20) = "skinWhite"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(bicor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(bicor(datExpr, brix20, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(brix20), sep="");
names(GSPvalue) = paste("p.GS.", names(brix20), sep="");

module = "turquoise"
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

#===============================================================================
# Visualizing of networks within R
#===============================================================================
## Visualizing the gene network
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 11);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^12; # color is faint at power=7


# discard the unassigned genes, and focus on the rest
# from: http://pklab.med.harvard.edu/scw2014/WGCNA.html
restGenes= (dynamicColors != "grey")
dissTOM=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower)

# look at the network w/out "grey" modules
colnames(dissTOM) = rownames(dissTOM) = geneNames[restGenes]
hier1=flashClust(as.dist(dissTOM), method = "average" )
plotDendroAndColors(hier1, colors = data.frame(dynamicColors[restGenes],moduleColors[restGenes]),
                    c("Dynamic-merged Tree Cut","Merged-module Tree Cut"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "Gene dendrogram and module colors")

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
diag(dissTOM) = NA

# Call the plot function
sizeGrWindow(12,9)
TOMplot(plotTOM^, hier1, as.character(moduleColors[restGenes]))
plotTOM = dissTOM^2; # color is faint at power=7
diag(plotTOM) = NA;
#TOMplot(plotTOM, hier1, moduleColors, main = "Network heatmap plot, all connected genes")

#===============================================================================
# based on the chpt 12 - all genes
#===============================================================================
A = adjacency(datExpr, power = softPower, type = "signed")
plotDendroAndColors(geneTree, colors = data.frame(dynamicColors,moduleColors),
                    c("Dynamic-merged Tree Cut","Merged-module Tree Cut"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "Gene dendrogram and module colors")
dissTOM = TOMdist(A)
diag(dissTOM) = NA
png("./Plots/TOMplot.brix.png", width = 600, height = 600)
TOMplot(dissim = dissTOM^7, geneTree, 
        moduleColors, 
        main = "Network heatmap plot, all genes")
dev.off()
#===============================================================================
# based on the chpt 12 - reomved grey genes
#===============================================================================
A2 = adjacency(datExpr[restGenes], power = softPower, type = "signed")
dissTOM2 = TOMdist(A2)
## generate a clustered gene tree
geneTree2 <- flashClust(as.dist(dissTOM2), method="average")

plotDendroAndColors(geneTree2, colors = data.frame(dynamicColors[restGenes],moduleColors[restGenes]),
                    c("Dynamic-merged Tree Cut","Merged-module Tree Cut"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors\nremoved unassigned gene modules (MEgrey)")

diag(dissTOM2) = NA
png("./Plots/TOMplot.brix.png", width = 600, height = 600)
TOMplot(dissim = dissTOM2^7, geneTree, 
        as.character(moduleColors[restGenes]), 
        main = "Network heatmap plot, removed Grey module genes")
dev.off()