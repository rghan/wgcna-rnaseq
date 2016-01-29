#===============================================================================
## Block-wise network construction and module detection, attempting to look at
## network connectivity of all genes at once.
#===============================================================================
bwnet = blockwiseModules(datExpr, maxBlockSize = 4000, # half of the # of genes
                         power = 10, TOMType = "signed", minModuleSize = 30,
                         corType = "bicor",
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "berrybrixTOM-blockwiseSigned.20counts",
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

