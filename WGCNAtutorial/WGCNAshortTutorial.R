###
# tutorial from https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA
###

source("http://bioconductor.org/biocLite.R")
biocLite("WGCNA")
install.packages("flashClust") 

library(WGCNA)
library(flashClust)

#STEP 1: uploading data into R and formatting it for WGCNA

#set your current working directory (where all your files are)
setwd("")

#this reads in the normalized counts file from DESeq2
time=read.csv("SampleTimeSeriesRLD.csv")

#Adjust file so it matches format WGCNA needs 
time=as.data.frame(time)
rownames(time) <-time$X
time$X=NULL
datExpr= as.data.frame(t(time[,]))
names(datExpr)= row.names(time)
rownames(datExpr)=names(time)
dim(datExpr)

#run this to check if there are gene outliers
gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK 
#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:

#if (!gsg$allOK)
#	{if (sum(!gsg$goodGenes)>0)
#		printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse= ", ")));
#		if (sum(!gsg$goodSamples)>0)
#			printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
#		datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
#		}

#upload your trait data
datTraits= read.csv("Traits_23May2015.csv")
head(datTraits)

#form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits) <- datTraits$Sample
datTraits$Sample <- NULL
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order

#expression data is in datExpr, corresponding clinical traits are datTraits


sampleTree2=flashClust(dist(datExpr), method="average")
traitColors= numbers2colors(datTraits, signed= FALSE)
plotDendroAndColors(sampleTree2, traitColors, groupLabels= names(datTraits), main="Sample Dendrogram and Trait heatmap")

save(datExpr, datTraits, file="SamplesAndTraits.RData")


### At this point you will need to identify sample outliers and chooose a soft threshold power. These are easy to do and are well documented in the online tutorials. Scripts for choosing a soft threshold are commented out below. It's important to choose the correct soft threhold for your dataset. 

#powers= c(c(1:10), seq(from =12, to=20, by=1)) #choosing a set of soft-thresholding powers
#sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="unsigned") #call network topology analysis function

#sizeGrWindow(9,5)
#par(mfrow= c(1,2))
#cex1=0.9
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
#abline(h=0.90, col="red")
#plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

#from this plot, we would choose a power of 12 becuase it's the lowest power for which the scale free topology index reaches 0.90

enableWGCNAThreads()
softPower=12
adjacency=adjacency(datExpr, power=softPower, type="unsigned")

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency, TOMType="unsigned")
dissTOM= 1-TOM

#generate a clustered gene tree
geneTree= flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module
minModuleSize=30 
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)

dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = 12)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")

#plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres= 0.0
merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors= merge$colors
mergedMEs= merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_unsigned_nomerge_RLDfiltered.RData")

#Define number of genes and samples
nGenes= ncol(datExpr)
nSamples= nrow(datExpr)

#Recalculate MEs with color labels
MEs0= moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs= orderMEs(MEs0)
moduleTraitCor= cor(MEs, datTraits, use= "p")
moduleTraitPvalue= corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(", 
						signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))
#display the corelation values with a heatmap plot
labeledHeatmap(Matrix= moduleTraitCor, 
			xLabels= names(datTraits), 
			yLabels= names(MEs), 
			ySymbols= names(MEs), 
			colorLabels= FALSE, 
			colors= blueWhiteRed(50), 
			textMatrix= textMatrix, 
			setStdMargins= FALSE, 
			cex.text= 0.5, 
			zlim= c(-1,1), 
			main= paste("Module-trait relationships"))
