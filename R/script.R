# script to follow for Cramer's lab to construct WGCNA networks for microarray analysis
# firt load the necessary libraries - already installed here
library(WGCNA)
library(cluster)
options(stringAsFactors=F)

# always a good idea to set your directory
# but if you are using R studio and create a new project, the directory will be set by default to your project directory
# to check
getwd()
# to set
setwd(file.choose())

# read in the microarray data 
# the data needs to be stored as a tab delimited file
mic.data <- read.delim(file.choose())

# check for dimension of the dataset
dim(mic.data)
# note that dataset may contain information on genes or different database identifier
# we need only one and the replicates - therefore, remove others 
# I will use the GeneId. Also, I will transpose the dataset
mic.data.reps <- as.data.frame(t(mic.data[,-c(2:3)]))
# set colnames to first col of mic.data
names(mic.data.reps) <- mic.data$V1_GeneID
# set rownames to colnames of mic.data without 2nd and 3rd colnames
rownames(mic.data.reps) <- names(mic.data)[-c(2:3)]
# remove first row
mic.data.reps <- mic.data.reps[-1,]
# copy the datafile for correlations based on the average
# its questionable if the average will be good enough, since we will only get about 10 datapoints, which is on the edge to get
# a strong correlation
mic.data.avg <- mic.data.reps  # we will calculate the average later, once we will have all the data

# now read in the metabolite data
met.data <- read.delim(file.choose())
dim(met.data)
names(met.data)
# we don't need cols 2 and 3 and transpose again
met.data.reps <- as.data.frame(t(met.data[,-c(2:3)]))
# set row and colnames
names(met.data.reps) <- met.data$metabolite.markers
rownames(met.data.reps) <- names(met.data)[-c(2:3)]
# remove first row
met.data.reps <- met.data.reps[-1,]

# rownames of the met.data.rep contain underscore "_"
# in order to synchronize the datasets we need to get rid of it
# loop over rownames
}
# look at rownames of datasets
rownames(met.data.reps)
rownames(mic.data.reps)

# match and sync rownames of datasets
match.rows <- match(tolower(rownames(met.data.reps)),tolower(rownames(mic.data.reps)))
# this shows that the variety name do not match
# Chardonnay is abbreviated CD in the microarray dataset and CH in the metabolite dataset
# thus fix names in the metabolite dataset
str <- c("CDW1","CDW2","CDW3","CDW4","CDW5","CDW6","CDD1","CDD2","CDD3","CDD4","CDD5","CDD6")
rownames(met.data.reps)[49:60] <- str
# repeat
match.rows <- match(tolower(rownames(mic.data.reps)),tolower(rownames(met.data.reps)))
if(length(match.rows)<1) match.rows <- match.rows[-which(is.na(match.rows))]

met.data.reps <- met.data.reps[match.rows,]
# show that rownames concur - should be all true
table(tolower(rownames(met.data.reps))==tolower(rownames(mic.data.reps)))

# find outlier in the microarray data set
# sample network based on squared Euclidean distance
# note that the data is transposed
A.mic.data <- adjacency(t(mic.data.reps),type="distance")
# this calculates the whole network connectivity
k <- as.numeric(apply(A.mic.data,2,sum))-1
# standardized connectivity
Z.k <- scale(k)

# designate samples as outlying
# don't do more than two rounds of this
# if the Z.value is below the threshold -5 or -2.5
# the color vector indicates outlyingness (red)
outlier.color <- ifelse(Z.k < -5, "red","black")

# calculat the cluster tree using flashClust or hclust
sample.tree <- flashClust(as.dist(1-A.mic.data),method="average")
# convert metabolites/traits to a color representation where red indicates high values
met.colors <- data.frame(numbers2colors(convert_to_num(as.matrix(met.data.reps)),signed=F))
dimnames(met.colors)[[2]] <- paste(names(met.data.reps),"C",sep="")
dat.colors <- data.frame(outlierC=outlier.color,met.colors)

# generate a dendrogram revealing the outliers
plotDendroAndColors(sample.tree,groupLables=names(dat.colors),
                    colors=dat.colors, main="Outlier dendrogram and trait heatmap")
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
# the first line in the dendrogram shows which sample is an outlier by highlighting the outlier with a red rectangle
# remove outliers in both datasets and recompute the sample network among the remaining samples
# as said above don't do more than 2 cycles. you may now find again new outliers

# no outliers detected - go on to calculate averages

# copy the datafile for correlations based on the average
# its questionable if the average will be good enough, since we will only get about 10 datapoints, which is on the edge to get
# a strong correlation
mic.data.avg <- mic.data.reps  # we will calculate the average later, once we will have all the data
# copy metabolite dataset for average calculations
met.data.avg <- met.data.reps  

# now we compute the average for the datasets
# first we need to remove the replicate number at the end of the replicate name
# we need to convert the dataset to a matrix first
met.data.avg <- as.matrix(met.data.avg)
for(i in 1:nrow(met.data.avg)){
  rownames(met.data.avg)[i] <- substr(rownames(met.data.avg)[i],1,nchar(rownames(met.data.avg)[i])-1)
}
# we need to make sure that the data is numeric - use function convert_to_num
met.data.avg <- convert_to_num(met.data.avg)
# do the same for mircroarray data
mic.data.avg <- as.matrix(mic.data.avg)
for(i in 1:nrow(mic.data.avg)){
  rownames(mic.data.avg)[i] <- substr(rownames(mic.data.avg)[i],1,nchar(rownames(mic.data.avg)[i])-1)
}
# this may take a few seconds
mic.data.avg <- convert_to_num(mic.data.avg)
# now compute averages - for that you will need to load the function calc_mean - ignore the warnings
met.data.avg <- calc_mean(met.data.avg)
# do for microarray data - this may take a while
mic.data.avg <- calc_mean(mic.data.avg)
  
# now go on to identify co-expressin modules
# choose a set of thresholding powers
powers <- c(1:20)
# choose power based on SFT criterion
sft.reps <- pickSoftThreshold(convert_to_num(as.matrix(mic.data.reps)),powerVector=powers, networkType="signed")
sft.avg <- pickSoftThreshold(mic.data.avg,powerVector=powers, networkType="signed")

# plot the results 
#with reps first

par(mfrow=c(1,2))
# SFT index as a functin of different powers
plot(sft.reps$fitIndices[,1],-sign(sft.reps$fitIndices[,3])*sft.reps$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="SFT, signed R^2", type="n", main=paste("Scale independence with reps"))
text(sft.reps$fitIndices[,1],-sign(sft.reps$fitIndices[,3])*sft.reps$fitIndices[,2], labels=powers, col="red")
# this line corresponds to using on R^2 cut-off of h
abline(h=0.9, col="red")
# mean connectivity as function of different powers
plot(sft.reps$fitIndices[,1], sft.reps$fitIndices[,5], type="n",
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", main=paste("Mean connectivity with reps"))
text(sft.reps$fitIndices[,1],sft.reps$fitIndices[,5],labels=powers,col="red")


# now for average
par(mfrow=c(1,2))
# SFT index as a functin of different powers
plot(sft.avg$fitIndices[,1],-sign(sft.avg$fitIndices[,3])*sft.avg$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="SFT, signed R^2", type="n", main=paste("Scale independence with average"))
text(sft.avg$fitIndices[,1],-sign(sft.avg$fitIndices[,3])*sft.avg$fitIndices[,2], labels=powers, col="red")
# this line corresponds to using on R^2 cut-off of h
abline(h=0.9, col="red")
# mean connectivity as function of different powers
plot(sft.avg$fitIndices[,1], sft.avg$fitIndices[,5], type="n",
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", main=paste("Mean connectivity with average"))
text(sft.avg$fitIndices[,1],sft.avg$fitIndices[,5],labels=powers,col="red")

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
beta=20     # select beta based on suggested powerEstimate and the plot (look at both sft.reps & sft.avg plots)
module.size <- 30      #30 for microarrays, 10 for proteins, and 5 for metabolites
# automatic module detection via dynamic tree cutting
mergingThresh <- 0.25


# this is automatic - skip {
net <- blockwiseModules(convert_to_num(as.matrix(mic.data.reps)), corType="pearson",
                        maxBlockSize=5000, networkType="signed",power=beta,minModuleSize=module.size,
                        mergeCutHeight=mergingThresh,numericLabels=T, saveTOMs=T,
                        pamRespectsDendro=F, saveTOMFileBase="mic2010TOM")      ### expect to get a weird error while its calculating
moduleLabelsAutomatic=net$colors
# convert labels to colors for plotting
moduleColorsAutomatic=labels2colors(moduleLabelsAutomatic)

# data frame of eigengenes
MEsAutomatic <- net$MEs
# }
# this is to first metabolite - tartarate
# include this if you want to correlate to one trait only
#*********************************************************************************
tartarate <- as.data.frame(met.data.reps[,1])
names(tartarate) <- "tartarate"
# now use this metabolite to define a gene significance variable
GS.tartarate <- as.numeric(cor(mic.data.reps,tartarate,use="p"))  # p = pearson
# this translates the numeric values into colors
GS.tartarateColor <- numbers2colors(GS.tartarate,signed=T)
#*********************************************************************************

# optional {
blocknumber <- 1
datColors <- data.frame(moduleColorsAutomatic, GS.tartarateColor)[net$blockGenes[[blocknumber]],]
# plot the dendrogram
plotDendroAndColors(net$dendrograms[[blocknumber]], colors=datColors,
                    groupLabels=c("Module colors","Tartaric acid"), dendroLabels=F,
                    hang=0.03, addGuide=T,guideHang=0.05)

# since we have quite a big dataset, we need to do split the module detection
# this may take a while

 # }
# this is the manual one - do this!!!!

bw.net <- blockwiseModules(convert_to_num(as.matrix(mic.data.reps)), corType="pearson",
                           maxBlockSize=5000, networkType="signed",power=beta,
                           minModuleSize=module.size, mergeCutHeight=mergingThresh, numericLabels=T,
                           saveTOMs=T,pamRespectsDendro=F, saveTOMFileBase="grape_mic_data_blockwise",verbose=3)

# relabel blockwise modules so that their labels match those from our pervious analysis
#module.labels.blockwise <- matchLabels(bw.net$colors, moduleLabelsAutomatic)
module.labels.blockwise <- bw.net$colors
# convert labels to colors for plotting
module.colors.blockwise=labels2colors(module.labels.blockwise)
datColors <- data.frame(moduleColorsAutomatic, GS.tartarateColor)[bw.net$blockGenes[[blocknumber]],]
# measuring agreement
mean(module.labels.blockwise==moduleLabelsAutomatic)
bn <- 2   # blocknumber - we have 7
plotDendroAndColors(bw.net$dendrograms[[bn]], 
                    module.colors.blockwise[bw.net$blockGenes[[bn]]],"Module colors",
                    main=paste("Dendogram and module colors in block", bn),
                    dendroLabels=F,hang=0.03, addGuide=T,guideHang=0.05)

MEs.blockwise <- bw.net$MEs

# remember that if you restart a session you have to reload the libraries

################################################################################################################
# let's see if we can get better results with manual module detection
# note that all these steps may take a while to process
A.man <- adjacency(convert_to_num(as.matrix(mic.data.reps)),power=beta, type="signed")

# define (dis)similarity based on the topological overlap
dissTOM <- TOMdist(A.man)
# hierarchical clustering
geneTree <- flashClust(as.dist(dissTOM), method="average")
# define modules by cutting branches
module.labels.manual.1 <- cutreeDynamic(dendro=geneTree, distM=dissTOM,
                                        method="hybrid", deepSplit=2,pamRespectsDendro=F,minClusterSize=module.size)

# relabel the manual modules so that their labels match those of the automatic clustering
module.labels.manual.2 <- matchLabels(module.labels.manual.1,module.labels.blockwise)
# convert labels to colors for plotting
# crucial step
module.colors.manual.2 <- labels2colors(module.labels.manual.2)
# calculate eigengenes
MEList.man <- moduleEigengenes(mic.data.reps,colors=module.colors.manual.2)
MEs.man <- MEList.man$eigengenes

# optional{
##################################################################################################################################
# add tartarate to the existing module eigengenes
ME.tar <- orderMEs(cbind(MEs.man,tartarate))
ME.tar.aut <- orderMEs(cbind(MEs.blockwise,tartarate))
# plat the relationships among the eigengenes and tartarate - automatic
plotEigengeneNetworks(ME.tar.aut,"",marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),cex.lab=0.5,xLabelsAngle=90)
# plat the relationships among the eigengenes and tartarate - manual
plotEigengeneNetworks_I(ME.tar,"",marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),cex.lab=0.5,xLabelsAngle=90)
##################################################################################################################################

# merge highly correlated modules in manual eigengenes
merge.man <- mergeCloseModules(mic.data.reps,module.colors.manual.2,cutHeight=mergingThresh)
# resulting merged module colors
module.colors.merge <- merge.man$colors
# eigengenes of the merged modules
MEs.man.merge <- merge.man$newMEs

# plat the relationships among the eigengenes and tartarate - manual
ME.tar.merge <- orderMEs(cbind(MEs.man.merge,tartarate))
plotEigengeneNetworks(ME.tar.merge,"",marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),cex.lab=0.5,xLabelsAngle=90)


# showing the effect of moudle mergig by plotting the original and merged module colors below the tree
datColors=data.frame(module.colors.blockwise, module.colors.merge, module.colors.manual.2,GS.tartarateColor)
plotDendroAndColors(geneTree,colors=datColors,
                   groupLabels=c("Blockwise","Manual hybrid","Manual","Tartaric acid"),
                   dendroLabels=F,hang=0.03,addGuide=T,guideHang=0.05)
# check agreement between manual and automatic module labels
# optional
# optional{
##################################################################################################################################
mean(module.colors.merge==module.colors.blockwise)
mean(module.colors.manual.2==module.colors.blockwise)
mean(module.colors.merge==moduleColorsAutomatic)
mean(module.colors.merge==module.colors.manual.2)

# 0.106973 - not a good match
# optional
##################################################################################################################################
plotDendroAndColors(geneTree,colors=moduleColorsAutomatic,
                    groupLabels="Automatic",
                    dendroLabels=F,hang=0.03,addGuide=T,guideHang=0.05)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# continuing with manual hybrid
# now we will try to detect the relationship between the eigengenes and the metabolites

# choosing the manual hybrid as our model
# define numbers of genes and samples
num.genes <- ncol(mic.data.reps)
num.samples <- nrow(mic.data.reps)
# recalculate the MEs with color labels
MEs.hybrid <- moduleEigengenes(mic.data.reps,module.colors.merge)$eigengenes
MEs.hybrid <- orderMEs(MEs.hybrid)
# correlating to the metabolites
mod.met.cor <- cor(MEs.hybrid,met.data.reps,use="p")
mod.met.pvalue <- corPvalueStudent(mod.met.cor,num.samples)

# display correlations and their p-values as heatmap
str.matrix <- paste(signif(mod.met.cor,2),"\n(",signif(mod.met.pvalue,1),")",sep="")
dim(str.matrix)=dim(mod.met.cor)
par(mfrow=c(1,1),mar=c(9,8,3,3))
labeledHeatmap(Matrix=mod.met.cor,xLabels=colnames(met.data.avg),yLabels=names(MEs.hybrid),
               ySymbols=names(MEs.hybrid), colorLabels=F, colors=blueWhiteRed(50),
               setStdMargins=F,cex.text=0.35,zlim=c(-1,1), main=paste("Module-metabolite relationships"),cex.lab=0.5)
# calculate the module membership values = module eigengene based connectivity kME
dat.KME <- signedKME(mic.data.reps,MEs.hybrid)


# as tartaric acid is the first metabolite in the dataset, let's look at these relationships specifically, i.e. identifying genes 
# significance for tartaric acid
# looking at the data we see strong positive correlation to indianred1 = row 9 and negative to khaki4 = row 25
# using GS and MM data
color.columns <- substring(names(dat.KME),4)

select.modules=c("indianred1","khaki4")
par(mfrow=c(length(select.modules)/2,length(select.modules)/2))
par(mfrow=c(1,2),mar=c(5,4,4,2))
for(module in select.modules){
  column <- match(module,color.columns)
  rest.module <- module.colors.merge==module
  verboseScatterplot(dat.KME[rest.module,column],GS.tartarate[rest.module],   # this is als your data
                     xlab=paste("Module Membership", module,"module"), ylab="Tartaric acid",
                     main=paste("kME.", module,"vs GS"),col=module,cex=0.6, cex.text=0.5, cex.main=0.7,
                     cex.lab=0.7, cex.axis=0.5)
}

# the relationship as indicated at the top of the plot shows that there is either a very strong positive or negative relationship 
# between the modules and tartaric acid, respectively

#*************************************************************************************************************************************************
#=================================================================================================================================================
#*************************************************************************************************************************************************
#=================================================================================================================================================
#*************************************************************************************************************************************************
# important script
# determining automatically high correlation modules between eigengenes and metabolites
# generate correlation between metabolites and genes - this may take a while

GS.mets  <- NULL
for(i in 1:ncol(met.data.reps))
  GS.mets <- cbind(GS.mets, as.numeric(cor(mic.data.reps,met.data.reps[,i],use="p")))

colnames(GS.mets) <- colnames(met.data.reps)
rownames(GS.mets) <- colnames(mic.data.reps)    

# 1. determine r coefficient and p-value threshold for eigengene-metabolite network
# make histogram of absolute r coefficient
hist.var <- vector("integer",10)
for(i in 1:length(hist.var)){hist.var[i] <- length(which(abs(mod.met.cor)>=i/10-0.1 & abs(mod.met.cor)<i/10))}
# we go for 0.5
fdr.eigen.met <- fdr.fun(p.value=mod.met.pvalue) 
# qvalue of 0.01 = p 0.01
# create data matrix containing only cells surpassing threshold
eigen.met.sig <- correlation.network.thresh(cor = mod.met.cor,p.cor = mod.met.pvalue,r=0.5,p=0.01)
# make scatter plot for each metabolite and its significant correlation to the eigengene
# columns are metabolites, rows are eigengenes
for(i in 1:ncol(eigen.met.sig)){
  print(colnames(eigen.met.sig)[i])
  rows <- which(eigen.met.sig[,i]!=0)
  if(length(rows)==0) next # skip
  select.modules <- rownames(eigen.met.sig)[rows]
  # omit first two characters
  select.modules <- substr(select.modules,3,nchar(select.modules))
  # preparing plot
  pdf(file=paste(colnames(eigen.met.sig)[i],".pdf",sep=""))
  if(length(select.modules)==1)
    par(mfrow=c(1,1), mar=c(5,4,4,2))
  else if(length(select.modules)==2) 
    par(mfrow=c(1,2), mar=c(5,4,4,2))
  else if(round(sqrt(length(select.modules)))==sqrt(length(select.modules))) 
    par(mfrow=c(sqrt(length(select.modules)),sqrt(length(select.modules))), mar=c(5,4,4,2))
  else if(round(sqrt(length(select.modules)))<sqrt(length(select.modules))) 
    par(mfrow=c(ceiling(sqrt(length(select.modules))),floor(sqrt(length(select.modules)))), mar=c(5,4,4,2))
  else if(round(sqrt(length(select.modules)))>sqrt(length(select.modules))) 
    par(mfrow=c(ceiling(sqrt(length(select.modules))),ceiling(sqrt(length(select.modules)))), mar=c(5,4,4,2))
  #print("ok1")
  for(module in select.modules){
    column <- match(module,color.columns)
   # print("ok2")
    rest.module <- module.colors.merge==module
    print("")
    verboseScatterplot(dat.KME[rest.module,column],GS.mets[rest.module,i],
                       xlab=paste("Module Membership", module,"module"), ylab=colnames(eigen.met.sig)[i],
                       main=paste("kME.", module,"vs GS"),col=module,cex=0.6, cex.text=0.5, cex.main=0.7,
                       cex.lab=0.7, cex.axis=0.5)
  }
  dev.off()
}

#*************************************************************************************************************************************************
#=================================================================================================================================================
#*************************************************************************************************************************************************
#=================================================================================================================================================
#*************************************************************************************************************************************************

####################################################################################################################################
####################################################################################################################################
# creating networks for cytoscape
# merge the eigengene network with the metabolite network
eigengene.met.matrix <- cbind(MEs.hybrid,met.data.reps)
eigengene.met.netw <- correlation.thresholds(eigengene.met.matrix,0,1)
eigengene.met.prop <- eigengene.met.netw@prop
# writing the property variable containing different p value and r coefficient thresholds
write.csv(eigengene.met.netw@prop,"Eigengen_metabolite_network_properties.csv")
# look at file and choose a r coefficient where network stabilizes
# here we choose 0.4
r.coef <- 0.4
# now calculate p-value threshold by applying an fdr
p.fdr <- fdr.fun(p.value=triangle(eigengene.met.netw@p.value.mat))
# store fdr values
write.csv(p.fdr,"qvalues.csv")
# look at file and search for qvalue of 0.01 = p value is 0.013
p.value <- 0.013
# make adjacency matrix with calculated r-coef and p-value
eigengene.met.adj <- correlation.network.thresh(eigengene.met.netw@cor.matrix, eigengene.met.netw@p.value.mat,r.coef,p.value)
# write this adj matrix for future usage
write.csv(eigengene.met.adj,"Weighted_correlation_adj_matrix_eigengene_metabolites.csv")
# make igraph object
eigengene.met.igraph <- get.igraph(eigengene.met.adj)
# calculate membership
eigengene.met.memb <- netw.membership.igraph(eigengene.met.igraph)
# load compound class file
cc <- read.delim(file.choose(),check.names=F)
# run cytoscape preparation funcion
cytoscape_preparation(data_matrix=eigengene.met.adj,comp_class_file=cc,membership=eigengene.met.memb,"Eigenge_metabolite")


####################################################################################################################################
####################################################################################################################################
# multidimensional scaling
cmd1 <- cmdscale(as.dist(dissTOM),2)
par(mfrow=c(1,1))
plot(cmd1,col=module.colors.merge,main="MDS plot",xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

# connectivity plot
# setting the diagonal
diag(dissTOM) <- NA
# transform dissTOM with beta
TOMplot(dissim=dissTOM^beta,dendro=geneTree, colors=module.colors.merge, main="Network heatmap plot of all genes")

