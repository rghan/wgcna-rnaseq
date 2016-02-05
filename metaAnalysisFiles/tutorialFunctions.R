userListEnrichment <- function (geneR, labelR, fnIn, catNmIn=1:1000, 
                                nameOut = "enrichment.csv"){
# This function measures list enrichment between inputted lists of genes
#  and files containing user-defined lists of genes.  Significant 
#  enrichment is measured using a hypergeometric test.
#
# geneR   = a vector of gene (or other) identifiers
# labelR  = a vector of labels (for example, module colors) corresponding
#           to the geneR list
# fnIn    = a vector of file names containing user-defined lists in one
#           of two formats: either gene lists ("Gene", "Var1") or module
#           membership tables ("PSID","Gene","Module",<other columns>).
# catNmIn = a vector of category names corresponding to each fnIn
# nameOut = name of the file to write the output enrichment information

  # Read in the data
  glIn = NULL
  for (i in 1:length(fnIn)){
    datIn = read.csv(fnIn[i])
    if (colnames(datIn)[2]=="Gene") datIn = datIn[,2:3]
    if (colnames(datIn)[2]=="Var1") datIn = datIn[,1:2]
    colnames(datIn)=c("Gene","Category")
    datIn[,2] = paste(catNmIn[i],datIn[,2],sep="_")
    glIn = rbind(glIn,datIn)
  }
  geneIn  = as.character(glIn[,1])
  labelIn = as.character(glIn[,2])
  
  # Format the data, excluding anything not in the inputted gene list 
  geneAll = sort(unique(geneR))
  keep    = is.element(geneIn,geneAll)
  geneIn  = geneIn[keep]
  labelIn = labelIn[keep]
  catsR   = sort(unique(labelR))
  catsIn  = sort(unique(labelIn))
  lenAll  = length(geneAll)
  
  # Determine the hypergeometric enrichment values 
  results = list(pValues = NULL, ovGenes = list(), sigOverlaps=NULL)
  namesOv = NULL
  for(r in 1:length(catsR)) for (i in 1:length(catsIn)) {
  	gr = geneR[labelR==catsR[r]]
  	gi = geneIn[labelIn==catsIn[i]]
  	go = intersect(gr,gi)
  	lr = length(gr)
  	li = length(gi)
  	lo = length(go)
  	pv = phyper2(lenAll,lr,li,lo,FALSE)
  	if (lo==0) pv=1
  	if (pv<0.0001) 	pv = phyper2(lenAll,lr,li,lo,TRUE)
  	pOut = c(catsR[r],catsIn[i],lr,li,lo,lenAll,pv)
  	results$pValues = rbind(results$pValues,pOut)
  	namesOv = c(namesOv,paste(catsR[r],"--",catsIn[i]))
  	results$ovGenes[[length(namesOv)]] = go
  }
  results$pValues = cbind(results$pValues, 
    apply(cbind(1,as.numeric(results$pValues[,7])*length(namesOv)),1,min))
  colnames(results$pValues) = c("InputCategories", "UserDefinedCategories", 
	"Input","UserList","Overlap","Total","Pvalues", "CorrectedPvalues")
  names(results$ovGenes) = namesOv
  results$sigOverlaps = results$pValues[as.numeric(results$pValues[,8])<0.05,c(1,2,8)]
  results$sigOverlaps = results$sigOverlaps[order(as.numeric(results$sigOverlaps[,3])),]
  
  write.csv(results$sigOverlaps,nameOut, row.names=FALSE)
  return(results)
}

##################################################

phyper2 <- function (total, group1, group2, overlap, verySig=TRUE ,lt=TRUE){
## This function is the same is phyper, just allows for more sensible input values
  q = overlap
  m = group1
  n = total-group1
  k = group2
  prob = phyper(q, m, n, k, log.p = verySig, lower.tail=lt)
  if (verySig) return(-prob)
  return(1-prob)
}

##################################################

matchModules <- function (gn1, mod1, gn2, mod2, omit="grey", allColors = standardColors()){
## This function converts the module colors from network 2 to the module colors
##  of best match in network 1, given the gene names (gn1,gn2) and corresponding
##  module assignments (mod1,mod2).  Non-overlapping modules will have unique labels.
## omit      = colors that should not be changed
## allColors = color set to choose from (default is all standardColors)
	
# Write out data from network 2 then read it back in to get module overlaps
	out2 = cbind(gn2, mod2);
	colnames(out2) = c("Gene","Var1")
	write.csv(out2,"eraseMe.csv",row.names=FALSE)
	overlaps = userListEnrichment(gn1, mod1, "eraseMe.csv", "X", "eraseMe.csv")
	c1 = as.character(overlaps$sigOverlaps[,1])
	c2 = as.character(overlaps$sigOverlaps[,2]);
	c2 = substr(c2,3,nchar(c2))
	kp = (!is.element(c1,omit))&(!is.element(c2,omit))
	c1 = c1[kp];  c2 = c2[kp]
	
# Change the labels in network 2 to the best overlapping module in network 1
	cOld <- cNew <- changed <- sort(unique(mod2))
	cUnused = setdiff(allColors,union(mod1,mod2))
	while(length(c1)>0){
		cNew[cOld==c2[1]]=c1[1]
		changed[cOld==c2[1]]="YES"
		kp = (c1!=c1[1])&(c2!=c2[1])
		c1 = c1[kp];  c2 = c2[kp]
	}
	changed[is.element(cOld,omit)] = "YES"
	cNew[changed!="YES"] = cUnused[1:sum(changed!="YES")]
	modOut = mod2
	for (i in 1:length(cNew)) modOut[mod2==cOld[i]] = cNew[i]
	write(paste("Old - New labels:",cOld,"-",cNew),"")
	return(modOut)
}

##################################################

t.test.l = function(x){
## Performs a t-test on a vector of genes
	tt = t.test(x[var[[1]]],x[var[[2]]])
	return(c(tt$est,sd(x[var[[1]]]),sd(x[var[[2]]]),tt$p.val))
}

##################################################

cor.test.l = function(x){
## Performs a Pearson correlation on a vector of genes
 ct = cor.test(x,var)
 return(c(ct$est,ct$p.val))
}

##################################################

visantPrepOverall <- function(colorFinal2, moduleColor, datExprrest, genes, numInt, power, signed=FALSE)
{
## This file returns the overall strongest connections within a module for one group of subjects

## USER inputs
# colorFinal2 = vector with module association for each probe
# moduleColor = color of the module... should be a member of colorFinal2
# datExprrest = expression data for the genes corresponding to cIndex
# genes = list of genes that correspond to the probes
# numInt = number of interactions to output to the visant file

cIndex = (colorFinal2==moduleColor)
datExpr=datExprrest[,cIndex]
if (signed){AdjMat =((1+cor(datExpr, use="p"))/2)^power
	} else {AdjMat =abs(cor(datExpr, use="p"))^power}
diag(AdjMat)=0
Degree <- apply(AdjMat,1,sum)
Degree = Degree/max(Degree)
meanExpr=apply(datExpr,2,mean)
ProbeSet=colnames(datExpr)
GeneSet=genes[cIndex];
Module=rep(moduleColor,length(meanExpr)) 

## This file summarizes intramodular connectivity and expression for each gene in each group:
fn=paste(moduleColor,"_connectivityOverall.csv",sep="")
datConn = cbind(ProbeSet,GeneSet,meanExpr,Degree,Module)
datConn = datConn[order(Degree, decreasing=TRUE),]
write.table(datConn,file=fn,sep=",",
            row.names=F, col.names= c("probes","genes","meanExpr","kin","module"))
write(paste(fn, "written."),"")

## TOM matrix
distTOM <- TOMdist1(AdjMat)
simTOM = 1-distTOM
diag(simTOM)=0
simTOMcutoff = simTOM

## Correlation matrix
Pearson = cor(datExpr ,use="p")
diag(Pearson)=0

## Dynamically determine the appropriate cutoff
cutoff = 0.24
len    = 10000
dir    = "increase"
loops  = 0
split  = 0.01
numInt = numInt+100
while(len>100){
 loops = loops+1
 if (dir == "increase") { cutoff = cutoff+split; }
 if (dir == "decrease") { cutoff = cutoff-split; }
 indices = (simTOMcutoff>cutoff)
 len = sum(sum(indices))
 if (len < numInt) {dir = "decrease";}
 if (len >=numInt) {dir = "increase";}
 len = abs(len-numInt)
 if (loops>500){ len=0;}  
 if ((loops%%100)==0){ split = split/2; }
}
write(c(loops,cutoff,len),"")

## Output using cutoffs:
indices = (simTOMcutoff[1,]>cutoff)
datout=cbind(
       rep(GeneSet[1], length(ProbeSet[indices])),
       GeneSet[indices],rep(0, length(ProbeSet[indices])),
       rep("M1002", length(ProbeSet[indices])),
       simTOM[1,][indices], Pearson[1,][indices])
colnames(datout) = c("gene1","gene2","zero","color","TO","Correlation")
for(i in seq(2,length(ProbeSet),by=1)){
 indices = (simTOMcutoff[i,]>cutoff)
 datout=rbind(datout,cbind(
        rep(GeneSet[i], length(ProbeSet[indices])),
        GeneSet[indices],rep(0, length(ProbeSet[indices])),
        rep("M1002", length(ProbeSet[indices])),
        simTOM[i,][indices], Pearson[i,][indices]))
}
datout = datout[order(datout[,5],decreasing=TRUE),]
fn = paste(moduleColor,"_visantOverall.csv",sep="")
write.table(datout,file=fn,sep=",",row.names=F, col.names=c("gene1","gene2",
            "zero","color","TO","Correlation"))
write(paste(fn, "written."),"")
}

##################################################

#The function TOMdist1 computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
if(exists("TOMdist1")) rm(TOMdist1);
TOMdist1=function(adjmat1, maxADJ=FALSE) {
diag(adjmat1)=0;
adjmat1[is.na(adjmat1)]=0;
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>10^(-12)   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
adjmat1= (adjmat1+ t(adjmat1) )/2
kk=apply(adjmat1,2,sum)
maxADJconst=1
if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1); 
gc();gc();
numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
#TOMmatrix=numTOM/denomTOM
# this turns the TOM matrix into a dissimilarity 
out1=1-as.matrix(numTOM/denomTOM) 
diag(out1)=1 
# setting the diagonal to 1 is unconventional (it should be 0)
# but it leads to nicer looking TOM plots... 
out1
}}}

#####################################################

visantPrep <- function(colorFinal2, moduleColor, ciIndex, msIndex, datExprrest, genes, cutoff, power=power, signed=FALSE)
{
## This file compares two groups of subjects and returns the connections that are differential between the two groups

## USER inputs
# colorFinal2 = vector with module association for each probe
# moduleColor = color of the module... should be a member of colorFinal2
# ciIndex = index of control subjects
# msIndex = index of AD subjects
# datExprrest = expression data for the genes corresponding to cIndex
# genes = list of genes that correspond to the probes
# cutoff = NOT A CUTOFF!!! NUMBER OF INTEGERS TO INCLUDE IN THE ANALYSIS!!!
# power = power of the network
# signed = is this a signed network TRUE/FALSE

cIndex = (colorFinal2==moduleColor)
datExprCI=datExprrest[ciIndex,cIndex]
datExprMS=datExprrest[msIndex,cIndex]
if (signed){AdjMatCI =((1+cor(datExprCI, use="p"))/2)^power
} else {AdjMatCI =abs(cor(datExprCI, use="p"))^power}
if (signed){AdjMatMS =((1+cor(datExprMS, use="p"))/2)^power
} else {AdjMatMS =abs(cor(datExprMS, use="p"))^power}
diag(AdjMatCI)=0
diag(AdjMatMS)=0
DegreeCI <- apply(AdjMatCI,1,sum)
DegreeMS <- apply(AdjMatMS,1,sum)
DegreeCI = DegreeCI/max(DegreeCI)
DegreeMS = DegreeMS/max(DegreeMS)
meanExprCI=apply(datExprCI,2,mean)
meanExprMS=apply(datExprMS,2,mean)
ProbeSet1=colnames(datExprMS)
GeneSet1=genes[cIndex];
Module=rep(moduleColor,length(meanExprMS)) 

if (cutoff>=0){
## This file summarizes intramodular connectivity and expression for each gene in each group:

write.table( cbind(ProbeSet1,GeneSet1,meanExprCI,meanExprMS,DegreeCI,DegreeMS,Module), 
  file= paste(moduleColor,"_connectivity.csv",sep=""),sep=",",row.names=F, 
  col.names= c("Probes","genes","meanExprCI","meanExprMS","kin_CI","kin_MS","Module"))

## Output correlation between kIn_CI and kIn_MS
write("kIn_CI vs. kIN_MS correlation","")
corr = cor(DegreeCI, DegreeMS, method="s")
write(corr,"")
}

## CI and MS TOM matrices

distTOMci <- TOMdist1(AdjMatCI)
distTOMms = TOMdist1(AdjMatMS)
simTOMms=1-distTOMms
simTOMci = 1-distTOMci
diag(simTOMci)=0
diag(simTOMms)=0

## Creating a matrix to describe relative similarity between human and chimp
# A value of 0.5 means that there is no difference between the two groups, 0 = higher similarity in MS group, 1 = higher similarity in CI group

simRatio=simTOMci/mean(simTOMci)/(simTOMci/mean(simTOMci)+ simTOMms/mean(simTOMms)+0.00001)

## Correlation matrices

PearsonCI = cor(datExprCI ,use="p")
PearsonMS = cor(datExprMS ,use="p")
diag(PearsonCI)=0
diag(PearsonMS)=0

## The code below lets you choose a cutoff within one of the TOM matrices.  You may want to experiment choosing cutoffs in both the low NFT and high NFT groups.

## Cutoffs for controls (cutoff>0) or AD (cutoff<0) group
if (cutoff>=0)
  simTOMcutoff = simTOMci
if (cutoff<0){
  simTOMcutoff = simTOMms
  cutoff = -cutoff
}

## Dynamically determine the appropriate cutoff
numInt = cutoff
cutoff = 0.24
len    = 10000
dir    = "increase"
loops  = 0
split  = 0.01
numInt = numInt+100
while(len>100){
 loops = loops+1
 if (dir == "increase") { cutoff = cutoff+split; }
 if (dir == "decrease") { cutoff = cutoff-split; }
 indices = (simTOMcutoff>cutoff)
 len = sum(sum(indices))
 if (len < numInt) {dir = "decrease";}
 if (len >=numInt) {dir = "increase";}
 len = abs(len-numInt)
 if (loops>500){ len=0;}  
 if ((loops%%100)==0){ split = split/2; }
}
write(c(loops,cutoff,len),"")

## Output using cutoffs:

for(i in c(1) ) {
indices = (simTOMcutoff[i,]>cutoff)
datout=data.frame(
rep(ProbeSet1[i], length(ProbeSet1[indices])),
rep(GeneSet1[i], length(ProbeSet1[indices])),
ProbeSet1[indices], GeneSet1[indices],
simTOMci[i,][indices], simTOMms[i,][indices],
simRatio[i,][indices], PearsonCI[i,][indices],
PearsonMS[i,][indices])
}

for(i in seq(2,length(ProbeSet1),by=1)) {
indices = (simTOMcutoff[i,]>cutoff)
datout=data.frame(rbind(datout,data.frame(
rep(ProbeSet1[i], length(ProbeSet1[indices])),
rep(GeneSet1[i], length(ProbeSet1[indices])),
ProbeSet1[indices], GeneSet1[indices],
simTOMci[i,][indices], simTOMms[i,][indices],
simRatio[i,][indices], PearsonCI[i,][indices],
PearsonMS[i,][indices])))
}

write.table(datout,file=paste(moduleColor,"_visant.csv",sep=""), sep=",",row.names=F, col.names=c("probe1","gene1","probe2","gene2","TO_CI","TO_MS","TO_Ratio","Correlation_CI","Correlation_MS"))
}

#######################################################

probe2Gene <- function(probes, probeList, geneList)
{
## This function converts a group of probenames to genenames, taking into
## account the fact that sometimes probelists have X added to the front

## USER inputs
# probes = the probes you want converted
# probeList = the list of probe names
# geneList = the list of gene names IN THE SAME ORDER as the probes

probes<-sub("X", "", probes)
names(geneList) = probeList
return(geneList[probes])
}

########################################################

write.geneList <- function(PG, filename, allProbes=0, allGenes=0, probe="g")
{
## These functions write a genelist / probelist to a file of geneNames

## USER inputs
# PG = the probe/gene you want written to a gene list
# allProbes = the list of probe names for the above probes
# allGenes = the list of gene names for the corresponding probes
# filename = the filename (can include folder)
# probe = the default ("g") says PG is a gene and doesn''t need to be converted
#         to a gene.  Otherwise PG is assumed to be a probe and converted

gene = PG
if (probe!="g") {
  gene = probe2Gene(PG,allProbes,allGenes)
}
write(gene,filename,sep="\n")

}

