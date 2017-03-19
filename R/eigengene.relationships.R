#===============================================================================
## Grant asked for this pg 12 of ch12, 05 Apr 16
# Calculate eigengenes 
# MEList=moduleEigengenes(datExprFemale,colors=moduleColorsManual2) 
# MEs = MEList$eigengenes

# Plot the relationships among the eigengenes and all traits
pdf(file = "Plots/eigengene.relationships.All_08Apr16.pdf", wi = 15, he = 12)
plotEigengeneNetworks(MEs,"",marDendro=c(0,3,1,2), marHeatmap=c(3,5,0.5,0),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## Red-skin
colorRed = as.data.frame(datTraits$skinRed)
# rename colname
names(colorRed) = "Red_skin"
# Add the weight to existing module eigengenes 
MEsRed <- orderMEs(cbind(MEs,colorRed))
# Plot the relationships among the eigengenes and the trait
pdf(file = "Plots/eigengene.relationships.Redskin_08Apr16.pdf", wi = 15, he = 12)
#sizeGrWindow(15, 15) 
plotEigengeneNetworks(MEsRed,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## White-skin
colorWhite = as.data.frame(datTraits$skinWhite)
# rename colname
names(colorWhite) = "White_skin"
# Add the weight to existing module eigengenes 
MEsWhite <- orderMEs(cbind(MEs,colorWhite))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.Whiteskin_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEsWhite,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## Brix20
brix20 = as.data.frame(datTraits$brix20)
# rename colname
names(brix20) = "Brix_20"
# Add the weight to existing module eigengenes 
MEsbrix20 <- orderMEs(cbind(MEs,brix20))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.Brix20_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEsbrix20,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## Brix22
brix22 = as.data.frame(datTraits$brix22)
# rename colname
names(brix22) = "Brix_22"
# Add the weight to existing module eigengenes 
MEsbrix22 <- orderMEs(cbind(MEs,brix22))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.Brix22_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEsbrix22,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## Brix24
brix24 = as.data.frame(datTraits$brix24)
# rename colname
names(brix24) = "Brix_24"
# Add the weight to existing module eigengenes 
MEsbrix24 <- orderMEs(cbind(MEs,brix24))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.Brix24_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEsbrix24,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## Brix26
brix26 = as.data.frame(datTraits$brix26)
# rename colname
names(brix26) = "Brix_26"
# Add the weight to existing module eigengenes 
MEsbrix26 <- orderMEs(cbind(MEs,brix26))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.Brix26_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEsbrix26,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## cultivarCS
cultCS = as.data.frame(datTraits$cultivarCS)
# rename colname
names(cultCS) = "cultCS"
# Add the weight to existing module eigengenes 
MEscultCS <- orderMEs(cbind(MEs,cultCS))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.cultivarCS_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEscultCS,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## cultivarCF
cultCF = as.data.frame(datTraits$cultivarCF)
# rename colname
names(cultCF) = "cultCF"
# Add the weight to existing module eigengenes 
MEscultCF <- orderMEs(cbind(MEs,cultCF))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.cultivarCF_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEscultCF,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## cultivarME
cultME = as.data.frame(datTraits$cultivarME)
# rename colname
names(cultME) = "cultME"
# Add the weight to existing module eigengenes 
MEscultME <- orderMEs(cbind(MEs,cultME))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.cultivarME_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEscultME,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## cultivarPN
cultPN = as.data.frame(datTraits$cultivarPN)
# rename colname
names(cultPN) = "cultPN"
# Add the weight to existing module eigengenes 
MEscultPN <- orderMEs(cbind(MEs,cultPN))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.cultivarPN_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEscultPN,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## cultivarCD
cultCD = as.data.frame(datTraits$cultivarCD)
# rename colname
names(cultCD) = "cultCD"
# Add the weight to existing module eigengenes 
MEscultCD <- orderMEs(cbind(MEs,cultCD))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.cultivarCD_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEscultCD,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## cultivarSB
cultSB = as.data.frame(datTraits$cultivarSB)
# rename colname
names(cultSB) = "cultSB"
# Add the weight to existing module eigengenes 
MEscultSB <- orderMEs(cbind(MEs,cultSB))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.cultivarSB_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEscultCF,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#======================================
## cultivarSM
cultSM = as.data.frame(datTraits$cultivarSM)
# rename colname
names(cultSM) = "cultSM"
# Add the weight to existing module eigengenes 
MEscultSM <- orderMEs(cbind(MEs,cultSM))
# Plot the relationships among the eigengenes and the trait 
pdf(file = "Plots/eigengene.relationships.cultivarSM_08Apr16.pdf", wi = 15, he = 12) 
plotEigengeneNetworks(MEscultSM,"",marDendro=c(0,5.5,0.5,2), marHeatmap=c(4.75,7,0.1,0.1),cex.lab=0.8,xLabelsAngle=90)
dev.off()
# save eigengen trait data
save(colorWhite, colorRed, 
	brix20, brix22, brix24, brix26, 
	cultCS, cultCF, cultME, cultPN, cultCD, cultSB, cultSM,
	MEsRed, MEsWhite,
	MEsbrix20, MEsbrix22, MEsbrix24, MEsbrix26, 
	MEscultCS, MEscultCF, MEscultME, MEscultPN, MEscultCD, MEscultSB, MEscultSM,
	file = "eigengene.relationship.traits_08Apr16.RData");

#===============================================================================
# Use the trait to define a gene significance variable 
# see Quantifying module-trait associations in wgcna_analysis.R section
#===============================================================================
GS.RedColor = numbers2colors(as.numeric(bicor(datExpr, colorRed, use = "p")), signed = TRUE);
GS.WhiteColor = numbers2colors(as.numeric(bicor(datExpr, colorWhite, use = "p")), signed = TRUE);
GS.brix20Color = numbers2colors(as.numeric(bicor(datExpr, brix20, use = "p")), signed = TRUE);
GS.brix22Color = numbers2colors(as.numeric(bicor(datExpr, brix22, use = "p")), signed = TRUE);
GS.brix24Color = numbers2colors(as.numeric(bicor(datExpr, brix24, use = "p")), signed = TRUE);
GS.brix26Color = numbers2colors(as.numeric(bicor(datExpr, brix26, use = "p")), signed = TRUE);
GS.cultCSColor = numbers2colors(as.numeric(bicor(datExpr, cultCS, use = "p")), signed = TRUE);
GS.cultCFColor = numbers2colors(as.numeric(bicor(datExpr, cultCF, use = "p")), signed = TRUE);
GS.cultMEColor = numbers2colors(as.numeric(bicor(datExpr, cultME, use = "p")), signed = TRUE);
GS.cultPNColor = numbers2colors(as.numeric(bicor(datExpr, cultPN, use = "p")), signed = TRUE);
GS.cultCDColor = numbers2colors(as.numeric(bicor(datExpr, cultCD, use = "p")), signed = TRUE);
GS.cultSBColor = numbers2colors(as.numeric(bicor(datExpr, cultSB, use = "p")), signed = TRUE);
GS.cultSMColor = numbers2colors(as.numeric(bicor(datExpr, cultSM, use = "p")), signed = TRUE)

pdf(file = "Plots/geneDendro.mod-trait.relationships_08Apr16.pdf", wi = 15, he = 9)
plotDendroAndColors(geneTree, cbind(mergedColors,
	GS.RedColor, GS.WhiteColor,
	GS.brix20Color, GS.brix22Color, GS.brix24Color, GS.brix26Color,
	GS.cultCSColor, GS.cultCFColor, GS.cultMEColor, GS.cultPNColor, GS.cultCDColor, GS.cultSBColor, GS.cultSMColor),
                    c("Merged dynamic", 
                    	"Red skin", "White skin",
						"20 Brix", "22 Brix", "24 Brix", "26 Brix",
						"Cabernet Sauvignon", "Cabernet Franc", "Merlot", "Pinot Noir", 
						"Chardonnay", "Sauvignon Blanc", "Semillon"),
                    dendroLabels = FALSE, hang = 0.03, marAll = c(0.5, 5.5, 3, 0),,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

