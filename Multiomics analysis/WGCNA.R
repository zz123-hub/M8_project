library(WGCNA)
library(dynamicTreeCut)
library(fastcluster)
options(stringsAsFactors = FALSE)
femData = read.csv("abundance.csv") 
dim(femData)
names(femData)
datExpr0 = as.data.frame(t(femData[, -c(1)])) 
names(datExpr0) = femData$ID
rownames(datExpr0) = names(femData)[-c(1)]
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

datExpr0 <- log2(datExpr0+1)
datExpr <- datExpr0
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,48) 
par(cex = 0.6);
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

allTraits = read.csv("clinic.csv");
dim(allTraits)
names(allTraits)
femaleSamples = rownames(datExpr);
femaleSamples
traitRows = match(femaleSamples, allTraits$IDX);
datTraits = allTraits[traitRows, -1];
datTraits
rownames(datTraits) = allTraits[traitRows, 1];
datTraits
collectGarbage()
sampleTree2 = hclust(dist(datExpr), method ="average")

sizeGrWindow(48,20)
par(cex = 0.5);
par(mar = c(0,0,0,0))
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels =names(datTraits),
                    main ="Sample heatmap")# 
enableWGCNAThreads() 
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow =c(1,2));	
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="SoftThreshold(power)",ylab="Scale,signedR^2",type="n",
     main =paste("Scaleindependence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="SoftThreshold(power)",ylab="MeanConnectivity", type="n",
     main =paste("Mean"))
text(sft$fitIndices[,1], sft$fitIndices[,5],labels=powers, cex=cex1,col="red")
sft$powerEstimate

adjacency = adjacency (datExpr, power = 18, type = "signed", 
                       corFnc = "bicor", 
                       corOptions = "use = 'pairwise.complete.obs'")

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
min.ModuleSize = 5

geneTree = hclust(as.dist(dissTOM), method = "complete");

moduleLabels1 = cutreeDynamic (dendro = geneTree, distM = dissTOM, method = "hybrid", 
                               deepSplit = 2, pamRespectsDendro = F, minClusterSize = min.ModuleSize)
table(moduleLabels1)

merge = mergeCloseModules (datExpr, moduleLabels1, corFnc = "cor", 
                           corOptions = list (use = 'p', method = 'spearman'), cutHeight = 0.1)
moduleColorsMeta = merge$colors
names (moduleColorsMeta) = colnames (datExpr)
MEsMeta = merge$newMEs
rownames (MEsMeta) = rownames (datExpr)
write.csv(MEsMeta,file = "indexMEs.csv")
write.csv(moduleColorsMeta,file = "indexcolors.csv")
nGenes =ncol(datExpr);
nSamples =nrow(datExpr)
moduleTraitCor =cor(MEsMeta, datTraits, method="spearman");	
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
write.csv(moduleTraitCor,file="indexcolorr.csv")
write.csv(moduleTraitPvalue,file = "indexcolorp.csv"