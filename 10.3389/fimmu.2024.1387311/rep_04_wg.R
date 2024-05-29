# ====== Tutorials ======

# https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA
# https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/wgcna.html


# ====== WGCNA analysis ======

data_path = './10.3389/fimmu.2024.1387311/data/'
load(paste0(data_path, 'rep_01_3_dp.RData'))
load(paste0(data_path, 'rep_02_deg.RData'))

# exp = merged_exp

library(WGCNA)
allowWGCNAThreads()

library(ggplotify)

datExpr = t(merged_exp)

# ====== A: soft-thresholding powers ======

# Analysis of the mean connectivity and scale-free fit index for different soft-thresholding powers
# Where the correlation coefficient is 0.9 and the matching soft-thresholding power is 8, the red line represents this location

# datExpr = t(merged_exp)
powers = c(1:30)

if (!file.exists(paste0(data_path, 'rep_04_wg_sft.RData'))) {

  sft = pickSoftThreshold(
    datExpr,
    #blockSize = 30,
    powerVector = powers,
    verbose = 5
  )

  save(sft, file = paste0(data_path, 'rep_04_wg_sft.RData'))

} else {
  load(paste0(data_path, 'rep_04_wg_sft.RData'))
}

# https://stackoverflow.com/questions/49602032/how-to-convert-plot-to-ggplot-in-r

plot_a = as.ggplot(
  function() {

par(mfrow = c(1,2))  # Set up the page as two plots side by side
cex1 = 0.9

plot(
  sft$fitIndices[, 1],
  # -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1, col = "red"
)

  }
)
plot_a

# ====== B: Gene dendrogram ======

# The cluster dendrogram of the top 25% of genes median absolute deviations.
# Each hue in the graphic below corresponds to a co-expression module,
# and each branch in the figure represents a single gene.

if (!file.exists(paste0(data_path, 'rep_04_wg_den.RData'))) {

softPower = 18
adjacency = adjacency(datExpr, power = softPower, type = "signed")  # specify network type
# head(adjacency)

TOM = TOMsimilarity(adjacency, TOMType="signed")  # specify network type
dissTOM = 1-TOM

# Generate a clustered gene tree
library(flashClust)
geneTree = flashClust(as.dist(dissTOM), method="average")

save(
  adjacency, TOM, dissTOM, geneTree,
  file = paste0(data_path, 'rep_04_wg_den.RData')
)

} else {
  load(paste0(data_path, 'rep_04_wg_den.RData'))
}



# load here

library(flashClust)

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

# plot_b = as.ggplot(
#   function() {
#     plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
#   }
# )
# plot_b


# ====== C: Module-trait relationships ======

#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#Recalculate MEs with color labels

if (!file.exists(paste0(data_path, 'rep_04_wg_cl.RData'))) {

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEDissThres = 0.0
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
moduleColors = mergedColors

save(
  moduleColors,
  file = paste0(data_path, 'rep_04_wg_cl.RData')
)

} else {

load(paste0(data_path, 'rep_04_wg_cl.RData'))
  
}

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

load(paste0(data_path, 'rep_01_3_dp.RData'))

datTraits = data.frame(group = group)
rownames(datTraits) = rownames(datExpr)
datTraits$RA <- ifelse(datTraits$group == "RA", 1, 0)
datTraits$NC <- ifelse(datTraits$group == "Control", 1, 0)
datTraits$group <- NULL

moduleTraitCor = cor(MEs, datTraits, use= "p")
nSamples = nrow(datExpr)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)


library(ggplotify)

plot_c = as.ggplot(
  function() {

par(mar= c(6, 8.5, 3, 3))

labeledHeatmap(
  Matrix=moduleTraitCor,
  xLabels=names(datTraits),
  yLabels=names(MEs),
  ySymbols=names(MEs),
  colorLabels=FALSE,
  colors=blueWhiteRed(50),
  textMatrix=textMatrix,
  setStdMargins=FALSE,
  cex.text=0.5,
  zlim=c(-1,1),
  main="Module-trait relationships"
)

  }
)

plot_c

# ====== D: Module membership vs. gene significance ======

# Do not restart R session
# Follow part C

# https://github.com/jmzeng1314/my_WGCNA

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


# Luminal = as.data.frame(design[,3]);
# names(Luminal) = "Luminal"
geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(datTraits), sep="")
names(GSPvalue) = paste("p.GS.", names(datTraits), sep="")

# module = "salmon"
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;

plot_d = as.ggplot(
  function() {

par(mfrow = c(1,1));
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for Luminal",
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
)

# add a line to the plot
abline(a = 1.6e-46, b = 0.73, col = "red")

  }
)
plot_d


# ====== E: Extract Hub Genes ======

hub_genes = chooseTopHubInEachModule(
  datExpr,
  moduleColors,
)

if (!file.exists(paste0(data_path, 'rep_04_wg_hub.RData'))) {
  
  save(
    hub_genes,
    file = paste0(data_path, 'rep_04_wg_hub.RData')
  )
  
} else {
  load(paste0(data_path, 'rep_04_wg_hub.RData'))
}
