data_path = './10.3389/fimmu.2024.1387311/data/'
load(paste0(data_path, 'rep_01_3_dp.RData'))
load(paste0(data_path, 'rep_02_deg.RData'))

exp = merged_exp


gene_up = deg$ENTREZID[deg$change == 'up']
gene_down = deg$ENTREZID[deg$change == 'down']
gene_diff = c(gene_up, gene_down)

library(clusterProfiler)
library(stringr)
library(ggplot2)
library(ggplotify)

# BiocManager::install('pathfindR')
library(pathfindR)

# BiocManager::install('DOSE')
library(DOSE)

# BiocManager::install('enrichplot')
library(enrichplot)


# ====== A: GO Enrichment Data ======

ego_file = paste0(data_path, 'rep_03_EGO.Rdata')
if (!file.exists(ego_file)) {
  print(Sys.time())
  print('ego')
  ego = enrichGO(
    gene = gene_diff,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    readable = T
  )
  print(Sys.time())
  print('ego_BP')
  ego_BP = enrichGO(
    gene = gene_diff,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    readable = T
  )
  print(Sys.time())
  print('Done.')
  save(ego, ego_BP, file = ego_file)
} else {
  load(ego_file)
}

# barplot(ego)
# dotplot(ego)



# ====== A: GO Enrichment Plot ======

ego_plot = dotplot(
  ego,
  split = "ONTOLOGY", font.size = 6, showCategory = 5,
  x = "Count",
  title = "GO Enrichment",
  size = "Count"
) +
  
  # geom_point(aes(size = p.adjust), alpha = 0.6) +
  geom_point(aes(size=Count,color=-log10(p.adjust))) +
  # https://github.com/YuLab-SMU/enrichplot/issues/130#issuecomment-1236270194

  # facet_grid(ONTOLOGY~., space = "free_y")
  # "space" is wrong
  facet_grid(ONTOLOGY ~ ., scales = "free_y") +
  
  # scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
  scale_y_discrete() +
  
  theme(axis.text.y = element_text(size = 6))  # +

  # scale_size_continuous(range = c(1, 10)) +

  # theme(legend.position = "none")
  # theme_bw()

ego_plot

plot_a = as.ggplot(ego_plot)

# ====== B: KEGG Enrichment Data ======

# COMMENT OUT HERE
if (FALSE) {

kegg_file = paste0(data_path, 'rep_03_KEGG.Rdata')
if (!file.exists(kegg_file)) {
  kegg_up = enrichKEGG(gene = gene_up, organism = 'hsa')
  kegg_down = enrichKEGG(gene = gene_down, organism = 'hsa')
  kegg_diff = enrichKEGG(gene = gene_diff, organism = 'hsa')
  save(kegg_up, kegg_down, kegg_diff, file = kegg_file)
} else {
  load(kegg_file)
}

table(kegg_up@result$p.adjust<0.05)
table(kegg_down@result$p.adjust<0.05)
table(kegg_diff@result$p.adjust<0.05)


kegg_up =
  kegg_up@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(group = 1)

kegg_down =
  kegg_down@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(group = -1)

kegg_diff =
  kegg_diff@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(group = 0)

}


pfr_file = paste0(data_path, 'rep_03_pfr.Rdata')
if (!file.exists(pfr_file)) {
  deg_for_pfr = data.frame(deg$symbol, deg$logFC, deg$adj.P.Val)
  colnames(deg_for_pfr) = c('symbol', 'logFC', 'adj.P.Val')
  
  pfr_df = run_pathfindR(deg_for_pfr)
  save(pfr_df, file = pfr_file)
} else {
  load(pfr_file)
}

# ====== B: KEGG Enrichment Plot ======

# Bubble diagram showing the KEGG enrichment analysis of DEGs
# x: enrichment factor
# y: KEGG pathway

# COMMENT OUT HERE

if (FALSE) {

kegg_plot = kegg_diff %>%
  ggplot(aes(x = Count, y = Description, size = Count, color = -log10(p.adjust))) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +
  theme_minimal() +
  # theme(legend.position = "none") +
  labs(
    x = "Count",
    y = "Description",
    title = "KEGG Enrichment"
  )
kegg_plot

}

pfr_plot = enrichment_chart(pfr_df)
pfr_plot

plot_b = as.ggplot(pfr_plot)


# ====== C: GSEA Analysis ======

# data(geneList, package="DOSE")
library(DOSE)
# data(geneList)

geneList = deg$logFC
names(geneList) = deg$ENTREZID
geneList = sort(geneList, decreasing = T)

# 共同基因
cnetplot(ego, categorySize = "pvalue", foldChange = geneList)
cnetplot(ego, showCategory = 3, foldChange = geneList)

# kk_gse = gseKEGG(geneList, organism = 'hsa', verbose = F)
# down_kegg = kk_gse[kk_gse$pvalue < 0.05 & kk_gse$enrichmentScore < 0,]
# up_kegg = kk_gse[kk_gse$pvalue < 0.05 & kk_gse$enrichmentScore > 0,]
# 
# g2 = kegg_plot(up_kegg, down_kegg)
# g2

# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html

gse_file = paste0(data_path, 'rep_03_gse.Rdata')
if (!file.exists(gse_file)) {
  gene_set = gseDO(geneList)
  save(gene_set, file = gse_file)
} else {
  load(gse_file)
}

gsea_plot = gseaplot2(gene_set, geneSetID = 1:5)
gsea_plot

plot_c = as.ggplot(gsea_plot)


# ====== D: GSVA ======

# A gene set variation analysis (GSVA) was performed
# using the GSVA R package to calculate a normalized enrichment score
# under the background of the hallmark gene set (c2.cp.kegg.v7.2)

# BiocManager::install('GSVA')
library(GSVA)
library(GSVAdata)

# first construct the gene set for GSVA
# gset_raw = clusterProfiler::read.gmt(paste0(data_path, '/c2.cp.kegg.v7.5.1.entrez.gmt'))

# https://www.biostars.org/p/206323/
gset = getGmt(paste0(data_path, 'c2.cp.kegg.v7.5.1.entrez.gmt'))

exp2id = exp
rownames(exp2id) = deg$ENTREZID[match(rownames(exp), deg$symbol)]
# exp2id = as.matrix(exp2id)

deg_param = gsvaParam(exp2id, gset, minSize=10, maxSize=500)
deg_gsva = gsva(deg_param)

library(limma)

mod <- model.matrix(~ factor(group))
colnames(mod) <- c("Control", "RA")
fit <- lmFit(deg_gsva, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.01)

tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]

# DEpwys_es <- exprs(deg_gsva[DEpwys, ])
DEpwys_es = deg_gsva[DEpwys, ]
colorLegend <- c("darkgreen", "forestgreen")
names(colorLegend) <- c("Control", "RA")
# sample.color.map <- colorLegend[pData(deg_gsva)[, "subtype"]]
sample.color.map <- colorLegend[group]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")
heatmap(
  DEpwys_es, ColSideColors=sample.color.map,
  # xlab="samples",
  # ylab="Pathways",
  margins=c(2, 20),
  # labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "", rownames(DEpwys_es))), 1, 35),
  labRow = rownames(DEpwys_es),
  labCol="",
  scale="row",
  # Colv=as.dendrogram(sampleClustering),
  Colv=NA,
  Rowv=as.dendrogram(geneSetClustering)
)
legend("right", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")

dev.off()

library(pheatmap)

annotation_col = data.frame(group = group)
rownames(annotation_col) = colnames(DEpwys_es)

gsva_heat_map = pheatmap(
  DEpwys_es,
  color = colorRampPalette(c("darkslateblue", "white", "firebrick"))(100),
  annotation_col = annotation_col,
  show_rownames = T,
  show_colnames = F,
  cluster_rows = T,
  cluster_cols = F,
  scale = "row",
  breaks = seq(-3, 3, length.out = 100),
  fontsize = 6
)
gsva_heat_map

plot_d = as.ggplot(gsva_heat_map)


# ====== FINAL PLOT ======

plot_comb = (plot_a + plot_c) / (plot_b + plot_d)
plot_comb
