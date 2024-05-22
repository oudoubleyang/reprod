data_path = './10.3389/fimmu.2024.1387311/data'
load(paste0(data_path, '/rep_02_dp.RData'))

exp = merged_exp
# remove(merged_exp)
exp_data = as.data.frame(t(exp))


# ====== DEG ======

library(limma)

# design = model.matrix(~0 + group)
design = model.matrix(~group)
fit = lmFit(exp, design)
fit = eBayes(fit)
deg = topTable(fit, coef = 2, number = Inf)

library(dplyr)
deg = mutate(deg, symbol = rownames(deg))
head(deg)

logFC_t = 1
P_Value_t = 0.05
down_cond = (deg$P.Value < P_Value_t) & (deg$logFC < -logFC_t)
up_cond   = (deg$P.Value < P_Value_t) & (deg$logFC >  logFC_t)
deg = mutate(
  deg,
  change = ifelse(
    down_cond, "down",
    ifelse(
      up_cond, "up",
      "normal"
    )
  )
)
table(deg$change)


library(clusterProfiler)
library(org.Hs.eg.db)
s2e = bitr(
  deg$symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
dim(deg)

deg = inner_join(deg, s2e, by = c("symbol" = "SYMBOL"))


# library(dplyr)
library(ggplot2)
# library(pheatmap)

exp_data = deg[!duplicated(deg$symbol),]
p = 
  ggplot(exp_data, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, aes(color = change)) +
  ylab("-log10(P.Value)") +
  # scale_color_manual(values = c("blue", "grey", "red")) +
  scale_color_manual(values = c("slategray", "sandybrown", "dodgerblue")) +
  geom_vline(xintercept = c(-logFC_t, logFC_t), lty=4, col="black", lwd=0.8) +
  geom_hline(yintercept = -log10(P_Value_t), lty=4, col="black", lwd=0.8) +
  theme_bw()
p

# label 20 most significant genes
# top_20 = head(exp_data[order(exp_data$P.Value),], 20)

# no, 10 for up, 10 for down
up_10 = head(exp_data[exp_data$change == 'up',][order(exp_data$P.Value),], 10)
down_10 = head(exp_data[exp_data$change == 'down',][order(exp_data$P.Value),], 10)
top_20 = rbind(up_10, down_10)

volcano_plot = 
  p +
  ggrepel::geom_text_repel(
    data = top_20,
    aes(label = symbol),
    fill = "white"
  )
volcano_plot

volcano_plot = 
  p +
  ggrepel::geom_label_repel(
    data = top_20,
    aes(label = symbol),
    fill = "white"
  )
volcano_plot

# ====== PCA Heat Map ======

# pick 50 most variable genes
# cg = names(tail(sort(apply(exp, 1, sd)), 50))

# pick 25 for most up, 25 for most down
up_25 = head(exp_data[exp_data$change == 'up',][order(exp_data$P.Value),], 25)
down_25 = head(exp_data[exp_data$change == 'down',][order(exp_data$P.Value),], 25)
cg = c(up_25$symbol, down_25$symbol)


library(FactoMineR)
library(factoextra)

exp_data = as.data.frame(t(merged_exp))

data_pca = PCA(exp_data, graph = F)
pca_plot = fviz_pca_ind(
  data_pca,
  geom.ind = "point",
  col.ind = group,
  # palette = c("#00AFBB", "#E7B800"),
  addEllipses = T,
  legend.title = "Group"
)
pca_plot
# bad


# ====== Heat Map ======

library(pheatmap)

annotation_col = data.frame(group = group)
n = merged_exp[cg,]
rownames(annotation_col) = colnames(n)

pca_heat_map = pheatmap(
  n,
  annotation_col = annotation_col,
  show_rownames = T,
  show_colnames = F,
  cluster_rows = F,
  cluster_cols = T,
  scale = "row",
  breaks = seq(-3, 3, length.out = 100)
)
pca_heat_map

# combine volcano plot and heat map into one figure

library(ggplotify)

volcano_plot + as.ggplot(pca_heat_map)


# save

save(
  deg,
  file = paste0(data_path, '/rep_03_deg.RData')
)
