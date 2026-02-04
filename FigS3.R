############## FigS3_CD8_T_cells ###############################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(Seurat)
  library(tidydr)
  library(tidyverse)
  library(scRNAtoolVis)
  library(cowplot)
  library(RColorBrewer)
  library(reshape2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggpubr)
  library(monocle)
  library(kernlab)
  library(ComplexHeatmap)
  library(circlize)
})

fig_dir <- "./FigS3"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path("utils", "plot_GP_pseudotime.R"))
source(file.path("utils", "plot_GP_pseudotime_group.R"))

# Input (not tracked; place under ./data/)
# - scRNA: `HCC_T_cells.rds` (preferred) or `HCC_scRNA.rds` (will subset CD4/CD8 T cells)
# - Optional: `HCC_CD8T_SCENIC.rds` (Seurat object with a SCENIC regulon AUC assay)
# - Optional (for FigS3i): `CD8T_top5_TFs_data.csv` (column: gene) and `TF_GeneNum_score.csv` (columns: TF, TF_genenum)
# - Optional (for FigS3j): `CD8T_top5_TF_Genes_score.csv` (columns: TF, Gene, Score)

################################################################################
# Load T cells and subset CD8 T cells
################################################################################

Hep <- readRDS(file.path(data_dir, "HCC_T_cells.rds"))

cd8_subtypes <- c(
  "CD8T_C1_GZMK", "CD8T_C2_PDCD1", "CD8T_C3_APOA2", "CD8T_C4_C1QB", "CD8T_C5_IFIT3", "MAIT"
)

CD8 <- subset(Hep, subcelltype %in% cd8_subtypes)
CD8@meta.data$subcelltype <- factor(CD8@meta.data$subcelltype, levels = cd8_subtypes)

cols_cd8 <- c("#36489E", "#DA2917", "#C6EE8F", "#B2DCEE", "#A680B9", "#28B461")
cols_cd8 <- setNames(cols_cd8, cd8_subtypes)

################################################################################
# FigS3a: Differentiation-associated genes dot plot across CD8 subsets
################################################################################

diff_genes <- c("IL7R", "TCF7", "CCR7", "FCGR3A", "GZMH", "GZMK", "PDCD1", "LAYN", "HAVCR2")

pS3a <- DotPlot(
  CD8,
  group.by = "subcelltype",
  features = diff_genes,
  cols = c("white", "#980201")
) +
  RotatedAxis() +
  ggtitle("") +
  theme(
    panel.border = element_rect(color = "black"),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
ggsave(file.path(fig_dir, "FigS3a_CD8_differentiation_genes_dotplot.pdf"), plot = pS3a, width = 9.5, height = 5.5)

################################################################################
# FigS3b: Downregulated GO BP pathways in AFP_Pos vs AFP_Neg (total CD8)
################################################################################

Idents(CD8) <- "AFP_status_group"
deg_cd8 <- FindMarkers(
  CD8,
  group.by = "AFP_status_group",
  ident.1 = "AFP_Pos",
  ident.2 = "AFP_Neg",
  min.pct = 0.1,
  logfc.threshold = 0,
  only.pos = FALSE
)

deg_cd8$gene <- rownames(deg_cd8)
write.csv(deg_cd8, file = file.path(fig_dir, "CD8_all_AFP_Pos_vs_AFP_Neg_FindMarkers.csv"), row.names = FALSE)

down_genes <- deg_cd8 %>%
  filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  pull(gene) %>%
  unique()

ego_down <- enrichGO(
  gene = down_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_df <- as.data.frame(ego_down@result)
if (nrow(ego_df) == 0) {
  stop("No GO BP enrichment results for downregulated genes in CD8 (AFP_Pos vs AFP_Neg).")
}

ego_df <- ego_df %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

pS3b <- ggplot(ego_df, aes(x = -log10(p.adjust), y = Description)) +
  geom_col(fill = "#1965B0") +
  labs(x = "-log10(adj.P)", y = "", title = "Downregulated GO BP in AFP_Pos (CD8 total)") +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(fig_dir, "FigS3b_CD8_down_GO_BP_barplot.pdf"), plot = pS3b, width = 10, height = 5.5)

################################################################################
# FigS3c: Cytotoxicity and exhaustion scores in total CD8 by AFP group
################################################################################

score_cols_S3c <- c("Cytotoxicity", "Exhaustion")

df_scores <- CD8@meta.data %>%
  dplyr::select(AFP_status_group, dplyr::all_of(score_cols_S3c)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(score_cols_S3c), names_to = "Signature", values_to = "Score")

pS3c <- ggplot(df_scores, aes(x = AFP_status_group, y = Score, fill = AFP_status_group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
  facet_wrap(~Signature, scales = "free_y") +
  scale_fill_manual(values = c(AFP_Pos = "#db5024", AFP_Neg = "#25acda")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(x = "", y = "Score") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none", strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "FigS3c_CD8_cytotoxicity_exhaustion_violin.pdf"), plot = pS3c, width = 8.5, height = 4.5)

################################################################################
# FigS3d: UMAP of CD8 T cell subsets
################################################################################

pS3d <- DimPlot(CD8, group.by = "subcelltype", reduction = "umap", cols = cols_cd8, raster = FALSE) +
  ggtitle("CD8 T cell subsets") +
  theme(panel.grid = element_blank())
ggsave(file.path(fig_dir, "FigS3d_CD8_subsets_umap.pdf"), plot = pS3d, width = 6.5, height = 4.8)

################################################################################
# FigS3e–f: Monocle2 pseudotime trajectory (by subtype and by pseudotime)
################################################################################

# Optional: downsample for trajectory to keep runtime reasonable
max_cells_traj <- 8000
if (ncol(CD8) > max_cells_traj) {
  set.seed(1)
  CD8_traj <- CD8[, sample(colnames(CD8), max_cells_traj)]
} else {
  CD8_traj <- CD8
}

CD8_traj <- FindVariableFeatures(CD8_traj, nfeatures = 2000)
ordering_genes <- VariableFeatures(CD8_traj)

expr_matrix <- as(as.matrix(CD8_traj@assays$RNA@counts), "sparseMatrix")
p_data <- CD8_traj@meta.data
f_data <- data.frame(gene_short_name = rownames(CD8_traj), row.names = rownames(CD8_traj))

pd <- new("AnnotatedDataFrame", data = p_data)
fd <- new("AnnotatedDataFrame", data = f_data)

cds <- newCellDataSet(
  expr_matrix,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)

pdf(file.path(fig_dir, "FigS3e_CD8_pseudotime_by_subtype.pdf"), width = 6.5, height = 5.8, useDingbats = FALSE)
plot_cell_trajectory(cds, color_by = "subcelltype", size = 0.8, show_backbone = TRUE) +
  scale_color_manual(values = cols_cd8) +
  theme(legend.title = element_blank())
dev.off()

pdf(file.path(fig_dir, "FigS3f_CD8_pseudotime_by_value.pdf"), width = 6.5, height = 5.8, useDingbats = FALSE)
plot_cell_trajectory(cds, color_by = "Pseudotime", size = 0.8, show_backbone = TRUE) +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd"))
dev.off()

################################################################################
# FigS3g–h: Gaussian process regression trends along pseudotime
################################################################################

plot_data <- pData(cds)

sigma <- 0.8

score_cols_S3gh <- c("Naive", "Cytotoxicity", "Activation/Effector function", "Exhaustion")

pS3g <- plot_GP_pseudotime(
  plot_data,
  pseudotime_col = "Pseudotime",
  sigma = sigma,
  y_col = score_cols_S3gh,
  line_colors = c("#33A02C", "#2c7bb6", "#55B1B1", "#d7191c")
)
ggsave(file.path(fig_dir, "FigS3g_GP_pseudotime_4_signatures.pdf"), plot = pS3g, width = 7.5, height = 4.6)

pS3h <- plot_GP_pseudotime_group(
  plot_data,
  pseudotime_col = "Pseudotime",
  y_col = "Cytotoxicity",
  group_col = "AFP_status_group",
  sigma = sigma,
  title = "",
  y_title = "Cytotoxicity score (GP mean ± 95% CI)",
  line_colors = c(AFP_Pos = "#db5024", AFP_Neg = "#25acda"),
  ribbon_alpha = 0.18
)
ggsave(file.path(fig_dir, "FigS3h_GP_pseudotime_cytotoxicity_by_AFP.pdf"), plot = pS3h, width = 7.2, height = 4.6)

################################################################################
# FigS3i: SCENIC top 5 TF activity bubble heatmap per CD8 subset
################################################################################

scenic_file <- file.path(data_dir, "HCC_CD8T_SCENIC.rds")
seurat.data <- readRDS(scenic_file)
DefaultAssay(seurat.data) <- "scenic_pos_genenum"

cd8T_top5_tf_file <- file.path(data_dir, "CD8T_top5_TFs_data.csv")
cd8T_top5_tf_data <- read.csv(cd8T_top5_tf_file, check.names = FALSE)
sel_TFs <- cd8T_top5_tf_data$gene

tf_genenum_file <- file.path(data_dir, "TF_GeneNum_score.csv")
TF_GeneNum <- read.csv(tf_genenum_file, check.names = FALSE)
rownames(TF_GeneNum) <- TF_GeneNum$TF
markergenes <- TF_GeneNum[sel_TFs, "TF_genenum"]
markergenes <- unique(markergenes[!is.na(markergenes)])

plot_order <- c("CD8T_C1_GZMK", "CD8T_C2_PDCD1", "CD8T_C3_APOA2", "CD8T_C4_C1QB", "CD8T_C5_IFIT3", "MAIT")
sub_seu <- subset(seurat.data, subcelltype %in% plot_order)
sub_seu$subcelltype <- factor(sub_seu$subcelltype, levels = plot_order)

DotPlot2 <- jjDotPlot(
  object = sub_seu,
  id = "subcelltype",
  assay = "scenic_pos_genenum",
  gene = markergenes,
  cluster.order = rev(plot_order),
  ytree = FALSE,
  legend.position = "top",
  base_size = 20,
  x.text.angle = 45,
  x.text.vjust = 1
) +
  theme(panel.grid = element_blank())

ggsave(file.path(fig_dir, "FigS3i_TFs_Clusters_Markers_DotPlot_gene_num.pdf"), plot = DotPlot2, width = 16, height = 15)

################################################################################
# FigS3j: Dynamic gene expression heatmap along pseudotime (Z-score)
################################################################################

tf_targets_file <- file.path(data_dir, "CD8T_top5_TF_Genes_score.csv")
tf_targets <- read.csv(tf_targets_file, check.names = FALSE)

genes_add <- c("GZMK", "GZMH", "PDCD1", "TIGIT", "TOX", "IFIT3", "CD69", "NKG7", "KLRD1", "FOS", "SLC4A10", "STMN1")

Filted_TFGenes_Score <- NULL
if (!is.null(tf_targets) && all(c("TF", "Gene", "Score") %in% colnames(tf_targets))) {
  Filted_TFGenes_Score <- tf_targets %>%
    dplyr::group_by(TF) %>%
    dplyr::mutate(any_above = any(Score > 3)) %>%
    dplyr::group_modify(~{
      if (.x$any_above[1]) {
        .x %>% dplyr::filter(Score > 3)
      } else {
        .x %>% dplyr::slice_max(order_by = Score, n = 1, with_ties = FALSE)
      }
    })
}

genes_show <- genes_add
if (!is.null(Filted_TFGenes_Score) && nrow(Filted_TFGenes_Score) > 0) {
  genes_show <- unique(c(genes_show, Filted_TFGenes_Score$Gene))
}
genes_tf <- character(0)
if (!is.null(tf_targets) && ("Gene" %in% colnames(tf_targets))) {
  genes_tf <- intersect(genes_show, unique(tf_targets$Gene))
}

genes_show <- intersect(genes_show, rownames(exprs(cds)))

pt <- as.numeric(pData(cds)$Pseudotime)
ord <- order(pt)
cell_order <- rownames(pData(cds))[ord]

expr_mat <- exprs(cds)[genes_show, cell_order, drop = FALSE]
expr_mat <- log2(expr_mat + 1)

max_cells_heat <- 3000
if (ncol(expr_mat) > max_cells_heat) {
  idx <- round(seq(1, ncol(expr_mat), length.out = max_cells_heat))
  expr_mat <- expr_mat[, idx, drop = FALSE]
  cell_order <- colnames(expr_mat)
  pt <- pt[ord][idx]
}

expr_z <- t(scale(t(as.matrix(expr_mat))))
expr_z[is.na(expr_z)] <- 0

anno <- data.frame(
  Pseudotime = pt,
  subcelltype = pData(cds)[cell_order, "subcelltype"],
  AFP_status_group = pData(cds)[cell_order, "AFP_status_group"]
)
anno$subcelltype <- factor(anno$subcelltype, levels = cd8_subtypes)

ha <- HeatmapAnnotation(
  df = anno[, c("subcelltype", "AFP_status_group")],
  col = list(
    subcelltype = cols_cd8,
    AFP_status_group = c(AFP_Pos = "#db5024", AFP_Neg = "#25acda")
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)

row_label_col <- rep("black", nrow(expr_z))
names(row_label_col) <- rownames(expr_z)
if (length(genes_tf) > 0) {
  row_label_col[names(row_label_col) %in% genes_tf] <- "red3"
}

ht_S3j <- Heatmap(
  expr_z,
  name = "Z",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(col = row_label_col, fontsize = 7),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  top_annotation = ha,
  column_title = "",
  row_title = ""
)

df_den <- data.frame(
  Pseudotime = as.numeric(pData(cds)$Pseudotime),
  subcelltype = pData(cds)$subcelltype
)
df_den$subcelltype <- factor(df_den$subcelltype, levels = cd8_subtypes)

p_den <- ggplot(df_den, aes(Pseudotime, colour = subcelltype, fill = subcelltype)) +
  geom_density(bw = 0.5, linewidth = 0.9, alpha = 0.35) +
  scale_color_manual(values = cols_cd8) +
  scale_fill_manual(values = cols_cd8) +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )

pdf(file.path(fig_dir, "FigS3j_pseudotime_density_and_heatmap.pdf"), width = 8.5, height = 8.5, useDingbats = FALSE)
print(p_den)
draw(ht_S3j, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
