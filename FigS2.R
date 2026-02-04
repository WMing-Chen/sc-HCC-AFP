############## FigS2_Tumor_SCENIC #############################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

fig_dir <- "./FigS2"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cluster_ids <- c(1, 3, 6)
cluster_colors <- setNames(c("#C6EE8F", "#B2DCEE", "#FA7E4D"), paste0("Cluster", cluster_ids))
group_colors <- setNames(c("#db5024", "#25acda"), c("AFP_Pos", "AFP_Neg"))

col_order <- unlist(lapply(cluster_ids, function(k) c(paste0("C", k, "_AFP_Pos"), paste0("C", k, "_AFP_Neg"))))

top5_tbl <- read.csv(file.path(data_dir, "Tumor_top5_TFs_data.csv"), check.names = FALSE)
top5_tbl <- top5_tbl %>%
  filter(subcelltype %in% col_order) %>%
  mutate(subcelltype = factor(subcelltype, levels = col_order, ordered = TRUE)) %>%
  arrange(subcelltype)

tf_order <- unique(gsub("\\(\\+\\)|\\(\\-\\)", "", top5_tbl$gene))

column_ha <- HeatmapAnnotation(
  AFP_Group = rep(c("AFP_Pos", "AFP_Neg"), length(cluster_ids)),
  Clusters = rep(paste0("Cluster", cluster_ids), each = 2),
  col = list(Clusters = cluster_colors, AFP_Group = group_colors),
  simple_anno_size = unit(3, "mm")
)

################################################################################
# FigS2a: TF activity heatmap (RSS Z-score) for C1/C3/C6 top5 TFs
################################################################################

rss_df <- read.csv(file.path(data_dir, "Tumor_RSS_results_Celltype_AFP_Group.csv"), row.names = 1, check.names = FALSE)
rownames(rss_df) <- gsub("\\(\\+\\)|\\(\\-\\)", "", rownames(rss_df))

mat_a <- as.matrix(rss_df[tf_order, col_order])
mat_a <- t(scale(t(mat_a)))

col_fun_a <- colorRamp2(c(-1, 0, 1), c("#1A5592", "white", "#B83D3D"))
ht_a <- Heatmap(
  mat_a,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = col_fun_a,
  top_annotation = column_ha,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 12),
  row_names_gp = grid::gpar(fontsize = 9),
  heatmap_legend_param = list(title = "Relative value")
)

pdf(file.path(fig_dir, "FigS2a_Cluster136_AFP_Group_Top5_TFs_activity_heatmap.pdf"), width = 6.5, height = 6, useDingbats = FALSE)
draw(ht_a, padding = unit(c(10, 20, 10, 10), "mm"), heatmap_legend_side = "right")
dev.off()

################################################################################
# FigS2b: Mean expression heatmap of TF genes for C1/C3/C6 top5 TFs
################################################################################

exp_df <- read.csv(
  file.path(data_dir, "Tumor_SCENIC_TF_genes_mean_exp_data_Cluster_AFP_group.csv"),
  row.names = 1,
  check.names = FALSE
)
rownames(exp_df) <- gsub("\\(\\+\\)|\\(\\-\\)", "", rownames(exp_df))

mat_b <- as.matrix(exp_df[tf_order, col_order])
mat_b <- t(scale(t(mat_b)))

col_fun_b <- colorRamp2(c(-2, 0, 2), c("#1A5592", "white", "#B83D3D"))
ht_b <- Heatmap(
  mat_b,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = col_fun_b,
  top_annotation = column_ha,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 12),
  row_names_gp = grid::gpar(fontsize = 9),
  heatmap_legend_param = list(title = "Exp (Z-score)")
)

pdf(file.path(fig_dir, "FigS2b_Cluster136_AFP_Group_Top5_TFs_gene_expression_heatmap.pdf"), width = 6.5, height = 6, useDingbats = FALSE)
draw(ht_b, padding = unit(c(10, 20, 10, 10), "mm"), heatmap_legend_side = "right")
dev.off()

