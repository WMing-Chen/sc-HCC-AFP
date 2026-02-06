############## Fig5_Treg_DC_CXCL10_cor_and_spatial ##############################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(corrplot)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(scRNAtoolVis)
  library(Seurat)
  library(semla)
  library(tidyr)
  library(readxl)
})

setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")

fig_dir <- "./00Figures/Fig5"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Fig5a: Co-occurrence correlation heatmap (Tumor C0/C3/C6 + T + DC subsets)
################################################################################

Tumor <- readRDS(file.path(data_dir, "HCC_Tumor.rds"))
Tcell <- readRDS(file.path(data_dir, "HCC_T_cells.rds"))
DC <- readRDS(file.path(data_dir, "HCC_DC.rds"))

tumor_subtypes <- c("Tumor_C0_FABP1", "Tumor_C3_STMN1", "Tumor_C6_AFP")
t_subtypes <- c(
  "CD4T_C1_GPR183", "CD4T_C2_PLCG2", "CD4T_C3_CXCL13",
  "CD8T_C1_GZMK", "CD8T_C2_PDCD1", "CD8T_C3_APOA2", "CD8T_C4_C1QB", "CD8T_C5_IFIT3",
  "Cycling_T_cells", "T_stress", "Treg", "NKT", "MAIT"
)
dc_subtypes <- c("DC_C1_CD1C", "DC_C2_STMN1", "DC_C3_CLEC9A", "DC_C4_LAMP3", "DC_C5_CXCL10", "DC_C6_APOC3")

tumor_per <- Tumor@meta.data %>%
  count(Sample = orig.ident, subcelltype) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  filter(subcelltype %in% tumor_subtypes) %>%
  select(Sample, subcelltype, prop) %>%
  tidyr::pivot_wider(names_from = subcelltype, values_from = prop, values_fill = 0)

t_per <- Tcell@meta.data %>%
  count(Sample = orig.ident, subcelltype) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  filter(subcelltype %in% t_subtypes) %>%
  select(Sample, subcelltype, prop) %>%
  tidyr::pivot_wider(names_from = subcelltype, values_from = prop, values_fill = 0)

dc_per <- DC@meta.data %>%
  count(Sample = orig.ident, subcelltype) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  filter(subcelltype %in% dc_subtypes) %>%
  select(Sample, subcelltype, prop) %>%
  tidyr::pivot_wider(names_from = subcelltype, values_from = prop, values_fill = 0)

subct_cluster_per <- tumor_per %>%
  full_join(t_per, by = "Sample") %>%
  full_join(dc_per, by = "Sample")
subct_cluster_per[is.na(subct_cluster_per)] <- 0

col_order <- c(tumor_subtypes, t_subtypes, dc_subtypes)
mat_prop <- as.matrix(subct_cluster_per[, col_order])

cor_mat <- cor(mat_prop, method = "spearman")

p_mat <- matrix(NA_real_, ncol = ncol(mat_prop), nrow = ncol(mat_prop))
colnames(p_mat) <- colnames(mat_prop)
rownames(p_mat) <- colnames(mat_prop)
diag(p_mat) <- 0
for (i in seq_len(ncol(mat_prop) - 1)) {
  for (j in (i + 1):ncol(mat_prop)) {
    tmp <- suppressWarnings(cor.test(mat_prop[, i], mat_prop[, j], method = "spearman"))
    p_mat[i, j] <- tmp$p.value
    p_mat[j, i] <- tmp$p.value
  }
}

colnames(cor_mat) <- gsub("^Tumor_C0_FABP1$", "Tumor_C1_FABP1", colnames(cor_mat))
rownames(cor_mat) <- gsub("^Tumor_C0_FABP1$", "Tumor_C1_FABP1", rownames(cor_mat))
colnames(p_mat) <- gsub("^Tumor_C0_FABP1$", "Tumor_C1_FABP1", colnames(p_mat))
rownames(p_mat) <- gsub("^Tumor_C0_FABP1$", "Tumor_C1_FABP1", rownames(p_mat))

col_pal <- c("#2166AC", "#4393C3", "#92C5DE", "white", "#F4A582", "#D6604D", "#B2182B")

pdf(file.path(fig_dir, "Fig5a_Tumor_T_DC_subtype_spearman_correlation_heatmap.pdf"), width = 8.2, height = 7.6, useDingbats = FALSE)
corrplot(
  cor_mat,
  method = "color",
  col = col_pal,
  tl.col = "black",
  tl.cex = 0.7,
  order = "hclust",
  addrect = 4,
  p.mat = p_mat,
  sig.level = 0.05,
  insig = "label_sig",
  pch.cex = 1.1
)
dev.off()

################################################################################
# Fig5b: Cell proportion vs serum AFP (Treg and DC_C5_CXCL10)
################################################################################

surival_data <- read.csv(file.path(data_dir, "surival_data.csv"))
rownames(surival_data) <- surival_data$Sample

df_afp <- subct_cluster_per %>%
  select(Sample, Treg, DC_C5_CXCL10) %>%
  mutate(
    AFP_ug_L = surival_data[Sample, "AFP_μg_L"],
    Log10_AFP = log10(AFP_ug_L)
  )

p5b_1 <- ggscatter(
  df_afp,
  x = "Log10_AFP",
  y = "Treg",
  size = 1.6,
  add = "reg.line",
  add.params = list(color = "#0AA1FF", fill = "#a5dff9", size = 1),
  conf.int = TRUE
) +
  stat_cor(method = "spearman", label.sep = "\n") +
  xlab("Serum AFP (log10)") +
  ylab("Treg proportion (of T cells)") +
  theme_bw() +
  theme(panel.grid = element_blank())

p5b_2 <- ggscatter(
  df_afp,
  x = "Log10_AFP",
  y = "DC_C5_CXCL10",
  size = 1.6,
  add = "reg.line",
  add.params = list(color = "#0AA1FF", fill = "#a5dff9", size = 1),
  conf.int = TRUE
) +
  stat_cor(method = "spearman", label.sep = "\n") +
  xlab("Serum AFP (log10)") +
  ylab("DC_C5_CXCL10 proportion (of DCs)") +
  theme_bw() +
  theme(panel.grid = element_blank())

p5b <- p5b_1 + p5b_2 + plot_layout(ncol = 2)
ggsave(file.path(fig_dir, "Fig5b_Treg_DC_CXCL10_proportion_vs_AFP_scatter.pdf"), plot = p5b, width = 10.2, height = 4.6, useDingbats = FALSE)

################################################################################
# Fig5c: Lollipop plot of bulk correlations (25 datasets)
################################################################################

bulk_df <- readxl::read_xlsx(file.path(data_dir, "Fig5c_HCCDB_2_Treg-DC_CXCL10.xlsx"))
bulk_df$R <- as.numeric(bulk_df$R)
bulk_df$p <- as.numeric(bulk_df$p)
bulk_df$Sample <- factor(bulk_df$Sample, levels = bulk_df$Sample[order(bulk_df$R)])

p_label_x <- max(bulk_df$R, na.rm = TRUE) + 0.05

p5c <- ggplot(bulk_df, aes(y = Sample, x = R)) +
  geom_segment(aes(y = Sample, yend = Sample, x = 0, xend = R), linewidth = 1.2, color = "#D6604D") +
  geom_point(shape = 21, size = 7.5, stroke = 1, color = "#D6604D", fill = "#D6604D") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text(aes(x = p_label_x, label = paste0("p<", signif(p, 2))), color = "black", size = 4, hjust = 0) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black")
  ) +
  labs(
    title = "Correlation between Treg and DC_CXCL10 score",
    x = "Correlation coefficient (R)",
    y = ""
  ) +
  expand_limits(x = p_label_x + 0.18)

ggsave(file.path(fig_dir, "Fig5c_Bulk_Treg_DC_CXCL10_score_correlation_lollipop.pdf"), plot = p5c, device = "pdf", width = 7.5, height = 11.5, units = "in", useDingbats = FALSE)

################################################################################
# Fig5d–f: scRNA-seq validation dataset (GSE151530)
################################################################################

Hep_gse151530 <- readRDS(file.path(data_dir, "GSE151530_Hep_majorcelltype.rds"))
Hep_gse151530$orig.ident <- Hep_gse151530$S_ID

# Fig5d: UMAP of major cell types
p5d <- DimPlot(Hep_gse151530, group.by = "Type", reduction = "umap", raster = FALSE) +
  theme_dr(
    xlength = 0.22,
    ylength = 0.22,
    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  theme(panel.grid = element_blank()) +
  ggtitle("GSE151530 major cell types")
ggsave(file.path(fig_dir, "Fig5d_GSE151530_major_celltypes_umap.pdf"), plot = p5d, width = 7.2, height = 5.0, useDingbats = FALSE)

# Fig5e: Signature scores in T and myeloid cells
dc_cxcl10_genes <- c(
  "CXCL10", "CXCL9", "GBP1", "CXCL11", "GBP5", "ISG15", "TNFSF10", "WARS", "GBP4", "STAT1",
  "IL4I1", "VAMP5", "GBP2", "TNFAIP2", "LAP3", "IFITM3", "IFI6", "MX1", "MT2A", "ANKRD22",
  "C15orf48", "IDO1", "IFIT3", "CD40", "IFIT2", "TMEM176B", "RSAD2", "SERPING1", "TYMP",
  "PARP14", "TNFSF13B", "NINJ1", "LY6E", "SLAMF7", "CTSC", "ATF5", "SAMD9L", "C1QB", "IRF1",
  "PSME2", "ISG20", "HAPLN3", "PPA1", "RNF213", "G0S2", "PLEK", "TAP1", "FCGR1A", "UBE2L6",
  "TFEC"
)

treg_genes <- c(
  "TNFRSF18", "TNFRSF4", "FOXP3", "BATF", "TIGIT", "IL2RA", "CTLA4", "CARD16", "LAYN", "LAIR2",
  "RTKN2", "TBC1D4", "CTSC", "IKZF2", "STAM", "GADD45A", "CORO1B", "GLRX", "DNPH1", "BEX3",
  "ICOS", "ARID5B", "TNFRSF9", "CD27", "UGP2", "GK", "NAMPT", "PMAIP1", "DUSP4", "GBP2", "IL32",
  "MAST4", "BACH1", "CCNG2", "PKM", "BTG3", "SAT1", "RAB11FIP1", "PELI1", "PHTF2", "TNFRSF1B",
  "RAB9A", "CD4", "HTATIP2", "PHLDA1", "DUSP16", "ACP5", "SELL", "TYMP", "ZNF292"
)

Mye_gse151530 <- subset(Hep_gse151530, Type %in% c("Myeloid Cells", "Myeloid cells", "Myeloid"))
T_gse151530 <- subset(Hep_gse151530, Type %in% c("T cells", "T cell", "T"))

Mye_gse151530 <- AddModuleScore(Mye_gse151530, features = list(dc_cxcl10_genes), ctrl = 100, name = "DC_CXCL10")
Mye_gse151530$DC_CXCL10_Score <- Mye_gse151530$DC_CXCL101

T_gse151530 <- AddModuleScore(T_gse151530, features = list(treg_genes), ctrl = 100, name = "Treg")
T_gse151530$Treg_Score <- T_gse151530$Treg1

p5e_left <- FeaturePlot(
  Mye_gse151530,
  features = "DC_CXCL10_Score",
  min.cutoff = -0.5,
  cols = c("lightgrey", "orange", "orange", "red", "red"),
  raster = FALSE
) +
  ggtitle("GSE151530 - DC_CXCL10 score")

p5e_right <- FeaturePlot(
  T_gse151530,
  features = "Treg_Score",
  min.cutoff = -0.3,
  cols = c("lightgrey", "orange", "orange", "red", "red"),
  raster = FALSE
) +
  ggtitle("GSE151530 - Treg score")

p5e <- p5e_left | p5e_right
ggsave(file.path(fig_dir, "Fig5e_GSE151530_Treg_DC_CXCL10_signature_scores_umap.pdf"), plot = p5e, width = 11.6, height = 5.0, useDingbats = FALSE)

# Fig5f: Per-sample mean score correlation
df_t <- T_gse151530@meta.data %>%
  group_by(orig.ident) %>%
  summarise(Treg = mean(Treg_Score, na.rm = TRUE), .groups = "drop")

df_m <- Mye_gse151530@meta.data %>%
  group_by(orig.ident) %>%
  summarise(DC_CXCL10 = mean(DC_CXCL10_Score, na.rm = TRUE), .groups = "drop")

merged_df <- inner_join(df_t, df_m, by = "orig.ident")

p5f <- ggscatter(
  merged_df,
  x = "Treg",
  y = "DC_CXCL10",
  size = 2.2,
  add = "reg.line",
  add.params = list(color = "#0AA1FF", fill = "#a5dff9", size = 1),
  conf.int = TRUE
) +
  stat_cor(method = "spearman", label.sep = "\n") +
  xlab("Mean Treg score per sample") +
  ylab("Mean DC_CXCL10 score per sample") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(file.path(fig_dir, "Fig5f_GSE151530_Treg_vs_DC_CXCL10_mean_score_scatter.pdf"), plot = p5f, width = 4.8, height = 4.6, useDingbats = FALSE)

################################################################################
# Fig5g: Spatial co-localization (HCC4L)
################################################################################

ST_list <- readRDS(file.path(data_dir, "ST_list_CellTrek.rds"))
HCC4L <- ST_list[["HCC4L"]]

cols_region <- c(Mal = "#DC050C", nMal = "#1965B0", Bdy = "#4DAF4A")
p5g_i <- SpatialPlot(HCC4L, group.by = "region", image.alpha = 0) +
  ggtitle("HCC4L") +
  scale_fill_manual(values = cols_region)
ggsave(file.path(fig_dir, "Fig5g_i_HCC4L_region_annotation.pdf"), plot = p5g_i, width = 5.6, height = 5.2, useDingbats = FALSE)

p5g_ii <- SpatialPlot(HCC4L, features = "AFP", image.alpha = 0, stroke = 0) +
  ggtitle("HCC4L - AFP")
ggsave(file.path(fig_dir, "Fig5g_ii_HCC4L_AFP_spatial.pdf"), plot = p5g_ii, width = 5.6, height = 5.2, useDingbats = FALSE)

HCC4L_semla <- UpdateSeuratForSemla(HCC4L)
cols_score <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#D6604D", "#B2182B")

p5g_iii <- MapMultipleFeatures(
  HCC4L_semla,
  pt_size = 2,
  max_cutoff = 0.997,
  override_plot_dims = TRUE,
  features = c("Treg_AUCell", "DC_CXCL10_AUCell"),
  colors = cols_score
)
ggsave(file.path(fig_dir, "Fig5g_iii_HCC4L_Treg_DC_CXCL10_signature_scores.pdf"), plot = p5g_iii, width = 6.5, height = 10.6, useDingbats = FALSE)

celltrek_obj <- HCC4L@misc$CellTrek

plot_cells <- c("Treg", dc_subtypes)
celltrek_sel <- subset(celltrek_obj, CellTrek_Cell_type %in% plot_cells)
celltrek_sel$CellTrek_Cell_type <- factor(celltrek_sel$CellTrek_Cell_type, levels = plot_cells)

cols_celltrek <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF", "#A6CEE3")
names(cols_celltrek) <- plot_cells

p5g_iv <- SpatialPlot(celltrek_sel, group.by = "CellTrek_Cell_type") +
  scale_fill_manual(values = cols_celltrek) +
  ggtitle("HCC4L - CellTrek mapping")
ggsave(file.path(fig_dir, "Fig5g_iv_HCC4L_CellTrek_mapping_Treg_DCs.pdf"), plot = p5g_iv, width = 6.5, height = 6.5, useDingbats = FALSE)

kdist_res <- HCC4L@misc$Treg_kdist
kdist_res <- kdist_res %>% filter(cell_names %in% dc_subtypes)
kdist_res$cell_names <- factor(kdist_res$cell_names, levels = dc_subtypes)

med_order <- kdist_res %>%
  group_by(cell_names) %>%
  summarise(med = median(Treg.Others, na.rm = TRUE), .groups = "drop") %>%
  arrange(med)
kdist_res$cell_names <- factor(kdist_res$cell_names, levels = med_order$cell_names)

p5g_v <- ggboxplot(
  data = kdist_res,
  x = "cell_names",
  y = "Treg.Others",
  fill = "cell_names",
  title = "K-distance to Treg"
) +
  scale_fill_manual(values = cols_celltrek[names(cols_celltrek) %in% levels(kdist_res$cell_names)]) +
  stat_compare_means(method = "kruskal.test", label.x = 4.5, label.y = max(kdist_res$Treg.Others, na.rm = TRUE) * 1.08) +
  theme(
    plot.title = element_text(color = "black", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(y = "Distance")

ggsave(file.path(fig_dir, "Fig5g_v_HCC4L_Treg_kdist_DC_subsets_boxplot.pdf"), plot = p5g_v, width = 6.5, height = 5.0, useDingbats = FALSE)
