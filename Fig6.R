############## Fig6_Treg_DC_CXCL10_and_immunotherapy ############################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(GSVA)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(patchwork)
  library(scRNAtoolVis)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(semla)
  library(survival)
  library(survminer)
  library(tidyr)
})

setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")

fig_dir <- "./00Figures/Fig6"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Input (not tracked; place under ./data/)
# - GSE202069 bulk: `gse202069.RData` with objects `gse202069` (expr: gene x sample) and `gse202069_cli` (clinical)
# - GSE235863 scRNA: `GSE235863/sce.all.rds` (processed Seurat object with `major_cluster`, `sample`)
# - ST + CellTrek: `ST_list_CellTrek.rds` (Seurat object list; must contain `region`, `AFP`, `Treg_AUCell`, `DC_CXCL10_AUCell`,
#   and `@misc$CellTrek`, `@misc$Treg_kdist` produced by `sc-HCC-AFP/CellTrek.R`)

################################################################################
# Fig6a: Kaplan–Meier survival (GSE202069 anti-PD-1)
################################################################################

load(file.path(data_dir, "gse202069.RData"))

Tumor_PD1_clin <- gse202069_cli %>% filter(group == "Tumor", PD1.therapy == "Yes")
Tumor_PD1_clin <- Tumor_PD1_clin[!is.na(Tumor_PD1_clin$sample), ]

Tumor_PD1 <- gse202069[, Tumor_PD1_clin$sample, drop = FALSE]

sig_genesets <- list(
  Tumor_C1_CYP3A4 = c(
    "CYP3A4", "SERPINC1", "CYP2A6", "A1BG", "ANG", "APOA5", "AZGP1", "HP", "CES2", "TTR", "HPR", "PLG",
    "TAT", "HSD11B1", "CD14", "HRG", "CPS1", "KNG1", "HPX", "SLC22A1", "SLC27A5", "APOC3", "HSD17B6",
    "TF", "CFHR3", "PCK1", "CYP8B1", "C3", "EPHX1", "AMBP", "CYP2E1", "GLUL", "CP", "SERPING1", "ITIH1",
    "CYP2D6", "ADH1B", "CFHR1", "CYP2C8", "APOH", "APOB", "G6PC", "SPP2", "UGT2B7", "RBP4", "PON3",
    "APCS", "SLC2A2", "ACSM2B", "CPB2"
  ),
  Tumor_C3_STMN1 = c(
    "STMN1", "UBE2C", "TUBA1B", "PTTG1", "H2AFZ", "TOP2A", "HMGB2", "NUSAP1", "HIST1H4C", "CDKN3", "CCNB1",
    "UBE2S", "CENPF", "BIRC5", "CKS2", "PCLAF", "HMGN2", "KPNA2", "TPX2", "CDK1", "TYMS", "ASPM", "TUBB",
    "CCNB2", "JPT1", "CDC20", "HMMR", "RRM2", "TUBB4B", "CENPW", "UBE2T", "CKS1B", "AURKA", "SMC4", "ZWINT",
    "MKI67", "TK1", "H2AFV", "TUBA1C", "HMGB1", "DEK", "HMGB3", "RAN", "PCNA", "RANBP1", "DTYMK", "DUT",
    "MAD2L1", "NQO1", "TMEM106C"
  ),
  Tumor_C6_AFP = c(
    "AFP", "GPC3", "PEG10", "A2M", "SCD", "SNORC", "ACSL4", "APOB", "ITIH2", "NPW", "EPCAM", "SPINT2", "MTRNR2L12",
    "C3", "MBNL3", "RRBP1", "LPGAT1", "C5", "CPD", "SERPINA5", "DLK1", "FDPS", "SERPINA1", "CYP3A7", "ANGPTL8",
    "GC", "AHNAK", "COL27A1", "ZKSCAN1", "TGOLN2", "FN1", "GLUD1", "SQLE", "FARP1", "DSP", "RELN", "S100A14",
    "FABP1", "FGFR4", "FGFR3", "APOA1", "FGFR2", "F10", "ATP2B1", "SPTBN1", "VTN", "HSD17B2", "DDX17", "ACSL3",
    "APLP2"
  ),
  Treg = c(
    "TNFRSF18", "TNFRSF4", "FOXP3", "BATF", "TIGIT", "IL2RA", "CTLA4", "CARD16", "LAYN", "LAIR2", "RTKN2", "TBC1D4",
    "CTSC", "IKZF2", "STAM", "GADD45A", "CORO1B", "GLRX", "DNPH1", "BEX3", "ICOS", "ARID5B", "TNFRSF9", "CD27",
    "UGP2", "GK", "NAMPT", "PMAIP1", "DUSP4", "GBP2", "IL32", "MAST4", "BACH1", "CCNG2", "PKM", "BTG3", "SAT1",
    "RAB11FIP1", "PELI1", "PHTF2", "TNFRSF1B", "RAB9A", "CD4", "HTATIP2", "PHLDA1", "DUSP16", "ACP5", "SELL", "TYMP",
    "ZNF292"
  )
)

ES <- gsva(as.matrix(Tumor_PD1), sig_genesets)
ES <- as.data.frame(t(ES))

surv_df <- Tumor_PD1_clin
if ("OS" %in% colnames(surv_df) && !("OS.state" %in% colnames(surv_df))) {
  colnames(surv_df)[colnames(surv_df) == "OS"] <- "OS.state"
}
if (!("OS.time_months" %in% colnames(surv_df)) && "OS.Time" %in% colnames(surv_df)) {
  surv_df$OS.time_months <- round(as.numeric(surv_df$OS.Time) / 30, 2)
}
surv_df$OS.state <- as.numeric(surv_df$OS.state)
surv_df$OS.time_months <- as.numeric(surv_df$OS.time_months)

plot_list <- list()
for (sig in names(sig_genesets)) {
  ES$Group <- ifelse(ES[[sig]] > median(ES[[sig]], na.rm = TRUE), "High", "Low")
  surv_df$Group <- ES[surv_df$sample, "Group"]

  fit <- survfit(Surv(as.numeric(OS.time_months), as.numeric(OS.state)) ~ Group, data = surv_df)
  p <- ggsurvplot(
    fit,
    data = surv_df,
    title = sig,
    pval = TRUE,
    pval.size = 4.0,
    size = 1.1,
    linetype = "solid",
    palette = c("#DC0000FF", "#1f78b4"),
    legend = c(0.25, 0.5),
    legend.title = "",
    xlab = "Time (months)"
  )$plot +
    theme_bw() +
    theme(panel.grid = element_blank())

  plot_list[[sig]] <- p
}

p6a <- (plot_list$Tumor_C1_CYP3A4 | plot_list$Tumor_C3_STMN1) / (plot_list$Tumor_C6_AFP | plot_list$Treg)
ggsave(file.path(fig_dir, "Fig6a_GSE202069_KM_TumorC1_C3_C6_Treg.pdf"), plot = p6a, width = 11.5, height = 8.2, useDingbats = FALSE)

################################################################################
# Fig6b–e: scRNA validation dataset (GSE235863; anti-PD-1 + Lenvatinib)
################################################################################

gse235863_dir <- file.path(data_dir, "GSE235863")
sce.all <- readRDS(file.path(gse235863_dir, "sce.all.rds"))

cols_major <- c(
  "#36489E", "#C6EE8F", "#28B461", "#B2DCEE", "#A680B9", "#DA2917", "#FA7E4D", "#F0C674",
  "#E7298A", "#E78AC3", "#aa8282", "#d4b7b7", "#999999", "#666666"
)

p6b <- DimPlot(sce.all, group.by = "major_cluster", reduction = "umap", label = TRUE, cols = cols_major, raster = FALSE) +
  theme_dr(
    xlength = 0.22,
    ylength = 0.22,
    arrow = grid::arrow(length = grid::unit(0.15, "inches"), type = "closed")
  ) +
  theme(panel.grid = element_blank()) +
  ggtitle("GSE235863 major cell types")
ggsave(file.path(fig_dir, "Fig6b_GSE235863_major_celltypes_umap.pdf"), plot = p6b, width = 7.2, height = 5.2, useDingbats = FALSE)

dc_cxcl10_genes <- c(
  "CXCL10", "CXCL9", "GBP1", "CXCL11", "GBP5", "ISG15", "TNFSF10", "WARS", "GBP4", "STAT1", "IL4I1", "VAMP5",
  "GBP2", "TNFAIP2", "LAP3", "IFITM3", "IFI6", "MX1", "MT2A", "ANKRD22", "C15orf48", "IDO1", "IFIT3", "CD40",
  "IFIT2", "TMEM176B", "RSAD2", "SERPING1", "TYMP", "PARP14", "TNFSF13B", "NINJ1", "LY6E", "SLAMF7", "CTSC",
  "ATF5", "SAMD9L", "C1QB", "IRF1", "PSME2", "ISG20", "HAPLN3", "PPA1", "RNF213", "G0S2", "PLEK", "TAP1",
  "FCGR1A", "UBE2L6", "TFEC"
)

treg_genes <- sig_genesets$Treg

Mye <- subset(sce.all, subset = major_cluster %in% c("Myeloid", "Myeloid Cells", "Myeloid cells"))
T_cells <- subset(sce.all, subset = major_cluster %in% c("CD8T", "CD4T", "gdT", "T cells", "T cell", "T"))

Mye <- AddModuleScore(Mye, features = list(dc_cxcl10_genes), ctrl = 100, name = "DC_CXCL10")
Mye$DC_CXCL10_Score <- Mye$DC_CXCL101

T_cells <- AddModuleScore(T_cells, features = list(treg_genes), ctrl = 100, name = "Treg")
T_cells$Treg_Score <- T_cells$Treg1

p6c_left <- FeaturePlot(
  Mye,
  features = "DC_CXCL10_Score",
  cols = c("lightgrey", "orange", "orange", "red", "red"),
  min.cutoff = -0.2,
  max.cutoff = 1.4,
  raster = FALSE
) +
  ggtitle("DC_C5_CXCL10 score (Myeloid)")

p6c_right <- FeaturePlot(
  T_cells,
  features = "Treg_Score",
  cols = c("lightgrey", "orange", "orange", "red", "red"),
  min.cutoff = -0.1,
  max.cutoff = 1.4,
  raster = FALSE
) +
  ggtitle("Treg score (T cells)")

p6c <- p6c_left | p6c_right
ggsave(file.path(fig_dir, "Fig6c_GSE235863_Treg_DC_CXCL10_signature_scores_umap.pdf"), plot = p6c, width = 11.6, height = 5.0, useDingbats = FALSE)

df_t <- T_cells@meta.data %>%
  group_by(sample) %>%
  summarise(Treg = mean(Treg_Score, na.rm = TRUE), .groups = "drop")

df_m <- Mye@meta.data %>%
  group_by(sample) %>%
  summarise(DC_CXCL10 = mean(DC_CXCL10_Score, na.rm = TRUE), .groups = "drop")

merged_df <- inner_join(df_t, df_m, by = "sample")
merged_df$treatment <- ifelse(grepl("pre", merged_df$sample, ignore.case = TRUE), "pre",
  ifelse(grepl("post", merged_df$sample, ignore.case = TRUE), "post", NA)
)
pre_df <- merged_df %>% filter(treatment == "pre")

p6d <- ggscatter(
  pre_df,
  x = "Treg",
  y = "DC_CXCL10",
  size = 2.0,
  add = "reg.line",
  add.params = list(color = "#0AA1FF", fill = "#a5dff9", size = 1),
  conf.int = TRUE
) +
  stat_cor(method = "spearman", label.sep = "\n") +
  ggtitle("Correlation of  Mean Scores in Pre-Samples") +
  xlab("Mean Treg score per sample") +
  ylab("Mean DC_C5_CXCL10 score per sample") +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(file.path(fig_dir, "Fig6d_GSE235863_pre_Treg_vs_DC_CXCL10_mean_score_scatter.pdf"), plot = p6d, width = 4.9, height = 4.6, useDingbats = FALSE)

Mye$Pre_Pos <- ifelse(grepl("pre", Mye$sample, ignore.case = TRUE), "Pre",
  ifelse(grepl("post|pos", Mye$sample, ignore.case = TRUE), "Pos", NA)
)
Mye$respond_Pre_Pos <- paste0(Mye$Pre_Pos, "_", Mye$respond)

p6e <- VlnPlot(
  Mye,
  group.by = "respond_Pre_Pos",
  features = "DC_CXCL10_Score",
  stack = FALSE,
  flip = TRUE,
  fill.by = "ident",
  cols = c("Pre_NR" = "#5999cc", "Pos_NR" = "#fb7e39", "Pre_R" = "#5999cc", "Pos_R" = "#fb7e39"),
  ncol = 1,
  pt.size = 0
) +
  guides(fill = "none") +
  labs(title = "DC_C5_CXCL10 Scores in Myeloids", x = "", y = "DC_C5_CXCL10 Scores") +
  theme(
    axis.text.x = element_text(angle = 45, size = 14, vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  ) +
  scale_x_discrete(limits = c("Pre_NR", "Pos_NR", "Pre_R", "Pos_R")) +
  scale_y_continuous(limits = c(-0.5, 2.15)) +
  geom_boxplot(width = 0.18, col = "black", fill = "white", outlier.shape = NA) +
  geom_signif(
    comparisons = list(c("Pre_NR", "Pos_NR"), c("Pre_R", "Pos_R"), c("Pos_NR", "Pos_R")),
    textsize = 4.2,
    y_position = c(1.5, 1.7, 1.93),
    map_signif_level = TRUE
  )
ggsave(file.path(fig_dir, "Fig6e_GSE235863_DC_CXCL10_scores_PrePost_NR_R_violin.pdf"), plot = p6e, width = 4.8, height = 5.4, useDingbats = FALSE)

################################################################################
# Fig6f–h: Spatial transcriptomics validation (Xun et al.; 8 samples)
################################################################################

ST_list <- readRDS(file.path(data_dir, "ST_list_CellTrek.rds"))
dist_summary_scaled <- attr(ST_list, "kdist_summary_scaled")
dist_cells_scaled <- attr(ST_list, "kdist_dc_cells_scaled")

dc_subtypes <- c("DC_C1_CD1C", "DC_C2_STMN1", "DC_C3_CLEC9A", "DC_C4_LAMP3", "DC_C5_CXCL10", "DC_C6_APOC3")
sample_order <- c("P1T", "P3T", "P5T", "P8T", "P11T", "P7T", "P9T", "P10T")
nr_samples <- c("P1T", "P3T", "P5T", "P8T", "P11T")
r_samples <- c("P7T", "P9T", "P10T")

mat_f <- dist_summary_scaled %>%
  filter(cell_names %in% dc_subtypes, sample %in% sample_order) %>%
  mutate(
    cell_names = factor(cell_names, levels = dc_subtypes),
    sample = factor(sample, levels = sample_order)
  ) %>%
  select(cell_names, sample, mean_dist_zscore) %>%
  tidyr::pivot_wider(names_from = sample, values_from = mean_dist_zscore) %>%
  tibble::column_to_rownames("cell_names") %>%
  as.matrix()
mat_f <- mat_f[, sample_order, drop = FALSE]

column_ha <- HeatmapAnnotation(
  Group = ifelse(sample_order %in% r_samples, "R", "NR"),
  col = list(Group = c(NR = "#DA2917", R = "#28B461")),
  simple_anno_size = grid::unit(3, "mm")
)

col_fun <- colorRamp2(c(-1, 0, 1), c("#1A5592", "white", "#B83D3D"))
ht_f <- Heatmap(
  mat_f,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = col_fun,
  top_annotation = column_ha,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 12),
  row_names_gp = grid::gpar(fontsize = 12),
  column_title = "Average K-distance to Treg",
  heatmap_legend_param = list(title = "Distance (scaled)")
)

pdf(file.path(fig_dir, "Fig6f_CellTrek_DCsub_to_Treg_distance_scaled_heatmap.pdf"), width = 8.5, height = 5.2, useDingbats = FALSE)
draw(ht_f, padding = grid::unit(c(10, 20, 10, 10), "mm"), heatmap_legend_side = "right")
dev.off()

df_g <- dist_cells_scaled %>%
  filter(sample %in% sample_order, cell_names == "DC_C5_CXCL10") %>%
  mutate(Group = ifelse(sample %in% r_samples, "R", "NR"))

my_comparisons <- list(c("R", "NR"))
p6g <- ggplot(df_g, aes(x = factor(Group, levels = c("NR", "R")), y = Treg.Others_zscore, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.2, color = "black", alpha = 0.6) +
  labs(title = "DC_C5_CXCL10 to Treg", x = "", y = "Distance (scaled)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = c(NR = "#DC050C", R = "#386cb0")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 5) +
  guides(fill = "none")
ggsave(file.path(fig_dir, "Fig6g_CellTrek_DC_C5_CXCL10_to_Treg_distance_scaled_boxplot.pdf"), plot = p6g, width = 4.6, height = 4.6, useDingbats = FALSE)

rep_samples <- list(NR = "P1T", R = "P7T")
cols_region <- c("#7fc97f", "#beaed4", "#fdc086", "#386cb0", "#f0027f", "#a34e3b", "#666666", "#1b9e77", "#d95f02", "#6c67a5", "#d01b2a", "#43acde", "#efbd25")
cols_score <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D", "#D6604D", "#B2182B")

for (grp in names(rep_samples)) {
  sample_sel <- rep_samples[[grp]]
  st_obj <- ST_list[[sample_sel]]

  p_i <- SpatialPlot(st_obj, group.by = "region", image.alpha = 0) +
    ggtitle(sample_sel)
  ggsave(file.path(fig_dir, paste0("Fig6h_", grp, "_i_", sample_sel, "_region.pdf")), plot = p_i, width = 5.6, height = 5.2, useDingbats = FALSE)

  p_ii <- SpatialPlot(st_obj, features = "AFP", image.alpha = 0, stroke = 0) +
    ggtitle(paste0(sample_sel, " - AFP"))
  ggsave(file.path(fig_dir, paste0("Fig6h_", grp, "_ii_", sample_sel, "_AFP_spatial.pdf")), plot = p_ii, width = 5.6, height = 5.2, useDingbats = FALSE)

  st_semla <- UpdateSeuratForSemla(st_obj)
  p_iii <- MapMultipleFeatures(
    st_semla,
    pt_size = 2,
    max_cutoff = 0.997,
    override_plot_dims = TRUE,
    features = c("DC_CXCL10_AUCell", "Treg_AUCell"),
    colors = cols_score
  ) +
    ggtitle(sample_sel)
  ggsave(file.path(fig_dir, paste0("Fig6h_", grp, "_iii_", sample_sel, "_Treg_DC_CXCL10_scores.pdf")), plot = p_iii, width = 6.5, height = 10.6, useDingbats = FALSE)

  celltrek_obj <- st_obj@misc$CellTrek
  plot_cells <- c("Treg", dc_subtypes)
  celltrek_sel <- subset(celltrek_obj, CellTrek_Cell_type %in% plot_cells)
  celltrek_sel$CellTrek_Cell_type <- factor(celltrek_sel$CellTrek_Cell_type, levels = plot_cells)

  cols_celltrek <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF", "#A6CEE3")
  names(cols_celltrek) <- plot_cells
  p_iv <- SpatialPlot(celltrek_sel, group.by = "CellTrek_Cell_type") +
    scale_fill_manual(values = cols_celltrek) +
    theme(legend.key = element_blank(), legend.position = "top") +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    ggtitle(sample_sel)
  ggsave(file.path(fig_dir, paste0("Fig6h_", grp, "_iv_", sample_sel, "_CellTrek_mapping.pdf")), plot = p_iv, width = 7, height = 6.5, useDingbats = FALSE)

  kdist_res <- st_obj@misc$Treg_kdist
  kdist_res <- kdist_res %>% filter(cell_names %in% dc_subtypes)
  med_order <- kdist_res %>%
    group_by(cell_names) %>%
    summarise(med = median(Treg.Others, na.rm = TRUE), .groups = "drop") %>%
    arrange(med)
  kdist_res$cell_names <- factor(kdist_res$cell_names, levels = med_order$cell_names)

  x_pos <- length(levels(kdist_res$cell_names)) / 2
  p_v <- ggboxplot(
    data = kdist_res,
    x = "cell_names",
    y = "Treg.Others",
    fill = "cell_names",
    title = "K-distance to Treg"
  ) +
    scale_fill_manual(values = cols_celltrek[names(cols_celltrek) %in% levels(kdist_res$cell_names)]) +
    stat_compare_means(method = "kruskal.test", label.x = x_pos, label.y = max(kdist_res$Treg.Others, na.rm = TRUE) * 1.1) +
    theme(
      plot.title = element_text(color = "black", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      legend.position = "none"
    ) +
    labs(y = "Distance")
  ggsave(file.path(fig_dir, paste0("Fig6h_", grp, "_v_", sample_sel, "_Treg_kdist_DC_subsets_boxplot.pdf")), plot = p_v, width = 6.5, height = 5.0, useDingbats = FALSE)
}
