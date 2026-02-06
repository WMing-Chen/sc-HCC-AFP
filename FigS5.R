############## FigS5_Treg_immunotherapy_and_spatial #############################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(patchwork)
  library(Seurat)
})

fig_dir <- "./FigS5"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# FigS5a: Treg scores in CD4+ T cells (GSE235863; pre vs post, NR vs R)
################################################################################

gse235863_dir <- file.path(data_dir, "GSE235863")
sce.all <- readRDS(file.path(gse235863_dir, "sce.all.rds"))

treg_genes <- c(
  "TNFRSF18", "TNFRSF4", "FOXP3", "BATF", "TIGIT", "IL2RA", "CTLA4", "CARD16", "LAYN", "LAIR2", "RTKN2", "TBC1D4",
  "CTSC", "IKZF2", "STAM", "GADD45A", "CORO1B", "GLRX", "DNPH1", "BEX3", "ICOS", "ARID5B", "TNFRSF9", "CD27",
  "UGP2", "GK", "NAMPT", "PMAIP1", "DUSP4", "GBP2", "IL32", "MAST4", "BACH1", "CCNG2", "PKM", "BTG3", "SAT1",
  "RAB11FIP1", "PELI1", "PHTF2", "TNFRSF1B", "RAB9A", "CD4", "HTATIP2", "PHLDA1", "DUSP16", "ACP5", "SELL", "TYMP",
  "ZNF292"
)

cd4 <- subset(sce.all, subset = major_cluster %in% c("CD4T"))
cd4 <- AddModuleScore(cd4, features = list(treg_genes), ctrl = 100, name = "Treg")
cd4$Treg_Score <- cd4$Treg1

cd4$Pre_Pos <- ifelse(grepl("pre", cd4$sample, ignore.case = TRUE), "Pre",
  ifelse(grepl("post|pos", cd4$sample, ignore.case = TRUE), "Pos", NA)
)
cd4$respond_Pre_Pos <- paste0(cd4$Pre_Pos, "_", cd4$respond)

p_s5a <- VlnPlot(
  cd4,
  group.by = "respond_Pre_Pos",
  features = "Treg_Score",
  stack = FALSE,
  flip = TRUE,
  fill.by = "ident",
  cols = c("Pre_NR" = "#5999cc", "Pos_NR" = "#fb7e39", "Pre_R" = "#5999cc", "Pos_R" = "#fb7e39"),
  ncol = 1,
  pt.size = 0
) +
  guides(fill = "none") +
  labs(title = "Treg Scores in CD4+ T cells", x = "", y = "Treg Scores") +
  theme(
    axis.text.x = element_text(angle = 45, size = 16, vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13),
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

ggsave(file.path(fig_dir, "FigS5a_GSE235863_CD4T_Treg_scores_violin.pdf"), plot = p_s5a, width = 4.8, height = 5.4, useDingbats = FALSE)

################################################################################
# FigS5b–c: CellTrek spatial mapping + DC-to-Treg K-distance
################################################################################

dc_subtypes <- c("DC_C1_CD1C", "DC_C2_STMN1", "DC_C3_CLEC9A", "DC_C4_LAMP3", "DC_C5_CXCL10", "DC_C6_APOC3")
plot_cells <- c("Treg", dc_subtypes)
cols_celltrek <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF", "#A6CEE3")
names(cols_celltrek) <- plot_cells

ST_list <- readRDS(file.path(data_dir, "ST_list_CellTrek.rds"))

rep_samples <- list(
  NR = c("P11T", "P3T", "P5T", "P8T"),
  R = c("P9T", "P10T")
)

for (grp in names(rep_samples)) {
  for (sample_sel in rep_samples[[grp]]) {
    st_obj <- ST_list[[sample_sel]]
    celltrek_obj <- st_obj@misc$CellTrek

    celltrek_sel <- subset(celltrek_obj, CellTrek_Cell_type %in% plot_cells)
    celltrek_sel$CellTrek_Cell_type <- factor(celltrek_sel$CellTrek_Cell_type, levels = plot_cells)

    p_left <- SpatialPlot(celltrek_sel, group.by = "CellTrek_Cell_type") +
      scale_fill_manual(values = cols_celltrek) +
      theme(legend.key = element_blank(), legend.position = "top") +
      guides(fill = guide_legend(override.aes = list(size = 3))) +
      ggtitle(paste0(sample_sel, " (", grp, ")"))

    kdist_res <- st_obj@misc$Treg_kdist
    kdist_res <- kdist_res %>% filter(cell_names %in% dc_subtypes)
    med_order <- kdist_res %>%
      group_by(cell_names) %>%
      summarise(med = median(Treg.Others, na.rm = TRUE), .groups = "drop") %>%
      arrange(med)
    kdist_res$cell_names <- factor(kdist_res$cell_names, levels = med_order$cell_names)

    x_pos <- length(levels(kdist_res$cell_names)) / 2
    p_right <- ggboxplot(
      data = kdist_res,
      x = "cell_names",
      y = "Treg.Others",
      fill = "cell_names",
      title = "K-distance to Treg"
    ) +
      scale_fill_manual(values = cols_celltrek[names(cols_celltrek) %in% levels(kdist_res$cell_names)]) +
      stat_compare_means(
        method = "kruskal.test",
        label.x = x_pos,
        label.y = max(kdist_res$Treg.Others, na.rm = TRUE) * 1.1
      ) +
      theme(
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = "none"
      ) +
      labs(y = "Distance")

    p_s5bc <- p_left | p_right
    out_name <- ifelse(grp == "NR", "FigS5b", "FigS5c")
    ggsave(
      file.path(fig_dir, paste0(out_name, "_", sample_sel, "_CellTrek_mapping_and_kdist.pdf")),
      plot = p_s5bc,
      width = 13.0,
      height = 6.0,
      useDingbats = FALSE
    )
  }
}
