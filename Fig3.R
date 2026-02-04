############## Fig3_T_cells #####################################################

rm(list = ls())
gc()

library(dplyr)
library(Matrix)
library(Seurat)
library(tidydr)
library(tidyverse)
library(harmony)
library(scRNAtoolVis)
library(cowplot)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(ggbeeswarm)
library(GSVA)
library(survival)
library(survminer)
library(org.Hs.eg.db)

setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")

fig_dir <- "./00Figures/Fig3"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Input (not tracked in repo; place under ./data/)
# - scRNA: `HCC_T_cells.rds` (preferred) or `HCC_scRNA.rds` (will subset CD4/CD8 T cells)
# - Gene programs: `T_cell_19_gene_programs.csv` (19 columns; each column is a program, values are gene symbols)
# - TCGA: `TCGA_matrix.txt` (gene x sample, gene symbols as rownames)
# - TCGA phenotype: `TCGA_phenotype_AFP_value.txt` (must include sample_id, OS.time_months, OS.state)

################################################################################
# Load T cells
################################################################################
Hep <- readRDS(file.path(data_dir, "HCC_T_cells.rds"))
table(Hep@meta.data$subcelltype)

################################################################################
# Fig3a–b: UMAP projections of T cell clusters and AFP group distribution
################################################################################
Hep@meta.data$Seu_Clusters <- Idents(Hep)
Hep@meta.data$Seu_Clusters <- as.character(Hep@meta.data$Seu_Clusters)
Hep@meta.data$Seu_Clusters <- factor(Hep@meta.data$Seu_Clusters, levels = sort(unique(as.integer(Hep@meta.data$Seu_Clusters))))

# Manual subcluster annotation (14 clusters; update mapping if cluster IDs differ)
Hep@meta.data$subcelltype <- as.character(Hep@meta.data$Seu_Clusters)
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(0)] <- "CD4T_C1_GPR183"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(6)] <- "CD4T_C2_PLCG2"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(10)] <- "CD4T_C3_CXCL13"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(2)] <- "Treg"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(8, 11, 13)] <- "Cycling_T_cells"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(5)] <- "MAIT"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(14)] <- "NKT"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(4)] <- "T_stress"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(1)] <- "CD8T_C1_GZMK"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(3)] <- "CD8T_C2_PDCD1"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(7)] <- "CD8T_C3_APOA2"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(9)] <- "CD8T_C4_CST3"
Hep@meta.data$subcelltype[Hep@meta.data$subcelltype %in% c(12)] <- "CD8T_C5_IFIT3"
table(Hep@meta.data$subcelltype)


plot_order <- c(
  "CD8T_C1_GZMK", "CD8T_C2_PDCD1", "CD8T_C3_APOA2", "CD8T_C4_CST3", "CD8T_C5_IFIT3",
  "MAIT", "T_stress",
  "CD4T_C1_GPR183", "CD4T_C2_PLCG2", "CD4T_C3_CXCL13",
  "Treg", "NKT", "Cycling_T_cells"
)
Hep@meta.data$subcelltype <- factor(Hep@meta.data$subcelltype, levels = plot_order)

cols <- c(
  "#28B461", "#8DD3C7", "#F0C674", "#FA7E4D", "#E7298A",
  "#A680B9", "#E78AC3",
  "#36489E", "#C6EE8F", "#d4b7b7",
  "#DA2917", "#B2DCEE", "#aa8282"
)

p3a <- DimPlot(Hep, group.by = "subcelltype", reduction = "umap", label = FALSE, cols = cols, raster = FALSE) +
  ggtitle("HCC T cells") +
  theme_dr(
    xlength = 0.22,
    ylength = 0.22,
    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 10)  
        )
p3a
ggsave(file.path(fig_dir, "Fig3a_T_cell_subsets_umap.pdf"), plot = p3a, width = 7, height = 4.5)

p3b <- DimPlot(
  Hep,
  reduction = "umap",
  group.by = "AFP_status_group",
  raster = FALSE,
  cols = c("#25acda", "#db5024" ),
  label = F)+
  theme_dr(
    xlength = 0.22,
    ylength = 0.22,
    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  theme(panel.grid = element_blank())
p3b
ggsave(file.path(fig_dir, "Fig3b_T_cell_subsets_AFP_status.pdf"), plot = p3b, width = 10, height = 4.5)

################################################################################
# Fig3c: Ro/e heatmap of T cell subsets (AFP_Pos vs AFP_Neg)
################################################################################

source("./utils/ROE_Heatmap_Function.R")
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:7])

plot.data <- ROIE(table(Hep@meta.data[, c("subcelltype", "AFP_status_group")])) %>%
  reshape2::melt() %>%
  mutate(value = pmin(value, 3))

plot.data <- plot.data %>%
  mutate(
    Var1 = factor(
      Var1,
      levels = as.character(unique(Var1[Var2 == "AFP_Pos"][order(-value[Var2 == "AFP_Pos"])]))
    )
  )

p3c <- ggplot(plot.data, aes(Var2, forcats::fct_rev(Var1), fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_gradientn(colours = myPalette(100)) +
  scale_x_discrete(limits = c("AFP_Pos", "AFP_Neg")) +
  labs(x = "", y = "", fill = "Ro/e") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12))
p3c
ggsave(file.path(fig_dir, "Fig3c_T_cell_ROE_heatmap.pdf"), plot = p3c, device = "pdf", width = 4, height = 6, units = "in")

################################################################################
# Fig3d: Canonical marker genes dot plot
################################################################################

Idents(Hep) <- "subcelltype"
markergenes <- c(
  "CD4", "CD8A", "GZMK", "GZMA", "PDCD1", "TIGIT",
  "APOA2", "APOC3", "C1QB", "C1QA", "IFIT3", "IFI6",
  "SLC4A10", "KLRB1", "FOS", "HSPA1B",
  "GPR183", "CCR6", "PLCG2", "GIMAP7", "CXCL13", "NMB",
  "FOXP3", "TNFRSF4", "GNLY", "FCGR3A", "STMN1", "MKI67"
)

p3d <- jjDotPlot(
  object = Hep,
  id = "subcelltype",
  gene = markergenes,
  cluster.order = rev(plot_order),
  ytree = FALSE,
  legend.position = "top",
  base_size = 20,
  x.text.angle = 45,
  x.text.vjust = 1
) +
  theme(panel.grid = element_blank())
p3d
ggsave(file.path(fig_dir, "Fig3d_T_cell_markers_dotplot.pdf"), plot = p3d, width = 11, height = 13)

################################################################################
# Fig3e: Heatmap of scaled signature scores (19 curated gene programs)
################################################################################

gene_program_file <- file.path(data_dir, "T_cell_19_gene_programs.csv")
program_table <- read.csv(gene_program_file, check.names = FALSE)
program_list <- lapply(seq_len(ncol(program_table)), function(i) na.omit(program_table[[i]]))
names(program_list) <- colnames(program_table)

names(program_list)[names(program_list) == "Naïve"] = "Naive"
names(program_list) <- gsub("[^A-Za-z0-9 ]", " ", names(program_list))
names(program_list)[names(program_list) == "Activation Effector function"] = "Activation/Effector function"

for (nm in names(program_list)) {
  Hep <- AddModuleScore(Hep, features = list(program_list[[nm]]), name = nm)
  colnames(Hep@meta.data)[colnames(Hep@meta.data) == paste0(nm, "1")] <- nm
}

score_mat <- as.matrix(Hep@meta.data[, names(program_list), drop = FALSE])
avg_scores <- sapply(levels(Hep@meta.data$subcelltype), function(ct) {
  colMeans(score_mat[Hep@meta.data$subcelltype == ct, , drop = FALSE], na.rm = TRUE)
})
avg_scores <- t(scale(t(avg_scores)))

column_ha <- HeatmapAnnotation(
  Clusters = levels(Hep@meta.data$subcelltype),
  col = list(Clusters = setNames(cols, levels(Hep@meta.data$subcelltype))),
  simple_anno_size = unit(3, "mm")
)
col_fun <- colorRamp2(c(-2, 0, 2), c("#1A5592", "white", "#B83D3D"))

ht <- ComplexHeatmap::Heatmap(
  as.matrix(avg_scores),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  top_annotation = column_ha,
  col = col_fun,
  row_title = "",
  column_title = "",
  row_names_side = "left"
)
pdf(file.path(fig_dir, "Fig3e_T_cell_gene_program_scores_heatmap.pdf"), width = 10, height = 6, useDingbats = FALSE)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

################################################################################
# Fig3f: Box plots of CD8T_C1 / CD8T_C5 / Treg abundance by AFP group
################################################################################

meta <- Hep@meta.data %>%
  dplyr::select(orig.ident, AFP_status_group, subcelltype)
colnames(meta) <- c("orig.ident", "Group", "celltype")

meta_percent <- meta %>%
  group_by(orig.ident, Group, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(values = count / sum(count)) %>%
  dplyr::select(orig.ident, values, Group, celltype)
meta_percent <- unique(meta_percent)

sel_celltypes <- c("CD8T_C1_GZMK", "CD8T_C5_IFIT3", "Treg")
meta_percent <- meta_percent %>% filter(celltype %in% sel_celltypes)

source("./utils/CellType_Group_BoxPlot_Function.R")
cols_group <- rev(c("#327db7", "#f18c8d"))

p3f <- Groups_box_plot(
  data = meta_percent,
  Group_box_colors = cols_group,
  Group_top = "celltype",
  Group_top_order = sel_celltypes,
  Group_box = "Group",
  Group_box_order = c("AFP_Pos", "AFP_Neg"),
  title_name = "Cell percentages across samples",
  sig_lable = TRUE
)
ggsave(file.path(fig_dir, "Fig3f_CD8T_C1_CD8T_C5_Treg_abundance_boxplot.pdf"), plot = p3f, device = "pdf", width = 10, height = 4.5)

################################################################################
# Fig3g: TCGA-LIHC survival based on scRNA marker signature scores (GSVA)
################################################################################

Idents(Hep) <- "subcelltype"
mk <- FindAllMarkers(Hep, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
mk$cluster <- as.character(mk$cluster)
mk_sig <- mk %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE)

geneset_list <- list(
  CD8T_C1_GZMK = unique(mk_sig$gene[mk_sig$cluster == "CD8T_C1_GZMK"]),
  Treg = unique(mk_sig$gene[mk_sig$cluster == "Treg"])
)

tcga_expr <- read.table(file.path(data_dir, "TCGA_matrix.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
tcga_expr <- as.matrix(tcga_expr)

tcga_pheno <- read.table(file.path(data_dir, "TCGA_phenotype_AFP_value.txt"), header = TRUE, sep = "\t", check.names = FALSE)
rownames(tcga_pheno) <- tcga_pheno$sample_id

ES <- gsva(tcga_expr, geneset_list)
ES <- as.data.frame(t(ES))

common_samples <- intersect(rownames(tcga_pheno), rownames(ES))
surival_data <- tcga_pheno[common_samples, , drop = FALSE]
for (sig in names(geneset_list)) {
  surival_data$Group <- ifelse(ES[rownames(surival_data), sig] > median(ES[[sig]], na.rm = TRUE), "High", "Low")
  fit <- survfit(Surv(as.numeric(OS.time_months), as.numeric(OS.state)) ~ Group, data = surival_data)
  p <- ggsurvplot(
    fit,
    data = surival_data,
    title = sig,
    pval = TRUE,
    size = 1.2,
    linetype = "solid",
    palette = c("#DC0000FF", "#1f78b4"),
    legend.title = " ",
    xlab = "Time(Months)"
  )
  ggsave(filename = file.path(fig_dir, paste0("Fig3g_TCGA_survival_", sig, ".pdf")), plot = p$plot, width = 4.5, height = 4.3)
}

################################################################################
# Fig3h: Violin plots of inhibitory immune checkpoint genes per subset
################################################################################

checkpoint_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2")
Idents(Hep) <- "subcelltype"
cols_subtype <- setNames(cols, levels(Hep@meta.data$subcelltype))

p3h <- VlnPlot(
  Hep,
  features = checkpoint_genes,
  group.by = "subcelltype",
  idents = levels(Hep@meta.data$subcelltype),
  stack = TRUE,
  flip = TRUE,
  fill.by = "ident",
  cols = cols_subtype,
  pt.size = 0
) +
  guides(fill = "none") +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave(file.path(fig_dir, "Fig3h_inhibitory_checkpoint_genes_violin.pdf"), plot = p3h, width = 12, height = 6)

################################################################################
# Fig3i: Differential expression (AFP_Pos vs AFP_Neg) within each subset
################################################################################

deg_list <- list()
for (ct in levels(Hep@meta.data$subcelltype)) {
  sub_obj <- subset(Hep, subcelltype == ct)

  deg <- FindMarkers(
    sub_obj,
    group.by = "AFP_status_group",
    ident.1 = "AFP_Pos",
    ident.2 = "AFP_Neg",
    min.pct = 0.1,
    logfc.threshold = 0,
    only.pos = FALSE
  )
  deg$gene <- rownames(deg)
  deg$cluster <- ct
  deg_list[[ct]] <- deg
}

deg_df <- bind_rows(deg_list)
write.csv(deg_df, file = file.path(fig_dir, "Fig3i_AFP_Pos_vs_AFP_Neg_DEGs_by_subset.csv"), row.names = FALSE)

cut_log2fc <- 0.1
deg_sig <- deg_df %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) >= cut_log2fc) %>%
  mutate(cluster = factor(cluster, levels = plot_order))

gene_counts <- deg_sig %>%
  group_by(cluster) %>%
  summarise(
    upregulated_count = sum(avg_log2FC > 0),
    downregulated_count = sum(avg_log2FC < 0),
    .groups = "drop"
  )

deg_sig <- left_join(deg_sig, gene_counts, by = "cluster")

color <- setNames(cols, plot_order)
new_labels <- c(
  "CD8T_C1", "CD8T_C2", "CD8T_C3", "CD8T_C4", "CD8T_C5",
  "MAIT", "T_stress",
  "CD4T_C1", "CD4T_C2", "CD4T_C3",
  "Treg", "NKT", "Cycling_T"
)

y_max <- max(deg_sig$avg_log2FC, na.rm = TRUE)
p3i <- ggplot(deg_sig, aes(x = cluster, y = avg_log2FC, color = cluster)) +
  geom_quasirandom(size = 0.08, method = "smiley") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  ggtitle("AFP_Pos vs AFP_Neg") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13.5),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 12),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = expression("Log"[2]~"(FC)")) +
  scale_color_manual(values = color) +
  scale_x_discrete(labels = new_labels) +
  geom_text(
    data = gene_counts,
    aes(x = cluster, y = y_max + 0.6, label = upregulated_count),
    color = "#CC3333",
    size = 3.8,
    fontface = "bold",
    vjust = 1.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = gene_counts,
    aes(x = cluster, y = y_max + 0.2, label = downregulated_count),
    color = "#2b2b95",
    size = 3.8,
    fontface = "bold",
    vjust = 1.2,
    inherit.aes = FALSE
  )
ggsave(file.path(fig_dir, "Fig3i_DE_log2FC_by_subset_counts.pdf"), plot = p3i, width = 12, height = 6)
