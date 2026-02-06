############## Fig4_DC_cells ###################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
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
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggpubr)
})

setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")

fig_dir <- "./00Figures/Fig4"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

Hep <- readRDS(file.path(data_dir, "HCC_DC.rds"))

cols_subtype <- c("#36489E", "#C6EE8F", "#28B461", "#B2DCEE", "#DA2917", "#A680B9")
plot_order <- c("DC_C1_CD1C", "DC_C2_STMN1", "DC_C3_CLEC9A", "DC_C4_LAMP3", "DC_C5_CXCL10", "DC_C6_APOC3")
cols_subtype <- setNames(cols_subtype, plot_order)

Hep@meta.data$subcelltype <- factor(Hep@meta.data$subcelltype, levels = plot_order)

################################################################################
# Fig4a: UMAP of DC subsets
################################################################################

p4a <- DimPlot(Hep, group.by = "subcelltype", reduction = "umap", cols = cols_subtype, raster = FALSE) +
  ggtitle("DC subsets") +
  theme_dr(
    xlength = 0.22,
    ylength = 0.22,
    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  theme(panel.grid = element_blank())
ggsave(file.path(fig_dir, "Fig4a_DC_subsets_umap.pdf"), plot = p4a, width = 6, height = 4.5)

################################################################################
# Fig4b: Canonical marker genes dot plot across DC subsets
################################################################################

cluster_gene_pairs <- list(
  DC_C1_CD1C = c("CD1C", "FCER1A", "CLEC10A"),
  DC_C2_STMN1 = c("STMN1", "MKI67", "TOP2A"),
  DC_C3_CLEC9A = c("CLEC9A", "XCR1", "CADM1"),
  DC_C4_LAMP3 = c("LAMP3", "CCR7", "FSCN1"),
  DC_C5_CXCL10 = c("CXCL10", "CXCL9", "CXCL11"),
  DC_C6_APOC3 = c("APOC3", "APOC1", "TTR")
)
marker_genes <- unique(unlist(cluster_gene_pairs))

p4b <- jjDotPlot(
  object = Hep,
  id = "subcelltype",
  gene = marker_genes,
  cluster.order = rev(plot_order),
  gene.order = marker_genes,
  ytree = FALSE,
  legend.position = "right",
  x.text.angle = 45,
  x.text.vjust = 1
) +
  theme(panel.grid = element_blank())
ggsave(file.path(fig_dir, "Fig4b_DC_markers_dotplot.pdf"), plot = p4b, width = 10, height = 12)

################################################################################
# Fig4c: DC subset proportions by AFP group (stacked bars)
################################################################################

p4c <- cellRatioPlot(
  object = Hep,
  sample.name = "AFP_status_group",
  celltype.name = "subcelltype",
  fill.col = unname(cols_subtype),
  flow.curve = 0.5
)
ggsave(file.path(fig_dir, "Fig4c_DC_subtype_proportions_bar.pdf"), plot = p4c, device = "pdf", width = 12, height = 14, units = "cm")

################################################################################
# Fig4d: Ro/e heatmap of DC subsets (AFP_Pos vs AFP_Neg)
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

p4d <- ggplot(plot.data, aes(Var2, forcats::fct_rev(Var1), fill = value)) +
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
ggsave(file.path(fig_dir, "Fig4d_DC_ROE_heatmap.pdf"), plot = p4d, device = "pdf", width = 4.5, height = 6.5, units = "in")

################################################################################
# Fig4e: Box plots of DC_C1 and DC_C5 abundance by AFP group
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

sel_celltypes <- c("DC_C1_CD1C", "DC_C5_CXCL10")
meta_percent <- meta_percent %>% filter(celltype %in% sel_celltypes)

source("./utils/CellType_Group_BoxPlot_Function.R")
cols_group <- rev(c("#327db7", "#f18c8d"))

p4e <- Groups_box_plot(
  data = meta_percent,
  Group_box_colors = cols_group,
  Group_top = "celltype",
  Group_top_order = sel_celltypes,
  Group_box = "Group",
  Group_box_order = c("AFP_Pos", "AFP_Neg"),
  title_name = "DC_C1 and DC_C5 abundance across samples",
  sig_lable = TRUE
)
ggsave(file.path(fig_dir, "Fig4e_DC_C1_C5_abundance_boxplot.pdf"), plot = p4e, device = "pdf", width = 8, height = 4, units = "in")

################################################################################
# Fig4f: DEGs between AFP_Pos and AFP_Neg tumors in each DC subset (volcano)
################################################################################

all_markers_list <- lapply(levels(Hep@meta.data$subcelltype), function(ct) {
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
  deg
})
all_markers <- bind_rows(all_markers_list)
all_markers$cluster <- as.character(all_markers$cluster)

p4f <- jjVolcano(diffData = all_markers, log2FC.cutoff = 0.15)
ggsave(file.path(fig_dir, "Fig4f_AFP_Pos_vs_AFP_Neg_jjVolcano_log2FC_0.15.pdf"), plot = p4f, width = 8, height = 6.5)

################################################################################
# Fig4g: KEGG enrichment of downregulated DEGs in DC_C1 and DC_C5 (AFP_Pos vs AFP_Neg)
################################################################################

deg_kegg_list <- list()
for (ct in c("DC_C1_CD1C", "DC_C5_CXCL10")) {
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
  down_genes <- deg %>% filter(p_val_adj < 0.05, avg_log2FC < 0) %>% pull(gene) %>% unique()

  eg <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ek <- enrichKEGG(gene = eg$ENTREZID, organism = "hsa")
  ek_df <- as.data.frame(ek@result)
  ek_df$Subset <- ct
  deg_kegg_list[[ct]] <- ek_df
}

kegg_df <- bind_rows(deg_kegg_list) %>%
  group_by(Subset) %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

p4g <- ggplot(kegg_df, aes(x = -log10(p.adjust), y = Description)) +
  geom_col(fill = "#1965B0") +
  facet_wrap(~Subset, scales = "free_y") +
  labs(x = "-log10(adj.P)", y = "", title = "KEGG enrichment (downregulated DEGs in AFP_Pos)") +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(fig_dir, "Fig4g_KEGG_downregulated_DEGs_DC_C1_C5.pdf"), plot = p4g, width = 10, height = 5.5)

################################################################################
# Fig4h: Antigen processing and presentation genes in DC_C1 and DC_C5 by AFP group
################################################################################

antigen_genes <- c("HLA-DQA1", "HLA-DRB1", "HLA-B", "HSPA1A", "HSPA6", "HSP90AA1", "HSPA1B", "LGMN")

sub_Hep <- subset(Hep, subcelltype %in% c("DC_C1_CD1C", "DC_C5_CXCL10"))
sub_Hep@meta.data$subcelltype_short <- substr(as.character(sub_Hep@meta.data$subcelltype), 1, 5)
sub_Hep@meta.data$subcelltype_AFPgroup <- paste0(sub_Hep@meta.data$subcelltype_short, "_", sub_Hep@meta.data$AFP_status_group)
plot_order_h <- c("DC_C1_AFP_Pos", "DC_C1_AFP_Neg", "DC_C5_AFP_Pos", "DC_C5_AFP_Neg")
sub_Hep@meta.data$subcelltype_AFPgroup <- factor(sub_Hep@meta.data$subcelltype_AFPgroup, levels = plot_order_h)

p4h <- jjDotPlot(
  object = sub_Hep,
  id = "subcelltype_AFPgroup",
  gene = antigen_genes,
  cluster.order = rev(plot_order_h),
  gene.order = antigen_genes,
  ytree = FALSE,
  legend.position = "right",
  x.text.angle = 45,
  x.text.vjust = 1
) +
  theme(panel.grid = element_blank()) +
  ggtitle("Antigen processing and presentation")
ggsave(file.path(fig_dir, "Fig4h_Antigen_processing_presentation_dotplot.pdf"), plot = p4h, width = 10, height = 5)

