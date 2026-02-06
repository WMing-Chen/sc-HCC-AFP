############## FigS2_Tumor_SCENIC_and_scMetabolism #############################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(ggsignif)
  library(scMetabolism)
  library(Seurat)
})

fig_dir <- "./FigS2"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cluster_ids <- c(1, 3, 6)
cluster_colors <- setNames(c("#C6EE8F", "#B2DCEE", "#FA7E4D"), paste0("Cluster", cluster_ids))
group_colors <- setNames(c("#db5024", "#25acda"), c("AFP_Pos", "AFP_Neg"))

col_order <- c()
for (k in cluster_ids) {
  col_order <- c(col_order, paste0("C", k, "_AFP_Pos"), paste0("C", k, "_AFP_Neg"))
}

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

################################################################################
# FigS2c: scMetabolism-based activity scores for 11 KEGG metabolic categories
################################################################################

# Input (not tracked; place under ./data/)
# - Tumor Seurat object with metadata `Seu_Clusters` and `AFP_status` (0/1): `HCC_Tumor_10per_sample_definited.rds`

tumor <- readRDS(file.path(data_dir, "HCC_Tumor_10per_sample_definited.rds"))
tumor@meta.data$AFP_status_group <- ifelse(tumor@meta.data$AFP_status == 0, "AFP_Neg", "AFP_Pos")

tumor <- sc.metabolism.Seurat(
  obj = tumor,
  method = "AUCell",
  imputation = FALSE,
  ncores = 2,
  metabolism.type = "KEGG"
)

metabolism_matrix <- as.matrix(tumor@assays$METABOLISM$score)

tumor[["scMetabolism"]] <- SeuratObject::CreateAssayObject(counts = metabolism_matrix)
tumor <- SeuratObject::SetAssayData(
  tumor,
  slot = "scale.data",
  new.data = metabolism_matrix,
  assay = "scMetabolism"
)

tumor_sub <- subset(tumor, subset = Seu_Clusters %in% c(1, 3, 6))

# Map KEGG pathways into 11 metabolism categories (class2) and sum the pathway scores.
kegg_class <- list(
  Metabolism = list(
    `Carbohydrate metabolism` = c(
      "Glycolysis / Gluconeogenesis",
      "Citrate cycle (TCA cycle)",
      "Pentose phosphate pathway",
      "Pentose and glucuronate interconversions",
      "Fructose and mannose metabolism",
      "Galactose metabolism",
      "Ascorbate and aldarate metabolism",
      "Starch and sucrose metabolism",
      "Pyruvate metabolism",
      "Glyoxylate and dicarboxylate metabolism",
      "Propanoate metabolism",
      "Butanoate metabolism",
      "Inositol phosphate metabolism",
      "Amino sugar and nucleotide sugar metabolism"
    ),
    `Energy metabolism` = c(
      "Oxidative phosphorylation",
      "Nitrogen metabolism",
      "Sulfur metabolism"
    ),
    `Lipid metabolism` = c(
      "Fatty acid biosynthesis",
      "Fatty acid elongation",
      "Fatty acid degradation",
      "Synthesis and degradation of ketone bodies",
      "Steroid biosynthesis",
      "Primary bile acid biosynthesis",
      "Steroid hormone biosynthesis",
      "Glycerolipid metabolism",
      "Glycerophospholipid metabolism",
      "Ether lipid metabolism",
      "Sphingolipid metabolism",
      "Arachidonic acid metabolism",
      "Linoleic acid metabolism",
      "alpha-Linolenic acid metabolism",
      "Biosynthesis of unsaturated fatty acids"
    ),
    `Nucleotide metabolism` = c(
      "Purine metabolism",
      "Pyrimidine metabolism"
    ),
    `Amino acid metabolism` = c(
      "Alanine, aspartate and glutamate metabolism",
      "Glycine, serine and threonine metabolism",
      "Cysteine and methionine metabolism",
      "Valine, leucine and isoleucine degradation",
      "Valine, leucine and isoleucine biosynthesis",
      "Lysine degradation",
      "Arginine biosynthesis",
      "Arginine and proline metabolism",
      "Histidine metabolism",
      "Tyrosine metabolism",
      "Phenylalanine metabolism",
      "Tryptophan metabolism",
      "Phenylalanine, tyrosine and tryptophan biosynthesis",
      "D-Glutamine and D-glutamate metabolism",
      "D-Arginine and D-ornithine metabolism"
    ),
    `Metabolism of other amino acids` = c(
      "beta-Alanine metabolism",
      "Taurine and hypotaurine metabolism",
      "Phosphonate and phosphinate metabolism",
      "Selenocompound metabolism",
      "Glutathione metabolism"
    ),
    `Glycan biosynthesis and metabolism` = c(
      "N-Glycan biosynthesis",
      "Mucin type O-glycan biosynthesis",
      "Mannose type O-glycan biosynthesis",
      "Other types of O-glycan biosynthesis",
      "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
      "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
      "Glycosaminoglycan biosynthesis - keratan sulfate",
      "Glycosaminoglycan degradation",
      "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
      "Glycosphingolipid biosynthesis - lacto and neolacto series",
      "Glycosphingolipid biosynthesis - globo and isoglobo series",
      "Glycosphingolipid biosynthesis - ganglio series",
      "Other glycan degradation"
    ),
    `Metabolism of cofactors and vitamins` = c(
      "Thiamine metabolism",
      "Riboflavin metabolism",
      "Vitamin B6 metabolism",
      "Nicotinate and nicotinamide metabolism",
      "Pantothenate and CoA biosynthesis",
      "Biotin metabolism",
      "Lipoic acid metabolism",
      "Folate biosynthesis",
      "One carbon pool by folate",
      "Retinol metabolism",
      "Porphyrin and chlorophyll metabolism",
      "Ubiquinone and other terpenoid-quinone biosynthesis"
    ),
    `Metabolism of terpenoids and polyketides` = c("Terpenoid backbone biosynthesis"),
    `Biosynthesis of other secondary metabolites` = c(
      "Caffeine metabolism",
      "Neomycin, kanamycin and gentamicin biosynthesis"
    ),
    `Xenobiotics biodegradation and metabolism` = c(
      "Metabolism of xenobiotics by cytochrome P450",
      "Drug metabolism - cytochrome P450",
      "Drug metabolism - other enzymes"
    )
  )
)

pathway2class <- data.frame(
  pathway = character(0),
  class1 = character(0),
  class2 = character(0),
  stringsAsFactors = FALSE
)
for (subcat in names(kegg_class$Metabolism)) {
  x <- kegg_class$Metabolism[[subcat]]
  if (length(x) == 0 || all(is.null(x))) next
  pathway2class <- rbind(
    pathway2class,
    data.frame(pathway = x, class1 = "Metabolism", class2 = subcat, stringsAsFactors = FALSE)
  )
}

common_pathway <- intersect(rownames(metabolism_matrix), pathway2class$pathway)
metabolism_matrix_sub <- metabolism_matrix[, Cells(tumor_sub), drop = FALSE]
mat_sub <- metabolism_matrix_sub[common_pathway, , drop = FALSE]
group <- pathway2class$class2[match(common_pathway, pathway2class$pathway)]

metabolism_class_sum <- rowsum(mat_sub, group = group, reorder = TRUE)

tumor_sub[["scMetabolism_Class_sum"]] <- SeuratObject::CreateAssayObject(counts = metabolism_class_sum)
tumor_sub <- SeuratObject::SetAssayData(
  tumor_sub,
  slot = "scale.data",
  new.data = metabolism_class_sum,
  assay = "scMetabolism_Class_sum"
)

DefaultAssay(tumor_sub) <- "scMetabolism_Class_sum"
tumor_sub@meta.data$AFP_Clus <- paste0("C", tumor_sub@meta.data$Seu_Clusters, "_", tumor_sub@meta.data$AFP_status_group)
tumor_sub@meta.data$AFP_Clus <- factor(tumor_sub@meta.data$AFP_Clus, levels = col_order, ordered = TRUE)
p2c <- DotPlot(
  tumor_sub,
  features = rownames(metabolism_class_sum),
  cols = c("lightgrey", "#DC050C"),
  group.by = "AFP_Clus",
  scale = TRUE
) +
  RotatedAxis() +
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_flip()
ggsave(file.path(fig_dir, "FigS2c_scMetabolism_11_categories_AFP_group_dotplot.pdf"), plot = p2c, width = 7.0, height = 6.0, useDingbats = FALSE)

################################################################################
# FigS2d–e: Representative pathway scores (AFP_Neg vs AFP_Pos)
################################################################################

DefaultAssay(tumor_sub) <- "scMetabolism"

p2d <- VlnPlot(
  tumor_sub,
  features = c("Oxidative phosphorylation"),
  group.by = "AFP_status_group",
  stack = FALSE,
  flip = TRUE,
  fill.by = "ident",
  cols = c("#25acda", "#db5024"),
  ncol = 1,
  pt.size = 0
) +
  guides(fill = "none") +
  labs(title = "OXPHOS Scores", x = "", y = "OXPHOS Scores") +
  theme(
    axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  ) +
  scale_x_discrete(limits = c("AFP_Pos", "AFP_Neg")) +
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  geom_boxplot(width = 0.25, col = "black", fill = "white", outlier.shape = NA) +
  geom_signif(
    comparisons = list(c("AFP_Neg", "AFP_Pos")),
    textsize = 5.2,
    y_position = 0.6,
    map_signif_level = TRUE
  )

p2e_left <- VlnPlot(
  tumor_sub,
  features = c("Glycolysis / Gluconeogenesis"),
  group.by = "AFP_status_group",
  stack = FALSE,
  flip = TRUE,
  fill.by = "ident",
  cols = c("#25acda", "#db5024"),
  ncol = 1,
  pt.size = 0
) +
  guides(fill = "none") +
  labs(title = "Glycolysis Scores", x = "", y = "Glycolysis Scores") +
  theme(
    axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  ) +
  scale_x_discrete(limits = c("AFP_Pos", "AFP_Neg")) +
  scale_y_continuous(limits = c(-0.05, 0.4)) +
  geom_boxplot(width = 0.25, col = "black", fill = "white", outlier.shape = NA) +
  geom_signif(
    comparisons = list(c("AFP_Neg", "AFP_Pos")),
    textsize = 5.2,
    y_position = 0.32,
    map_signif_level = TRUE
  )

p2e_right <- VlnPlot(
  tumor_sub,
  features = c("Fatty acid degradation"),
  group.by = "AFP_status_group",
  stack = FALSE,
  flip = TRUE,
  fill.by = "ident",
  cols = c("#25acda", "#db5024"),
  ncol = 1,
  pt.size = 0
) +
  guides(fill = "none") +
  labs(title = "Fatty acid degradation Scores", x = "", y = "Fatty acid degradation Scores") +
  theme(
    axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  ) +
  scale_x_discrete(limits = c("AFP_Pos", "AFP_Neg")) +
  scale_y_continuous(limits = c(-0.1, 0.7)) +
  geom_boxplot(width = 0.25, col = "black", fill = "white", outlier.shape = NA) +
  geom_signif(
    comparisons = list(c("AFP_Neg", "AFP_Pos")),
    textsize = 5.2,
    y_position = 0.6,
    map_signif_level = TRUE
  )

p2de <- p2d | p2e_left | p2e_right
ggsave(file.path(fig_dir, "FigS2d-e_scMetabolism_representative_pathway_scores_violin.pdf"), plot = p2de, width = 12.5, height = 4.9, useDingbats = FALSE)
