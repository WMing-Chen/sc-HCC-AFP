############## Fig2_Tumor_cells ###############################################

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
library(GSVA)
library(survival)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")

fig_dir <- "./00Figures/Fig2"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Fig2a–b: Tumor clustering and AFP group distribution
################################################################################
Hep <- readRDS(file.path(data_dir, "HCC_Tumor.rds"))
table(Hep@meta.data$AFP_status_group)

Hep <- NormalizeData(Hep) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(verbose = TRUE)

Hep <- RunHarmony(Hep, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
Hep <- RunUMAP(Hep, dims = 1:30, reduction = "harmony", reduction.name = "umap")
Hep <- RunTSNE(Hep, reduction = "harmony", dims = 1:30)

Hep <- FindNeighbors(Hep, reduction = "harmony", dims = 1:30)
Hep <- FindClusters(Hep, resolution = 0.2)
Hep@meta.data$Seu_Clusters <- Idents(Hep)
Hep@meta.data$Seu_Clusters_str <- paste0("Cluster", Hep@meta.data$Seu_Clusters)
Hep@meta.data$Seu_Clusters <- factor(Hep@meta.data$Seu_Clusters, levels = sort(unique(as.integer(as.character(Hep@meta.data$Seu_Clusters)))))

table(Hep@meta.data$subcelltype)

cols <- c(
  "#36489E", "#C6EE8F", "#28B461", "#B2DCEE", "#A680B9", "#DA2917", "#FA7E4D", "#F0C674",
  "#E7298A", "#E78AC3", "#aa8282", "#d4b7b7", "#999999", "#666666", "#A6761D", "#E6AB02"
)

p2a <- DimPlot(Hep, group.by = "subcelltype", reduction = "umap", label = F, cols = cols, raster = FALSE) +
  ggtitle("HCC Tumor") +
  theme_dr(
    xlength = 0.22,
    ylength = 0.22,
    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  theme(panel.grid = element_blank())
p2a
ggsave(file.path(fig_dir, "Fig2a_Tumor_clusters_umap.pdf"), plot = p2a, width = 7, height = 4.5)

p2b <- DimPlot(
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
p2b
ggsave(file.path(fig_dir, "Fig2b_Tumor_clusters_AFP_status.pdf"), plot = p2b, width=6, height = 4.5)

################################################################################
# Fig2c: Canonical marker expression across tumor subsets (dot plot)
################################################################################

Idents(Hep) <- "subcelltype"
Markers_Findall <- FindAllMarkers(
  object = Hep,
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(Markers_Findall, file = file.path(fig_dir, "Tumor_Clusters_Markers_Findall.csv"), row.names = FALSE)


Markers_Findall = read.csv("/home/chenweiming/Project/HCC_scRNAseq/luo/Tumor/Tumor_Clusters_Markers_Findall.csv")
head(Markers_Findall)

Markers_Findall$cluster <- as.character(Markers_Findall$cluster)
top3 <- Markers_Findall %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE)

plot_order <- c("Tumor_C8_PAGE1","Tumor_C7_PCSK1N","Tumor_C9_CD164","Tumor_C3_STMN1","Tumor_C0_FABP1","Tumor_C5_HLA-DRA","Tumor_C1_CYP3A4","Tumor_C2_SOX4","Tumor_C4_MLXIPL","Tumor_C6_AFP")

top3 <- top3 %>%
  mutate(cluster = factor(cluster, levels = plot_order, ordered = TRUE))
top3 <- top3 %>%
  arrange(cluster)
levels(top3$cluster)

Idents(Hep) <- "subcelltype"

DotPlot2 <- jjDotPlot(
  object = Hep,
  id = "subcelltype",
  gene = top3$gene,
  cluster.order = plot_order,
  gene.order = rev(top3$gene),
  ytree = FALSE,
  legend.position = "right",
  x.text.angle = 45,
  x.text.vjust = 1
) +
  theme(panel.grid = element_blank()) +
  coord_flip()
DotPlot2
ggsave(file.path(fig_dir, "Fig2c_Tumor_markers_dotplot.pdf"), plot = DotPlot2, width = 7.5, height = 7.5)

################################################################################
# Tumor subcelltype labels (Tumor_C{cluster}_{top1_marker})
################################################################################

top1 <- Markers_Findall %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1, with_ties = FALSE)

sub_name <- paste0("Tumor_C", top1$cluster, "_", top1$gene)
cluster_mapping <- setNames(sub_name, top1$cluster)
Hep@meta.data$tumor_subcelltype <- as.character(Hep@meta.data$Seu_Clusters)
Hep@meta.data$tumor_subcelltype <- cluster_mapping[Hep@meta.data$tumor_subcelltype]
Hep@meta.data$tumor_subcelltype <- factor(Hep@meta.data$tumor_subcelltype, levels = sub_name)

################################################################################
# Fig2d: Ro/e heatmap (AFP_Pos vs AFP_Neg)
################################################################################

source("./utils/ROE_Heatmap_Function.R")
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:7])

plot.data <- ROIE(table(Hep@meta.data[, c("tumor_subcelltype", "AFP_status_group")])) %>%
  reshape2::melt() %>%
  mutate(value = pmin(value, 3))

plot.data <- plot.data %>%
  mutate(
    Var1 = factor(
      Var1,
      levels = as.character(unique(Var1[Var2 == "AFP_Pos"][order(-value[Var2 == "AFP_Pos"])]))
    )
  )

p2d <- ggplot(plot.data, aes(Var2, forcats::fct_rev(Var1), fill = value)) +
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
ggsave(file.path(fig_dir, "Fig2d_Tumor_ROE_heatmap.pdf"), plot = p2d, device = "pdf", width = 7.5, height = 9, units = "in")

################################################################################
# Fig2e: Box plots of Tumor_C1/C3/C6 abundance by AFP group
################################################################################

meta <- Hep@meta.data %>%
  dplyr::select(orig.ident, AFP_status_group, tumor_subcelltype)
colnames(meta) <- c("orig.ident", "Group", "celltype")

meta_percent <- meta %>%
  group_by(orig.ident, Group, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(values = count / sum(count)) %>%
  dplyr::select(orig.ident, values, Group, celltype)
meta_percent <- unique(meta_percent)

sel_celltypes <- c(
  levels(Hep@meta.data$tumor_subcelltype)[grepl("^Tumor_C1_", levels(Hep@meta.data$tumor_subcelltype))],
  levels(Hep@meta.data$tumor_subcelltype)[grepl("^Tumor_C3_", levels(Hep@meta.data$tumor_subcelltype))],
  levels(Hep@meta.data$tumor_subcelltype)[grepl("^Tumor_C6_", levels(Hep@meta.data$tumor_subcelltype))]
)
meta_percent <- meta_percent %>% filter(celltype %in% sel_celltypes)

source("./utils/CellType_Group_BoxPlot_Function.R")
cols_group <- rev(c("#327db7", "#f18c8d"))

p2e <- Groups_box_plot(
  data = meta_percent,
  Group_box_colors = cols_group,
  Group_top = "celltype",
  Group_top_order = sel_celltypes,
  Group_box = "Group",
  Group_box_order = c("AFP_Pos", "AFP_Neg"),
  title_name = "Tumor subcelltype percentages across samples",
  sig_lable = TRUE
)
ggsave(file.path(fig_dir, "Fig2e_Tumor_C1_C3_C6_abundance_boxplot.pdf"), plot = p2e, device = "pdf", width = 12, height = 5, units = "in")

################################################################################
# Fig2f–g: TCGA-LIHC (GSVA signature scores + survival)
################################################################################

tcga_expr <- read.table(file.path(data_dir, "TCGA_matrix.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
tcga_expr <- as.matrix(tcga_expr)

tcga_pheno <- read.table(file.path(data_dir, "TCGA_phenotype_AFP_value.txt"), header = TRUE, sep = "\t", check.names = FALSE)
tcga_pheno$AFP_Group <- factor(tcga_pheno$AFP_Group, levels = c("AFP+", "AFP-"))
rownames(tcga_pheno) <- tcga_pheno$sample_id

top50 <- Markers_Findall %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE)

geneset_list <- list(
  Tumor_C1 = unique(top50$gene[top50$cluster %in% c(1)]),
  Tumor_C3 = unique(top50$gene[top50$cluster %in% c(3)]),
  Tumor_C6 = unique(top50$gene[top50$cluster %in% c(6)])
)

ES <- gsva(tcga_expr, geneset_list)
ES <- as.data.frame(t(ES))
ES$AFP_Group <- tcga_pheno[rownames(ES), "AFP_Group"]

ES_long <- ES %>%
  tibble::rownames_to_column("sample_id") %>%
  pivot_longer(cols = c("Tumor_C1", "Tumor_C3", "Tumor_C6"), names_to = "Signature", values_to = "Score")

p2f <- ggplot(ES_long, aes(x = AFP_Group, y = Score, fill = AFP_Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.2, color = "black", alpha = 0.6) +
  facet_wrap(~Signature, scales = "free_y") +
  scale_fill_manual(values = c("AFP+" = "#DC050C", "AFP-" = "#386cb0")) +
  labs(x = "", y = "GSVA score") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(fig_dir, "Fig2f_TCGA_Tumor_C1_C3_C6_signature_boxplots.pdf"), plot = p2f, width = 9, height = 4.5)

surival_data <- tcga_pheno
for (sig in c("Tumor_C1", "Tumor_C3", "Tumor_C6")) {
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
  ggsave(filename = file.path(fig_dir, paste0("Fig2g_TCGA_survival_", sig, ".pdf")), plot = p$plot, width = 4.5, height = 4.3)
}

################################################################################
# Fig2h: GO BP enrichment of Tumor_C1/C3/C6 marker genes (scRNA-seq)
################################################################################

go_list <- list(
  Tumor_C1 = geneset_list$Tumor_C1,
  Tumor_C3 = geneset_list$Tumor_C3,
  Tumor_C6 = geneset_list$Tumor_C6
)

go_res <- list()
for (nm in names(go_list)) {
  ego <- enrichGO(
    gene = go_list[[nm]],
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  res <- as.data.frame(ego@result)
  if (nrow(res) > 0) {
    res$Subset <- nm
    go_res[[nm]] <- res
  }
}

go_df <- bind_rows(go_res) %>%
  group_by(Subset) %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup()

p2h_list <- list()
for (nm in unique(go_df$Subset)) {
  tmp <- go_df %>% filter(Subset == nm) %>% arrange(p.adjust)
  tmp$Description <- factor(tmp$Description, levels = rev(tmp$Description))
  p_tmp <- ggplot(tmp, aes(x = -log10(p.adjust), y = Description)) +
    geom_col(fill = "#C50B05") +
    labs(x = "-log10(adj.P)", y = "", title = nm) +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
  p2h_list[[nm]] <- p_tmp
  ggsave(file.path(fig_dir, paste0("Fig2h_GO_BP_top_markers_", nm, ".pdf")), plot = p_tmp, width = 9, height = 4)
}

if (length(p2h_list) > 0) {
  p2h_all <- plot_grid(plotlist = p2h_list, ncol = 1, align = "v")
  ggsave(file.path(fig_dir, "Fig2h_GO_BP_top_markers_combined.pdf"), plot = p2h_all, width = 9, height = 4 * length(p2h_list))
}

################################################################################
# Fig2i: GO BP enrichment of TCGA DEGs (AFP+ vs AFP-)
################################################################################

library(limma)
common_samples <- intersect(colnames(tcga_expr), rownames(tcga_pheno))
expr_use <- tcga_expr[, common_samples, drop = FALSE]
pheno_use <- tcga_pheno[common_samples, , drop = FALSE]

design <- model.matrix(~ AFP_Group, data = pheno_use)
fit <- lmFit(expr_use, design)
fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
tt$gene <- rownames(tt)

deg <- tt %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)
if (nrow(deg) == 0) {
  stop("No TCGA DEGs found (adj.P.Val < 0.05 and abs(logFC) > 0.25).")
}
ego_deg <- enrichGO(
  gene = deg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
)
deg_go <- as.data.frame(ego_deg@result)
if (nrow(deg_go) == 0) {
  stop("No GO BP enrichment results for TCGA DEGs.")
}
deg_go <- deg_go %>%
  slice_min(order_by = p.adjust, n = 15, with_ties = FALSE) %>%
  mutate(Description = reorder(Description, -log10(p.adjust)))

p2i <- ggplot(deg_go, aes(x = -log10(p.adjust), y = Description)) +
  geom_col(fill = "#1965B0") +
  labs(x = "-log10(adj.P)", y = "", title = "TCGA-LIHC AFP+ vs AFP- DEGs (GO BP)") +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(fig_dir, "Fig2i_TCGA_DEG_GO_BP.pdf"), plot = p2i, width = 9, height = 5.5)
