############## CellTrek_processing ##############################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(CellTrek)
  library(dplyr)
  library(Seurat)
})

setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")

data_dir <- "./data"

# Input (not tracked in repo; place under ./data/)
# - scRNA reference Seurat object: `HCC_10per_cells_definited.rds`
#   - expected meta columns: `Major_Celltype`, `sub_cell_type`
# - Spatial transcriptomics Seurat list: `ST_list.rds`
#
# Output (written to ./data/)
# - `ST_list_CellTrek.rds`: ST Seurat list with CellTrek results embedded
#   - for each sample object:
#     - `@misc$CellTrek`: CellTrek-mapped Seurat object
#     - `@misc$Treg_kdist`: per-cell k-distance table (DC subsets -> Treg)
#     - `@misc$Treg_kdist_summary`: per-celltype summary (median/mean/n)
#   - list attributes (combined tables across samples):
#     - `kdist_dc_cells`, `kdist_dc_cells_scaled`
#     - `kdist_summary`, `kdist_summary_scaled`

################################################################################
# Prepare scRNA reference for CellTrek
################################################################################

sc_ref <- readRDS(file.path(data_dir, "HCC_10per_cells_definited.rds"))

sc_ref@meta.data$CellTrek_Cell_type <- as.character(sc_ref@meta.data$Major_Celltype)

idx <- sc_ref@meta.data$CellTrek_Cell_type == "CD4+ T cells"
sc_ref@meta.data$CellTrek_Cell_type[idx] <- sc_ref@meta.data$sub_cell_type[idx]

idx <- sc_ref@meta.data$CellTrek_Cell_type == "CD8+ T cells"
sc_ref@meta.data$CellTrek_Cell_type[idx] <- sc_ref@meta.data$sub_cell_type[idx]

idx <- sc_ref@meta.data$CellTrek_Cell_type == "DC"
sc_ref@meta.data$CellTrek_Cell_type[idx] <- sc_ref@meta.data$sub_cell_type[idx]

sc_ref@meta.data$CellTrek_Cell_type[sc_ref@meta.data$CellTrek_Cell_type == "Mesenchymal cells"] <- "MSCs"
sc_ref@meta.data$CellTrek_Cell_type[sc_ref@meta.data$CellTrek_Cell_type %in% c("Cycling_T_cells")] <- "Cycling_T"

sc_ref <- RenameCells(sc_ref, new.names = make.names(Cells(sc_ref)))

################################################################################
# Load ST objects
################################################################################

ST_list <- readRDS(file.path(data_dir, "ST_list.rds"))

samples_ST <- names(ST_list)

dc_subtypes <- c("DC_C1_CD1C", "DC_C2_STMN1", "DC_C3_CLEC9A", "DC_C4_LAMP3", "DC_C5_CXCL10", "DC_C6_APOC3")
k_nn <- 10

################################################################################
# CellTrek mapping + k-distance (DC subsets -> Treg)
################################################################################

all_results <- list()
all_med_results <- list()

for (sample_sel in samples_ST) {
  st_obj <- ST_list[[sample_sel]]
  st_obj <- RenameCells(st_obj, new.names = make.names(Cells(st_obj)))

  st_sc_int <- CellTrek::traint(
    st_data = st_obj,
    sc_data = sc_ref,
    sc_assay = "RNA",
    cell_names = "CellTrek_Cell_type"
  )

  celltrek_obj <- CellTrek::celltrek(
    st_sc_int = st_sc_int,
    int_assay = "traint",
    sc_data = sc_ref,
    sc_assay = "RNA",
    reduction = "pca",
    intp = TRUE,
    intp_pnt = 5000,
    intp_lin = FALSE,
    nPCs = 30,
    ntree = 1000,
    dist_thresh = 0.55,
    top_spot = 2,
    spot_n = 10,
    repel_r = 20,
    repel_iter = 20,
    keep_model = TRUE
  )$celltrek

  inp_df <- celltrek_obj@meta.data %>%
    dplyr::select(cell_names = dplyr::one_of("CellTrek_Cell_type"), coord_x, coord_y)

  inp_df$coord_x <- max(inp_df$coord_x, na.rm = TRUE) - inp_df$coord_x

  ref_celltype <- "Treg"
  other_cell <- dc_subtypes

  kdist_out <- CellTrek::kdist(
    inp_df = inp_df,
    ref = ref_celltype,
    ref_type = "all",
    que = other_cell,
    k = k_nn,
    new_name = paste0(ref_celltype, ".Others"),
    keep_nn = FALSE
  )

  res <- kdist_out$kdist_df
  res$barcode <- rownames(res)
  inp_df$barcode <- rownames(inp_df)
  res <- dplyr::left_join(res, inp_df, by = "barcode")

  summary_tbl <- res %>%
    group_by(cell_names) %>%
    summarise(
      n = dplyr::n(),
      median_dist = median(get(paste0(ref_celltype, ".Others")), na.rm = TRUE),
      mean_dist = mean(get(paste0(ref_celltype, ".Others")), na.rm = TRUE),
      .groups = "drop"
    )

  ST_list[[sample_sel]]@misc$CellTrek <- celltrek_obj
  ST_list[[sample_sel]]@misc$Treg_kdist <- res
  ST_list[[sample_sel]]@misc$Treg_kdist_summary <- summary_tbl

  all_results[[sample_sel]] <- res
  all_med_results[[sample_sel]] <- summary_tbl %>% mutate(sample = sample_sel)
}

combined_results <- bind_rows(all_results, .id = "sample")
combined_results_scaled <- combined_results %>%
  group_by(sample) %>%
  mutate(Treg.Others_zscore = scale(Treg.Others)) %>%
  ungroup()
colnames(combined_results_scaled)[colnames(combined_results_scaled) == "Treg.Others_zscore[,1]"] <- "Treg.Others_zscore"

combined_summary <- bind_rows(all_med_results)
combined_summary_scaled <- combined_summary %>%
  group_by(sample) %>%
  mutate(
    mean_dist_zscore = scale(mean_dist),
    median_dist_zscore = scale(median_dist)
  ) %>%
  ungroup()
colnames(combined_summary_scaled)[colnames(combined_summary_scaled) == "mean_dist_zscore[,1]"] <- "mean_dist_zscore"
colnames(combined_summary_scaled)[colnames(combined_summary_scaled) == "median_dist_zscore[,1]"] <- "median_dist_zscore"

attr(ST_list, "kdist_dc_cells") <- combined_results
attr(ST_list, "kdist_dc_cells_scaled") <- combined_results_scaled
attr(ST_list, "kdist_summary") <- combined_summary
attr(ST_list, "kdist_summary_scaled") <- combined_summary_scaled

saveRDS(ST_list, file.path(data_dir, "ST_list_CellTrek.rds"))
