

# FigS1a -----------------------------
# TCGA-LIHC clinical characteristics

library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(grid)

dir.create("./FigS1", showWarnings = FALSE, recursive = TRUE)

data <- read.csv("./data/TCGA_LIHC.csv", check.names = FALSE)

data$AFP_status <- ifelse(as.numeric(data$AFP) >= 20, "AFP_POS", "AFP_NEG")

data <- data %>%
  mutate(
    Virus = case_when(
      is.na(Virus) ~ NA_character_,
      grepl("Hepatitis B", Virus) & grepl("Hepatitis C", Virus) ~ "BandC",
      grepl("Hepatitis B", Virus) ~ "HBV",
      grepl("Hepatitis C", Virus) ~ "HCV",
      grepl("Other", Virus, ignore.case = TRUE) ~ "NBNC",
      grepl("Alcohol consumption", Virus, ignore.case = TRUE) ~ "NBNC",
      TRUE ~ "NBNC"
    ),
    Cirrhosis = case_when(
      fibrosis_ishak_score == "0 - No Fibrosis" ~ "NO",
      !is.na(fibrosis_ishak_score) ~ "YES",
      TRUE ~ NA_character_
    ),
    AFP = log10(as.numeric(AFP) + 1),
    OS_state = case_when(
      OS_status == "LIVING" ~ "YES",
      OS_status == "DECEASED" ~ "NO",
      TRUE ~ as.character(OS_status)
    ),
    TNM_stage = case_when(
      TNM_stage %in% c("Stage I") ~ "I",
      TNM_stage %in% c("Stage II") ~ "II",
      TNM_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "III",
      TNM_stage %in% c("Stage IV", "Stage IVA", "Stage IVB") ~ "IV",
      TRUE ~ as.character(TNM_stage)
    ),
    T = case_when(
      T %in% c("T2", "T2a", "T2b") ~ "T2",
      T %in% c("T3", "T3a", "T3b") ~ "T3",
      TRUE ~ as.character(T)
    )
  )

colnames(data)[colnames(data) == "AFP"] <- "Log10_AFP+1"

new_data <- data %>%
  select(AFP_status, `Log10_AFP+1`, everything()) %>%
  filter(!is.na(`Log10_AFP+1`))
new_data <- new_data[order(new_data$AFP_status), ]

zero_row_mat <- matrix(nrow = 0, ncol = nrow(new_data))

median_value <- median(new_data[["Log10_AFP+1"]])
min_value <- min(new_data[["Log10_AFP+1"]])
max_value <- max(new_data[["Log10_AFP+1"]])
col_fun_afp <- colorRamp2(c(min_value, median_value, max_value), c("#1f78b4", "#FFFFFF", "#DC0000FF"))

col_fun_os_time <- colorRamp2(c(0, 3600), c("#FFFFFF", "#8B4789"))

ha <- HeatmapAnnotation(
  df = new_data,
  col = list(
    "Log10_AFP+1" = col_fun_afp,
    AFP_status = c(AFP_NEG = "#2874A6", AFP_POS = "#DD3F4E"),
    OS_state = c(YES = "#D1EEEE", NO = "#FFA07A"),
    Cirrhosis = c(NO = "#CDAF95", YES = "#FF9999"),
    T = c(T1 = "#8B636C", T2 = "#CD919E", T3 = "#EEA9B8", T4 = "#FFC0CB"),
    N = c(N0 = "#CD9B1D", N1 = "#EEB422", NX = "#FFC125"),
    M = c(M0 = "#607B8B", M1 = "#8DB6CD", MX = "#A4D3EE"),
    New_tumor_event = c(YES = "#53868B", NO = "#ADD8E6"),
    Vascular_invasion = c(Micro = "#FFEC8B", Macro = "#FFC125", None = "#EEEED1"),
    OS_time = col_fun_os_time,
    Virus = c(BandC = "#698B69", HBV = "#9BCD9B", HCV = "#B4EEB4", NBNC = "#E0EEE0"),
    TNM_stage = c(I = "#C1FFC1", II = "#B0C4DE", III = "#41ECC5", IV = "#00688B")
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(6, "mm")
)

Hm <- Heatmap(zero_row_mat, top_annotation = ha)
pdf("./FigS1/FigS1a_TCGA_LIHC_clinical_annotation.pdf", width = 7.2, height = 2.2, useDingbats = FALSE)
draw(
  Hm,
  merge_legend = TRUE,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  width = unit(16, "cm"),
  height = unit(1, "cm")
)
grid.text("TCGA-LIHC", x = 0.5, y = 0.9, gp = gpar(fontsize = 20, fontface = "bold"))
dev.off()


# FigS1b Major Celltype Components Bar -------
cols <- c("#E78AC3", "#C6EE8F", "#28B461", "#aa8282", "#d4b7b7", "#55B1B1", "#36489E","#B2DCEE","#FA7E4D","#A680B9","#F0C674","#DA2917")
plot_order = c('CD4+ T cells', 'CD8+ T cells', 'NK cells', 'B cells', "Plasma cells",
               'Macrophage', 'DC', 'Monocyte', 
               'Mast cells','Endothelial', 'Mesenchymal cells', 
               'Tumor')
seu@meta.data$Major_Celltype = factor(seu@meta.data$Major_Celltype, levels = plot_order)
# Grouping
table(seu$Major_Celltype)
library(scRNAtoolVis)
plot_group <- cellRatioPlot(object = seu,
                            sample.name = "AFP_status_group",
                            celltype.name = "Major_Celltype",
                            fill.col = cols,
                            flow.curve = 0.5 # curve strength
)
plot_group
ggsave("./FigS1/Major_Celltype Components Bar.pdf", plot = plot_group, device = "pdf",width = 12,height = 14,units = "cm") 
