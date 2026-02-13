
# Fig1b -----------------------------

library(ComplexHeatmap)
library(dplyr)
library(circlize)

setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")

fig_dir <- "./00Figures/Fig1"
data_dir <- "./data"

clinical <- read.csv("./data/Sample_information.csv", check.names = FALSE)

clinical[["AFP..g.L."]] <- log10(as.numeric(clinical[["AFP..g.L."]]) + 1)
clinical[["CA199..U.ml."]] <- log10(as.numeric(clinical[["CA199..U.ml."]]) + 1)
clinical[["CEA..g.L."]] <- log10(as.numeric(clinical[["CEA..g.L."]]) + 1)

colnames(clinical)[colnames(clinical) == "AFP..g.L."] <- "Log10_AFP+1"
colnames(clinical)[colnames(clinical) == "CA199..U.ml."] <- "Log10_CA199+1"
colnames(clinical)[colnames(clinical) == "CEA..g.L."] <- "Log10_CEA+1"

colnames(clinical)[colnames(clinical) == "Relapse_state..Yes.0."] <- "Relapse_state"
colnames(clinical)[colnames(clinical) == "OS_state..Yes.0."] <- "OS_state"

clinical$TNM_stage[clinical$TNM_stage %in% c("IVA", "IV", "IVB")] <- "IV"

clinical <- clinical %>%
  mutate(
    OS_state = recode(as.character(OS_state), `0` = "YES", `1` = "NO"),
    Relapse_state = recode(as.character(Relapse_state), `0` = "YES", `1` = "NO"),
    T = recode(as.character(T), `1` = "T1", `2` = "T2", `3` = "T3", `4` = "T4"),
    N = recode(as.character(N), `0` = "N0", `1` = "N1"),
    M = recode(as.character(M), `0` = "M0", `1` = "M1"),
    AFP_status = recode(as.character(AFP_status), `0` = "AFP_NEG", `1` = "AFP_POS")
  )

new_data <- clinical %>%
  select(AFP_status, `Log10_AFP+1`, everything()) %>%
  filter(!is.na(`Log10_AFP+1`))
new_data <- new_data[order(new_data$AFP_status), ]

zero_row_mat <- matrix(nrow = 0, ncol = nrow(new_data))

median_value <- median(new_data[["Log10_AFP+1"]])
min_value <- min(new_data[["Log10_AFP+1"]])
max_value <- max(new_data[["Log10_AFP+1"]])

col_fun_afp <- colorRamp2(c(min_value, median_value, max_value), c("#1f78b4", "#FFFFFF", "#DC0000FF"))
col_fun_ca199 <- colorRamp2(c(0, 3), c("#FFFFFF", "#EE7942"))
col_fun_cea <- colorRamp2(c(0, 3), c("#FFFFFF", "#FF0000"))
col_fun_pfs <- colorRamp2(c(0, 33), c("#FFFFFF", "#8B4789"))
col_fun_os <- colorRamp2(c(1, 33), c("#FDE725FF", "#9AFF9A"))

ha <- HeatmapAnnotation(
  df = new_data,
  col = list(
    "Log10_AFP+1" = col_fun_afp,
    AFP_status = c(AFP_NEG = "#2874A6", AFP_POS = "#DD3F4E"),
    Differentiation = c(
      Low = "#87CEFA",
      `Low-Moderate` = "#B0E2FF",
      Moderate = "#A4D3EE",
      `Moderate-High` = "#8DB6CD",
      High = "#607B8B"
    ),
    Venous_invasion = c(`0` = "#F0E68C", `1` = "#EE9A49"),
    T = c(T1 = "#8B636C", T2 = "#CD919E", T3 = "#EEA9B8", T4 = "#FFC0CB"),
    N = c(N0 = "#98FB98", N1 = "#87CEEB"),
    M = c(M0 = "#CCE2CB", M1 = "#D8B7DD"),
    TNM_stage = c(I = "#C1FFC1", II = "#B0C4DE", IIIA = "#41ECC5", IV = "#00688B"),
    BCLC_stage = c(D = "#698B69", C = "#9BCD9B", B = "#B4EEB4", A = "#E0EEE0"),
    Cirrhosis = c(`0` = "#FF9999", `1` = "#CDAF95"),
    "Log10_CA199+1" = col_fun_ca199,
    "Log10_CEA+1" = col_fun_cea,
    PFS_time = col_fun_pfs,
    OS_time = col_fun_os,
    OS_state = c(YES = "#FFA07A", NO = "#D1EEEE"),
    Relapse_state = c(YES = "#53868B", NO = "#ADD8E6"),
    CA199_status = c(`0` = "#CCE2CB", `1` = "#339933"),
    CEA_status = c(`1` = "#EBCC96", `0` = "#EEEED1"),
    Virus = c(HBV = "#CD9B1D", HCV = "#EEB422", NBNC = "#FFC125")
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(6, "mm")
)

Hm <- Heatmap(zero_row_mat, top_annotation = ha)
pdf("./Fig1/Fig1b_scRNA_clinical_annotation.pdf", width = 7.2, height = 2.2, useDingbats = FALSE)
draw(
  Hm,
  merge_legend = TRUE,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "left",
  width = unit(16, "cm"),
  height = unit(1, "cm")
)
dev.off()



# Fig1c -----------------------------
# HCC scRNA-seq Survival Analysis
library("survival")
library("survminer")


surival_data<-read.csv("./data/surival_data.csv")
head(surival_data)

surival_data$AFP_status_group[surival_data$AFP_status_group == "AFP_Neg"] <- "AFP-"
surival_data$AFP_status_group[surival_data$AFP_status_group == "AFP_Pos"] <- "AFP+"


surival_data$Group = surival_data$AFP_status_group
table(surival_data$Group)

fit <- survfit(Surv(as.numeric(OS_time), as.numeric(OS_state)) ~ Group, data = surival_data)

ggsurvplot(
  fit, 
  data = surival_data, 
  title = " ",     # title
  font.main = c(14, "bold"),  # title font
  font.x = 18, 
  font.y = 18,  # axis label font
  font.tickslab = 14,  # tick label font
  pval = TRUE,  # show p-value
  pval.size = 6,  # p-value font size
  size = 1.5,  # line width
  linetype = "solid",  # solid lines
  palette = c("#1f78b4", "#DC0000FF"),  # line colors
  legend = c(0.2, 0.45),  # legend position
  legend.title = " ",  # legend title
  legend.labs = c("AFP-(35)", "AFP+(38)"),  # keep counts consistent with sample size
  font.legend = 14, 
  xlab = "Time(Months)"  # x-axis label
)

ggsave(filename = "./00Figures/Fig1/surival_scRNA_HCC_AFP.pdf", width = 4.5,height = 4.3)


# TCGA-LIHC Survival Analysis
library("survival")
library("survminer")

surival_data<-read.table("./data/TCGA_phenotype_AFP_value.txt", sep = "\t", header = TRUE)
head(surival_data)

# Group by AFP_value: >=20 as AFP+ (high), <20 as AFP- (low)
surival_data$Group = ifelse(surival_data$AFP_value >= 20, "AFP+", "AFP-")

fit <- survfit(Surv(as.numeric(OS.time_months), as.numeric(OS.state)) ~ Group, data = surival_data)

ggsurvplot(
  fit, 
  data = surival_data, 
  title = " ",     # title
  font.main = c(14, "bold"),  # title font
  font.x = 18, 
  font.y = 18,  # axis label font
  font.tickslab = 14,  # tick label font
  pval = TRUE,  # show p-value
  pval.size = 6,  # p-value font size
  size = 1.5,  # line width
  linetype = "solid",  # solid lines
  palette = c("#DC0000FF", "#1f78b4"),  # line colors
  legend = c(0.8, 0.85),  # legend position
  legend.title = " ",  # legend title
  legend.labs = c("AFP+(147)", "AFP-(130)"),  # keep counts consistent with sample size
  font.legend = 14, 
  xlab = "Time(Months)"  # x-axis label
)

ggsave(filename = "./00Figures/Fig1/surival_TCGA_LIHC_AFP.pdf", width = 4.5,height = 4.3)
