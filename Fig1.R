
# Fig1b -----------------------------



# Fig1c -----------------------------
# HCC scRNA-seq Survival Analysis
library("survival")
library("survminer")


surival_data<-read.csv("/home/chenweiming/Project/HCC_scRNAseq/luo/data/surival_data.csv")
surival_data$Group = surival_data$AFP_status_group

fit <- survfit(Surv(as.numeric(OS.time_months), as.numeric(OS.state)) ~ Group, data = surival_data)

ggsurvplot(
  fit, 
  data = surival_data, 
  title = " ",     # 更改标题 
  font.main = c(14, "bold"),  # 标题字体
  font.x = 14, 
  font.y = 14,  # 坐标轴标题字体大小
  font.tickslab = 12,  # 刻度字体大小
  pval = TRUE,  # 添加 p 值
  pval.size = 4.5,  # 修改 p 值字体大小
  size = 1.2,  # 线条大小
  linetype = "solid",  # 设置所有线条为实线
  palette = c("#DC0000FF", "#1f78b4"),  # 线条颜色
  legend = c(0.85, 0.85),  # 改变 legend 位置
  legend.title = " ",  # 改变 legend 的标题
  legend.labs = c("AFP+(147)", "AFP-(130)"),  # !!这里数量要对应上面分组样本数
  xlab = "Time(Months)"  # 修改 x 轴标题
)

ggsave(filename = "./Fig1/surival_scRNA_HCC_AFP.pdf", width = 4.5,height = 4.3)


# TCGA-LIHC Survival Analysis
library("survival")
library("survminer")

surival_data<-read.table("/home/chenweiming/Project/HCC_scRNAseq/luo/data/TCGA_phenotype_AFP_value.txt", sep = "\t")

# 根据AFP_value列分组，>=20为AFP+，高表达组，<20为AFP-，低表达组
surival_data$Group = ifelse(surival_data$AFP_value >= 20, "AFP+", "AFP-")

fit <- survfit(Surv(as.numeric(OS.time_months), as.numeric(OS.state)) ~ Group, data = surival_data)

ggsurvplot(
  fit, 
  data = surival_data, 
  title = " ",     # 更改标题 
  font.main = c(14, "bold"),  # 标题字体
  font.x = 14, 
  font.y = 14,  # 坐标轴标题字体大小
  font.tickslab = 12,  # 刻度字体大小
  pval = TRUE,  # 添加 p 值
  pval.size = 4.5,  # 修改 p 值字体大小
  size = 1.2,  # 线条大小
  linetype = "solid",  # 设置所有线条为实线
  palette = c("#DC0000FF", "#1f78b4"),  # 线条颜色
  legend = c(0.85, 0.85),  # 改变 legend 位置
  legend.title = " ",  # 改变 legend 的标题
  legend.labs = c("AFP+(147)", "AFP-(130)"),  # !!这里数量要对应上面分组样本数
  xlab = "Time(Months)"  # 修改 x 轴标题
)

ggsave(filename = "./Fig1/surival_TCGA_LIHC_AFP.pdf", width = 4.5,height = 4.3)
