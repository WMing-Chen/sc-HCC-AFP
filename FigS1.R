

# FigS1b Major Celltype Components Bar -------
cols <- c("#E78AC3", "#C6EE8F", "#28B461", "#aa8282", "#d4b7b7", "#55B1B1", "#36489E","#B2DCEE","#FA7E4D","#A680B9","#F0C674","#DA2917")
plot_order = c('CD4+ T cells', 'CD8+ T cells', 'NK cells', 'B cells', "Plasma cells",
               'Macrophage', 'DC', 'Monocyte', 
               'Mast cells','Endothelial', 'Mesenchymal cells', 
               'Tumor')
seu@meta.data$Major_Celltype = factor(seu@meta.data$Major_Celltype, levels = plot_order)
# 分组
# seu@meta.data$CA199_status = paste0("statu_", seu@meta.data$CA199_status)
table(seu$Major_Celltype)
library(scRNAtoolVis)
plot_group <- cellRatioPlot(object = seu,
                            sample.name = "AFP_status_group",
                            celltype.name = "Major_Celltype",
                            fill.col = cols,
                            flow.curve = 0.5 # 连线弯曲程度
)
plot_group
ggsave("./FigS1/Major_Celltype Components Bar.pdf", plot = plot_group, device = "pdf",width = 12,height = 14,units = "cm") 