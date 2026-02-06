############## FigS4_DC_Treg_CellChat_CXCL #####################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(CellChat)
})

fig_dir <- "./FigS4"
data_dir <- "./data"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Input (not tracked; place under ./data/)
# - `cellchat_merged_DC_subcelltype_Treg.RData`
#   - expected object: `cellchat` (merged CellChat object with group comparison)

load("cellchat_merged_DC_subcelltype_Treg.RData")

pathways_show <- "CXCL"

################################################################################
# FigS4a: CXCL signaling role network (sender/receiver/mediator/influencer)
################################################################################

p_s4a <- netAnalysis_signalingRole_network(
  cellchat,
  signaling = pathways_show,
  width = 10,
  height = 5,
  font.size = 10,
  color.heatmap = "RdPu"
)

pdf(file.path(fig_dir, "FigS4a_CXCL_signaling_role_network.pdf"), width = 10, height = 10, useDingbats = FALSE)
print(p_s4a)
dev.off()

################################################################################
# FigS4b: Ligand–receptor comparison bubble (AFP_Pos vs AFP_Neg)
################################################################################

pdf(file.path(fig_dir, "FigS4b_CXCL_LR_comparison_bubble.pdf"), width = 5, height = 5, useDingbats = FALSE)
netVisual_bubble(
  cellchat,
  sources.use = 14:19,
  targets.use = 3,
  max.quantile = 0.75,
  comparison = c(1, 2),
  angle.x = 45,
  color.text = c("#25acda", "#db5024")
)
dev.off()

