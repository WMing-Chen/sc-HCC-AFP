library(GO.db)
library(GenomeInfoDbData)
library(HDO.db)
library(clusterProfiler)
library(AnnotationHub) 
library(org.Hs.eg.db) 
library(org.Mm.eg.db) 
library(ggplot2)
library(DOSE)
library(ggplot2)
library(tidyverse)
library(dplyr)

Run_GO_KEGG <- function(DEGs_TAMsubsets2, wz = "has", sel_clus, updown = F){

  if (wz == "hsa"){
    org_db = org.Hs.eg.db
    keytypes(org_db)
    # org.Hs.egENSEMBL, org.Mm.egENSEMBL 
    org_ENS = org.Hs.egENSEMBL
    # Kegg 'hsa' 'mmu' 
    organism = 'hsa'
  }else if (wz == "mmu"){
    org_db = org.Mm.eg.db
    keytypes(org_db)
    # org.Hs.egENSEMBL, org.Mm.egENSEMBL 
    org_ENS = org.Mm.egENSEMBL
    organism = 'mmu'
  }else{
    print("'wz'--Species parameter error")
  }
  
  
  symbol_to_Entrez_Gene_ID <- function(gene_symbol_vector, org_db, org_ENS) {
    gene_id_trans <- bitr(gene_symbol_vector, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org_db, drop = TRUE)
    
    ensembl_df <- as.data.frame(gene_id_trans[, 2])
    colnames(ensembl_df)[1] <- ''
    colnames(ensembl_df) <- "V1"
    
    ensembl_column <- ensembl_df$V1
    
    EG2Ensembl <- toTable(org_ENS)
    gene_lists <- data.frame(ensembl_id = ensembl_column)
    results <- merge(gene_lists, EG2Ensembl, by = 'ensembl_id', all.x = TRUE)
    
    gene_id <- na.omit(results$gene_id)
    
    return(gene_id)
  }
  
  if (!is.null(sel_clus)) {
    sel_clus_gene <- DEGs_TAMsubsets2[which(DEGs_TAMsubsets2$cluster %in% sel_clus), ]
  }else{
    sel_clus_gene <- DEGs_TAMsubsets2
  }
  
  if (updown){
    sel_clus_gene$group <- ifelse(sel_clus_gene$avg_log2FC > 0, "UP", "DOWN")
    sel_up = sel_clus_gene[which(sel_clus_gene$group == "UP"), ]
    sel_down = sel_clus_gene[which(sel_clus_gene$group == "DOWN"), ]
    
    go_ythdf2 <- sel_up$gene 
    gene_id = symbol_to_Entrez_Gene_ID(go_ythdf2, org_db, org_ENS)
    GO <- enrichGO(OrgDb=org_db, gene = gene_id, ont = "ALL",  readable= TRUE) 
    GO_UP = GO@result
    GO_UP$group = rep("UP", times = dim(GO_UP)[1])
    KEGG <- enrichKEGG(gene= gene_id, organism  = organism, pvalueCutoff = 0.05)

    KEGG = setReadable(KEGG, 
                       OrgDb = org_db, 
                       keyType = "ENTREZID") 
    KEGG@result$Description <- gsub("- Mus musculus \\(house mouse\\)", "", KEGG@result$Description)
    KEGG_UP = KEGG@result
    KEGG_UP$group = rep("UP", times = dim(KEGG_UP)[1])
    
    go_ythdf2 <- sel_down$gene 
    gene_id = symbol_to_Entrez_Gene_ID(go_ythdf2, org_db, org_ENS)
    GO <- enrichGO(OrgDb=org_db, gene = gene_id, ont = "ALL",  readable= TRUE) 
    GO_DOWN = GO@result
    GO_DOWN$group = rep("DOWN", times = dim(GO_DOWN)[1])
    KEGG <- enrichKEGG(gene= gene_id, organism  = organism, pvalueCutoff = 0.05)
    KEGG = setReadable(KEGG, 
                       OrgDb = org_db, 
                       keyType = "ENTREZID") 
    KEGG@result$Description <- gsub("- Mus musculus \\(house mouse\\)", "", KEGG@result$Description)
    KEGG_DOWN = KEGG@result
    KEGG_DOWN$group = rep("DOWN", times = dim(KEGG_DOWN)[1])
    
    GO_result = rbind(GO_UP, GO_DOWN)
    KEGG_result = rbind(KEGG_UP, KEGG_DOWN)
    
  }else{
    go_ythdf2 <- sel_clus_gene$gene
    gene_id = symbol_to_Entrez_Gene_ID(go_ythdf2, org_db, org_ENS)
    ####1 GO
    GO <- enrichGO(OrgDb=org_db, gene = gene_id, ont = "ALL",  readable= TRUE) 
    GO_result = GO@result
    
    ######KEGG -------------------------
    KEGG <- enrichKEGG(gene= gene_id, organism  = organism, pvalueCutoff = 0.05)
    KEGG = setReadable(KEGG, 
                       OrgDb = org_db, 
                       keyType = "ENTREZID") 
    KEGG@result$Description <- gsub("- Mus musculus \\(house mouse\\)", "", KEGG@result$Description)
    KEGG_result = KEGG@result
  }
  
  ONTOLOGY = rep("KEGG", times = dim(KEGG_result)[1])
  KEGG_result = cbind(ONTOLOGY, KEGG_result)
  GO_KEGG_result = rbind(GO_result, KEGG_result)
  
  return(GO_KEGG_result)
}



plot_GO_bilateral_bar <- function(GOresult1, tile_name, ONTOLOGY, x_value = "Conut") {
  # Function to create a bilateral bar plot for Gene Ontology (GO) enrichment results.
  
  # Arguments:
  # - GOresult1: Data frame containing the GO enrichment results. (Output from Run_GO_KEGG())
  # - tile_name: Character, part of the title indicating the cluster name.
  # - ONTOLOGY: Character, the type of ontology to display ("BP", "MF", "CC", "KEGG").
  # - x_value: Character, the variable to be plotted on the x-axis ("Count", "GeneRatio", "pvalue", "qvalue", or "p.adjust").
  
  GOresult1$p.adjust <- as.numeric(GOresult1$p.adjust)
  GOresult1$p.adjust <- sprintf("%.2e", GOresult1$p.adjust)
  
  filter_condition <- GOresult1$ONTOLOGY == ONTOLOGY
  GOresult <- GOresult1[filter_condition, ]
  
  Descri <- as.data.frame(table(GOresult$Description))
  Des <- Descri[Descri$Freq > 1, ]$Var1
  filter_condition <- which(GOresult$Description %in% Des & GOresult$group == "DOWN")
  if (length(filter_condition) != 0) {
    GOresult <- GOresult[-filter_condition, ]
  }
  
  GOresult$xvalue <- switch(
    x_value,
    Count = GOresult$Count,
    GeneRatio = sapply(strsplit(as.character(GOresult$GeneRatio), "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
    pvalue = -log10(GOresult$pvalue + 1e-1000),
    qvalue = -log10(GOresult$qvalue + 1e-1000),
    p.adjust = -log10(as.numeric(GOresult$p.adjust) + 1e-1000),
    stop("Invalid x_value. Supported values: 'Count', 'GeneRatio', 'pvalue', 'qvalue', 'p.adjust'")
  )
  if (x_value %in% c("pvalue", "qvalue", "p.adjust")){
    x_value = paste0("-log10(",x_value,")")
  }
  
  max_xvalue <- max(GOresult$xvalue)
  max_10 <- ceiling(max_xvalue/10) * 10
  
  GOresult[GOresult$group == "DOWN", ]$xvalue <- -GOresult[GOresult$group == "DOWN", ]$xvalue
  GOresult <- GOresult[order(GOresult$xvalue), ]
  GOresult <- rbind(head(GOresult, 10), tail(GOresult, 10))
  
  gg_title <- if (ONTOLOGY == "KEGG") {
    paste0("KEGG | Cluster ", tile_name)
  } else {
    paste0("GO ", ONTOLOGY, " | Cluster ", tile_name)
  }
  
  p <- ggplot(GOresult, aes(reorder(Description, xvalue), xvalue, fill = group)) +
    geom_col() +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.title = element_blank(),
      legend.position = 'none',
      axis.text = element_text(color = "black", size = 10),
      axis.line.x = element_line(color = 'black'),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    coord_flip() +
    geom_segment(aes(y = 0, yend = 0, x = 0, xend = 20.5)) +
    geom_text(data = GOresult[GOresult$group == "UP", ], aes(x = Description, y = -0.01, label = Description), hjust = 1.005, size = 4.5) +
    geom_text(data = GOresult[GOresult$group == "DOWN", ], aes(x = Description, y = 0.01, label = Description), hjust = -0.01, size = 4.5) +
    scale_fill_manual(values = c("#1965B0", "#C50B05")) +
    scale_y_continuous(
      breaks = c(-max_10, -max_10 * 0.5, 0, max_10 * 0.5, max_10),
      labels = c(max_10, max_10 * 0.5, 0, max_10 * 0.5, max_10),
      limits = c(-max_10 * 1.1, max_10 * 1.1)
    ) +
    labs(x = '', y = x_value) +
    ggtitle(gg_title)
  
  return(p)
}



plot_GO_bar <- function(GOresult1, tile_name, ONTOLOGY, x_value = "Conut", col = "#C50B05", top_n = 10) {
  # Function to create a ggplot bar plot for the top genes based on specified criteria.
  
  # Arguments:
  # - GOresult1: Data frame containing the GO enrichment results.
  # - tile_name: Character, name of the cluster for which the plot is generated.
  # - ONTOLOGY: Character, the ontology type ("BP", "MF", "CC", or "KEGG").
  # - x_value: Character, the variable to be plotted on the x-axis ("Count", "GeneRatio", "pvalue", "qvalue", or "p.adjust").
  # - top_n: Numeric, number of top genes to display in the plot (default is 10).
  
  # Filter data based on ONTOLOGY
  filter_condition <- GOresult1$ONTOLOGY == ONTOLOGY
  GOresult <- GOresult1[filter_condition, ]
  
  # Set x-axis variable
  GOresult$xvalue <- switch(
    x_value,
    Count = GOresult$Count,
    GeneRatio = sapply(strsplit(as.character(GOresult$GeneRatio), "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
    pvalue = -log10(GOresult$pvalue + 1e-1000),
    qvalue = -log10(GOresult$qvalue + 1e-1000),
    p.adjust = -log10(as.numeric(GOresult$p.adjust) + 1e-1000),
    stop("Invalid x_value. Supported values: 'Count', 'GeneRatio', 'pvalue', 'qvalue', 'p.adjust'")
  )
  if (x_value %in% c("pvalue", "qvalue", "p.adjust")){
    x_value = paste0("-log10(",x_value,")")
  }
  
  # Take the top_n rows based on xvalue
  GOresult <- GOresult[order(GOresult$xvalue, decreasing = TRUE), ][1:top_n,]
  
  gg_title <- if (ONTOLOGY == "KEGG") {
    paste0("KEGG | Cluster ", tile_name)
  } else {
    paste0("GO ", ONTOLOGY," | Cluster ", tile_name)
  }
  
  p <- ggplot(GOresult, aes(reorder(Description, xvalue), xvalue)) +
    geom_col(fill = col) +   # #C50B05  #1965B0
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.position = 'none',
      axis.text = element_text(color = "black", size = 11),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 1, size = 13)
    ) +
    coord_flip() +
    labs(x = '', y = x_value) +
    ggtitle(gg_title)
  
  return(p)
}



