create_group_column <- function(data, Group_top, Group_box, celltype_order = NULL, group_order = NULL) {
  
  if (!is.null(celltype_order)) {
    data[[Group_top]] <- factor(data[[Group_top]], levels = celltype_order)
  } else if (is.factor(data[[Group_top]])) {
    celltype_order <- levels(data[[Group_top]])
  } else {
    celltype_order <- unique(data[[Group_top]])
  }
  
  if (!is.null(group_order)) {
    data[[Group_box]] <- factor(data[[Group_box]], levels = group_order)
  } else if (is.factor(data[[Group_box]])) {
    group_order <- levels(data[[Group_box]])
  } else {
    group_order <- unique(data[[Group_box]])
  }
  
  data$ct_group <- paste0(data[[Group_top]], "_lianjiefuhao_", data[[Group_box]])
  level_order <- paste0(rep(celltype_order, each = length(group_order)), "_lianjiefuhao_", rep(group_order, length(celltype_order)))
  
  missing_combinations <- setdiff(level_order, data$ct_group)
  if (length(missing_combinations) > 0) {
    missing_rows <- data.frame(
      orig.ident = "miss_value",                
      values = 0,
      Group = sub(".*_lianjiefuhao_", "", missing_combinations),       
      celltype = sub("_lianjiefuhao_.*", "", missing_combinations),
      ct_group = missing_combinations
    )
    data <- rbind(data, missing_rows)
  }
  
  data$ct_group <- paste0(data[[Group_top]], "_", data[[Group_box]])
  level_order <- paste0(rep(celltype_order, each = length(group_order)), "_", rep(group_order, length(celltype_order)))
  
  data[[Group_top]] <- factor(data[[Group_top]], levels = celltype_order)
  data[[Group_box]] <- factor(data[[Group_box]], levels = group_order)
  data$ct_group <- factor(data$ct_group, levels = level_order)
  
  return(data)
}

create_comparison_list <- function(celltypes, groups) {
  comparisons <- list()
  
  for (celltype in celltypes) {
    celltype_groups <- paste0(celltype, "_", groups)
    group_combinations <- combn(celltype_groups, 2, simplify = FALSE)
    comparisons <- append(comparisons, group_combinations)
  }
  
  return(comparisons)
}

generate_median_line_data <- function(data, ct_group_col_name, value_col) {
  n_groups <- length(unique(data[[ct_group_col_name]]))  # 计算分组数量
  x_value <- c(0.7, 1.3)
  for (i in 1:(n_groups - 1)) {
    x_value <- append(x_value, x_value[((i - 1) * 2 + 1):(i * 2)] + 1)
  }
  
  medians <- data %>%
    group_by(!!sym(ct_group_col_name)) %>%
    summarise(median_value = median(!!sym(value_col)))
  
  med_value <- rep(medians$median_value, each = 2)
  
  tmp_data <- data.frame(
    x_value = x_value,
    med_value = med_value,
    ct_group = rep(unique(data[[ct_group_col_name]]), each = 2)
  )
  
  return(tmp_data)
}



library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggsci)

plot_levels_box <- function(data,
                            ct_group_col_name = "ct_group", 
                            value_col = "values", 
                            Group_top = "celltype", 
                            Group_col = "Group",
                            Group_box_colors = NULL,
                            comparisons, 
                            title_name,
                            sig_lable = TRUE,
                            test_method = "wilcox.test") {
  
  tmp_data <- generate_median_line_data(data, ct_group_col_name, value_col)
  
  y_min <- pmin(min(data[[value_col]]), 0) 
  values_max <- max(data[[value_col]])
  y_max <- values_max * 1.2 
  Group_levels <- levels(data[[Group_col]])
  num_Group = length(Group_levels)
  if (is.null(Group_box_colors)){
    Group_box_colors = pal_npg("nrc", alpha = 0.75)(num_Group) 
  }
    
  num_ct_group = length(unique(data[[ct_group_col_name]]))
  rect_loc_ymin <- y_max * 0.9  
  rect_loc_ymax <- y_max  
  rect_loc_xmax <- num_ct_group + 1
  
  num_celltype = length(unique(data[[Group_top]]))
  x_labels <- rep(Group_levels, num_celltype)
  
  celltype_labels <- levels(data[[Group_top]])
  celltype_labels_loc_y <- y_max * 0.94
  celltype_labels_loc_x <- seq((length(Group_levels) + 1) / 2, rect_loc_xmax, length(Group_levels))

  line_pos = seq(length(Group_levels) + 0.5,  by = length(Group_levels), length.out = num_celltype-1)
  
  psp = rev(seq(from = 0.88, by = -0.05, length.out = choose(num_Group, 2)))
  sig_text_loc <- rep(rect_loc_ymin * psp, num_celltype)
  
  p <- ggplot(data) +
    stat_boxplot(aes(x = .data[[ct_group_col_name]], y = .data[[value_col]]), geom = "errorbar", 
                 position = position_dodge(width = 0.2), 
                 width = 0.1) +
    geom_boxplot(aes(x = .data[[ct_group_col_name]], y = .data[[value_col]], fill = .data[[Group_col]], color = .data[[Group_col]]), 
                 width = 0.5, outlier.shape = NA) +
    
    geom_rect(aes(xmin = 0, xmax = rect_loc_xmax, ymin = rect_loc_ymin, ymax = rect_loc_ymax), 
              fill = "#eaeae0") +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    
    scale_fill_manual(values = Group_box_colors) +
    scale_color_manual(values =Group_box_colors) +
    theme_bw() +
    
    geom_line(data = tmp_data, aes(x_value, med_value, group = ct_group), color = "#ffffff") +
    
    theme(panel.grid = element_blank(), 
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black")) +
    
    ggtitle(title_name) +
    
    annotate("text", x = celltype_labels_loc_x, y = celltype_labels_loc_y, 
             label = celltype_labels) 
  
  if (sig_lable){
    p= p +
      geom_signif(test = test_method,
                  aes(x = .data[[ct_group_col_name]], y = .data[[value_col]]),
                  comparisons = comparisons,
                  map_signif_level = TRUE,
                  vjust = 0, 
                  tip_length = rep(0, length(comparisons)), 
                  textsize = 2.7,
                  y_position = sig_text_loc)
  }
  
  return(p)
}

Groups_box_plot <- function(data, Group_top, Group_top_order=NULL, Group_box, Group_box_colors=NULL, Group_box_order=NULL, title_name=NULL, sig_lable = TRUE, test_method = "wilcox.test"){
  
  if (!is.null(Group_top_order)){
    data[[Group_top]] <- factor(data[[Group_top]], levels = Group_top_order)
  }
  if (!is.null(Group_top_order)){
    data[[Group_box]] <- factor(data[[Group_box]], levels = Group_box_order)
  }
  
  data <- create_group_column(data, Group_top = Group_top, Group_box = Group_box)

  celltypes = levels(data[[Group_top]])
  groups = levels(data[[Group_box]])
  my_comparisons <- create_comparison_list(celltypes, groups)
  
  p <- plot_levels_box(
    data = data,
    comparisons = my_comparisons,
    Group_box_colors = Group_box_colors,
    title_name = title_name,
    sig_lable = sig_lable
  )
  return(p)
}

  



