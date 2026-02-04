#' Plot grouped Gaussian process regression curves with confidence intervals (multi-group supported)
#'
#' @param data Input data.frame. Must contain pseudotime, observation, and grouping columns.
#' @param pseudotime_col Pseudotime column name (string).
#' @param y_col Observation column name (single string).
#' @param group_col Grouping column name (string).
#' @param sigma Gaussian process sigma parameter controlling smoothness (default: 0.5).
#' @param title Plot title (default: auto-generated: "{y_col} along the pseudotime").
#' @param line_colors Color vector. Supports a named vector whose names match group levels.
#' @param ribbon_alpha Confidence interval transparency (default: 0.3).
#' @param y_title Y-axis title (default: "{y_col} Expression level").
#' @return A ggplot2 object.
#' @examples 
#' # Generate example data
#' set.seed(2)
#' df <- data.frame(
#'   pseudotime = sort(runif(100, 0, 10)),
#'   expression = c(sin(seq(0, 10, length=50)) + rnorm(50, sd=0.3),
#'                  cos(seq(0, 10, length=50)) + rnorm(50, sd=0.3)),
#'   group = rep(c("A", "B"), each = 50)
#' )
#' # Auto-generate title when not provided
#' plot_GP_pseudotime_group(df, "pseudotime", "expression", "group")
#' 
#' # Use a named color vector
#' plot_GP_pseudotime_group(df, "pseudotime", "expression", "group",
#'                          line_colors = c(A = "#2c7bb6", B = "#d7191c"))

plot_GP_pseudotime_group <- function(data, pseudotime_col, y_col, group_col, 
                                     sigma = 0.5, 
                                     title = paste(y_col, "along the pseudotime"),
                                     y_title = paste(y_col, "Expression level"),
                                     line_colors = NULL, ribbon_alpha = 0.3) {
  # Validate inputs
  stopifnot(
    is.data.frame(data),
    pseudotime_col %in% colnames(data),
    y_col %in% colnames(data),
    group_col %in% colnames(data),
    is.null(line_colors) || (length(line_colors) >= length(unique(data[[group_col]])))
  )
  
  # Group information
  groups <- unique(data[[group_col]])
  n_groups <- length(groups)
  
  # Create color mapping (supports named vectors)
  if (is.null(line_colors)) {
    line_colors <- scales::hue_pal()(n_groups)
    names(line_colors) <- groups
  } else {
    # Validate that named colors match group levels
    if (!is.null(names(line_colors))) {
      missing_groups <- setdiff(groups, names(line_colors))
      if (length(missing_groups) > 0) {
        stop("The color vector is missing names for groups: ", paste(missing_groups, collapse = ", "))
      }
    } else {
      names(line_colors) <- groups[1:length(line_colors)]
    }
  }
  
  # Build a common prediction grid
  x_range <- range(data[[pseudotime_col]])
  new_data <- data.frame(
    seq(x_range[1], x_range[2], length.out = 200)
  )
  colnames(new_data) <- pseudotime_col
  
  # Fit models per group
  predictions <- lapply(groups, function(grp) {
    subset_data <- data[data[[group_col]] == grp, ]
    
    formula <- as.formula(paste(y_col, "~", pseudotime_col))
    gp_model <- kernlab::gausspr(
      formula,
      data = subset_data,
      kernel = "rbfdot",
      kpar = list(sigma = sigma),
      variance.model = TRUE
    )
    
    pred <- predict(gp_model, new_data, type = "response")
    pred_se <- predict(gp_model, new_data, type = "sdev")
    
    data.frame(
      x = new_data[[pseudotime_col]],
      fit = pred,
      se = pred_se,
      group = factor(grp, levels = groups)
    )
  })
  
  plot_df <- do.call(rbind, predictions)
  colnames(plot_df)[1] <- pseudotime_col
  
  # Build plot object
  p <- ggplot(plot_df, aes(x = !!sym(pseudotime_col), 
                           color = !!sym("group"), 
                           fill = !!sym("group"))) +
    geom_ribbon(
      aes(ymin = .data$fit - 1.96*.data$se, ymax = .data$fit + 1.96*.data$se),
      alpha = ribbon_alpha,
      colour = NA
    ) +
    geom_line(
      aes(y = .data$fit),
      linewidth = 1
    ) +
    scale_color_manual(values = line_colors, name = "Group") +
    scale_fill_manual(values = line_colors, name = "Group") +
    labs(
      x = pseudotime_col,
      y = y_title,
      title = title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      panel.background = element_rect(fill = "white"),
      legend.position = "right",
      axis.ticks = element_line()
    )
  
  return(p)
}
