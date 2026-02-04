plot_GP_pseudotime <- function(
  plot_data,
  pseudotime_col = "Pseudotime",
  sigma = 0.8,
  y_col,
  line_colors = NULL,
  ribbon_alpha = 0.18,
  grid_n = 200
) {
  if (!pseudotime_col %in% colnames(plot_data)) {
    stop("Missing column: ", pseudotime_col)
  }
  missing_y <- setdiff(y_col, colnames(plot_data))
  if (length(missing_y) > 0) {
    stop("Missing y columns: ", paste(missing_y, collapse = ", "))
  }

  x <- as.numeric(plot_data[[pseudotime_col]])
  x_grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = grid_n)

  gp_df_list <- vector("list", length(y_col))
  names(gp_df_list) <- y_col

  for (nm in y_col) {
    y <- as.numeric(plot_data[[nm]])
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 50) next

    gp <- kernlab::gausspr(as.matrix(x[keep]), y[keep], kernel = "rbfdot", kpar = list(sigma = sigma))
    mu <- as.numeric(predict(gp, as.matrix(x_grid)))
    vv <- as.numeric(predict(gp, as.matrix(x_grid), type = "variances"))
    se <- sqrt(pmax(vv, 0))

    gp_df_list[[nm]] <- data.frame(
      Pseudotime = x_grid,
      Feature = nm,
      mean = mu,
      lower = mu - 1.96 * se,
      upper = mu + 1.96 * se
    )
  }

  gp_df <- dplyr::bind_rows(gp_df_list)
  if (nrow(gp_df) == 0) stop("No valid GP fits were generated (insufficient finite points).")

  if (is.null(line_colors)) {
    line_colors <- scales::hue_pal()(length(y_col))
  }
  if (is.null(names(line_colors)) && length(line_colors) == length(y_col)) {
    line_colors <- setNames(line_colors, y_col)
  }

  ggplot2::ggplot(gp_df, ggplot2::aes(x = Pseudotime, y = mean, color = Feature, fill = Feature)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = ribbon_alpha, color = NA) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::scale_color_manual(values = line_colors) +
    ggplot2::scale_fill_manual(values = line_colors) +
    ggplot2::labs(x = "Pseudotime", y = "Score (GP mean ± 95% CI)") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.title = ggplot2::element_blank())
}

