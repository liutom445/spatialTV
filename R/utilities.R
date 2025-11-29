#' Evaluate Objective Function
#'
#' Computes the objective function value and its components without optimization.
#'
#' @param S Data matrix (N x K)
#' @param A Sparse incidence matrix (m x N)
#' @param X Mixture matrix (N x D)
#' @param B Dictionary matrix (K x D)
#' @param w Edge weights (length m)
#' @param lambda Regularization parameter
#' @param delta Smoothing parameter. Default: 1e-6
#'
#' @return A list with components:
#' \item{obj}{Total objective value}
#' \item{nll}{Negative log-likelihood}
#' \item{tv}{Total variation penalty}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- evaluate_objective(S, A, X, B, w, lambda = 2.0)
#' print(obj$obj)
#' }
evaluate_objective <- function(S, A, X, B, w, lambda, delta = 1e-6) {
  
  if (!is.matrix(S) || !is.matrix(X) || !is.matrix(B)) {
    stop("S, X, and B must be matrices")
  }
  
  if (!inherits(A, "sparseMatrix")) {
    stop("A must be a sparse matrix")
  }
  
  result <- eval_softL1_w_cpp(S, A, X, B, w, lambda, delta)
  return(result)
}


#' Convert Component Mixtures to Cell Type Proportions
#'
#' Transforms the low-rank component representation (X) back to cell type
#' proportions using the dictionary (B).
#'
#' @param X Mixture matrix (N x D)
#' @param B Dictionary matrix (K x D)
#' @param normalize If TRUE, normalize rows to sum to 1. Default: TRUE
#'
#' @return Matrix (N x K) of cell type proportions
#'
#' @details
#' The transformation is: \eqn{P = X B^T} followed by row normalization.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After optimization
#' result <- fit_spatial_mixture(S, A, B, lambda = 2)
#' 
#' # Convert to cell type proportions
#' cell_type_props <- components_to_celltypes(result$X, result$B)
#' 
#' # Check dimensions
#' dim(cell_type_props)  # N x K
#' rowSums(cell_type_props)  # All 1s
#' }
components_to_celltypes <- function(X, B, normalize = TRUE) {
  
  if (ncol(X) != ncol(B)) {
    stop("Number of columns in X must equal number of columns in B")
  }
  
  # Compute X * B^T
  P <- X %*% t(B)
  
  if (normalize) {
    # Row normalize
    row_sums <- rowSums(P)
    row_sums[row_sums < 1e-12] <- 1e-12  # Prevent division by zero
    P <- P / row_sums
  }
  
  return(P)
}


#' Get Dominant Cell Type per Spot
#'
#' Extracts the dominant (argmax) cell type for each spatial spot.
#'
#' @param X Mixture matrix (N x D) or cell type matrix (N x K)
#' @param B Optional dictionary matrix (K x D). If provided, converts X to
#'   cell types first.
#' @param cell_type_names Optional character vector of cell type names
#'
#' @return Factor vector of length N with dominant cell type per spot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get dominant cell types
#' dominant <- get_dominant_celltype(result$X, result$B, 
#'                                   cell_type_names = colnames(S))
#' table(dominant)
#' }
get_dominant_celltype <- function(X, B = NULL, cell_type_names = NULL) {
  
  if (!is.null(B)) {
    # Convert to cell type proportions
    X <- components_to_celltypes(X, B)
  }
  
  # Get argmax
  dominant_idx <- apply(X, 1, which.max)
  
  if (!is.null(cell_type_names)) {
    if (length(cell_type_names) != ncol(X)) {
      stop("Length of cell_type_names must equal ncol(X)")
    }
    dominant <- factor(cell_type_names[dominant_idx], levels = cell_type_names)
  } else {
    dominant <- factor(dominant_idx)
  }
  
  return(dominant)
}


#' Extract Top K Cell Types per Spot
#'
#' For each spot, returns the top K cell types and their proportions.
#'
#' @param X Mixture matrix (N x D) or cell type matrix (N x K)
#' @param B Optional dictionary matrix (K x D)
#' @param K Number of top cell types to return. Default: 3
#' @param cell_type_names Optional character vector of cell type names
#'
#' @return A list with:
#' \item{types}{Matrix (N x K) of cell type indices/names}
#' \item{proportions}{Matrix (N x K) of proportions}
#'
#' @export
get_top_celltypes <- function(X, B = NULL, K = 3, cell_type_names = NULL) {
  
  if (!is.null(B)) {
    X <- components_to_celltypes(X, B)
  }
  
  N <- nrow(X)
  n_types <- ncol(X)
  
  if (K > n_types) {
    warning("K larger than number of cell types. Using K = ", n_types)
    K <- n_types
  }
  
  top_types <- matrix(0, nrow = N, ncol = K)
  top_props <- matrix(0, nrow = N, ncol = K)
  
  for (i in seq_len(N)) {
    ord <- order(X[i, ], decreasing = TRUE)[seq_len(K)]
    top_types[i, ] <- ord
    top_props[i, ] <- X[i, ord]
  }
  
  if (!is.null(cell_type_names)) {
    type_names <- matrix(cell_type_names[top_types], nrow = N, ncol = K)
    list(types = type_names, proportions = top_props)
  } else {
    list(types = top_types, proportions = top_props)
  }
}


#' Print Method for spatialTV Objects
#'
#' @param x A spatialTV object
#' @param ... Additional arguments (unused)
#'
#' @export
print.spatialTV <- function(x, ...) {
  cat("spatialTV model fit\n")
  cat("===================\n\n")
  cat("Dimensions:\n")
  cat("  Spots (N):      ", nrow(x$X), "\n")
  cat("  Cell types (K): ", nrow(x$B), "\n")
  cat("  Components (D): ", ncol(x$X), "\n\n")
  cat("Objective:        ", sprintf("%.6e", x$obj), "\n")
  cat("  NLL:            ", sprintf("%.6e", x$nll), "\n")
  cat("  TV penalty:     ", sprintf("%.6e", x$tv), "\n\n")
  cat("Step sizes:\n")
  cat("  X-step:         ", sprintf("%.2e", x$stepX), "\n")
  cat("  B-step:         ", sprintf("%.2e", x$stepB), "\n")
  invisible(x)
}


#' Print Method for spatialTV_adaptive Objects
#'
#' @param x A spatialTV_adaptive object
#' @param ... Additional arguments (unused)
#'
#' @export
print.spatialTV_adaptive <- function(x, ...) {
  cat("spatialTV model fit (adaptive weights)\n")
  cat("======================================\n\n")
  cat("Dimensions:\n")
  cat("  Spots (N):      ", nrow(x$X), "\n")
  cat("  Cell types (K): ", nrow(x$B), "\n")
  cat("  Components (D): ", ncol(x$X), "\n")
  cat("  Edges (m):      ", length(x$w), "\n\n")
  cat("Convergence:\n")
  cat("  Outer iters:    ", nrow(x$history), "\n")
  cat("  Final obj:      ", sprintf("%.6e", tail(x$history$obj, 1)), "\n\n")
  cat("Edge weights:\n")
  cat("  Mean:           ", sprintf("%.3f", mean(x$w)), "\n")
  cat("  Median:         ", sprintf("%.3f", median(x$w)), "\n")
  cat("  Range:          [", sprintf("%.3f", min(x$w)), ", ",
      sprintf("%.3f", max(x$w)), "]\n\n")
  cat("Final metrics:\n")
  for (name in names(x$final_metrics)) {
    cat("  ", name, ": ", sprintf("%.4f", x$final_metrics[[name]]), "\n", sep = "")
  }
  invisible(x)
}


#' Plot Component Proportions Across Space
#'
#' Creates high-quality spatial visualizations of mixture component proportions.
#'
#' @param X Mixture matrix (N x D)
#' @param coords Spatial coordinates (N x 2 matrix)
#' @param point_size Point size for plotting. Default: 2.5
#' @param color_palette Color palette to use. Options: "viridis" (default),
#'   "magma", "plasma", "inferno", "cividis", "rocket", "mako", "turbo"
#' @param ncol Number of columns in the grid layout. Default: 2
#' @param alpha Transparency level (0-1). Default: 0.9
#' @param add_contours Add contour lines. Default: FALSE
#' @param contour_bins Number of contour bins. Default: 10
#' @param show_colorbar Show color bar legend. Default: TRUE
#' @param title_prefix Prefix for subplot titles. Default: "Component"
#' @param high_res Use high-resolution settings. Default: TRUE
#'
#' @return A ggplot object (if single component) or grid of plots (if multiple)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After fitting model
#' result <- fit_spatial_mixture(S, A, B, lambda = 2)
#'
#' # Basic plot
#' plot_components_spatial(result$X, coords)
#'
#' # Enhanced plot with custom settings
#' plot_components_spatial(
#'   result$X,
#'   coords,
#'   point_size = 3,
#'   color_palette = "magma",
#'   ncol = 3,
#'   add_contours = TRUE
#' )
#'
#' # Save high-quality PDF
#' p <- plot_components_spatial(result$X, coords)
#' ggsave("components.pdf", p, width = 10, height = 8, dpi = 300)
#' }
plot_components_spatial <- function(X,
                                   coords,
                                   point_size = 2.5,
                                   color_palette = "viridis",
                                   ncol = 2,
                                   alpha = 0.9,
                                   add_contours = FALSE,
                                   contour_bins = 10,
                                   show_colorbar = TRUE,
                                   title_prefix = "Component",
                                   high_res = TRUE) {

  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  # Validate inputs
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("coords must be an N x 2 matrix")
  }
  if (nrow(X) != nrow(coords)) {
    stop("Number of rows in X must match number of rows in coords")
  }

  D <- ncol(X)
  N <- nrow(X)

  # Prepare data frame
  plot_df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    X
  )

  # Set color scale function
  color_scale <- switch(color_palette,
    "viridis" = ggplot2::scale_color_viridis_c,
    "magma" = function(...) ggplot2::scale_color_viridis_c(option = "magma", ...),
    "plasma" = function(...) ggplot2::scale_color_viridis_c(option = "plasma", ...),
    "inferno" = function(...) ggplot2::scale_color_viridis_c(option = "inferno", ...),
    "cividis" = function(...) ggplot2::scale_color_viridis_c(option = "cividis", ...),
    "rocket" = function(...) ggplot2::scale_color_viridis_c(option = "rocket", ...),
    "mako" = function(...) ggplot2::scale_color_viridis_c(option = "mako", ...),
    "turbo" = function(...) ggplot2::scale_color_viridis_c(option = "turbo", ...),
    ggplot2::scale_color_viridis_c  # default
  )

  # Create individual plots
  plots <- list()
  for (d in 1:D) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = X[, d])) +
      ggplot2::geom_point(size = point_size, alpha = alpha) +
      color_scale(
        limits = c(0, 1),
        name = "Proportion",
        guide = if (show_colorbar) ggplot2::guide_colorbar(
          barwidth = 1,
          barheight = 8,
          title.position = "top",
          title.hjust = 0.5
        ) else "none"
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal(base_size = if (high_res) 14 else 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = if (high_res) 16 else 14),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
        axis.title = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      ) +
      ggplot2::labs(
        title = paste0(title_prefix, " ", d),
        x = "X coordinate",
        y = "Y coordinate"
      )

    # Add contours if requested
    if (add_contours) {
      p <- p + ggplot2::geom_contour(
        ggplot2::aes(z = X[, d]),
        bins = contour_bins,
        color = "white",
        alpha = 0.3,
        linewidth = 0.3
      )
    }

    plots[[d]] <- p
  }

  # Return single plot or grid
  if (D == 1) {
    return(plots[[1]])
  } else {
    # Check if gridExtra is available
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      warning("Package 'gridExtra' not found. Returning list of plots instead of grid.\n",
              "Install with: install.packages('gridExtra')")
      return(plots)
    }

    return(gridExtra::grid.arrange(grobs = plots, ncol = ncol))
  }
}
