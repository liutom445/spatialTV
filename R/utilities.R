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
