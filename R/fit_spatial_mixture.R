#' Optimize Spatial Mixture Model with FISTA
#'
#' Performs alternating optimization of cell type mixtures (X) and dictionary (B)
#' using Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) with optional
#' weighted Total Variation regularization.
#'
#' @param S Data matrix (N x K) containing cell type likelihood scores from
#'   deconvolution methods like RCTD. Each row corresponds to a spatial spot,
#'   each column to a cell type.
#' @param A Sparse incidence matrix (m x N) for the spatial graph, typically
#'   from a Delaunay triangulation. Each row represents an edge with +1/-1
#'   entries at the two incident spots.
#' @param B Initial dictionary matrix (K x D) where K is the number of cell types
#'   and D is the number of mixture components. Columns must sum to 1.
#' @param w Edge weights vector (length m) for weighted TV. Default is uniform
#'   weights of 1.
#' @param lambda Regularization parameter controlling spatial smoothness.
#'   Larger values enforce stronger spatial coherence. Default: 1.0.
#' @param delta Smoothing parameter for soft-L1 norm. Default: 1e-6.
#' @param it_X Maximum iterations for X optimization. Default: 1500.
#' @param eps_proj_X Floor constraint for X simplex projection. Default: 1e-6.
#' @param eps_step_X Step size denominator bound for X. Default: 1e-3.
#' @param normalize_for_bound Use normalized bound for Lipschitz constant.
#'   Default: TRUE.
#' @param safety_X Safety factor for X step size (0 < safety < 1). Default: 0.95.
#' @param log_every_X Logging frequency for X optimization (0 = no logging).
#'   Default: 250.
#' @param it_B Maximum iterations for B optimization. Default: 500.
#' @param eps_proj_B Floor constraint for B simplex projection. Default: 1e-6.
#' @param eps_step_B Step size denominator bound for B. Default: 1e-3.
#' @param safety_B Safety factor for B step size. Default: 0.95.
#' @param log_every_B Logging frequency for B optimization. Default: 250.
#' @param update_B Whether to update dictionary B. Set to FALSE to keep B fixed.
#'   Default: TRUE.
#' @param X0 Optional initial X matrix (N x D). If NULL, initialized uniformly.
#'
#' @return A list with components:
#' \item{X}{Optimized mixture matrix (N x D) with rows on simplex}
#' \item{B}{Optimized dictionary matrix (K x D) with columns on simplex}
#' \item{obj}{Final objective value}
#' \item{nll}{Negative log-likelihood component}
#' \item{tv}{Total variation component}
#' \item{stepX}{Step size used for X optimization}
#' \item{stepB}{Step size used for B optimization}
#'
#' @details
#' The objective function minimized is:
#' \deqn{-\sum_{i=1}^N \log z_i + \lambda \sum_{e=1}^m w_e \sum_{d=1}^D \sqrt{(AX)_{ed}^2 + \delta}}
#' where \eqn{z_i = \sum_k S_{ik} (XB^T)_{ik}}.
#'
#' The algorithm alternates between:
#' \enumerate{
#'   \item X-step: Optimize mixture proportions with fixed dictionary
#'   \item B-step: Optimize dictionary with fixed mixtures
#' }
#'
#' @references
#' Beck, A. and Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences,
#' 2(1), 183-202.
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(example_data)
#' 
#' # Create spatial graph
#' A <- build_delaunay_graph(coords)
#' 
#' # Initialize dictionary
#' B <- initialize_dictionary(S, D = 3)
#' 
#' # Run optimization
#' result <- fit_spatial_mixture(S, A, B, lambda = 2.0)
#' 
#' # Extract results
#' X <- result$X  # Mixture proportions
#' B <- result$B  # Cell type dictionary
#' }
#'
#' @export
fit_spatial_mixture <- function(S, A, B, w = NULL,
                               lambda = 1.0, delta = 1e-6,
                               it_X = 1500, eps_proj_X = 1e-6, eps_step_X = 1e-3,
                               normalize_for_bound = TRUE,
                               safety_X = 0.95, log_every_X = 250,
                               it_B = 500, eps_proj_B = 1e-6, eps_step_B = 1e-3,
                               safety_B = 0.95, log_every_B = 250,
                               update_B = TRUE, X0 = NULL) {
  
  # Input validation
  if (!is.matrix(S)) {
    stop("S must be a matrix")
  }
  if (!inherits(A, "dgCMatrix") && !inherits(A, "sparseMatrix")) {
    stop("A must be a sparse matrix (Matrix package)")
  }
  if (!is.matrix(B)) {
    stop("B must be a matrix")
  }
  
  # Dimension checks
  N <- nrow(S)
  K <- ncol(S)
  m <- nrow(A)
  
  if (ncol(A) != N) {
    stop("Number of columns in A must equal number of rows in S")
  }
  if (nrow(B) != K) {
    stop("Number of rows in B must equal number of columns in S")
  }
  
  # Set default weights if not provided
  if (is.null(w)) {
    w <- rep(1, m)
  }
  
  if (length(w) != m) {
    stop("Length of w must equal number of rows in A")
  }
  
  # Convert to proper sparse matrix format if needed
  if (!inherits(A, "dgCMatrix")) {
    A <- as(A, "dgCMatrix")
  }
  
  # Call C++ implementation
  result <- xb_fista_wfixed(
    S_mat = S,
    A = A,
    B = B,
    w = w,
    lambda = lambda,
    delta = delta,
    it_X = it_X,
    eps_proj_X = eps_proj_X,
    eps_step_X = eps_step_X,
    normalize_for_bound = normalize_for_bound,
    safety_X = safety_X,
    log_every_X = log_every_X,
    it_B = it_B,
    eps_proj_B = eps_proj_B,
    eps_step_B = eps_step_B,
    safety_B = safety_B,
    log_every_B = log_every_B,
    update_B = update_B,
    X0_opt = X0
  )
  
  class(result) <- c("spatialTV", "list")
  return(result)
}
