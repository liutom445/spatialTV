#' Fit Spatial Mixture Model with Learned Edge Weights
#'
#' Performs alternating optimization of mixture proportions (X), dictionary (B),
#' and edge weights (w) using FISTA for X and B, and linear programming for w.
#'
#' @param S Data matrix (N x K) containing cell type likelihood scores
#' @param A Sparse incidence matrix (m x N) for the spatial graph
#' @param B_init Initial dictionary matrix (K x D)
#' @param lambda Regularization parameter. Default: 1.0
#' @param delta Smoothing parameter for soft-L1 norm. Default: 1e-6
#' @param outer_iter Number of outer alternating iterations. Default: 3
#' @param it_X Maximum iterations for X optimization per outer iteration. Default: 1500
#' @param it_B Maximum iterations for B optimization per outer iteration. Default: 500
#' @param eps_proj Floor constraint for simplex projection. Default: 1e-6
#' @param eps_step Step size denominator bound. Default: 1e-3
#' @param safety Safety factor for step size. Default: 0.95
#' @param log_every Logging frequency. Default: 250
#' @param wmin Minimum edge weight. Default: 0.01
#' @param wmax Maximum edge weight. Default: 0.50
#' @param learn_weights If TRUE, optimize edge weights; if FALSE, keep uniform.
#'   Default: TRUE
#' @param X0 Optional initial X matrix
#'
#' @return A list with components:
#' \item{X}{Optimized mixture matrix (N x D)}
#' \item{B}{Optimized dictionary matrix (K x D)}
#' \item{w}{Optimized edge weights (length m)}
#' \item{history}{Data frame with objective values and weight statistics per iteration}
#' \item{final_metrics}{List of final performance metrics}
#'
#' @details
#' The algorithm alternates between three steps:
#' \enumerate{
#'   \item X, B-step: Optimize using FISTA with fixed weights
#'   \item w-step: Solve linear program to optimize edge weights
#' }
#'
#' The weight optimization problem is:
#' \deqn{\min_w w^T c \quad \text{s.t.} \quad w_{min} \le w \le w_{max}, \quad |A|^T w = 1}
#' where \eqn{c_e = \sum_d |(AX)_{ed}|}.
#'
#' @note Requires the CVXR package for weight optimization. If CVXR is not
#'   available, set \code{learn_weights = FALSE}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(example_data)
#' 
#' # Build spatial graph
#' A <- build_spatial_graph(coords)
#' 
#' # Initialize dictionary
#' B_init <- initialize_dictionary(S, D = 3)
#' 
#' # Run alternating optimization
#' result <- fit_spatial_mixture_adaptive(
#'   S, A, B_init,
#'   lambda = 2.0,
#'   outer_iter = 5,
#'   wmin = 0.01,
#'   wmax = 0.50
#' )
#' 
#' # Check convergence
#' plot(result$history$iter, result$history$obj, type = "b")
#' 
#' # Edge weight distribution
#' hist(result$w, breaks = 30)
#' }
fit_spatial_mixture_adaptive <- function(S, A, B_init,
                                        lambda = 1.0, delta = 1e-6,
                                        outer_iter = 3,
                                        it_X = 1500, it_B = 500,
                                        eps_proj = 1e-6, eps_step = 1e-3,
                                        safety = 0.95, log_every = 250,
                                        wmin = 0.01, wmax = 0.50,
                                        learn_weights = TRUE,
                                        X0 = NULL) {
  
  # Input validation
  if (!is.matrix(S)) stop("S must be a matrix")
  if (!inherits(A, "sparseMatrix")) stop("A must be a sparse matrix")
  if (!is.matrix(B_init)) stop("B_init must be a matrix")
  
  N <- nrow(S)
  m <- nrow(A)
  D <- ncol(B_init)
  
  # Check CVXR availability if learning weights
  if (learn_weights && wmin <= 1) {
    if (!requireNamespace("CVXR", quietly = TRUE)) {
      warning("CVXR package not available. Disabling weight learning.")
      learn_weights <- FALSE
    }
  }
  
  # Initialize
  B <- B_init
  w <- rep(1, m)
  X <- X0
  
  history <- data.frame(
    iter = integer(0),
    obj = numeric(0),
    nll = numeric(0),
    tv = numeric(0),
    mean_w = numeric(0),
    median_w = numeric(0),
    sd_w = numeric(0)
  )
  
  # Alternating optimization
  for (iter in seq_len(outer_iter)) {
    cat(sprintf("\n========== Outer Iteration %d / %d ==========\n", 
                iter, outer_iter))
    
    # Step 1: Optimize X and B with fixed w
    res <- fit_spatial_mixture(
      S = S,
      A = A,
      B = B,
      w = w,
      lambda = lambda,
      delta = delta,
      it_X = it_X,
      eps_proj_X = eps_proj,
      eps_step_X = eps_step,
      safety_X = safety,
      log_every_X = log_every,
      it_B = it_B,
      eps_proj_B = eps_proj,
      eps_step_B = eps_step,
      safety_B = safety,
      log_every_B = log_every,
      X0 = X
    )
    
    X <- res$X
    B <- res$B
    
    # Step 2: Optimize w with fixed X (if enabled)
    if (learn_weights && wmin <= 1) {
      w_new <- optimize_edge_weights(X, A, wmin, wmax)
      if (!is.null(w_new)) {
        w <- w_new
      } else {
        warning(sprintf("Weight optimization failed at iteration %d. Keeping previous weights.", iter))
      }
    }
    
    # Record history
    history <- rbind(history, data.frame(
      iter = iter,
      obj = res$obj,
      nll = res$nll,
      tv = res$tv,
      mean_w = mean(w),
      median_w = median(w),
      sd_w = sd(w)
    ))
    
    cat(sprintf("Outer %d done | obj: %.6e | mean_w: %.3f | median_w: %.3f\n",
                iter, res$obj, mean(w), median(w)))
  }
  
  # Compute final metrics
  final_metrics <- compute_metrics(X, B, S, A, w, lambda, delta)
  
  result <- list(
    X = X,
    B = B,
    w = w,
    history = history,
    final_metrics = final_metrics,
    lambda = lambda,
    delta = delta
  )
  
  class(result) <- c("spatialTV_adaptive", "list")
  return(result)
}


#' Optimize Edge Weights via Linear Programming
#'
#' Solves a linear program to find edge weights that minimize the total
#' variation while satisfying normalization constraints.
#'
#' @param X Current mixture matrix (N x D)
#' @param A Sparse incidence matrix (m x N)
#' @param wmin Minimum weight value
#' @param wmax Maximum weight value
#'
#' @return Vector of optimized weights (length m), or NULL if optimization fails
#'
#' @keywords internal
optimize_edge_weights <- function(X, A, wmin, wmax) {
  
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    return(NULL)
  }
  
  m <- nrow(A)
  
  # Compute cost vector: sum of absolute edge differences across all components
  AX <- abs(as.matrix(A %*% X))
  cvec <- rowSums(AX)
  
  # Define optimization problem
  w_var <- CVXR::Variable(m)
  
  constraints <- list(
    w_var >= wmin,
    w_var <= wmax,
    Matrix::t(abs(A)) %*% w_var == 1
  )
  
  objective <- CVXR::Minimize(t(cvec) %*% w_var)
  problem <- CVXR::Problem(objective, constraints)
  
  # Solve
  solution <- tryCatch(
    CVXR::solve(problem, solver = "ECOS"),
    error = function(e) {
      warning("CVXR solver failed: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(solution)) {
    return(NULL)
  }
  
  if (solution$status %in% c("infeasible", "unbounded")) {
    return(NULL)
  }
  
  w_opt <- as.vector(solution$getValue(w_var))
  return(w_opt)
}


#' Compute Performance Metrics
#'
#' Calculates various metrics to assess solution quality.
#'
#' @param X Mixture matrix (N x D)
#' @param B Dictionary matrix (K x D)
#' @param S Data matrix (N x K)
#' @param A Incidence matrix (m x N)
#' @param w Edge weights (length m)
#' @param lambda Regularization parameter
#' @param delta Smoothing parameter
#'
#' @return List of metrics
#'
#' @keywords internal
compute_metrics <- function(X, B, S, A, w, lambda, delta) {
  
  # Compute objective components
  eval_result <- eval_softL1_w_cpp(S, A, X, B, w, lambda, delta)
  
  # Neighbor agreement (fraction of edges with same dominant component)
  cls <- apply(X, 1, which.max)
  edges <- Matrix::summary(abs(A))[, 1:2, drop = FALSE]
  neigh_agree <- mean(cls[edges[, 1]] == cls[edges[, 2]])
  
  # Entropy (normalized, 0 = deterministic, 1 = uniform)
  row_entropy <- function(p) {
    p <- pmax(p, 1e-12)
    -sum(p * log(p)) / log(length(p))
  }
  entropy <- mean(apply(X, 1, row_entropy))
  
  # Sharpness (mean of maximum component per row)
  sharpness <- mean(apply(X, 1, max))
  
  list(
    objective = eval_result$obj,
    nll = eval_result$nll,
    tv = eval_result$tv,
    neighbor_agreement = neigh_agree,
    entropy = entropy,
    sharpness = sharpness
  )
}
