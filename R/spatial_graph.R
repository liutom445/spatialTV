#' Build Delaunay Triangulation Graph
#'
#' Constructs a sparse incidence matrix representing a Delaunay triangulation
#' of spatial coordinates. Each edge in the triangulation becomes a row in the
#' incidence matrix with +1 and -1 entries at the two incident nodes.
#'
#' @param coords Matrix (N x 2) of spatial coordinates (x, y) for each spot
#' @param method Method for triangulation. Currently only "delaunay" is supported.
#'
#' @return A sparse matrix (dgCMatrix) of size m x N where m is the number of
#'   edges. Each row has exactly two non-zero entries: +1 and -1.
#'
#' @details
#' The Delaunay triangulation maximizes the minimum angle of all triangles,
#' which tends to produce a well-conditioned graph for spatial regularization.
#'
#' @importFrom Matrix Matrix sparseMatrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate random spatial coordinates
#' coords <- cbind(x = runif(100, 0, 10), y = runif(100, 0, 10))
#' 
#' # Build graph
#' A <- build_spatial_graph(coords)
#' 
#' # Check properties
#' dim(A)  # m x 100 where m is number of edges
#' sum(A != 0) / length(A)  # Sparsity
#' }
build_spatial_graph <- function(coords, method = "delaunay") {
  
  if (!is.matrix(coords) && !is.data.frame(coords)) {
    stop("coords must be a matrix or data.frame")
  }
  
  coords <- as.matrix(coords)
  
  if (ncol(coords) != 2) {
    stop("coords must have exactly 2 columns (x, y)")
  }
  
  N <- nrow(coords)
  
  if (N < 3) {
    stop("Need at least 3 points for triangulation")
  }
  
  if (method == "delaunay") {
    if (!requireNamespace("interp", quietly = TRUE)) {
      stop("Package 'interp' is required for Delaunay triangulation. Please install it.")
    }

    # Build triangulation (tri.mesh expects separate x and y vectors)
    tri <- interp::tri.mesh(coords[, 1], coords[, 2])
    m <- nrow(tri$arcs)
    
    # Create sparse incidence matrix
    A <- Matrix::Matrix(0, nrow = m, ncol = N, sparse = TRUE)
    A[cbind(seq_len(m), tri$arcs[, "from"])] <- 1
    A[cbind(seq_len(m), tri$arcs[, "to"])] <- -1
    
    return(A)
    
  } else {
    stop("Unknown method. Currently only 'delaunay' is supported.")
  }
}


#' Initialize Dictionary Matrix
#'
#' Creates an initial dictionary matrix B using k-means clustering on the
#' cell type score matrix S. This provides a data-driven initialization
#' for the optimization.
#'
#' @param S Data matrix (N x K) of cell type scores
#' @param D Number of mixture components (dictionary atoms)
#' @param method Initialization method. Options: "kmeans", "random", "uniform"
#' @param eps Floor value for simplex constraint. Default: 1e-6
#' @param seed Random seed for reproducibility. Default: NULL
#'
#' @return Dictionary matrix (K x D) with columns on the simplex
#'
#' @details
#' For \code{method = "kmeans"}: Performs k-means clustering on S and uses
#' cluster centers as initial dictionary atoms after normalizing columns.
#'
#' For \code{method = "random"}: Samples from Dirichlet distribution.
#'
#' For \code{method = "uniform"}: Uses uniform distribution on simplex.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' S <- matrix(rexp(100 * 10), nrow = 100, ncol = 10)
#' 
#' # Initialize with k-means
#' B <- initialize_dictionary(S, D = 3, method = "kmeans")
#' 
#' # Check constraints
#' colSums(B)  # Should be all 1s
#' min(B)      # Should be >= eps
#' }
initialize_dictionary <- function(S, D, method = "kmeans",
                                 eps = 1e-6, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  K <- ncol(S)
  
  if (D < 1 || D > K) {
    stop("D must be between 1 and ncol(S)")
  }
  
  if (method == "kmeans") {
    # Use k-means on S to get initial centers
    km <- stats::kmeans(S, centers = D, nstart = 10)
    centers <- km$centers  # D x K
    
    # Transpose and project to column simplex
    B <- t(centers)  # K x D
    B <- pmax(B, eps)
    B <- sweep(B, 2, colSums(B), "/")
    B <- project_columns_to_simplex(B, eps = eps)
    
  } else if (method == "random") {
    # Random initialization from Dirichlet
    B <- matrix(0, nrow = K, ncol = D)
    for (d in seq_len(D)) {
      alpha <- rep(1, K)
      B[, d] <- stats::rgamma(K, shape = alpha)
      B[, d] <- B[, d] / sum(B[, d])
    }
    B <- pmax(B, eps)
    B <- sweep(B, 2, colSums(B), "/")
    
  } else if (method == "uniform") {
    # Uniform initialization
    B <- matrix(1 / K, nrow = K, ncol = D)
    
  } else {
    stop("Unknown method. Options: 'kmeans', 'random', 'uniform'")
  }
  
  return(B)
}


#' Project Columns to Simplex with Floor Constraint
#'
#' Projects each column of a matrix onto the simplex with a minimum value
#' constraint: sum(col) = 1 and all elements >= eps.
#'
#' @param B Matrix to project
#' @param eps Minimum value for each element. Default: 1e-6
#'
#' @return Matrix with same dimensions as B, columns projected to simplex
#'
#' @keywords internal
project_columns_to_simplex <- function(B, eps = 1e-6) {
  K <- nrow(B)
  D <- ncol(B)
  target_sum <- 1 - K * eps
  
  for (d in seq_len(D)) {
    v <- B[, d]
    # Shift by floor
    v_shifted <- pmax(v - eps, 0)
    # Project to simplex
    v_proj <- project_to_simplex(v_shifted, target_sum)
    # Add floor back
    B[, d] <- v_proj + eps
  }
  
  return(B)
}


#' Project Vector to Simplex
#'
#' Projects a vector onto the probability simplex using efficient sorting.
#'
#' @param v Vector to project
#' @param target_sum Target sum (default: 1)
#'
#' @return Projected vector
#'
#' @keywords internal
project_to_simplex <- function(v, target_sum = 1) {
  n <- length(v)
  v_sorted <- sort(v, decreasing = TRUE)
  cssv <- cumsum(v_sorted) - target_sum
  
  rho <- 0
  for (j in seq_along(v_sorted)) {
    if (v_sorted[j] > cssv[j] / j) {
      rho <- j
    }
  }
  
  theta <- if (rho > 0) cssv[rho] / rho else 0
  pmax(v - theta, 0)
}
