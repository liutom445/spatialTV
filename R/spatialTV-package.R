#' spatialTV: Spatial Total Variation Regularization for Transcriptomics
#'
#' Implements spatially-aware cell type deconvolution for spatial transcriptomics
#' data using Total Variation (TV) regularization. The package uses the Fast
#' Iterative Shrinkage-Thresholding Algorithm (FISTA) for efficient optimization.
#'
#' @description
#' The spatialTV package provides tools for:
#' \itemize{
#'   \item Spatial cell type deconvolution with TV regularization
#'   \item Low-rank dictionary learning for cell type mixtures
#'   \item Adaptive edge weight learning via linear programming
#'   \item Efficient C++ implementation via Rcpp/RcppArmadillo
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{fit_spatial_mixture}}: Basic FISTA optimization with fixed weights
#'   \item \code{\link{fit_spatial_mixture_adaptive}}: Alternating optimization with learned weights
#'   \item \code{\link{build_spatial_graph}}: Construct Delaunay triangulation
#'   \item \code{\link{initialize_dictionary}}: Initialize dictionary matrix
#'   \item \code{\link{evaluate_objective}}: Evaluate objective function
#' }
#'
#' @section Mathematical Formulation:
#' The optimization problem is:
#' \deqn{\min_{X,B} -\sum_{i=1}^N \log z_i + \lambda \sum_{e=1}^m w_e \sum_{d=1}^D \sqrt{(AX)_{ed}^2 + \delta}}
#' subject to simplex constraints on X and B, where:
#' \itemize{
#'   \item X: per-spot mixture proportions (N x D)
#'   \item B: cell type dictionary (K x D)
#'   \item A: spatial graph incidence matrix (m x N)
#'   \item w: edge weights
#'   \item \eqn{z_i = \sum_k S_{ik}(XB^T)_{ik}}: per-spot evidence
#' }
#'
#' @docType package
#' @name spatialTV-package
#' @aliases spatialTV
#' @useDynLib spatialTV, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix Matrix sparseMatrix
#' @importFrom methods as is
#' @importFrom stats kmeans rgamma
#'
#' @references
#' Beck, A. and Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
