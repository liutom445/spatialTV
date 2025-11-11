test_that("initialize_dictionary works with different methods", {
  set.seed(123)
  S <- matrix(rexp(100 * 10), nrow = 100, ncol = 10)
  
  # K-means initialization
  B_km <- initialize_dictionary(S, D = 3, method = "kmeans", seed = 123)
  expect_equal(dim(B_km), c(10, 3))
  expect_true(all(colSums(B_km) > 0.99 & colSums(B_km) < 1.01))
  expect_true(all(B_km >= 1e-7))
  
  # Random initialization
  B_rand <- initialize_dictionary(S, D = 3, method = "random", seed = 123)
  expect_equal(dim(B_rand), c(10, 3))
  expect_true(all(colSums(B_rand) > 0.99 & colSums(B_rand) < 1.01))
  
  # Uniform initialization
  B_unif <- initialize_dictionary(S, D = 3, method = "uniform")
  expect_equal(dim(B_unif), c(10, 3))
  expect_true(all(abs(B_unif - 0.1) < 1e-6))
})

test_that("initialize_dictionary validates inputs", {
  S <- matrix(rexp(100 * 10), nrow = 100, ncol = 10)
  
  expect_error(initialize_dictionary(S, D = 0), "D must be between")
  expect_error(initialize_dictionary(S, D = 15), "D must be between")
  expect_error(initialize_dictionary(S, D = 3, method = "invalid"), "Unknown method")
})

test_that("project_to_simplex works correctly", {
  # Test basic projection
  v <- c(1, 2, 3)
  v_proj <- spatialTV:::project_to_simplex(v, target_sum = 1)
  
  expect_equal(sum(v_proj), 1, tolerance = 1e-10)
  expect_true(all(v_proj >= 0))
  
  # Already on simplex
  v_simplex <- c(0.2, 0.3, 0.5)
  v_proj2 <- spatialTV:::project_to_simplex(v_simplex, target_sum = 1)
  expect_equal(v_proj2, v_simplex, tolerance = 1e-10)
  
  # Negative values
  v_neg <- c(-1, 2, 3)
  v_proj3 <- spatialTV:::project_to_simplex(v_neg, target_sum = 1)
  expect_equal(sum(v_proj3), 1, tolerance = 1e-10)
  expect_true(all(v_proj3 >= 0))
})

test_that("components_to_celltypes works", {
  X <- matrix(c(0.5, 0.3, 0.2,
                0.1, 0.7, 0.2), byrow = TRUE, nrow = 2, ncol = 3)
  B <- matrix(runif(5 * 3), nrow = 5, ncol = 3)
  B <- sweep(B, 2, colSums(B), "/")  # Normalize columns
  
  P <- components_to_celltypes(X, B)
  
  expect_equal(dim(P), c(2, 5))
  expect_true(all(rowSums(P) > 0.99 & rowSums(P) < 1.01))
  expect_true(all(P >= 0))
})

test_that("get_dominant_celltype works", {
  X <- matrix(c(0.1, 0.7, 0.2,
                0.6, 0.3, 0.1,
                0.2, 0.2, 0.6), byrow = TRUE, nrow = 3, ncol = 3)
  
  dominant <- get_dominant_celltype(X)
  expect_equal(as.numeric(dominant), c(2, 1, 3))
  
  # With cell type names
  names <- c("TypeA", "TypeB", "TypeC")
  dominant_named <- get_dominant_celltype(X, cell_type_names = names)
  expect_equal(as.character(dominant_named), c("TypeB", "TypeA", "TypeC"))
})

test_that("get_top_celltypes works", {
  X <- matrix(c(0.1, 0.7, 0.2,
                0.6, 0.3, 0.1), byrow = TRUE, nrow = 2, ncol = 3)
  
  top <- get_top_celltypes(X, K = 2)
  
  expect_equal(dim(top$types), c(2, 2))
  expect_equal(dim(top$proportions), c(2, 2))
  expect_equal(top$types[1, ], c(2, 3))  # Second and third highest
  expect_equal(top$proportions[1, ], c(0.7, 0.2))
})

test_that("Input validation works for fit_spatial_mixture", {
  library(Matrix)
  
  S <- matrix(rexp(50 * 5), nrow = 50, ncol = 5)
  A <- Matrix(0, nrow = 100, ncol = 50, sparse = TRUE)
  A[1:100, sample(50, 2)] <- c(1, -1)  # Random edges
  B <- matrix(runif(5 * 2), nrow = 5, ncol = 2)
  B <- sweep(B, 2, colSums(B), "/")
  
  # Should work
  expect_silent({
    result <- fit_spatial_mixture(S, A, B, lambda = 1, it_X = 10, it_B = 5, 
                                   log_every_X = 0, log_every_B = 0)
  })
  
  # Wrong S type
  expect_error(fit_spatial_mixture(as.data.frame(S), A, B), "S must be a matrix")
  
  # Wrong A type
  expect_error(fit_spatial_mixture(S, as.matrix(A), B), "A must be a sparse")
  
  # Dimension mismatch
  B_wrong <- matrix(runif(4 * 2), nrow = 4, ncol = 2)
  expect_error(fit_spatial_mixture(S, A, B_wrong), "Number of rows in B")
})
