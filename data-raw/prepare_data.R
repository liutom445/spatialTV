## Code to prepare example data for the package
## This script should be run manually to generate the data
##
## Design philosophy:
## - S and A are user inputs (loaded from external sources)
## - B is generated from S using initialize_dictionary()

library(Matrix)
devtools::load_all()  # Load package functions including initialize_dictionary()

# Load the CSV files (S and A are user inputs)
S_csv <- read.csv("../../../mnt/project/S.csv", row.names = 1)
A_csv <- read.csv("../../../mnt/project/A.csv", row.names = 1)

# Convert S to matrix
example_S <- as.matrix(S_csv)

# Convert A to sparse matrix
example_A <- as(as.matrix(A_csv), "dgCMatrix")

# Generate B from S using initialize_dictionary (as per package design)
# B should be generated based on S, not loaded from external files
example_B <- initialize_dictionary(example_S, D = 3, method = "kmeans", seed = 42)

# Save as internal data
usethis::use_data(example_S, example_A, example_B, 
                  internal = FALSE, overwrite = TRUE)

# Create a smaller example for quick demos (first 100 spots)
example_S_small <- example_S[1:100, ]
example_A_small <- example_A[, 1:100]
# Keep only edges that connect spots in the subset
edge_mask <- rowSums(abs(example_A_small)) == 2
example_A_small <- example_A_small[edge_mask, ]

usethis::use_data(example_S_small, example_A_small, 
                  internal = FALSE, overwrite = TRUE)
