## Code to prepare example data for the package
## This script should be run manually to generate the data

library(Matrix)

# Load the CSV files
S_csv <- read.csv("../../../mnt/project/S.csv", row.names = 1)
A_csv <- read.csv("../../../mnt/project/A.csv", row.names = 1)
B_csv <- read.csv("../../../mnt/project/B.csv", row.names = 1)

# Convert S to matrix
example_S <- as.matrix(S_csv)

# Convert A to sparse matrix
example_A <- as(as.matrix(A_csv), "dgCMatrix")

# Convert B to matrix
example_B <- as.matrix(B_csv)

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
