# spatialTV 0.1.0

## Initial Release

### New Features

* `fit_spatial_mixture()`: Core FISTA optimization with fixed edge weights
* `fit_spatial_mixture_adaptive()`: Alternating optimization with learned edge weights via LP
* `build_spatial_graph()`: Construct Delaunay triangulation from spatial coordinates
* `initialize_dictionary()`: Multiple initialization methods (k-means, random, uniform)
* `evaluate_objective()`: Compute objective value and components
* `components_to_celltypes()`: Convert component mixtures to cell type proportions
* `get_dominant_celltype()`: Extract dominant cell type per spot
* `get_top_celltypes()`: Get top K cell types per spot

### Implementation

* Efficient C++ backend using Rcpp and RcppArmadillo
* Sparse matrix support via Matrix package
* Optional CVXR integration for weight learning
* Comprehensive input validation and error messages

### Documentation

* Complete function documentation with examples
* Introductory vignette with detailed workflow
* Mathematical formulation and algorithm details
* Performance benchmarks and tuning guidelines

### Testing

* Unit tests for all major functions
* Input validation tests
* Simplex projection tests
* Integration tests

### Infrastructure

* GitHub Actions CI/CD workflow
* Code coverage tracking
* CRAN compliance checks
