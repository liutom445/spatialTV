# spatialTV Package: Installation and Usage Guide

## For Package Developers

### Prerequisites

Ensure you have the following installed:

```r
# Required R packages
install.packages(c(
  "Rcpp",
  "RcppArmadillo",
  "Matrix",
  "devtools",
  "roxygen2",
  "testthat",
  "knitr",
  "rmarkdown"
))

# Suggested packages
install.packages(c(
  "CVXR",       # For adaptive weight learning
  "interp",     # For Delaunay triangulation
  "ggplot2",    # For visualization
  "viridis"     # For color scales
))
```

### System Requirements

- **C++ Compiler**: C++11 compatible
  - Linux: g++ >= 4.8
  - macOS: clang with Xcode Command Line Tools
  - Windows: Rtools >= 4.0

- **Libraries**: BLAS/LAPACK (usually included with R)

### Building from Source

```bash
# Clone or navigate to package directory
cd spatialTV

# Generate Rcpp exports
Rscript -e "Rcpp::compileAttributes('.')"

# Generate documentation
Rscript -e "roxygen2::roxygenise()"

# Build package
R CMD build .

# Check package
R CMD check spatialTV_0.1.0.tar.gz --as-cran

# Install
R CMD INSTALL spatialTV_0.1.0.tar.gz
```

Or use devtools:

```r
setwd("spatialTV")
devtools::load_all()      # For development
devtools::document()      # Generate documentation
devtools::test()          # Run tests
devtools::check()         # Full package check
devtools::install()       # Install package
```

### Quick Development Workflow

```r
# 1. Make changes to R or C++ code

# 2. Recompile and reload
devtools::clean_dll()
devtools::load_all()

# 3. Test changes
devtools::test()

# 4. Update documentation if needed
devtools::document()

# 5. Full check before commit
devtools::check()
```

## For Package Users

### Installation from GitHub

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install spatialTV
devtools::install_github("yourusername/spatialTV")

# With vignettes (slower but includes tutorials)
devtools::install_github("yourusername/spatialTV", 
                        build_vignettes = TRUE)
```

### Installation from CRAN (when published)

```r
install.packages("spatialTV")
```

### Basic Usage

```r
library(spatialTV)

# View vignettes
vignette("introduction", package = "spatialTV")

# Get help
?fit_spatial_mixture
?build_spatial_graph

# Run examples
example(fit_spatial_mixture)
```

## Common Issues and Solutions

### Issue: Compilation errors

**Solution 1**: Update Rcpp and RcppArmadillo
```r
install.packages(c("Rcpp", "RcppArmadillo"))
```

**Solution 2** (macOS): Install Xcode Command Line Tools
```bash
xcode-select --install
```

**Solution 3** (Windows): Install/update Rtools
- Download from: https://cran.r-project.org/bin/windows/Rtools/

### Issue: "cannot find -lgfortran" (Linux)

**Solution**: Install gfortran
```bash
# Ubuntu/Debian
sudo apt-get install gfortran

# CentOS/RHEL
sudo yum install gcc-gfortran
```

### Issue: CVXR solver fails

**Solution**: Install additional solvers
```r
install.packages("ECOSolveR")  # Usually included with CVXR
```

Or disable adaptive weighting:
```r
result <- fit_spatial_mixture_adaptive(..., learn_weights = FALSE)
```

### Issue: Memory issues with large datasets

**Solution 1**: Reduce iterations
```r
result <- fit_spatial_mixture(S, A, B, 
                              it_X = 500,   # Default: 1500
                              it_B = 200)   # Default: 500
```

**Solution 2**: Use smaller D (number of components)
```r
B <- initialize_dictionary(S, D = 3)  # Instead of D = 8
```

**Solution 3**: Subset data spatially
```r
# Work with spatial subregions
subset_idx <- 1:500  # First 500 spots
S_sub <- S[subset_idx, ]
A_sub <- A[, subset_idx]
```

## Package Structure

```
spatialTV/
├── DESCRIPTION              # Package metadata
├── NAMESPACE               # Exported functions
├── LICENSE                 # License file
├── README.md              # Package overview
├── NEWS.md                # Version history
│
├── R/                     # R functions
│   ├── fit_spatial_mixture.R
│   ├── fit_adaptive.R
│   ├── spatial_graph.R
│   ├── utilities.R
│   ├── spatialTV-package.R
│   └── RcppExports.R
│
├── src/                   # C++ source
│   ├── FISTA_w.cpp
│   ├── RcppExports.cpp
│   ├── Makevars
│   └── Makevars.win
│
├── man/                   # Documentation (auto-generated)
│   └── *.Rd
│
├── vignettes/            # Long-form documentation
│   └── introduction.Rmd
│
├── tests/                # Unit tests
│   ├── testthat.R
│   └── testthat/
│       └── test-*.R
│
├── data/                 # Example datasets (optional)
│   └── *.rda
│
├── data-raw/            # Data preparation scripts
│   └── prepare_data.R
│
└── inst/                # Additional files
    └── extdata/         # External data files
```

## Development Guidelines

### Adding New Features

1. Write R wrapper function in `R/`
2. Add roxygen2 documentation
3. Write unit tests in `tests/testthat/`
4. Update NEWS.md
5. Run `devtools::check()`

### C++ Modifications

1. Edit `src/FISTA_w.cpp`
2. Ensure functions are exported with `// [[Rcpp::export]]`
3. Rebuild: `Rcpp::compileAttributes()`
4. Test: `devtools::test()`

### Documentation

- Use roxygen2 format in R files
- Run `devtools::document()` to generate `.Rd` files
- Update vignettes in `vignettes/`

### Testing

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-main-functions.R")

# Check code coverage
covr::package_coverage()
```

## Preparing for CRAN Submission

1. **Check package thoroughly**:
```r
devtools::check(cran = TRUE)
devtools::check_rhub()
devtools::check_win_devel()
```

2. **Update version** in DESCRIPTION

3. **Update NEWS.md** with changes

4. **Check all URLs** are valid

5. **Ensure examples run** quickly (< 5 seconds each)

6. **Submit**:
```r
devtools::submit_cran()
```

## Continuous Integration

The package includes GitHub Actions workflows for:

- R-CMD-check on multiple platforms (Linux, macOS, Windows)
- Multiple R versions (devel, release, oldrel)
- Code coverage tracking

Workflows are in `.github/workflows/`

## Getting Help

- **Documentation**: `?functionname` or `help(package = "spatialTV")`
- **Vignettes**: `browseVignettes("spatialTV")`
- **Issues**: https://github.com/yourusername/spatialTV/issues
- **Email**: your.email@example.com

## Citation

```r
citation("spatialTV")
```

## License

MIT License - see LICENSE file for details.
