# Testing Environment for microFGT

This document outlines the recommended testing environment setup for microFGT package development.

## Local Testing Environment

### Recommended Setup: R Package Development with `renv`

The [`renv`](https://rstudio.github.io/renv/) package provides project-local R dependency management, similar to virtual environments in Python. This ensures consistent, reproducible package testing across different development environments.

#### Setup Instructions

1. **Initialize renv in your project:**

```r
# Install renv if not already installed
install.packages("renv")

# Initialize renv in the project
renv::init()
```

2. **Install package dependencies:**

```r
# Install dependencies from DESCRIPTION
renv::install()

# Install additional development dependencies
renv::install(c("covr", "devtools", "lintr", "styler", "pkgdown"))
```

3. **Snapshot the environment:**

```r
# Save the current state of dependencies
renv::snapshot()
```

4. **Run tests in the isolated environment:**

```r
# Load the project environment
renv::load()

# Run tests
devtools::test()
```

### Bioconductor Dependencies

Since microFGT relies on several Bioconductor packages, ensure they're installed correctly:

```r
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor dependencies
BiocManager::install(c(
    "TreeSummarizedExperiment",
    "SummarizedExperiment",
    "S4Vectors",
    "BiocParallel",
    "Biostrings",
    "DESeq2"
))
```

## Docker-Based Testing

For completely isolated and reproducible testing, consider using Docker:

### Example Dockerfile

```dockerfile
FROM bioconductor/bioconductor_docker:RELEASE_3_17

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev

# Install R dependencies
RUN R -e "install.packages(c('devtools', 'testthat', 'covr', 'renv'))"
RUN R -e "BiocManager::install(c('TreeSummarizedExperiment', 'SummarizedExperiment', 'S4Vectors', 'BiocParallel', 'Biostrings', 'DESeq2'))"

# Copy the package
COPY . /microFGT
WORKDIR /microFGT

# Check and test the package
CMD R CMD build . && R CMD check microFGT_*.tar.gz && R -e "testthat::test_package('microFGT')"
```

### Building and running the Docker container:

```bash
# Build the Docker image
docker build -t microfgt-test .

# Run tests in the container
docker run microfgt-test
```

## GitHub Actions for CI/CD

For continuous integration, configure GitHub Actions to automatically test the package:

Create a file at `.github/workflows/R-CMD-check.yml`:

```yaml
on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main, develop]

name: R-CMD-check

jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}
    name: ${{ matrix.config.os }} (${{ matrix.config.r }})
    
    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: ubuntu-latest, r: 'release', bioc: '3.17'}
          - {os: macos-latest, r: 'release', bioc: '3.17'}

    env:
      R_REMOTES_NO_ERRORS_FROM_WARNINGS: true
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - uses: actions/checkout@v3

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}

      - uses: r-lib/actions/setup-bioc@v2
        with:
          bioc-version: ${{ matrix.config.bioc }}

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
          needs: check

      - uses: r-lib/actions/check-r-package@v2
        with:
          upload-snapshots: true
```

## Test Coverage Reporting

To track test coverage:

```r
# Install covr if needed
if (!requireNamespace("covr", quietly = TRUE))
    install.packages("covr")

# Generate coverage report
covr::report()

# Calculate coverage percentage
covr::package_coverage()
```

## Recommended Testing Workflow

1. **Local Development and Testing**:
   - Use `renv` for local dependency management
   - Run `devtools::test()` frequently during development
   - Check coverage with `covr::report()`

2. **CI/CD Pipeline**:
   - Push to GitHub to trigger GitHub Actions workflows
   - Review test results and coverage reports

3. **Comprehensive Testing**:
   - Periodically run full tests in Docker for complete isolation

## Potential Issues and Solutions

### Bioconductor Dependency Conflicts

If you encounter version conflicts with Bioconductor packages:

```r
# Use specific Bioconductor version
BiocManager::install(version = "3.17")
BiocManager::valid()  # Check for problems
```

### R Package Building Issues

For troubleshooting package building:

```r
# Check for common issues
devtools::check()

# More detailed diagnostics
rcmdcheck::rcmdcheck(args = c("--as-cran"))
```

### Testing with CRAN Policies

To ensure CRAN compliance:

```r
# Run with --as-cran flag
devtools::check(cran = TRUE)
```

### Memory Issues During Testing

If tests require significant memory:

```r
# Set memory limit in .Renviron
cat("R_MAX_VSIZE=32Gb\n", file = file.path(Sys.getenv("HOME"), ".Renviron"), append = TRUE)
```