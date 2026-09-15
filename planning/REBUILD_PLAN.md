# microFGT Package Rebuild Implementation Plan

## Overview

This document outlines a comprehensive plan for rebuilding the microFGT package with a focus on robust, maintainable architecture. The rebuild will prioritize the core data structure and essential functionality, deferring advanced analytical features until the foundation is solid.

## Core Design Principles

1. **Minimalism**: Focus on essential features first
2. **Correctness**: Ensure proper S4 implementation and registration
3. **Testability**: Design for comprehensive testing
4. **Maintainability**: Clear documentation and consistent patterns
5. **Reproducibility**: Reliable build and test process

## Phase 1: Core Architecture

### 1.1 File Organization

```
R/
├── AllClasses.R          # S4 class definitions only
├── AllGenerics.R         # Generic function definitions only
├── methods-FGTExperiment.R  # Methods for FGTExperiment class
├── constructor.R         # Constructor functions
├── validators.R          # Validation functions
├── accessors.R           # Getter/setter function implementations
├── coerce.R              # Coercion methods
├── transform.R           # Core transformation functions
├── import-export.R       # Import/export utilities
├── utils.R               # Internal utility functions
└── zzz.R                 # Package lifecycle hooks

tests/
├── testthat.R            # Test runner
└── testthat/
    ├── test-class-definition.R      # Class definition tests
    ├── test-class-validation.R      # Validation function tests
    ├── test-constructor.R           # Constructor tests
    ├── test-accessors.R             # Accessor function tests
    ├── test-coerce.R                # Coercion method tests
    ├── test-transform.R             # Transformation function tests
    └── test-import-export.R         # Import/export function tests

man/                      # Documentation (roxygen2 generated)
data-raw/                 # Data processing scripts
inst/                     # Installed files
vignettes/                # Package vignettes
src/                      # Compiled code (if needed)
```

### 1.2 S4 Class Design

#### FGTExperiment Class

```r
# In AllClasses.R
#' FGTExperiment Class
#'
#' @slot experimentType Character string indicating the experiment type
#' @slot fgtMetadata SimpleList containing additional metadata
#'
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors SimpleList
#'
#' @export
#' @exportClass FGTExperiment
setClass("FGTExperiment",
         contains = "TreeSummarizedExperiment",
         slots = list(
           experimentType = "character",
           fgtMetadata = "SimpleList"
         ),
         prototype = list(
           experimentType = "amplicon",
           fgtMetadata = S4Vectors::SimpleList()
         ))

# ValidityMethod
setValidity("FGTExperiment", function(object) {
  msg <- NULL
  
  # Check experimentType
  validTypes <- c("amplicon", "metagenomic", "integrated") 
  if (!(object@experimentType %in% validTypes)) {
    msg <- c(msg, paste("experimentType must be one of:", 
                       paste(validTypes, collapse=", ")))
  }
  
  # Return messages or TRUE
  if (is.null(msg)) TRUE else msg
})
```

### 1.3 Generic Functions

```r
# In AllGenerics.R
#' Get the experiment type
#'
#' @param object FGTExperiment object
#'
#' @return Character string indicating the experiment type
#' @export
setGeneric("experimentType", function(object) standardGeneric("experimentType"))

#' Set the experiment type
#'
#' @param object FGTExperiment object
#' @param value New experiment type value
#'
#' @return Updated FGTExperiment object
#' @export
setGeneric("experimentType<-", function(object, value) standardGeneric("experimentType<-"))

#' Get the FGT metadata
#'
#' @param object FGTExperiment object
#'
#' @return SimpleList containing the metadata
#' @export
setGeneric("fgtMetadata", function(object) standardGeneric("fgtMetadata"))

#' Set the FGT metadata
#'
#' @param object FGTExperiment object
#' @param value SimpleList containing the new metadata
#'
#' @return Updated FGTExperiment object
#' @export
setGeneric("fgtMetadata<-", function(object, value) standardGeneric("fgtMetadata<-"))
```

### 1.4 Method Implementations

```r
# In methods-FGTExperiment.R
#' @describeIn FGTExperiment Get experiment type
#' @export
setMethod("experimentType", "FGTExperiment", function(object) {
  object@experimentType
})

#' @describeIn FGTExperiment Set experiment type
#' @export
setMethod("experimentType<-", "FGTExperiment", function(object, value) {
  validTypes <- c("amplicon", "metagenomic", "integrated")
  if (!value %in% validTypes) {
    stop("experimentType must be one of: ", paste(validTypes, collapse=", "))
  }
  
  object@experimentType <- value
  validObject(object)
  return(object)
})

#' @describeIn FGTExperiment Get FGT metadata
#' @export
setMethod("fgtMetadata", "FGTExperiment", function(object) {
  object@fgtMetadata
})

#' @describeIn FGTExperiment Set FGT metadata
#' @export
setMethod("fgtMetadata<-", "FGTExperiment", function(object, value) {
  if (!methods::is(value, "SimpleList")) {
    stop("fgtMetadata must be a SimpleList")
  }
  
  object@fgtMetadata <- value
  validObject(object)
  return(object)
})
```

### 1.5 Constructor Implementation

```r
# In constructor.R
#' Create a new FGTExperiment object
#'
#' @param assays List of matrices or similar objects containing count data
#' @param rowData DataFrame of feature metadata
#' @param colData DataFrame of sample metadata
#' @param rowTree Phylogenetic tree of type phylo
#' @param experimentType Type of experiment ("amplicon", "metagenomic", "integrated")
#' @param fgtMetadata Additional metadata as a SimpleList
#' @param ... Additional arguments passed to TreeSummarizedExperiment
#'
#' @return A new FGTExperiment object
#' @export
FGTExperiment <- function(assays, rowData = NULL, colData = NULL, rowTree = NULL,
                        experimentType = "amplicon", fgtMetadata = S4Vectors::SimpleList(), ...) {
  # Input validation
  if (!is.list(assays) || length(assays) == 0) {
    stop("'assays' must be a non-empty list")
  }

  # Verify the first assay 
  first_assay <- assays[[1]]
  if (!is.matrix(first_assay) && !methods::is(first_assay, "Matrix")) {
    stop("The first element of 'assays' must be a matrix or Matrix object")
  }

  # Ensure row and column names are present
  if (is.null(rownames(first_assay))) {
    rownames(first_assay) <- paste0("Feature", seq_len(nrow(first_assay)))
    assays[[1]] <- first_assay
  }

  if (is.null(colnames(first_assay))) {
    colnames(first_assay) <- paste0("Sample", seq_len(ncol(first_assay)))
    assays[[1]] <- first_assay
  }

  # Prepare rowData if needed
  if (is.null(rowData)) {
    rowData <- S4Vectors::DataFrame(row.names = rownames(first_assay))
  } else if (!is.data.frame(rowData) && !methods::is(rowData, "DataFrame")) {
    rowData <- S4Vectors::DataFrame(rowData, row.names = rownames(first_assay))
  }

  # Prepare colData if needed
  if (is.null(colData)) {
    colData <- S4Vectors::DataFrame(row.names = colnames(first_assay))
  } else if (!is.data.frame(colData) && !methods::is(colData, "DataFrame")) {
    colData <- S4Vectors::DataFrame(colData, row.names = colnames(first_assay))
  }

  # Create SummarizedExperiment 
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = rowData,
    colData = colData,
    ...
  )

  # Convert to TreeSummarizedExperiment
  tse <- tryCatch({
    TreeSummarizedExperiment::TreeSummarizedExperiment(se, rowTree = rowTree)
  }, error = function(e) {
    warning("Error creating TreeSummarizedExperiment: ", conditionMessage(e),
            ". Using SummarizedExperiment instead.")
    se
  })

  # Create FGTExperiment object
  obj <- methods::new("FGTExperiment",
                     tse,
                     experimentType = experimentType,
                     fgtMetadata = fgtMetadata)

  validObject(obj)
  return(obj)
}
```

## Phase 2: Testing Framework

### 2.1 Unit Tests for Class Definition and Methods

```r
# In test-class-definition.R
test_that("FGTExperiment class is properly defined", {
  expect_true(isClass("FGTExperiment"))
  expect_true("TreeSummarizedExperiment" %in% getClass("FGTExperiment")@contains@names)
  expect_true(all(c("experimentType", "fgtMetadata") %in% slotNames("FGTExperiment")))
})

test_that("FGTExperiment validity method works", {
  # Create valid object
  counts <- matrix(1:12, nrow=3, ncol=4)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Sample", 1:4)
  valid_obj <- FGTExperiment(assays=list(counts=counts))
  
  # Test validation succeeds for valid object
  expect_true(validObject(valid_obj))
  
  # Create invalid object
  invalid_obj <- valid_obj
  invalid_obj@experimentType <- "invalid_type"
  
  # Test validation fails for invalid object
  expect_error(validObject(invalid_obj))
})

# In test-accessors.R
test_that("experimentType accessor methods work", {
  counts <- matrix(1:12, nrow=3, ncol=4)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Sample", 1:4)
  obj <- FGTExperiment(assays=list(counts=counts))
  
  # Test getter
  expect_equal(experimentType(obj), "amplicon")
  
  # Test setter
  experimentType(obj) <- "metagenomic"
  expect_equal(experimentType(obj), "metagenomic")
  
  # Test validation in setter
  expect_error(experimentType(obj) <- "invalid_type")
})
```

### 2.2 Constructor Tests

```r
# In test-constructor.R
test_that("FGTExperiment constructor works with minimal input", {
  counts <- matrix(1:12, nrow=3, ncol=4)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Sample", 1:4)
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  expect_s4_class(obj, "FGTExperiment")
  expect_equal(dim(obj), c(3, 4))
  expect_equal(experimentType(obj), "amplicon")
  expect_true(length(fgtMetadata(obj)) == 0)
})

test_that("FGTExperiment constructor handles missing row/colnames", {
  counts <- matrix(1:12, nrow=3, ncol=4)
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  expect_equal(rownames(obj), c("Feature1", "Feature2", "Feature3"))
  expect_equal(colnames(obj), c("Sample1", "Sample2", "Sample3", "Sample4"))
})

test_that("FGTExperiment constructor validates inputs", {
  expect_error(FGTExperiment(assays=list()))
  expect_error(FGTExperiment(assays=list(counts="not_a_matrix")))
})
```

### 2.3 Test Helper Functions

```r
# In helper-microFGT.R
#' Create a test FGTExperiment object
#'
#' @param rows Number of rows
#' @param cols Number of columns
#' @param withTree Include a phylogenetic tree
#' @param experimentType Type of experiment
#'
#' @return FGTExperiment object
create_test_object <- function(rows=10, cols=5, withTree=FALSE, experimentType="amplicon") {
  # Create count matrix
  counts <- matrix(rpois(rows*cols, lambda=10), nrow=rows, ncol=cols)
  rownames(counts) <- paste0("Feature", seq_len(rows))
  colnames(counts) <- paste0("Sample", seq_len(cols))
  
  # Create row data
  rowData <- S4Vectors::DataFrame(
    FeatureID = rownames(counts),
    row.names = rownames(counts)
  )
  
  # Create column data
  colData <- S4Vectors::DataFrame(
    SampleID = colnames(counts),
    Condition = rep(c("A", "B"), length.out=cols),
    row.names = colnames(counts)
  )
  
  # Create tree if requested
  rowTree <- NULL
  if (withTree) {
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Package 'ape' needed to create test tree")
    }
    rowTree <- ape::rtree(rows)
    rowTree$tip.label <- rownames(counts)
  }
  
  # Create FGTExperiment
  FGTExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData,
    rowTree = rowTree,
    experimentType = experimentType,
    fgtMetadata = S4Vectors::SimpleList(
      CreatedBy = "test_helper",
      Date = Sys.Date()
    )
  )
}
```

## Phase 3: Transformation Functions

### 3.1 Core Transformation Functions

```r
# In transform.R
#' Transform abundance values
#'
#' @param object FGTExperiment or SummarizedExperiment object
#' @param type Transformation type: "relative", "log", "clr", "presence"
#' @param assay_name Name of assay to transform
#' @param pseudocount Value to add before log transformations
#'
#' @return FGTExperiment or SummarizedExperiment with added assay containing the transformed values
#' @export
setGeneric("transformAbundance", function(object, type = "relative", assay_name = "counts", pseudocount = 1) 
  standardGeneric("transformAbundance"))

#' @describeIn transformAbundance Transform abundance for FGTExperiment objects
#' @export
setMethod("transformAbundance", "FGTExperiment", function(object, type = "relative", assay_name = "counts", pseudocount = 1) {
  # Validate inputs
  validTypes <- c("relative", "log", "clr", "presence")
  if (!type %in% validTypes) {
    stop("type must be one of: ", paste(validTypes, collapse=", "))
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
    stop("assay '", assay_name, "' not found")
  }
  
  # Get assay data
  counts <- SummarizedExperiment::assay(object, assay_name)
  
  # Perform transformation
  result <- switch(
    type,
    "relative" = {
      # Relative abundance (proportions)
      t(apply(counts, 2, function(x) {
        if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) {
          return(rep(NA, length(x)))
        }
        return(x / sum(x, na.rm = TRUE))
      }))
    },
    "log" = {
      # Log transformation (log(x + pseudocount))
      log1p(counts)
    },
    "clr" = {
      # Centered log-ratio transformation
      t(apply(counts + pseudocount, 2, function(x) {
        log_x <- log(x)
        log_x - mean(log_x)
      }))
    },
    "presence" = {
      # Presence/absence (logical)
      counts > 0
    }
  )
  
  # Ensure dimensions and names are preserved
  dimnames(result) <- dimnames(counts)
  
  # Create a copy with the new assay
  new_object <- object
  SummarizedExperiment::assays(new_object)[[paste0(type, "_", assay_name)]] <- result
  
  return(new_object)
})

#' @describeIn transformAbundance Transform abundance for SummarizedExperiment objects
#' @export
setMethod("transformAbundance", "SummarizedExperiment", function(object, type = "relative", assay_name = "counts", pseudocount = 1) {
  # Implementation identical to FGTExperiment method
  # This avoids code duplication through inheritance
  validTypes <- c("relative", "log", "clr", "presence")
  if (!type %in% validTypes) {
    stop("type must be one of: ", paste(validTypes, collapse=", "))
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
    stop("assay '", assay_name, "' not found")
  }
  
  # Get assay data
  counts <- SummarizedExperiment::assay(object, assay_name)
  
  # Perform transformation
  result <- switch(
    type,
    "relative" = {
      # Relative abundance (proportions)
      t(apply(counts, 2, function(x) {
        if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) {
          return(rep(NA, length(x)))
        }
        return(x / sum(x, na.rm = TRUE))
      }))
    },
    "log" = {
      # Log transformation (log(x + pseudocount))
      log1p(counts)
    },
    "clr" = {
      # Centered log-ratio transformation
      t(apply(counts + pseudocount, 2, function(x) {
        log_x <- log(x)
        log_x - mean(log_x)
      }))
    },
    "presence" = {
      # Presence/absence (logical)
      counts > 0
    }
  )
  
  # Ensure dimensions and names are preserved
  dimnames(result) <- dimnames(counts)
  
  # Create a copy with the new assay
  new_object <- object
  SummarizedExperiment::assays(new_object)[[paste0(type, "_", assay_name)]] <- result
  
  return(new_object)
})
```

## Phase 4: CI/CD Setup

### 4.1 GitHub Actions Workflows

#### R-CMD-check.yaml

```yaml
name: R-CMD-check

on:
  push:
    branches: [main, develop, rebuild-core]
    paths-ignore:
      - '**.md'
      - 'docs/**'
      - '.github/workflows/**'
  pull_request:
    branches: [main, develop]
    paths-ignore:
      - '**.md'
      - 'docs/**'
  workflow_dispatch:

concurrency:
  group: ${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true

jobs:
  R-CMD-check:
    runs-on: ${{ matrix.config.os }}
    name: ${{ matrix.config.os }} (${{ matrix.config.r }})

    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: ubuntu-22.04, r: 'release', bioc: '3.17'}
          - {os: ubuntu-22.04, r: 'devel', bioc: '3.18'}
          - {os: macOS-latest, r: 'release', bioc: '3.17'}
          - {os: windows-latest, r: 'release', bioc: '3.17'}

    env:
      R_REMOTES_NO_ERRORS_FROM_WARNINGS: true
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-pandoc@v2

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}
          http-user-agent: ${{ matrix.config.http-user-agent }}
          use-public-rspm: true
          bioc-version: ${{ matrix.config.bioc }}

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
          needs: check

      - name: Run namespace validation
        run: |
          Rscript -e "source('tests/scripts/validate_namespace.R'); validate_namespace()"

      - uses: r-lib/actions/check-r-package@v2
        with:
          upload-snapshots: true
```

#### s4-class-check.yaml

```yaml
name: S4 Class Validation

on:
  push:
    branches: [main, develop, rebuild-core]
    paths:
      - "R/**"
      - "NAMESPACE"
      - "DESCRIPTION"
  pull_request:
    branches: [main, develop]
    paths:
      - "R/**"
      - "NAMESPACE"
      - "DESCRIPTION"
  workflow_dispatch:

concurrency:
  group: ${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true

jobs:
  s4-class-check:
    name: Validate S4 Classes
    runs-on: ubuntu-22.04
    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - name: Check out repository
        uses: actions/checkout@v4

      - name: Set up R
        uses: r-lib/actions/setup-r@v2
        with:
          r-version: 'release'
          use-public-rspm: true

      - name: Install dependencies
        uses: r-lib/actions/setup-r-dependencies@v2
        with:
          dependencies: '"all"'
          extra-packages: |
            any::rcmdcheck
            any::roxygen2

      - name: Install package
        run: |
          R CMD build .
          PKG_FILE=$(ls -1t *.tar.gz | head -n 1)
          R CMD INSTALL ${PKG_FILE}

      - name: Check S4 class registration
        run: |
          Rscript -e '
          library(methods)
          library(microFGT)

          # Check class registration
          check_class <- function(class_name) {
            is_registered <- isClass(class_name)
            if (!is_registered) {
              cat("✗ Class", class_name, "is not properly registered\n")
              return(FALSE)
            }
            
            cat("✓ Class", class_name, "is properly registered\n")
            return(TRUE)
          }

          # Check method registration
          check_method <- function(generic_name, class_name) {
            method <- tryCatch({
              selectMethod(generic_name, class_name, optional = TRUE)
            }, error = function(e) NULL)
            
            if (is.null(method)) {
              cat("✗ Method", generic_name, "for", class_name, "not found\n")
              return(FALSE)
            }
            
            cat("✓ Method", generic_name, "for", class_name, "resolves correctly\n")
            return(TRUE)
          }

          # Check classes
          classes_ok <- check_class("FGTExperiment")

          # Check methods
          methods_ok <- c(
            check_method("experimentType", "FGTExperiment"),
            check_method("experimentType<-", "FGTExperiment"),
            check_method("fgtMetadata", "FGTExperiment"),
            check_method("fgtMetadata<-", "FGTExperiment"),
            check_method("transformAbundance", "FGTExperiment")
          )
          
          # Check if we can create an object
          obj_ok <- tryCatch({
            counts <- matrix(1:12, nrow=3, ncol=4)
            rownames(counts) <- paste0("Gene", 1:3)
            colnames(counts) <- paste0("Sample", 1:4)
            obj <- FGTExperiment(assays=list(counts=counts))
            cat("✓ Successfully created FGTExperiment object\n")
            TRUE
          }, error = function(e) {
            cat("✗ Failed to create FGTExperiment object:", conditionMessage(e), "\n")
            FALSE
          })
          
          # Overall result
          all_ok <- all(classes_ok, all(methods_ok), obj_ok)
          if (!all_ok) {
            cat("✗ S4 class validation failed\n")
            quit(status = 1)
          }
          
          cat("✓ All S4 class checks passed\n")
          '
```

### 4.2 Pre-commit Hook Configuration

```bash
#!/bin/bash
# Pre-commit hook for microFGT package

# Print colored output
green() { echo -e "\033[0;32m$1\033[0m"; }
yellow() { echo -e "\033[0;33m$1\033[0m"; }
red() { echo -e "\033[0;31m$1\033[0m"; }

# Capture the start time
start_time=$(date +%s)

# Show header
echo "=========================================================="
echo "              microFGT Package Pre-commit Hook            "
echo "=========================================================="
echo "Running pre-commit checks - this may take a minute..."

# Initialize error flag
errors_found=0

# Function to run a check and report status
run_check() {
  local name="$1"
  local command="$2"
  
  echo "Running check: $name"
  if eval "$command"; then
    green "✓ $name check passed"
    return 0
  else
    red "✗ $name check failed"
    errors_found=1
    return 1
  fi
}

# Regenerate NAMESPACE
run_check "Roxygen documentation update" "Rscript -e 'roxygen2::roxygenize()'" || {
  yellow "Tip: Roxygen documentation needs updating"
}

# Check NAMESPACE validity
run_check "NAMESPACE validation" "Rscript tests/scripts/validate_namespace.R" || {
  yellow "Tip: Fix NAMESPACE issues by updating roxygen documentation"
}

# Check S4 class definitions
run_check "S4 class validation" "Rscript tests/scripts/validate_s4_classes.R" || {
  yellow "Tip: Check for proper S4 class registration and method definitions"
}

# Run tests
run_check "Unit tests" "Rscript -e 'testthat::test_package(\"microFGT\")'" || {
  yellow "Tip: Fix failing tests before committing"
}

# Run R CMD check without building vignettes or manual
run_check "R CMD check" "R CMD check --no-manual --no-build-vignettes --no-examples ." || {
  yellow "Tip: Fix R CMD check issues before committing"
}

# Calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))

# Report result
echo "=========================================================="
if [ $errors_found -eq 0 ]; then
  green "✓ All pre-commit checks passed! (completed in ${duration}s)"
  echo "=========================================================="
  exit 0
else
  red "✗ Pre-commit checks failed (completed in ${duration}s)"
  yellow "Fix the issues above before committing your changes."
  echo "=========================================================="
  exit 1
fi
```

## Phase 5: Development Workflow

### 5.1 Package Documentation

The package will be fully documented using roxygen2 with appropriate tags for:
- Class descriptions (@slot, @exportClass)
- Method documentation (@describeIn, @export)
- Function parameters (@param, @return)
- Examples (@examples)

### 5.2 Development Guidelines

1. **NAMESPACE Management**:
   - Never manually edit NAMESPACE file
   - Always use roxygen2 tags for exports
   - Run roxygen2::roxygenize() before commits

2. **Class Extensions**:
   - Keep class hierarchy simple
   - Validate objects using validObject()
   - Follow Bioconductor S4 patterns

3. **Testing**:
   - Write tests first (TDD approach)
   - Test both normal and edge cases
   - Maintain 100% test coverage for core functionality

4. **File Organization**:
   - One class definition per file
   - Group related methods in the same file
   - Use snake_case for file names and function names
   - Use CamelCase for class names

### 5.3 Development Process

1. Create feature branch
2. Implement feature with tests
3. Validate with pre-commit hook
4. Submit pull request
5. Review and address CI failures
6. Merge to develop branch
7. Periodically merge develop to main for releases

## Implementation Timeline

1. **Week 1**: Core architecture setup
   - Class definitions
   - Method implementations
   - Constructor functions
   - Basic tests

2. **Week 2**: Core functionality
   - Transformation functions
   - Import/export utilities
   - CI/CD pipelines
   - Detailed tests

3. **Week 3**: Documentation and examples
   - Package vignettes
   - Example datasets
   - Developer documentation

## Success Criteria

- All S4 classes are properly registered
- All methods resolve correctly
- Package passes R CMD check without warnings
- Test coverage is >90%
- Documentation is complete and accurate
- CI/CD pipeline works reliably
- Pre-commit hook catches NAMESPACE issues

## Summary

This rebuild plan provides a structured approach to reimplementing the microFGT package with a focus on robust S4 class design and maintainability. By following this plan, we will create a solid foundation that can be extended with more advanced features in the future, while avoiding the NAMESPACE and S4 registration issues that have plagued the package so far.