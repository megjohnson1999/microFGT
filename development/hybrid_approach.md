# Hybrid Approach for microFGT

## Overview

We're adopting a hybrid approach for the microFGT package to balance Bioconductor compatibility with ease of development and testing:

1. **Minimal S4 Class Layer**: A thin layer of S4 classes for extending TreeSummarizedExperiment
2. **Function-Centric Design**: Most functionality implemented as standalone functions
3. **S3 Methods**: Using S3 for visualization and common operations

## Implementation Plan

### 1. Core S4 Class Structure

We'll simplify the S4 class structure to only include:

```r
# Core S4 class - minimal implementation
setClass("FGTExperiment",
         contains = "TreeSummarizedExperiment",
         slots = list(
           experimentType = "character",
           fgtMetadata = "SimpleList"
         ))

# Constructor - forwards to TreeSummarizedExperiment
FGTExperiment <- function(counts, ...) {
  # Create TSE
  tse <- TreeSummarizedExperiment(counts, ...)
  
  # Wrap in FGTExperiment
  new("FGTExperiment", 
      tse,
      experimentType = "amplicon",
      fgtMetadata = SimpleList())
}

# Basic accessor methods
setGeneric("experimentType", function(x) standardGeneric("experimentType"))
setMethod("experimentType", "FGTExperiment", function(x) x@experimentType)

setGeneric("fgtMetadata", function(x) standardGeneric("fgtMetadata"))
setMethod("fgtMetadata", "FGTExperiment", function(x) x@fgtMetadata)
```

### 2. Function-Based Design

Most functionality will be implemented as standalone functions:

```r
# Transformation function - works with any SummarizedExperiment-like object
transform_abundance <- function(se, method = "relative", assay_name = "counts") {
  # Extract the assay
  counts <- assays(se)[[assay_name]]
  
  # Transform
  transformed <- switch(method,
                        "relative" = { /* implementation */ },
                        "log" = { /* implementation */ })
  
  # Add as new assay
  assays(se)[[method]] <- transformed
  return(se)
}

# Filter function - works with any SummarizedExperiment-like object
filter_taxa <- function(se, min_prevalence = 0.1, ...) {
  # Implementation
  # Returns the filtered object
}
```

### 3. S3 Methods for Visualization

```r
# S3 generic
plot_composition <- function(x, ...) {
  UseMethod("plot_composition")
}

# Method for FGTExperiment
plot_composition.FGTExperiment <- function(x, assay_name = "counts", ...) {
  # Implementation for FGTExperiment objects
}

# Method for other SummarizedExperiment objects
plot_composition.SummarizedExperiment <- function(x, assay_name = "counts", ...) {
  # Implementation for general SummarizedExperiment objects
}

# Default method
plot_composition.default <- function(x, ...) {
  stop("Don't know how to plot composition for objects of class ", class(x))
}
```

## Benefits of the Hybrid Approach

1. **Simpler Code Structure**:
   - Clear separation between class definitions and functionality
   - Most functions are standalone and easier to test
   - Minimal S4 complexity

2. **Better Interoperability**:
   - Functions work with any SummarizedExperiment-like object
   - Users can opt-in to FGTExperiment or stay with TreeSummarizedExperiment
   - Gradual adoption path for users

3. **Easier Testing**:
   - Test S4 classes once, separately
   - Most tests can use simple matrices or SummarizedExperiment objects
   - Reduced complexity in test fixtures

4. **Maintainable Long-Term**:
   - New functionality can be added without touching the class system
   - Clear separation of concerns
   - More aligned with functional programming paradigms common in R

## Implementation Sequence

1. Create minimal S4 class layer
2. Implement core functions working with any SummarizedExperiment
3. Add S3 visualization methods
4. Create comprehensive test suite