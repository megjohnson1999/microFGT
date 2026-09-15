#!/usr/bin/env Rscript

# microFGT Validation Script
# 
# This script validates core functionality of the microFGT package.
# Run it after installation to check if the package is working correctly.

# Load required packages
library(microFGT)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)

# Helper function to report test results
report <- function(name, result) {
  status <- if (result) "\u2705 PASS" else "\u274C FAIL"
  cat(sprintf("[%s] %s\n", status, name))
  return(result)
}

# Test 1: Create FGTExperiment with original constructor
test_original_constructor <- function() {
  tryCatch({
    counts <- matrix(rpois(50, 20), nrow=10, ncol=5)
    rownames(counts) <- paste0("Feature", 1:10)
    colnames(counts) <- paste0("Sample", 1:5)
    
    fgt <- FGTExperiment(assays = list(counts = counts))
    
    return(is(fgt, "FGTExperiment"))
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Test 2: Create FGTExperiment with the workaround method
test_workaround <- function() {
  tryCatch({
    counts <- matrix(rpois(50, 20), nrow=10, ncol=5)
    rownames(counts) <- paste0("Feature", 1:10)
    colnames(counts) <- paste0("Sample", 1:5)
    
    # Create SummarizedExperiment
    se <- SummarizedExperiment(assays = list(counts = counts))
    
    # Convert to TreeSummarizedExperiment
    tse <- TreeSummarizedExperiment(se)
    
    # Create FGTExperiment
    fgt <- methods::new("FGTExperiment",
                      tse,
                      experimentType = "amplicon",
                      fgtMetadata = S4Vectors::SimpleList())
    
    # Ensure assay is named
    assayNames(fgt)[1] <- "counts"
    
    return(is(fgt, "FGTExperiment") && identical(assayNames(fgt)[1], "counts"))
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Test 3: Create FGTExperiment with fixed constructor (if available)
test_fixed_constructor <- function() {
  tryCatch({
    if (!exists("create_FGTExperiment")) {
      # Try to source the fixed constructor
      if (file.exists(file.path(system.file(package = "microFGT"), "R", "fixed_constructor.R"))) {
        source(file.path(system.file(package = "microFGT"), "R", "fixed_constructor.R"))
      } else if (file.exists("../R/fixed_constructor.R")) {
        source("../R/fixed_constructor.R")
      } else {
        cat("  Fixed constructor not found\n")
        return(FALSE)
      }
    }
    
    counts <- matrix(rpois(50, 20), nrow=10, ncol=5)
    rownames(counts) <- paste0("Feature", 1:10)
    colnames(counts) <- paste0("Sample", 1:5)
    
    fgt <- create_FGTExperiment(assays = list(counts = counts))
    
    return(is(fgt, "FGTExperiment") && identical(assayNames(fgt)[1], "counts"))
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Test 4: Test transform_abundance function
test_transform_abundance <- function() {
  tryCatch({
    # First create a valid FGTExperiment
    counts <- matrix(rpois(50, 20), nrow=10, ncol=5)
    rownames(counts) <- paste0("Feature", 1:10)
    colnames(counts) <- paste0("Sample", 1:5)
    
    # Create with workaround method
    se <- SummarizedExperiment(assays = list(counts = counts))
    tse <- TreeSummarizedExperiment(se)
    fgt <- methods::new("FGTExperiment", tse, experimentType = "amplicon", 
                      fgtMetadata = S4Vectors::SimpleList())
    assayNames(fgt)[1] <- "counts"
    
    # Check if function exists
    if (!exists("transform_abundance", mode = "function", where = "package:microFGT")) {
      cat("  transform_abundance function not found\n")
      return(FALSE)
    }
    
    # Try transformation
    tryCatch({
      # Note: using try() because the function might exist but have issues
      fgt_rel <- transform_abundance(fgt, type = "relative", assay_name = "counts")
      return(is(fgt_rel, "FGTExperiment") && "relative" %in% assayNames(fgt_rel))
    }, error = function(e) {
      cat("  Error in transform_abundance:", conditionMessage(e), "\n")
      return(FALSE)
    })
  }, error = function(e) {
    cat("  Error in test setup:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Test 5: Test filter_taxa function
test_filter_taxa <- function() {
  tryCatch({
    # First create a valid FGTExperiment
    counts <- matrix(rpois(50, 20), nrow=10, ncol=5)
    rownames(counts) <- paste0("Feature", 1:10)
    colnames(counts) <- paste0("Sample", 1:5)
    
    # Create with workaround method
    se <- SummarizedExperiment(assays = list(counts = counts))
    tse <- TreeSummarizedExperiment(se)
    fgt <- methods::new("FGTExperiment", tse, experimentType = "amplicon", 
                      fgtMetadata = S4Vectors::SimpleList())
    assayNames(fgt)[1] <- "counts"
    
    # Check if function exists
    if (!exists("filter_taxa", mode = "function", where = "package:microFGT")) {
      cat("  filter_taxa function not found\n")
      return(FALSE)
    }
    
    # Try filtering
    tryCatch({
      # Note: using try() because the function might exist but have issues
      fgt_filtered <- filter_taxa(fgt, min_prevalence = 0.5, min_abundance = 0.001, 
                                  assay_name = "counts")
      return(is(fgt_filtered, "FGTExperiment"))
    }, error = function(e) {
      cat("  Error in filter_taxa:", conditionMessage(e), "\n")
      return(FALSE)
    })
  }, error = function(e) {
    cat("  Error in test setup:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Run the tests
cat("\n=== microFGT Validation Results ===\n\n")

cat("Testing core object creation:\n")
report("Original constructor (FGTExperiment)", test_original_constructor())
report("Workaround method (new + SummarizedExperiment)", test_workaround())
report("Fixed constructor (create_FGTExperiment)", test_fixed_constructor())

cat("\nTesting core functions:\n")
report("transform_abundance", test_transform_abundance())
report("filter_taxa", test_filter_taxa())

cat("\n=== Summary of Installed Packages ===\n")
installed_pkgs <- c("microFGT", "TreeSummarizedExperiment", "SummarizedExperiment", 
                     "MultiAssayExperiment", "S4Vectors", "Biostrings", "BiocParallel")
                     
for (pkg in installed_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("- %s: Installed, version %s\n", pkg, 
                as.character(packageVersion(pkg))))
  } else {
    cat(sprintf("- %s: Not installed\n", pkg))
  }
}

cat("\n=== For help with issues, see the Development Guide ===\n")
cat("DEVELOPMENT.md contains known issues and workarounds\n")