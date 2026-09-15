# Comprehensive edge case tests for microFGT package

library(testthat)
library(microFGT)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)

# Source helper functions
source("helper-test-reporting.R")

# Skip if FGTExperiment is not available
skip_if_no_FGTExperiment <- function() {
  if (!exists("FGTExperiment", mode = "function")) {
    skip("FGTExperiment class not available")
  }
}

# Test suite for edge cases
context("Comprehensive Edge Case Testing")

# 1. Empty Data Structures
test_that("FGTExperiment handles completely empty data", {
  skip_if_no_FGTExperiment()
  
  # Test with zero rows and columns
  expect_error(
    FGTExperiment(
      assays = list(counts = matrix(numeric(0), nrow = 0, ncol = 0)),
      experimentType = "amplicon"
    ),
    "must have at least one row and one column"
  )
  
  # Test with empty but properly dimensioned matrix
  empty_counts <- matrix(0, nrow = 5, ncol = 3)
  rownames(empty_counts) <- paste0("Feature", 1:5)
  colnames(empty_counts) <- paste0("Sample", 1:3)
  
  fgt <- FGTExperiment(
    assays = list(counts = empty_counts),
    experimentType = "amplicon"
  )
  
  expect_equal(sum(assay(experimentData(fgt))), 0)
  expect_equal(dim(fgt), c(5, 3))
})

# 2. Single Element Data
test_that("FGTExperiment handles single sample/feature data", {
  skip_if_no_FGTExperiment()
  
  # Single sample
  single_sample <- matrix(10, nrow = 5, ncol = 1)
  rownames(single_sample) <- paste0("Feature", 1:5)
  colnames(single_sample) <- "Sample1"
  
  fgt_single <- FGTExperiment(
    assays = list(counts = single_sample),
    experimentType = "amplicon"
  )
  
  expect_equal(ncol(fgt_single), 1)
  expect_equal(nrow(fgt_single), 5)
  
  # Single feature
  single_feature <- matrix(10, nrow = 1, ncol = 5)
  rownames(single_feature) <- "Feature1"
  colnames(single_feature) <- paste0("Sample", 1:5)
  
  fgt_single_feature <- FGTExperiment(
    assays = list(counts = single_feature),
    experimentType = "amplicon"
  )
  
  expect_equal(nrow(fgt_single_feature), 1)
  expect_equal(ncol(fgt_single_feature), 5)
  
  # Single cell (1x1 matrix)
  single_cell <- matrix(42, nrow = 1, ncol = 1)
  rownames(single_cell) <- "Feature1"
  colnames(single_cell) <- "Sample1"
  
  fgt_single_cell <- FGTExperiment(
    assays = list(counts = single_cell),
    experimentType = "amplicon"
  )
  
  expect_equal(dim(fgt_single_cell), c(1, 1))
  expect_equal(assay(experimentData(fgt_single_cell))[1,1], 42)
})

# 3. Missing Data (NA values)
test_that("FGTExperiment handles NA values appropriately", {
  skip_if_no_FGTExperiment()
  
  # NA in counts (should error or convert to 0)
  na_counts <- matrix(c(10, NA, 20, 30, NA, 40), nrow = 3, ncol = 2)
  rownames(na_counts) <- paste0("Feature", 1:3)
  colnames(na_counts) <- paste0("Sample", 1:2)
  
  expect_warning(
    fgt_na <- FGTExperiment(
      assays = list(counts = na_counts),
      experimentType = "amplicon"
    ),
    "NA values in count matrix"
  )
  
  # NA in metadata (should be allowed)
  metadata_with_na <- DataFrame(
    Group = c("A", NA, "B"),
    Value = c(1, 2, NA)
  )
  rownames(metadata_with_na) <- paste0("Sample", 1:3)
  
  clean_counts <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(clean_counts) <- paste0("Feature", 1:3)
  colnames(clean_counts) <- paste0("Sample", 1:3)
  
  fgt_na_metadata <- FGTExperiment(
    assays = list(counts = clean_counts),
    colData = metadata_with_na,
    experimentType = "amplicon"
  )
  
  expect_true(any(is.na(colData(experimentData(fgt_na_metadata)))))
  
  # NA in taxonomy (should be allowed)
  taxonomy_with_na <- DataFrame(
    Genus = c("Lactobacillus", NA, "Gardnerella"),
    Species = c(NA, "vaginalis", NA)
  )
  rownames(taxonomy_with_na) <- paste0("Feature", 1:3)
  
  fgt_na_taxonomy <- FGTExperiment(
    assays = list(counts = clean_counts),
    rowData = taxonomy_with_na,
    experimentType = "amplicon"
  )
  
  expect_true(any(is.na(rowData(experimentData(fgt_na_taxonomy)))))
})

# 4. Invalid Data Types
test_that("FGTExperiment rejects invalid data types", {
  skip_if_no_FGTExperiment()
  
  # Character counts (should error)
  char_counts <- matrix(c("10", "20", "30", "40"), nrow = 2, ncol = 2)
  rownames(char_counts) <- paste0("Feature", 1:2)
  colnames(char_counts) <- paste0("Sample", 1:2)
  
  expect_error(
    FGTExperiment(
      assays = list(counts = char_counts),
      experimentType = "amplicon"
    ),
    "must be numeric"
  )
  
  # Negative counts (should warn or error)
  negative_counts <- matrix(c(10, -5, 20, -10), nrow = 2, ncol = 2)
  rownames(negative_counts) <- paste0("Feature", 1:2)
  colnames(negative_counts) <- paste0("Sample", 1:2)
  
  expect_warning(
    fgt_negative <- FGTExperiment(
      assays = list(counts = negative_counts),
      experimentType = "amplicon"
    ),
    "Negative values in count matrix"
  )
  
  # Non-integer counts (should be allowed but rounded)
  decimal_counts <- matrix(c(10.5, 20.3, 30.7, 40.9), nrow = 2, ncol = 2)
  rownames(decimal_counts) <- paste0("Feature", 1:2)
  colnames(decimal_counts) <- paste0("Sample", 1:2)
  
  fgt_decimal <- FGTExperiment(
    assays = list(counts = decimal_counts),
    experimentType = "amplicon"
  )
  
  # Check if values were handled appropriately
  result_counts <- assay(experimentData(fgt_decimal))
  expect_true(all(result_counts == round(result_counts)))
})

# 5. Mismatched Dimensions
test_that("FGTExperiment handles dimension mismatches", {
  skip_if_no_FGTExperiment()
  
  # Create base count matrix
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- paste0("Feature", 1:3)
  colnames(counts) <- paste0("Sample", 1:4)
  
  # Mismatched column metadata
  wrong_col_metadata <- DataFrame(
    Group = c("A", "B"),  # Only 2 rows, but 4 samples
    Value = c(1, 2)
  )
  rownames(wrong_col_metadata) <- paste0("Sample", 1:2)
  
  expect_error(
    FGTExperiment(
      assays = list(counts = counts),
      colData = wrong_col_metadata,
      experimentType = "amplicon"
    ),
    "colData.*does not match.*samples"
  )
  
  # Mismatched row metadata
  wrong_row_metadata <- DataFrame(
    Genus = c("Lactobacillus", "Gardnerella"),  # Only 2 rows, but 3 features
    Species = c("crispatus", "vaginalis")
  )
  rownames(wrong_row_metadata) <- paste0("Feature", 1:2)
  
  expect_error(
    FGTExperiment(
      assays = list(counts = counts),
      rowData = wrong_row_metadata,
      experimentType = "amplicon"
    ),
    "rowData.*does not match.*features"
  )
})

# 6. Very Large and Sparse Matrices
test_that("FGTExperiment handles large sparse matrices efficiently", {
  skip_if_no_FGTExperiment()
  
  # Generate large sparse matrix
  large_sparse <- generate_test_fixture("sparse", rows = 10000, cols = 100, sparsity = 0.95)
  
  # Time the creation
  create_time <- system.time({
    fgt_large <- FGTExperiment(
      assays = list(counts = large_sparse$counts),
      rowData = DataFrame(large_sparse$taxonomy),
      colData = DataFrame(large_sparse$metadata),
      experimentType = "amplicon"
    )
  })
  
  # Should create in reasonable time (< 5 seconds)
  expect_lt(create_time["elapsed"], 5)
  
  # Check object validity
  expect_equal(dim(fgt_large), c(10000, 100))
  
  # Test operations on large object
  transform_time <- system.time({
    fgt_rel <- transformAbundance(fgt_large, type = "relative")
  })
  
  # Transformation should be reasonably fast (< 10 seconds)
  expect_lt(transform_time["elapsed"], 10)
})

# 7. Unicode and Special Characters
test_that("FGTExperiment handles special characters in names", {
  skip_if_no_FGTExperiment()
  
  # Create count matrix with special characters
  special_counts <- matrix(1:12, nrow = 3, ncol = 4)
  
  # Unicode and special characters in row names
  rownames(special_counts) <- c(
    "Lactobacillus_crispatus_αβγ",
    "Gardnerella_vaginalis_<>",
    "Feature with spaces & symbols!@#"
  )
  
  # Special characters in column names
  colnames(special_counts) <- c(
    "Sample-1",
    "Sample_2",
    "Sample.3",
    "Sample 4 (special)"
  )
  
  # Create with special characters
  fgt_special <- FGTExperiment(
    assays = list(counts = special_counts),
    experimentType = "amplicon"
  )
  
  # Verify names are preserved
  expect_equal(rownames(fgt_special), rownames(special_counts))
  expect_equal(colnames(fgt_special), colnames(special_counts))
  
  # Test with metadata containing special characters
  special_metadata <- DataFrame(
    Group = c("Control", "Treatment (high dose)", "Treatment (low dose)", "Placebo"),
    Description = c("Normal μ", "Test β", "Test γ", "Control δ")
  )
  rownames(special_metadata) <- colnames(special_counts)
  
  fgt_special_meta <- FGTExperiment(
    assays = list(counts = special_counts),
    colData = special_metadata,
    experimentType = "amplicon"
  )
  
  # Verify metadata is preserved
  expect_equal(colData(experimentData(fgt_special_meta))$Description[1], "Normal μ")
})

# 8. Extreme Values
test_that("FGTExperiment handles extreme values", {
  skip_if_no_FGTExperiment()
  
  # Very large counts
  large_value_counts <- matrix(c(1e9, 1e10, 1e11, 1e12), nrow = 2, ncol = 2)
  rownames(large_value_counts) <- paste0("Feature", 1:2)
  colnames(large_value_counts) <- paste0("Sample", 1:2)
  
  fgt_large_values <- FGTExperiment(
    assays = list(counts = large_value_counts),
    experimentType = "amplicon"
  )
  
  expect_equal(max(assay(experimentData(fgt_large_values))), 1e12)
  
  # Mix of very large and very small values
  mixed_extreme <- matrix(c(0, 1e-10, 1e10, 0.0001), nrow = 2, ncol = 2)
  rownames(mixed_extreme) <- paste0("Feature", 1:2)
  colnames(mixed_extreme) <- paste0("Sample", 1:2)
  
  fgt_mixed <- FGTExperiment(
    assays = list(counts = mixed_extreme),
    experimentType = "amplicon"
  )
  
  # Test transformations with extreme values
  expect_no_error({
    rel_abund <- transformAbundance(fgt_mixed, type = "relative")
    log_trans <- transformAbundance(fgt_mixed, type = "log")
  })
})

# 9. Memory Stress Test
test_that("FGTExperiment memory usage is reasonable", {
  skip_if_no_FGTExperiment()
  skip_on_cran()  # Skip on CRAN due to resource constraints
  
  # Generate progressively larger datasets
  sizes <- c(1000, 5000, 10000)
  memory_usage <- numeric(length(sizes))
  
  for (i in seq_along(sizes)) {
    n_features <- sizes[i]
    n_samples <- 50
    
    # Generate test data
    test_data <- generate_test_fixture("large", rows = n_features, cols = n_samples)
    
    # Measure memory before
    gc()
    mem_before <- gc()[2, 2]  # Current memory usage in MB
    
    # Create object
    fgt_test <- FGTExperiment(
      assays = list(counts = test_data$counts),
      rowData = DataFrame(test_data$taxonomy),
      colData = DataFrame(test_data$metadata),
      experimentType = "amplicon"
    )
    
    # Force garbage collection and measure memory after
    gc()
    mem_after <- gc()[2, 2]
    
    memory_usage[i] <- mem_after - mem_before
    
    # Clean up
    rm(fgt_test, test_data)
    gc()
  }
  
  # Memory usage should scale linearly with data size
  # Check that doubling the features roughly doubles the memory
  ratio1 <- memory_usage[2] / memory_usage[1]
  ratio2 <- memory_usage[3] / memory_usage[2]
  
  # Allow some variance (between 1.5x and 3x)
  expect_gte(ratio1, 1.5)
  expect_lte(ratio1, 3.0)
  expect_gte(ratio2, 1.5)
  expect_lte(ratio2, 3.0)
})

# 10. File System Edge Cases
test_that("Import functions handle file system edge cases", {
  skip_if_no_FGTExperiment()
  
  # Test with non-existent file
  expect_error(
    import_speciateit("/path/that/does/not/exist"),
    "does not exist|cannot be found"
  )
  
  # Test with empty directory
  temp_dir <- tempdir()
  empty_dir <- file.path(temp_dir, "empty_test_dir")
  dir.create(empty_dir, showWarnings = FALSE)
  
  expect_error(
    import_speciateit(empty_dir),
    "No valid files found|Empty directory"
  )
  
  # Test with file instead of directory
  temp_file <- tempfile(fileext = ".txt")
  writeLines("Not a valid data file", temp_file)
  
  expect_error(
    import_speciateit(temp_file),
    "Expected directory|Not a directory"
  )
  
  # Clean up
  unlink(empty_dir, recursive = TRUE)
  unlink(temp_file)
})

# 11. Platform-Specific Edge Cases
test_that("Platform imports handle malformed data", {
  skip_if_no_FGTExperiment()
  
  # Create temporary directory for test files
  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "platform_tests")
  dir.create(test_dir, showWarnings = FALSE)
  
  # Test SpeciateIT with malformed count file
  speciateit_dir <- file.path(test_dir, "speciateit")
  dir.create(speciateit_dir, showWarnings = FALSE)
  
  # Write malformed count file
  malformed_counts <- "This is not a valid count matrix\nJust some text"
  writeLines(malformed_counts, file.path(speciateit_dir, "counts.txt"))
  
  expect_error(
    import_speciateit(speciateit_dir),
    "Invalid count matrix|Failed to parse"
  )
  
  # Test VALENCIA with missing CST assignments
  valencia_dir <- file.path(test_dir, "valencia")
  dir.create(valencia_dir, showWarnings = FALSE)
  
  # Write count file without CST info
  valid_counts <- matrix(1:12, nrow = 3, ncol = 4)
  write.table(valid_counts, file.path(valencia_dir, "counts.txt"), 
              sep = "\t", quote = FALSE)
  
  expect_warning(
    import_valencia(valencia_dir),
    "CST assignments not found|Missing CST"
  )
  
  # Clean up
  unlink(test_dir, recursive = TRUE)
})

# 12. Concurrent Access
test_that("Package handles concurrent operations safely", {
  skip_if_no_FGTExperiment()
  skip_on_cran()  # Skip on CRAN
  skip_if_not_installed("parallel")
  
  # Generate test data
  test_data <- generate_test_fixture("standard")
  
  # Create base FGTExperiment
  fgt_base <- FGTExperiment(
    assays = list(counts = test_data$counts),
    rowData = DataFrame(test_data$taxonomy),
    colData = DataFrame(test_data$metadata),
    experimentType = "amplicon"
  )
  
  # Test concurrent transformations
  library(parallel)
  cl <- makeCluster(2)
  
  # Export necessary objects to cluster
  clusterExport(cl, c("fgt_base", "transformAbundance"))
  clusterEvalQ(cl, library(microFGT))
  
  # Run transformations in parallel
  results <- parLapply(cl, c("relative", "log"), function(type) {
    transformAbundance(fgt_base, type = type)
  })
  
  stopCluster(cl)
  
  # Both transformations should succeed
  expect_length(results, 2)
  expect_s4_class(results[[1]], "FGTExperiment")
  expect_s4_class(results[[2]], "FGTExperiment")
})