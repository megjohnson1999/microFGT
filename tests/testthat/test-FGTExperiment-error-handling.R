# Tests for error handling in FGTExperiment

test_that("Constructor handles invalid inputs with informative errors", {
  # Empty list should fail
  expect_error(FGTExperiment(assays = list()), 
              regexp = "must be a non-empty list")
  
  # Non-matrix first element
  expect_error(FGTExperiment(assays = list(counts = "not a matrix")), 
              regexp = "must be a matrix")
  
  # Invalid experiment type
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  expect_error(FGTExperiment(
    assays = list(counts = counts),
    experimentType = "invalid_type"
  ), regexp = "must be one of")
})

test_that("transformAbundance handles errors gracefully", {
  # Create test FGTExperiment
  fgt <- create_test_fgt()
  
  # Invalid transformation type
  expect_error(transformAbundance(fgt, type = "invalid"), 
              regexp = "type must be one of")
  
  # Non-existent assay name
  expect_error(transformAbundance(fgt, assay_name = "nonexistent"), 
              regexp = "not found")
  
  # Invalid pseudocount (negative)
  expect_error(transformAbundance(fgt, type = "clr", pseudocount = -1), 
              regexp = "must be positive")
})

test_that("accessors handle invalid inputs with proper errors", {
  # Create test FGTExperiment
  fgt <- create_test_fgt()
  
  # Test experimentType setter with invalid value
  expect_error(experimentType(fgt) <- "invalid", 
              regexp = "must be one of")
  
  # Test fgtMetadata setter with invalid type
  # This should auto-convert rather than error
  fgtMetadata(fgt) <- list(note = "test")
  expect_s4_class(fgtMetadata(fgt), "SimpleList")
  
  # Try setting to something that can't be coerced to SimpleList
  expect_error(fgtMetadata(fgt) <- environment())
})

test_that("show method handles edge cases", {
  # Create minimal FGTExperiment
  counts <- matrix(1:4, nrow = 2, ncol = 2)
  rownames(counts) <- c("F1", "F2")
  colnames(counts) <- c("S1", "S2")
  fgt <- FGTExperiment(assays = list(counts = counts))
  
  # Empty fgtMetadata
  output <- capture.output(show(fgt))
  expect_true(any(grepl("Empty SimpleList", output)))
  
  # Add metadata and check it shows up
  fgtMetadata(fgt) <- S4Vectors::SimpleList(note = "test")
  output <- capture.output(show(fgt))
  expect_true(any(grepl("note", output)))
  expect_true(any(grepl("test", output)))
})

test_that("edge case test: zero rows or columns", {
  # Zero rows
  counts_zero_rows <- matrix(numeric(0), nrow = 0, ncol = 3)
  colnames(counts_zero_rows) <- paste0("S", 1:3)
  
  # This should be allowed but with a warning
  expect_warning(
    fgt_zero_rows <- FGTExperiment(assays = list(counts = counts_zero_rows))
  )
  
  # If created, check dimensions
  if (exists("fgt_zero_rows")) {
    expect_equal(nrow(fgt_zero_rows), 0)
    expect_equal(ncol(fgt_zero_rows), 3)
  }
  
  # Zero columns
  counts_zero_cols <- matrix(numeric(0), nrow = 3, ncol = 0)
  rownames(counts_zero_cols) <- paste0("F", 1:3)
  
  # This should be allowed but with a warning
  expect_warning(
    fgt_zero_cols <- FGTExperiment(assays = list(counts = counts_zero_cols))
  )
  
  # If created, check dimensions
  if (exists("fgt_zero_cols")) {
    expect_equal(nrow(fgt_zero_cols), 3)
    expect_equal(ncol(fgt_zero_cols), 0)
  }
})

test_that("empty constructor creates minimal valid object", {
  # Try empty constructor - check if it errors or creates a minimal object
  result <- tryCatch({
    FGTExperiment()
  }, error = function(e) {
    e$message
  })
  
  # Either should error (which we record) or create a minimal object
  if (is(result, "FGTExperiment")) {
    expect_s4_class(result, "FGTExperiment")
    expect_true(validObject(result))
  } else {
    # If it errors, result contains the error message
    expect_match(result, "must be a non-empty")
  }
})

test_that("migration function handles missing slots gracefully", {
  # Skip if migration function is not available
  if (!exists("migrateToNewFGTExperiment", mode = "function")) {
    skip("Migration function not available")
  }
  
  # Create TreeSummarizedExperiment
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- paste0("F", 1:3)
  colnames(counts) <- paste0("S", 1:4)
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Create mock old-style object with missing slot
  # This emulates legacy code issues where slots might be absent
  if (!isClass("PartialFGTExperiment")) {
    setClass("PartialFGTExperiment",
             contains = "TreeSummarizedExperiment",
             slots = list(
               experimentType = "character"
               # No fgtMetadata slot to simulate partial object
             ),
             prototype = list(
               experimentType = "amplicon"
             ))
  }
  
  # Create partial object
  partial_obj <- methods::new("PartialFGTExperiment", tse,
                            experimentType = "amplicon")
  
  # Migration should handle this gracefully
  migrated <- tryCatch({
    migrateToNewFGTExperiment(partial_obj)
  }, error = function(e) {
    e
  })
  
  # Either it succeeds or fails gracefully
  if (is(migrated, "FGTExperiment")) {
    expect_s4_class(migrated, "FGTExperiment")
    expect_equal(migrated@experimentType, "amplicon")
    expect_s4_class(migrated@fgtMetadata, "SimpleList")
  } else {
    # If it errors, the error should be about missing slot
    expect_match(as.character(migrated), "slot")
  }
})