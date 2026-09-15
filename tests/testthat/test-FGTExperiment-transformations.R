# Comprehensive tests for FGTExperiment transformation methods

test_that("transformAbundance generic exists", {
  expect_true(isGeneric("transformAbundance"))
})

test_that("relative abundance transformation works correctly", {
  # Create test object
  fgt <- create_test_fgt(rows = 5, cols = 3)
  
  # Apply relative abundance transformation
  fgt_rel <- transformAbundance(fgt, type = "relative")
  
  # Check that new assay was created
  expect_true("relative_counts" %in% assayNames(fgt_rel))
  
  # Check that columns sum to 1
  rel_counts <- assay(fgt_rel, "relative_counts")
  expect_equal(colSums(rel_counts), rep(1, ncol(rel_counts)))
  
  # Check that the original assay is unchanged
  expect_equal(assay(fgt_rel, "counts"), assay(fgt, "counts"))
  
  # Verify transformation using helper function
  expect_true(verify_transformation(fgt, fgt_rel, type = "relative"))
})

test_that("log transformation works correctly", {
  # Create test object
  fgt <- create_test_fgt(rows = 5, cols = 3)
  
  # Apply log transformation
  fgt_log <- transformAbundance(fgt, type = "log")
  
  # Check that new assay was created
  expect_true("log_counts" %in% assayNames(fgt_log))
  
  # Check that values match log1p of original
  log_counts <- assay(fgt_log, "log_counts")
  original_counts <- assay(fgt, "counts")
  expected_log <- log1p(original_counts)
  expect_equal(log_counts, expected_log)
  
  # Verify transformation using helper function
  expect_true(verify_transformation(fgt, fgt_log, type = "log"))
})

test_that("clr transformation works correctly", {
  # Create test object
  fgt <- create_test_fgt(rows = 5, cols = 3)
  
  # Apply CLR transformation
  fgt_clr <- transformAbundance(fgt, type = "clr")
  
  # Check that new assay was created
  expect_true("clr_counts" %in% assayNames(fgt_clr))
  
  # Check that column means are approximately zero
  clr_counts <- assay(fgt_clr, "clr_counts")
  col_means <- colMeans(clr_counts)
  expect_true(all(abs(col_means) < 1e-10))
  
  # Verify transformation using helper function
  expect_true(verify_transformation(fgt, fgt_clr, type = "clr"))
})

test_that("presence transformation works correctly", {
  # Create test object
  fgt <- create_test_fgt(rows = 5, cols = 3)
  
  # Ensure some zeros in the data for better testing
  counts <- assay(fgt, "counts")
  counts[sample(seq_len(length(counts)), size = 5)] <- 0
  assays(fgt@experimentData)$counts <- counts
  
  # Apply presence/absence transformation
  fgt_pres <- transformAbundance(fgt, type = "presence")
  
  # Check that new assay was created
  expect_true("presence_counts" %in% assayNames(fgt_pres))
  
  # Check that values are binary (0 or 1)
  pres_counts <- assay(fgt_pres, "presence_counts")
  expect_true(all(pres_counts %in% c(0, 1)))
  
  # Check that presence matches where original is greater than zero
  original_counts <- assay(fgt, "counts")
  expected_pres <- as.numeric(original_counts > 0)
  dim(expected_pres) <- dim(original_counts)
  rownames(expected_pres) <- rownames(original_counts)
  colnames(expected_pres) <- colnames(original_counts)
  expect_equal(pres_counts, expected_pres)
  
  # Verify transformation using helper function
  expect_true(verify_transformation(fgt, fgt_pres, type = "presence"))
})

test_that("transformAbundance handles edge cases", {
  # Create test object
  fgt <- create_test_fgt(rows = 5, cols = 3)
  
  # Test with invalid transformation type
  expect_error(transformAbundance(fgt, type = "invalid"))
  
  # Test with non-existent assay
  expect_error(transformAbundance(fgt, assay_name = "nonexistent"))
  
  # Test with empty columns (all zeros)
  counts <- assay(fgt, "counts")
  counts[, 1] <- 0  # Make first column all zeros
  assays(fgt@experimentData)$counts <- counts
  
  # Relative abundance should handle zero columns gracefully
  fgt_rel <- transformAbundance(fgt, type = "relative")
  rel_counts <- assay(fgt_rel, "relative_counts")
  expect_true(all(is.na(rel_counts[, 1])))
  
  # CLR should handle zeros with pseudocount
  fgt_clr <- transformAbundance(fgt, type = "clr", pseudocount = 1)
  clr_counts <- assay(fgt_clr, "clr_counts")
  expect_false(any(is.na(clr_counts)))
  
  # Test different pseudocount values
  fgt_clr2 <- transformAbundance(fgt, type = "clr", pseudocount = 0.5)
  clr_counts2 <- assay(fgt_clr2, "clr_counts")
  expect_false(identical(clr_counts, clr_counts2))
  
  # Presence/absence should work with all zeros
  fgt_pres <- transformAbundance(fgt, type = "presence")
  pres_counts <- assay(fgt_pres, "presence_counts")
  expect_equal(pres_counts[, 1], rep(0, nrow(pres_counts)))
})

test_that("multiple transformations can be applied sequentially", {
  # Create test object
  fgt <- create_test_fgt(rows = 5, cols = 3)
  
  # Apply relative abundance transformation
  fgt_rel <- transformAbundance(fgt, type = "relative")
  
  # Apply CLR to the relative abundance assay
  fgt_rel_clr <- transformAbundance(fgt_rel, type = "clr", assay_name = "relative_counts")
  
  # Check that both assays exist
  expect_true("relative_counts" %in% assayNames(fgt_rel_clr))
  expect_true("clr_relative_counts" %in% assayNames(fgt_rel_clr))
  
  # Original assay should still be unchanged
  expect_equal(assay(fgt_rel_clr, "counts"), assay(fgt, "counts"))
  
  # CLR should have been applied to relative abundance
  clr_rel_counts <- assay(fgt_rel_clr, "clr_relative_counts")
  col_means <- colMeans(clr_rel_counts)
  expect_true(all(abs(col_means) < 1e-10))
})

test_that("transformAbundance works with TreeSummarizedExperiment", {
  # Create a TreeSummarizedExperiment
  counts <- matrix(rpois(50, lambda = 20), nrow = 10, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:10)
  colnames(counts) <- paste0("Sample", 1:5)
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Apply transformation
  tse_rel <- transformAbundance(tse, type = "relative")
  
  # Check that new assay was created
  expect_true("relative_counts" %in% assayNames(tse_rel))
  
  # Check that columns sum to 1
  rel_counts <- assay(tse_rel, "relative_counts")
  expect_equal(colSums(rel_counts), rep(1, ncol(rel_counts)))
})

test_that("transformAbundance works with SummarizedExperiment", {
  # Create a SummarizedExperiment
  counts <- matrix(rpois(50, lambda = 20), nrow = 10, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:10)
  colnames(counts) <- paste0("Sample", 1:5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Apply transformation
  se_rel <- transformAbundance(se, type = "relative")
  
  # Check that new assay was created
  expect_true("relative_counts" %in% assayNames(se_rel))
  
  # Check that columns sum to 1
  rel_counts <- assay(se_rel, "relative_counts")
  expect_equal(colSums(rel_counts), rep(1, ncol(rel_counts)))
})