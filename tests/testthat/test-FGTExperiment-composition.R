# Test for the compositional FGTExperiment class implementation

test_that("FGTExperiment can be constructed correctly", {
  # Create test data
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  
  # Create FGTExperiment
  fgt <- FGTExperiment(
    assays = list(counts = counts),
    experimentType = "amplicon",
    fgtMetadata = S4Vectors::SimpleList(notes = "Test object")
  )
  
  # Check class
  expect_s4_class(fgt, "FGTExperiment")
  
  # Check slots
  expect_s4_class(fgt@experimentData, "TreeSummarizedExperiment")
  expect_equal(fgt@experimentType, "amplicon")
  expect_s4_class(fgt@fgtMetadata, "SimpleList")
  expect_equal(fgt@fgtMetadata$notes, "Test object")
  
  # Check dimensions
  expect_equal(dim(fgt), c(4, 5))
  expect_equal(rownames(fgt), paste0("Feature", 1:4))
  expect_equal(colnames(fgt), paste0("Sample", 1:5))
  
  # Check accessor methods
  expect_equal(experimentType(fgt), "amplicon")
  expect_equal(fgtMetadata(fgt)$notes, "Test object")
})

test_that("FGTExperiment can be created from SummarizedExperiment", {
  # Create SummarizedExperiment
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Convert to FGTExperiment
  fgt <- FGTExperiment(se)
  
  # Check class
  expect_s4_class(fgt, "FGTExperiment")
  
  # Check dimensions
  expect_equal(dim(fgt), c(4, 5))
  expect_equal(rownames(fgt), paste0("Feature", 1:4))
  expect_equal(colnames(fgt), paste0("Sample", 1:5))
  
  # Check that assay data is preserved
  expect_equal(assay(fgt), counts)
})

test_that("transformAbundance works correctly", {
  # Create test data
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  
  # Create FGTExperiment
  fgt <- FGTExperiment(
    assays = list(counts = counts)
  )
  
  # Apply transformations
  fgt_rel <- transformAbundance(fgt, type = "relative")
  fgt_log <- transformAbundance(fgt, type = "log")
  fgt_clr <- transformAbundance(fgt, type = "clr")
  fgt_presence <- transformAbundance(fgt, type = "presence")
  
  # Check that transformations were applied correctly
  expect_true("relative_counts" %in% assayNames(fgt_rel))
  expect_true("log_counts" %in% assayNames(fgt_log))
  expect_true("clr_counts" %in% assayNames(fgt_clr))
  expect_true("presence_counts" %in% assayNames(fgt_presence))
  
  # Verify relative transformation
  rel_counts <- assay(fgt_rel, "relative_counts")
  expect_equal(colSums(rel_counts), rep(1, 5))
  
  # Verify presence transformation (should be all TRUE since counts are all > 0)
  pres_counts <- assay(fgt_presence, "presence_counts")
  expect_true(all(pres_counts))
})

test_that("Migration function works for SummarizedExperiment", {
  # Create SummarizedExperiment
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Convert using migration function
  fgt <- migrateToNewFGTExperiment(se)
  
  # Check class
  expect_s4_class(fgt, "FGTExperiment")
  expect_s4_class(fgt@experimentData, "TreeSummarizedExperiment")
  
  # Check dimensions
  expect_equal(dim(fgt), c(4, 5))
  
  # Check assay data
  expect_equal(assay(fgt), counts)
})