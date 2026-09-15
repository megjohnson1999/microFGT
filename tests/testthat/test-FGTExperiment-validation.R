# Validation test for FGTExperiment class

test_that("FGTExperiment class exists", {
  expect_true(isClass("FGTExperiment"))
})

test_that("FGTExperiment constructor works", {
  # Create test data
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  
  # Test basic construction
  fgt <- FGTExperiment(
    assays = list(counts = counts),
    experimentType = "amplicon",
    fgtMetadata = S4Vectors::SimpleList(note = "test")
  )
  
  # Check if object is created successfully
  expect_s4_class(fgt, "FGTExperiment")
  expect_equal(dim(fgt), c(4, 5))
  expect_equal(experimentType(fgt), "amplicon")
  expect_equal(fgtMetadata(fgt)$note, "test")
  
  # Test accessor methods
  expect_equal(rownames(fgt), paste0("Feature", 1:4))
  expect_equal(colnames(fgt), paste0("Sample", 1:5))
  
  # Test creation from SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
  fgt2 <- FGTExperiment(se)
  expect_s4_class(fgt2, "FGTExperiment")
  expect_equal(dim(fgt2), c(4, 5))
})

test_that("transformAbundance method works", {
  # Skip if method not available
  if (!isGeneric("transformAbundance")) {
    skip("transformAbundance generic not defined")
  }

  # Create test data with non-zero values
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  
  # Create FGTExperiment
  fgt <- FGTExperiment(
    assays = list(counts = counts)
  )
  
  # Test relative transformation
  rel_fgt <- transformAbundance(fgt, type = "relative")
  expect_true("relative_counts" %in% SummarizedExperiment::assayNames(rel_fgt))
  
  # Check sum to 1
  rel_counts <- SummarizedExperiment::assay(rel_fgt, "relative_counts")
  expect_equal(colSums(rel_counts), rep(1, 5))
})
