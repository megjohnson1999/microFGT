# Test directly that FGTExperiment validation works as expected

test_that("FGTExperiment_direct creates a valid object", {
  # Skip if function is not available
  skip_if_not_installed("microFGT")
  skip_if_not(exists("FGTExperiment_direct", mode = "function"))
  
  # Create a test object
  fgt <- FGTExperiment_direct()
  
  # Check class
  expect_s4_class(fgt, "FGTExperiment")
  
  # Check slots
  expect_s4_class(fgt@experimentData, "TreeSummarizedExperiment")
  expect_type(fgt@experimentType, "character")
  expect_s4_class(fgt@fgtMetadata, "SimpleList")
  
  # Check accessors
  expect_equal(experimentType(fgt), "amplicon")
  expect_s4_class(fgtMetadata(fgt), "SimpleList")
  
  # Check dimensions
  expect_length(dim(fgt), 2)
  expect_true(all(dim(fgt) > 0))
  
  # Check that assay exists
  expect_true("counts" %in% SummarizedExperiment::assayNames(fgt))
  
  # Check row and column names
  expect_type(rownames(fgt), "character")
  expect_type(colnames(fgt), "character")
})

test_that("transformAbundance_direct creates correct transformations", {
  # Skip if functions not available
  skip_if_not_installed("microFGT")
  skip_if_not(exists("FGTExperiment_direct", mode = "function"))
  skip_if_not(exists("transformAbundance_direct", mode = "function"))
  
  # Create test object
  fgt <- FGTExperiment_direct()
  
  # Test relative abundance transformation
  fgt_rel <- transformAbundance_direct(fgt, type = "relative")
  expect_true("relative_counts" %in% SummarizedExperiment::assayNames(fgt_rel))
  rel_counts <- SummarizedExperiment::assay(fgt_rel, "relative_counts")
  expect_equal(colSums(rel_counts), rep(1, ncol(rel_counts)))
  
  # Test log transformation
  fgt_log <- transformAbundance_direct(fgt, type = "log")
  expect_true("log_counts" %in% SummarizedExperiment::assayNames(fgt_log))
  
  # Test CLR transformation
  fgt_clr <- transformAbundance_direct(fgt, type = "clr")
  expect_true("clr_counts" %in% SummarizedExperiment::assayNames(fgt_clr))
  
  # Test presence transformation
  fgt_pres <- transformAbundance_direct(fgt, type = "presence")
  expect_true("presence_counts" %in% SummarizedExperiment::assayNames(fgt_pres))
})

test_that("FGTExperiment validation prevents invalid objects", {
  # Skip if function is not available
  skip_if_not_installed("microFGT")
  skip_if_not(exists("FGTExperiment_direct", mode = "function"))
  
  # Create a valid object
  fgt <- FGTExperiment_direct()
  
  # Test invalid experimentType
  expect_error({
    fgt@experimentType <- "invalid"
    validObject(fgt)
  }, "must be one of")
  
  # Reset to valid value
  fgt@experimentType <- "amplicon"
  expect_true(validObject(fgt))
  
  # Test setting via accessor
  expect_error(
    experimentType(fgt) <- "invalid",
    "must be one of"
  )
  
  # Test valid experimentType values
  experimentType(fgt) <- "metagenomic"
  expect_equal(experimentType(fgt), "metagenomic")
  
  experimentType(fgt) <- "integrated"
  expect_equal(experimentType(fgt), "integrated")
})