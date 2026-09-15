# Test for the direct implementation of FGTExperiment class

test_that("FGTExperiment class exists", {
  source(file.path("../../R/AllGenerics.R"))
  source(file.path("../../R/FGTExperiment-complete.R"))
  expect_true(isClass("FGTExperiment"))
})

test_that("FGTExperiment constructor works", {
  source(file.path("../../R/AllGenerics.R"))
  source(file.path("../../R/FGTExperiment-complete.R"))
  
  # Create test data
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  
  # Create FGTExperiment
  fgt <- FGTExperiment(
    assays = list(counts = counts),
    experimentType = "amplicon",
    fgtMetadata = S4Vectors::SimpleList(note = "test object")
  )
  
  # Check basic properties
  expect_s4_class(fgt, "FGTExperiment")
  expect_equal(dim(fgt), c(4, 5))
  expect_equal(experimentType(fgt), "amplicon")
  expect_equal(fgtMetadata(fgt)$note, "test object")
  
  # Check accessors
  expect_equal(rownames(fgt), paste0("Feature", 1:4))
  expect_equal(colnames(fgt), paste0("Sample", 1:5))
  expect_equal(assayNames(fgt), "counts")
})

test_that("transformAbundance works", {
  source(file.path("../../R/AllGenerics.R"))
  source(file.path("../../R/FGTExperiment-complete.R"))
  
  # Create test data
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  
  # Create FGTExperiment
  fgt <- FGTExperiment(
    assays = list(counts = counts)
  )
  
  # Test transformations
  rel_fgt <- transformAbundance(fgt, type = "relative")
  log_fgt <- transformAbundance(fgt, type = "log")
  clr_fgt <- transformAbundance(fgt, type = "clr")
  pres_fgt <- transformAbundance(fgt, type = "presence")
  
  # Check that new assays were created
  expect_true("relative_counts" %in% assayNames(rel_fgt))
  expect_true("log_counts" %in% assayNames(log_fgt))
  expect_true("clr_counts" %in% assayNames(clr_fgt))
  expect_true("presence_counts" %in% assayNames(pres_fgt))
  
  # Check transformation results
  rel_counts <- assay(rel_fgt, "relative_counts")
  expect_equal(colSums(rel_counts), rep(1, 5))
})