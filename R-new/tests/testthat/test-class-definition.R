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
  
  # Create invalid object with direct slot access (unsafe, only for testing)
  invalid_obj <- valid_obj
  invalid_obj@experimentType <- "invalid_type"
  
  # Test validation fails for invalid object
  expect_error(validObject(invalid_obj))
  
  # Test that fgtMetadata must be SimpleList
  invalid_obj2 <- valid_obj
  # Use unsafe direct slot manipulation for testing purposes only
  invalid_obj2@fgtMetadata <- list(a = 1)
  expect_error(validObject(invalid_obj2))
})