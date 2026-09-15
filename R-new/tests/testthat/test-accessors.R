test_that("experimentType accessor methods work", {
  obj <- create_test_object(rows=3, cols=3)
  
  # Test getter
  expect_equal(experimentType(obj), "amplicon")
  
  # Test setter
  experimentType(obj) <- "metagenomic"
  expect_equal(experimentType(obj), "metagenomic")
  
  # Test validation in setter
  expect_error(experimentType(obj) <- "invalid_type", "must be one of")
})

test_that("fgtMetadata accessor methods work", {
  obj <- create_test_object(rows=3, cols=3)
  
  # Test getter (default metadata from create_test_object)
  expect_s4_class(fgtMetadata(obj), "SimpleList")
  expect_true("CreatedBy" %in% names(fgtMetadata(obj)))
  expect_equal(fgtMetadata(obj)$CreatedBy, "test_helper")
  
  # Test setter
  new_metadata <- S4Vectors::SimpleList(source = "test_accessors", version = "1.0")
  fgtMetadata(obj) <- new_metadata
  
  expect_equal(fgtMetadata(obj), new_metadata)
  expect_equal(fgtMetadata(obj)$source, "test_accessors")
  expect_equal(fgtMetadata(obj)$version, "1.0")
  
  # Test validation in setter
  expect_error(fgtMetadata(obj) <- list(a = 1), "must be a SimpleList")
})