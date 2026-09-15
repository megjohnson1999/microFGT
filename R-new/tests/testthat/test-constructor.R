test_that("FGTExperiment constructor works with minimal input", {
  counts <- matrix(1:12, nrow=3, ncol=4)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Sample", 1:4)
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  expect_s4_class(obj, "FGTExperiment")
  expect_equal(dim(obj), c(3, 4))
  expect_equal(experimentType(obj), "amplicon")
  expect_true(length(fgtMetadata(obj)) == 0)
})

test_that("FGTExperiment constructor handles missing row/colnames", {
  counts <- matrix(1:12, nrow=3, ncol=4)
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  expect_equal(rownames(obj), c("Feature1", "Feature2", "Feature3"))
  expect_equal(colnames(obj), c("Sample1", "Sample2", "Sample3", "Sample4"))
})

test_that("FGTExperiment constructor validates inputs", {
  expect_error(FGTExperiment(assays=list()), "'assays' must be a non-empty list")
  expect_error(FGTExperiment(assays=list(counts="not_a_matrix")), "must be a matrix or Matrix object")
})

test_that("FGTExperiment constructor accepts experimentType parameter", {
  counts <- matrix(1:12, nrow=3, ncol=4)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Sample", 1:4)
  
  # Test default
  obj1 <- FGTExperiment(assays=list(counts=counts))
  expect_equal(experimentType(obj1), "amplicon")
  
  # Test custom value
  obj2 <- FGTExperiment(assays=list(counts=counts), experimentType="metagenomic")
  expect_equal(experimentType(obj2), "metagenomic")
  
  # Test invalid value
  expect_error(FGTExperiment(assays=list(counts=counts), experimentType="invalid"))
})

test_that("FGTExperiment constructor accepts metadata", {
  counts <- matrix(1:12, nrow=3, ncol=4)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Sample", 1:4)
  
  metadata <- S4Vectors::SimpleList(
    source = "test",
    date = Sys.Date()
  )
  
  obj <- FGTExperiment(assays=list(counts=counts), fgtMetadata=metadata)
  
  expect_equal(fgtMetadata(obj), metadata)
  expect_equal(fgtMetadata(obj)$source, "test")
})