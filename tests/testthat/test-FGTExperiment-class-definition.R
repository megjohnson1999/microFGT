# Comprehensive tests for FGTExperiment class definition

test_that("FGTExperiment class exists and has correct structure", {
  # Check class exists
  expect_true(isClass("FGTExperiment"))
  
  # Check slot definitions
  slots <- getSlots("FGTExperiment")
  expect_true("experimentData" %in% names(slots))
  expect_true("experimentType" %in% names(slots))
  expect_true("fgtMetadata" %in% names(slots))
  
  # Check slot types
  expect_equal(slots["experimentData"], "TreeSummarizedExperiment")
  expect_equal(slots["experimentType"], "character")
  expect_equal(slots["fgtMetadata"], "SimpleList")
})

test_that("FGTExperiment constructor handles various input types", {
  # Basic constructor with count matrix
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  
  fgt1 <- FGTExperiment(assays = list(counts = counts))
  expect_s4_class(fgt1, "FGTExperiment")
  expect_s4_class(fgt1@experimentData, "TreeSummarizedExperiment")
  expect_equal(dim(fgt1), c(4, 5))
  
  # Constructor with SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
  fgt2 <- FGTExperiment(se)
  expect_s4_class(fgt2, "FGTExperiment")
  expect_equal(dim(fgt2), c(4, 5))
  
  # Constructor with TreeSummarizedExperiment
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(assays = list(counts = counts))
  fgt3 <- FGTExperiment(tse)
  expect_s4_class(fgt3, "FGTExperiment")
  expect_equal(dim(fgt3), c(4, 5))
  
  # Constructor with additional metadata
  fgt4 <- FGTExperiment(
    assays = list(counts = counts),
    experimentType = "metagenomic",
    fgtMetadata = S4Vectors::SimpleList(notes = "test data")
  )
  expect_s4_class(fgt4, "FGTExperiment")
  expect_equal(fgt4@experimentType, "metagenomic")
  expect_equal(fgt4@fgtMetadata$notes, "test data")
  
  # Check that constructor adds names if missing
  counts_no_names <- matrix(1:20, nrow = 4, ncol = 5)
  fgt5 <- FGTExperiment(assays = list(counts = counts_no_names))
  expect_false(is.null(rownames(fgt5)))
  expect_false(is.null(colnames(fgt5)))
})

test_that("FGTExperiment constructor handles edge cases and errors", {
  # Empty list should fail
  expect_error(FGTExperiment(assays = list()))
  
  # Non-matrix assay should fail
  expect_error(FGTExperiment(assays = list(counts = "not a matrix")))
  
  # Invalid experiment type should fail
  counts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:5)
  expect_error(FGTExperiment(
    assays = list(counts = counts),
    experimentType = "invalid"
  ))
  
  # Invalid fgtMetadata should be converted to SimpleList
  fgt <- FGTExperiment(
    assays = list(counts = counts),
    fgtMetadata = list(notes = "test data")
  )
  expect_s4_class(fgt@fgtMetadata, "SimpleList")
})

test_that("FGTExperiment validation works correctly", {
  # Create valid object
  fgt <- create_test_fgt()
  expect_true(validObject(fgt))
  
  # Test invalid experimentType
  fgt_invalid <- fgt
  fgt_invalid@experimentType <- "invalid"
  expect_error(validObject(fgt_invalid))
  
  # Test invalid fgtMetadata
  fgt_invalid2 <- fgt
  fgt_invalid2@fgtMetadata <- "not a SimpleList"
  expect_error(validObject(fgt_invalid2))
})

test_that("FGTExperiment accessors return correct values", {
  fgt <- create_test_fgt(rows = 5, cols = 3)
  
  # Basic accessors
  expect_equal(dim(fgt), c(5, 3))
  expect_equal(nrow(fgt), 5)
  expect_equal(ncol(fgt), 3)
  expect_equal(rownames(fgt), paste0("Feature", 1:5))
  expect_equal(colnames(fgt), paste0("Sample", 1:3))
  
  # Slot accessors
  expect_equal(experimentType(fgt), "amplicon")
  expect_s4_class(fgtMetadata(fgt), "SimpleList")
  
  # Check assay accessors
  expect_true("counts" %in% assayNames(fgt))
  expect_equal(dim(assay(fgt)), c(5, 3))
  expect_equal(dim(assays(fgt)$counts), c(5, 3))
  
  # Check rowData and colData
  expect_s4_class(rowData(fgt), "DataFrame")
  expect_s4_class(colData(fgt), "DataFrame")
  expect_equal(nrow(rowData(fgt)), 5)
  expect_equal(nrow(colData(fgt)), 3)
})

test_that("FGTExperiment setter methods work correctly", {
  fgt <- create_test_fgt()
  
  # Set experimentType
  experimentType(fgt) <- "metagenomic"
  expect_equal(experimentType(fgt), "metagenomic")
  
  # Try invalid experimentType
  expect_error(experimentType(fgt) <- "invalid")
  
  # Set fgtMetadata
  new_metadata <- S4Vectors::SimpleList(notes = "new notes", version = "2.0")
  fgtMetadata(fgt) <- new_metadata
  expect_equal(fgtMetadata(fgt)$notes, "new notes")
  expect_equal(fgtMetadata(fgt)$version, "2.0")
  
  # Set with list instead of SimpleList (should be converted)
  fgtMetadata(fgt) <- list(notes = "converted notes")
  expect_s4_class(fgtMetadata(fgt), "SimpleList")
  expect_equal(fgtMetadata(fgt)$notes, "converted notes")
})

test_that("show method displays the object correctly", {
  fgt <- create_test_fgt(rows = 3, cols = 2)
  
  # Capture the show method output
  output <- capture.output(show(fgt))
  
  # Check for key components in the output
  expect_true(any(grepl("FGTExperiment object with 3 features and 2 samples", output)))
  expect_true(any(grepl("experimentType: amplicon", output)))
  expect_true(any(grepl("fgtMetadata: SimpleList with", output)))
  expect_true(any(grepl("assays", output)))
  expect_true(any(grepl("rownames", output)))
  expect_true(any(grepl("colnames", output)))
})

test_that("coercion methods work correctly", {
  fgt <- create_test_fgt()
  
  # Test coercion to TreeSummarizedExperiment
  tse <- as(fgt, "TreeSummarizedExperiment")
  expect_s4_class(tse, "TreeSummarizedExperiment")
  expect_equal(dim(tse), dim(fgt))
  expect_true("counts" %in% assayNames(tse))
  
  # Test coercion to SummarizedExperiment
  se <- as(fgt, "SummarizedExperiment")
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), dim(fgt))
  expect_true("counts" %in% assayNames(se))
  
  # Test coercion from TreeSummarizedExperiment
  fgt2 <- as(tse, "FGTExperiment")
  expect_s4_class(fgt2, "FGTExperiment")
  expect_equal(dim(fgt2), dim(tse))
  
  # Test coercion from SummarizedExperiment
  fgt3 <- as(se, "FGTExperiment")
  expect_s4_class(fgt3, "FGTExperiment")
  expect_equal(dim(fgt3), dim(se))
})