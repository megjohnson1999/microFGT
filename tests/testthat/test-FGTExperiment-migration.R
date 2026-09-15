# Tests for migrating between old and new FGTExperiment structures

# Set up a function to mock an old-style FGTExperiment object
create_mock_old_fgt <- function() {
  # Skip if TreeSummarizedExperiment is not available
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    skip("TreeSummarizedExperiment package not available")
  }
  
  # Create a mock class that simulates the old inheritance-based structure
  if (!isClass("MockOldFGTExperiment")) {
    setClass("MockOldFGTExperiment",
             contains = "TreeSummarizedExperiment",
             slots = list(
               experimentType = "character",
               fgtMetadata = "SimpleList"
             ),
             prototype = list(
               experimentType = "amplicon",
               fgtMetadata = S4Vectors::SimpleList()
             ))
  }
  
  # Create TreeSummarizedExperiment as the base
  counts <- matrix(rpois(50, lambda = 20), nrow = 10, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:10)
  colnames(counts) <- paste0("Sample", 1:5)
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Create the mock old-style object
  old_fgt <- methods::new("MockOldFGTExperiment", tse,
                         experimentType = "metagenomic",
                         fgtMetadata = S4Vectors::SimpleList(note = "old-style"))
  
  return(old_fgt)
}

test_that("migrateToNewFGTExperiment migrates old-style objects", {
  # Skip if migration function is not available
  if (!exists("migrateToNewFGTExperiment", mode = "function")) {
    skip("Migration function not available")
  }
  
  # Create a mock old-style object
  old_fgt <- create_mock_old_fgt()
  
  # Apply the migration function
  new_fgt <- migrateToNewFGTExperiment(old_fgt)
  
  # Check that the result is a new-style FGTExperiment
  expect_s4_class(new_fgt, "FGTExperiment")
  expect_true(methods::existsSlot(new_fgt, "experimentData"))
  expect_s4_class(new_fgt@experimentData, "TreeSummarizedExperiment")
  
  # Check that the metadata was preserved
  expect_equal(new_fgt@experimentType, "metagenomic")
  expect_equal(new_fgt@fgtMetadata$note, "old-style")
  
  # Check that the dimensions are the same
  expect_equal(dim(new_fgt), c(10, 5))
  
  # Check that the assay data was preserved
  expect_true("counts" %in% assayNames(new_fgt))
  
  # Additional test: try migrating something that's already new-style
  # This should return the object unchanged
  new_fgt2 <- migrateToNewFGTExperiment(new_fgt)
  expect_identical(new_fgt2, new_fgt)
})

test_that("migrateToNewFGTExperiment handles TreeSummarizedExperiment", {
  # Skip if migration function is not available
  if (!exists("migrateToNewFGTExperiment", mode = "function")) {
    skip("Migration function not available")
  }
  
  # Create a TreeSummarizedExperiment
  counts <- matrix(rpois(50, lambda = 20), nrow = 10, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:10)
  colnames(counts) <- paste0("Sample", 1:5)
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Apply the migration function
  fgt <- migrateToNewFGTExperiment(tse)
  
  # Check that the result is a new-style FGTExperiment
  expect_s4_class(fgt, "FGTExperiment")
  expect_s4_class(fgt@experimentData, "TreeSummarizedExperiment")
  
  # Check that the dimensions are the same
  expect_equal(dim(fgt), dim(tse))
  
  # Check that the assay data was preserved
  expect_equal(assayNames(fgt), assayNames(tse))
})

test_that("migrateToNewFGTExperiment handles SummarizedExperiment", {
  # Skip if migration function is not available
  if (!exists("migrateToNewFGTExperiment", mode = "function")) {
    skip("Migration function not available")
  }
  
  # Create a SummarizedExperiment
  counts <- matrix(rpois(50, lambda = 20), nrow = 10, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:10)
  colnames(counts) <- paste0("Sample", 1:5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Apply the migration function
  fgt <- migrateToNewFGTExperiment(se)
  
  # Check that the result is a new-style FGTExperiment
  expect_s4_class(fgt, "FGTExperiment")
  expect_s4_class(fgt@experimentData, "TreeSummarizedExperiment")
  
  # Check that the dimensions are the same
  expect_equal(dim(fgt), dim(se))
  
  # Check that the assay data was preserved
  expect_equal(assayNames(fgt), assayNames(se))
})

test_that("migrateToNewFGTExperiment handles errors gracefully", {
  # Skip if migration function is not available
  if (!exists("migrateToNewFGTExperiment", mode = "function")) {
    skip("Migration function not available")
  }
  
  # Try with an invalid object type
  expect_error(migrateToNewFGTExperiment("not a valid object"))
  expect_error(migrateToNewFGTExperiment(data.frame()))
  
  # Try with NULL
  expect_error(migrateToNewFGTExperiment(NULL))
})