test_that("transformAbundance creates expected output type", {
  obj <- create_test_object(rows=3, cols=3)
  
  # Test relative abundance transformation
  rel_obj <- transformAbundance(obj, type="relative")
  
  # Check object class
  expect_s4_class(rel_obj, "FGTExperiment")
  
  # Check assay name
  expect_true("relative_counts" %in% assayNames(rel_obj))
})

test_that("relative abundance transformation calculates correctly", {
  # Create a simple count matrix
  counts <- matrix(c(1, 2, 3, 4, 5, 6), nrow=2, ncol=3)
  rownames(counts) <- c("A", "B")
  colnames(counts) <- c("S1", "S2", "S3")
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  # Calculate expected relative abundances
  expected_rel <- matrix(
    c(1/3, 2/3,   # S1: 1/(1+2) and 2/(1+2)
      3/7, 4/7,   # S2: 3/(3+4) and 4/(3+4)
      5/11, 6/11),# S3: 5/(5+6) and 6/(5+6)
    nrow=2, ncol=3
  )
  rownames(expected_rel) <- c("A", "B")
  colnames(expected_rel) <- c("S1", "S2", "S3")
  
  # Transform and check
  rel_obj <- transformAbundance(obj, type="relative")
  
  # Compare matrices
  expect_equal(assay(rel_obj, "relative_counts"), expected_rel)
  
  # Check column sums
  expect_equal(colSums(assay(rel_obj, "relative_counts")), c(S1=1, S2=1, S3=1))
})

test_that("log transformation calculates correctly", {
  # Create a simple count matrix
  counts <- matrix(c(0, 1, 2, 3), nrow=2, ncol=2)
  rownames(counts) <- c("A", "B")
  colnames(counts) <- c("S1", "S2")
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  # Expected: log(x + 1)
  expected_log <- log1p(counts)
  
  # Transform and check
  log_obj <- transformAbundance(obj, type="log")
  
  # Compare matrices
  expect_equal(assay(log_obj, "log_counts"), expected_log)
})

test_that("clr transformation calculates correctly", {
  # Create a simple count matrix with no zeros
  counts <- matrix(c(1, 2, 3, 4), nrow=2, ncol=2)
  rownames(counts) <- c("A", "B")
  colnames(counts) <- c("S1", "S2")
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  # Expected: log(x) - mean(log(x)) for each column
  s1_log <- log(c(1, 2))
  s1_expected <- s1_log - mean(s1_log)
  
  s2_log <- log(c(3, 4))
  s2_expected <- s2_log - mean(s2_log)
  
  expected_clr <- cbind(s1_expected, s2_expected)
  colnames(expected_clr) <- c("S1", "S2")
  rownames(expected_clr) <- c("A", "B")
  
  # Transform with pseudocount=0 since we have no zeros
  clr_obj <- transformAbundance(obj, type="clr", pseudocount=0)
  
  # Compare matrices
  expect_equal(assay(clr_obj, "clr_counts"), expected_clr)
})

test_that("presence transformation calculates correctly", {
  # Create a matrix with zeros and non-zeros
  counts <- matrix(c(0, 1, 2, 0), nrow=2, ncol=2)
  rownames(counts) <- c("A", "B")
  colnames(counts) <- c("S1", "S2")
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  # Expected: TRUE for non-zero, FALSE for zero
  expected_presence <- counts > 0
  
  # Transform
  presence_obj <- transformAbundance(obj, type="presence")
  
  # Compare matrices
  expect_equal(assay(presence_obj, "presence_counts"), expected_presence)
})

test_that("transformAbundance handles edge cases", {
  # Create a matrix with a column of all zeros
  counts <- matrix(c(0, 0, 1, 2), nrow=2, ncol=2)
  rownames(counts) <- c("A", "B")
  colnames(counts) <- c("S1", "S2")
  
  obj <- FGTExperiment(assays=list(counts=counts))
  
  # For relative abundance, a column of all zeros should result in NAs
  rel_obj <- transformAbundance(obj, type="relative")
  expect_true(all(is.na(assay(rel_obj, "relative_counts")[,1])))
  expect_false(any(is.na(assay(rel_obj, "relative_counts")[,2])))
  
  # Log transformation should handle zeros by adding pseudocount
  log_obj <- transformAbundance(obj, type="log")
  expect_equal(assay(log_obj, "log_counts")[1,1], log1p(0))
})

test_that("transformAbundance validates inputs", {
  obj <- create_test_object(rows=3, cols=3)
  
  # Test invalid type
  expect_error(transformAbundance(obj, type="invalid"), "type must be one of")
  
  # Test non-existent assay
  expect_error(transformAbundance(obj, assay_name="nonexistent"), "not found")
})

test_that("transformAbundance works with SummarizedExperiment input", {
  # Create a SummarizedExperiment
  counts <- matrix(1:4, nrow=2, ncol=2)
  rownames(counts) <- c("A", "B")
  colnames(counts) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts))
  
  # Transform
  rel_se <- transformAbundance(se, type="relative")
  
  # Check object class
  expect_s4_class(rel_se, "SummarizedExperiment")
  
  # Check assay and values
  expect_true("relative_counts" %in% assayNames(rel_se))
  expect_equal(colSums(assay(rel_se, "relative_counts")), c(S1=1, S2=1))
})