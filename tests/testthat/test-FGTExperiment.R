test_that("FGTExperiment constructor works", {
  # Create a simple count matrix
  counts <- matrix(1:12, nrow = 4, ncol = 3)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:3)
  
  # Create a basic FGTExperiment
  fgt_exp <- FGTExperiment(assays = list(counts = counts))
  
  # Check that it's the right class
  expect_s4_class(fgt_exp, "FGTExperiment")
  expect_s4_class(fgt_exp, "TreeSummarizedExperiment")
  
  # Check dimensions
  expect_equal(dim(fgt_exp), c(4, 3))
  
  # Check default experiment type
  expect_equal(experimentType(fgt_exp), "amplicon")
  
  # Check assay names
  expect_equal(assayNames(fgt_exp), "counts")
})

test_that("transform_abundance works", {
  # Create a simple count matrix
  counts <- matrix(1:12, nrow = 4, ncol = 3)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:3)
  
  # Create a basic FGTExperiment
  fgt_exp <- FGTExperiment(assays = list(counts = counts))
  
  # Transform to relative abundance
  fgt_rel <- transform_abundance(fgt_exp, type = "relative")
  
  # Check that a new assay was added
  expect_true("relative" %in% assayNames(fgt_rel))
  
  # Check that relative abundances sum to 1 per sample
  rel_sums <- colSums(assays(fgt_rel)$relative)
  expect_equal(rel_sums, c(1, 1, 1), tolerance = 1e-10)
})

test_that("filter_taxa works", {
  # Create a count matrix with some zeros
  counts <- matrix(c(5, 0, 10, 0, 0, 5, 0, 10, 5, 0, 0, 5), nrow = 4, ncol = 3)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:3)
  
  # Create a basic FGTExperiment
  fgt_exp <- FGTExperiment(assays = list(counts = counts))
  
  # Filter with default parameters
  fgt_filtered <- filter_taxa(fgt_exp)
  
  # Only features 1 and 4 should remain (with prevalence > 0.1)
  expect_equal(dim(fgt_filtered), c(2, 3))
  expect_equal(rownames(fgt_filtered), c("Feature1", "Feature4"))
})

test_that("load_example_data works", {
  skip_if_not_installed("microFGT")
  
  # Try to load example data
  example_dir <- system.file("extdata", package = "microFGT")
  
  # Skip if no example data is available
  example_files <- list.files(example_dir, pattern = "microFGT_example_.*\\.rds")
  skip_if(length(example_files) == 0, "No example data available")
  
  # Try loading as list
  example_list <- load_example_data(as_fgt_experiment = FALSE)
  
  # Check structure
  expect_type(example_list, "list")
  expect_true(all(c("counts", "taxonomy", "metadata") %in% names(example_list)))
})

test_that("generate_fgt_example_data works", {
  # Set seed for reproducibility
  set.seed(123)
  
  # Generate small dataset
  example_data <- generate_fgt_example_data(
    n_samples = 5,
    n_features = 10,
    format = "list"
  )
  
  # Check structure
  expect_type(example_data, "list")
  expect_true(all(c("counts", "taxonomy", "metadata") %in% names(example_data)))
  
  # Check dimensions
  expect_equal(dim(example_data$counts), c(10, 5))
  expect_equal(dim(example_data$taxonomy), c(10, 7))
  expect_equal(dim(example_data$metadata), c(5, 5))
  
  # Generate as FGTExperiment
  fgt_example <- generate_fgt_example_data(
    n_samples = 5,
    n_features = 10,
    format = "FGTExperiment"
  )
  
  # Check class
  expect_s4_class(fgt_example, "FGTExperiment")
})