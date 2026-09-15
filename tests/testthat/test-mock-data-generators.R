# Tests for mock data generators

# The functions should be available through the package
library(microFGT)

# SpeciateIT Tests
test_that("generate_mock_speciateit creates valid data", {
  # Test with default parameters
  mock_data <- generate_mock_speciateit(n_sequences = 100)
  
  # Check structure
  expect_s3_class(mock_data, "data.frame")
  expect_equal(nrow(mock_data), 100)
  expect_equal(ncol(mock_data), 4)
  expect_equal(
    colnames(mock_data), 
    c("Sequence ID", "Classification", "posterior probability", "number of Decisions")
  )
  
  # Check data types
  expect_type(mock_data[["Sequence ID"]], "character")
  expect_type(mock_data[["Classification"]], "character")
  expect_type(mock_data[["posterior probability"]], "double")
  expect_type(mock_data[["number of Decisions"]], "integer")
  
  # Check value ranges
  expect_true(all(mock_data[["posterior probability"]] >= 0.3))
  expect_true(all(mock_data[["posterior probability"]] <= 1.0))
  expect_true(all(mock_data[["number of Decisions"]] >= 5))
  expect_true(all(mock_data[["number of Decisions"]] <= 20))
})

test_that("generate_mock_speciateit respects custom parameters", {
  # Custom species distribution
  species_dist <- c(
    "Species1" = 0.7,
    "Species2" = 0.3
  )
  
  mock_data <- generate_mock_speciateit(
    n_sequences = 50,
    seq_id_prefix = "Custom",
    species_distribution = species_dist,
    confidence_range = c(0.5, 0.9),
    decisions_range = c(10, 15),
    include_unclassified = 0,
    seed = 123
  )
  
  # Check structure
  expect_equal(nrow(mock_data), 50)
  
  # Check sequence ID prefix
  expect_true(all(grepl("^Custom_", mock_data[["Sequence ID"]])))
  
  # Check species
  unique_species <- unique(mock_data[["Classification"]])
  expect_true(all(unique_species %in% c("Species1", "Species2")))
  
  # Check confidence range
  expect_true(all(mock_data[["posterior probability"]] >= 0.5))
  expect_true(all(mock_data[["posterior probability"]] <= 0.9))
  
  # Check decisions range
  expect_true(all(mock_data[["number of Decisions"]] >= 10))
  expect_true(all(mock_data[["number of Decisions"]] <= 15))
})

test_that("generate_mock_speciateit handles unclassified sequences", {
  # Set high unclassified percentage
  mock_data <- generate_mock_speciateit(
    n_sequences = 100,
    include_unclassified = 0.2,
    seed = 456
  )
  
  # Check unclassified count (allow for small random variation)
  unclassified_count <- sum(mock_data[["Classification"]] == "Unclassified")
  expect_true(unclassified_count >= 15 && unclassified_count <= 25)
  
  # Check unclassified confidence scores are lower
  unclassified_indices <- which(mock_data[["Classification"]] == "Unclassified")
  unclassified_confidence <- mock_data[["posterior probability"]][unclassified_indices]
  expect_true(all(unclassified_confidence >= 0.3 & unclassified_confidence <= 0.6))
})

test_that("convert_speciateit_to_fgt correctly formats data", {
  # Generate mock data
  mock_data <- generate_mock_speciateit(n_sequences = 20, seed = 789)
  
  # Convert to FGT format
  fgt_data <- convert_speciateit_to_fgt(mock_data)
  
  # Check structure
  expect_type(fgt_data, "list")
  expect_true("assignments" %in% names(fgt_data))
  expect_true("version" %in% names(fgt_data))
  expect_equal(fgt_data$version, "mock")
  
  # Check assignments
  expect_s3_class(fgt_data$assignments, "data.frame")
  expect_equal(nrow(fgt_data$assignments), 20)
  expect_true(all(c("Feature", "Species", "Confidence") %in% colnames(fgt_data$assignments)))
  
  # Check data correspondence
  expect_equal(fgt_data$assignments$Feature, mock_data[["Sequence ID"]])
  expect_equal(fgt_data$assignments$Species, mock_data[["Classification"]])
  expect_equal(fgt_data$assignments$Confidence, mock_data[["posterior probability"]])
})

test_that("write_mock_speciateit writes file correctly", {
  # Create temporary file
  temp_file <- tempfile(fileext = ".txt")
  
  # Write mock data
  result <- write_mock_speciateit(temp_file, n_sequences = 10, seed = 101)
  
  # Check file exists and is non-empty
  expect_true(file.exists(temp_file))
  expect_gt(file.size(temp_file), 0)
  
  # Read the file back and check structure
  file_data <- read.delim(temp_file, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(file_data), 10)
  expect_equal(
    colnames(file_data), 
    c("Sequence ID", "Classification", "posterior probability", "number of Decisions")
  )
  
  # Clean up
  unlink(temp_file)
})

# VIRGO Tests
test_that("generate_mock_virgo creates valid data", {
  # Test with minimal parameters for speed
  mock_data <- generate_mock_virgo(n_genes = 50, n_samples = 3, seed = 123)
  
  # Check structure
  expect_type(mock_data, "list")
  expect_true(all(c("counts", "genes", "metadata", "test_out") %in% names(mock_data)))
  
  # Check counts matrix
  expect_true(is.matrix(mock_data$counts))
  expect_equal(nrow(mock_data$counts), 50)
  expect_equal(ncol(mock_data$counts), 3)
  
  # Check genes data frame
  expect_s3_class(mock_data$genes, "data.frame")
  expect_equal(nrow(mock_data$genes), 50)
  expect_true(all(c("gene_id", "length") %in% colnames(mock_data$genes)))
  
  # Check metadata data frame
  expect_s3_class(mock_data$metadata, "data.frame")
  expect_equal(nrow(mock_data$metadata), 3)
  expect_true(all(c("Sample", "CST", "Total_reads") %in% colnames(mock_data$metadata)))
  
  # Check test_out data frame
  expect_s3_class(mock_data$test_out, "data.frame")
  expect_true(all(c("gene_id", "read_count", "gene_length", "Sample") %in% colnames(mock_data$test_out)))
})

test_that("generate_mock_virgo respects custom parameters", {
  # Custom CST distribution
  cst_dist <- c("CST I" = 0.7, "CST III" = 0.3)
  
  mock_data <- generate_mock_virgo(
    n_genes = 40,
    n_samples = 5,
    sample_id_prefix = "VirgoTest",
    include_zeros = 0.5,
    simulation_type = "realistic",
    cst_distribution = cst_dist,
    seed = 456
  )
  
  # Check structure
  expect_equal(ncol(mock_data$counts), 5)
  expect_equal(nrow(mock_data$counts), 40)
  
  # Check sample ID prefix
  expect_true(all(grepl("^VirgoTest_", colnames(mock_data$counts))))
  
  # Check CST distribution (should only contain our custom CSTs)
  unique_csts <- unique(mock_data$metadata$CST)
  expect_true(all(unique_csts %in% c("CST I", "CST III")))
  
  # Check gene IDs format (VIRGO uses "V" followed by 7 digits)
  expect_true(all(grepl("^V\\d{7}$", mock_data$genes$gene_id)))
})

test_that("convert_virgo_to_fgt correctly formats data", {
  # Generate mock data
  mock_data <- generate_mock_virgo(n_genes = 30, n_samples = 2, seed = 789)
  
  # Convert to FGT format
  fgt_data <- convert_virgo_to_fgt(mock_data)
  
  # Check structure
  expect_type(fgt_data, "list")
  expect_true("genes" %in% names(fgt_data))
  expect_true("version" %in% names(fgt_data))
  expect_equal(fgt_data$version, "mock")
  
  # Check genes matrix
  expect_true(is.matrix(fgt_data$genes))
  expect_equal(dim(fgt_data$genes), c(30, 2))
})

test_that("write_mock_virgo writes files correctly", {
  # Create temporary directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Write mock data
  result <- write_mock_virgo(temp_dir, "test_virgo", n_genes = 20, n_samples = 2, seed = 101)
  
  # Check files exist
  expect_true(file.exists(file.path(temp_dir, "test_virgo_test.out")))
  expect_true(file.exists(file.path(temp_dir, "test_virgo_genes.txt")))
  expect_true(file.exists(file.path(temp_dir, "test_virgo_metadata.txt")))
  expect_true(file.exists(file.path(temp_dir, "test_virgo_counts.tsv")))
  
  # Read files back and check structure
  test_out <- read.delim(file.path(temp_dir, "test_virgo_test.out"), stringsAsFactors = FALSE)
  genes <- read.delim(file.path(temp_dir, "test_virgo_genes.txt"), stringsAsFactors = FALSE)
  
  expect_true(all(c("gene_id", "read_count", "gene_length", "Sample") %in% colnames(test_out)))
  expect_true(all(c("gene_id", "length") %in% colnames(genes)))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

# VALENCIA Tests
test_that("generate_mock_valencia creates valid data", {
  # Test with default parameters
  mock_data <- generate_mock_valencia(n_samples = 5, seed = 123)
  
  # Check structure
  expect_type(mock_data, "list")
  expect_true("cst" %in% names(mock_data))
  
  # Check cst data frame
  expect_s3_class(mock_data$cst, "data.frame")
  expect_equal(nrow(mock_data$cst), 5)
  expect_true(all(c("Sample", "CST") %in% colnames(mock_data$cst)))
  
  # Check CST values are valid
  valid_csts <- c("I", "II", "III", "IV-A", "IV-B", "IV-C", "V")
  expect_true(all(mock_data$cst$CST %in% valid_csts))
})

test_that("generate_mock_valencia includes optional components", {
  # Test with scores and abundance data
  mock_data <- generate_mock_valencia(
    n_samples = 3,
    include_scores = TRUE,
    add_taxonomic_data = TRUE,
    seed = 456
  )
  
  # Check scores matrix
  expect_true("scores" %in% names(mock_data))
  expect_true(is.matrix(mock_data$scores))
  expect_equal(nrow(mock_data$scores), 3)
  expect_equal(ncol(mock_data$scores), 7)  # 7 CSTs
  
  # Check that scores are probabilities
  expect_true(all(mock_data$scores >= 0 & mock_data$scores <= 1))
  # Check that scores sum to approximately 1 for each sample
  expect_true(all(abs(rowSums(mock_data$scores) - 1) < 1e-10))
  
  # Check abundance data
  expect_true("abundance" %in% names(mock_data))
  expect_s3_class(mock_data$abundance, "data.frame")
  expect_true("Sample" %in% colnames(mock_data$abundance))
  expect_true("read_count" %in% colnames(mock_data$abundance))
  
  # Check that taxonomic columns exist (at least the common species)
  common_species <- c("Lactobacillus_crispatus", "Lactobacillus_iners", "Gardnerella_vaginalis")
  expect_true(all(common_species %in% colnames(mock_data$abundance)))
})

test_that("convert_valencia_to_fgt correctly formats data", {
  # Generate mock data
  mock_data <- generate_mock_valencia(n_samples = 4, include_scores = TRUE, seed = 789)
  
  # Convert to FGT format
  fgt_data <- convert_valencia_to_fgt(mock_data)
  
  # Check structure
  expect_type(fgt_data, "list")
  expect_true("cst" %in% names(fgt_data))
  expect_true("scores" %in% names(fgt_data))
  expect_true("version" %in% names(fgt_data))
  expect_equal(fgt_data$version, "mock")
  
  # Check cst data frame
  expect_s3_class(fgt_data$cst, "data.frame")
  expect_equal(nrow(fgt_data$cst), 4)
  expect_true(all(c("Sample", "CST") %in% colnames(fgt_data$cst)))
  
  # Check scores matrix
  expect_true(is.matrix(fgt_data$scores))
  expect_equal(nrow(fgt_data$scores), 4)
})

test_that("write_mock_valencia writes files correctly", {
  # Create temporary directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Write mock data
  result <- write_mock_valencia(
    temp_dir, 
    "test_valencia", 
    n_samples = 3, 
    include_scores = TRUE,
    add_taxonomic_data = TRUE,
    seed = 101
  )
  
  # Check files exist
  expect_true(file.exists(file.path(temp_dir, "test_valencia_cst.csv")))
  expect_true(file.exists(file.path(temp_dir, "test_valencia_scores.csv")))
  expect_true(file.exists(file.path(temp_dir, "test_valencia_abundance.csv")))
  
  # Read files back and check structure
  cst_data <- read.csv(file.path(temp_dir, "test_valencia_cst.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(cst_data), 3)
  expect_true(all(c("Sample", "CST") %in% colnames(cst_data)))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

# Combined Dataset Tests
test_that("generate_mock_fgt_dataset creates coordinated data", {
  # Test with minimal parameters for speed
  mock_data <- generate_mock_fgt_dataset(
    n_samples = 3,
    n_sequences = 30,
    n_genes = 40,
    seed = 123
  )
  
  # Check structure
  expect_type(mock_data, "list")
  expect_true(all(c("speciateit", "virgo", "valencia", "metadata") %in% names(mock_data)))
  
  # Check metadata
  expect_equal(nrow(mock_data$metadata), 3)
  expect_true(all(c("Sample", "CST") %in% colnames(mock_data$metadata)))
  
  # Check coordination - CST assignments should match across tools
  valencia_csts <- mock_data$valencia$raw$cst$CST
  metadata_csts <- mock_data$metadata$CST
  expect_equal(valencia_csts, metadata_csts)
  
  # Check FGT-compatible formats are present
  expect_true("fgt_compatible" %in% names(mock_data$speciateit))
  expect_true("fgt_compatible" %in% names(mock_data$virgo))
  expect_true("fgt_compatible" %in% names(mock_data$valencia))
})

test_that("generate_mock_fgt_dataset can create files", {
  # Create temporary directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Generate dataset with files
  mock_data <- generate_mock_fgt_dataset(
    n_samples = 2,
    n_sequences = 20,
    n_genes = 30,
    create_files = TRUE,
    output_dir = temp_dir,
    seed = 456
  )
  
  # Check that file paths are returned
  expect_true("file_paths" %in% names(mock_data))
  
  # Check that directories were created
  expect_true(dir.exists(file.path(temp_dir, "speciateit")))
  expect_true(dir.exists(file.path(temp_dir, "virgo")))
  expect_true(dir.exists(file.path(temp_dir, "valencia")))
  
  # Check that key files exist
  expect_true(file.exists(file.path(temp_dir, "speciateit", "MC_order7_results.txt")))
  expect_true(file.exists(file.path(temp_dir, "virgo", "mock_virgo_test.out")))
  expect_true(file.exists(file.path(temp_dir, "valencia", "mock_valencia_cst.csv")))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})