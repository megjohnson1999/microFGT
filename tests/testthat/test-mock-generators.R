# Tests for mock data generators

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