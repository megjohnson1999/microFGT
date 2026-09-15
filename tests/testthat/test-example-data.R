test_that("create_fgt_count_matrix works", {
  # Skip if not interactive or on CI/CD
  skip_if_not(interactive() || Sys.getenv("CI") == "true")
  
  # Generate a small count matrix
  count_matrix <- create_fgt_count_matrix(
    n_samples = 5,
    n_features = 20,
    community_types = c("CST-I", "CST-III", "CST-IV")
  )
  
  # Basic validation
  expect_is(count_matrix, "matrix")
  expect_equal(dim(count_matrix), c(20, 5))
  expect_true(all(count_matrix >= 0))
  expect_true(all(is.finite(count_matrix)))
  expect_true(all(count_matrix == floor(count_matrix)))  # Check they are integers
  
  # Check row and column names
  expect_equal(length(rownames(count_matrix)), 20)
  expect_equal(length(colnames(count_matrix)), 5)
})

test_that("create_fgt_taxonomy works", {
  # Skip if not interactive or on CI/CD
  skip_if_not(interactive() || Sys.getenv("CI") == "true")
  
  # Generate feature IDs
  feature_ids <- paste0("Feature", 1:20)
  
  # Generate taxonomy
  taxonomy <- create_fgt_taxonomy(
    feature_ids = feature_ids,
    community_types = c("CST-I", "CST-III")
  )
  
  # Basic validation
  expect_is(taxonomy, "data.frame")
  expect_equal(nrow(taxonomy), 20)
  expect_true(all(rownames(taxonomy) == feature_ids))
  
  # Check taxonomic columns
  expect_true(all(c("Kingdom", "Phylum", "Genus", "Species") %in% colnames(taxonomy)))
  
  # Check for specific taxa
  expect_true(any(grepl("Lactobacillus", taxonomy$Genus, fixed = TRUE)))
})

test_that("create_fgt_sample_metadata works", {
  # Skip if not interactive or on CI/CD
  skip_if_not(interactive() || Sys.getenv("CI") == "true")
  
  # Generate sample IDs
  sample_ids <- paste0("Sample", 1:10)
  
  # Generate metadata
  metadata <- create_fgt_sample_metadata(
    sample_ids = sample_ids,
    groups = c("Healthy", "BV"),
    group_proportions = c(0.7, 0.3)
  )
  
  # Basic validation
  expect_is(metadata, "data.frame")
  expect_equal(nrow(metadata), 10)
  expect_true(all(rownames(metadata) == sample_ids))
  
  # Check columns
  expect_true("condition" %in% colnames(metadata))
  expect_true("pH" %in% colnames(metadata))
  expect_true("Nugent_Score" %in% colnames(metadata))
  
  # Check distribution of groups
  expect_equal(sum(metadata$condition == "Healthy"), 7)
  expect_equal(sum(metadata$condition == "BV"), 3)
  
  # Check clinical parameters
  expect_true(all(metadata$pH >= 3.8 & metadata$pH <= 6.0))
  expect_true(all(metadata$Nugent_Score >= 0 & metadata$Nugent_Score <= 10))
})

test_that("generate_fgt_example_data works", {
  # Skip if not interactive or on CI/CD
  skip_if_not(interactive() || Sys.getenv("CI") == "true")
  
  # Generate a small example dataset
  example_data <- generate_fgt_example_data(
    n_samples = 5,
    n_features = 15,
    sample_groups = c("Healthy", "BV"),
    group_proportions = c(0.6, 0.4),
    community_types = c("CST-I", "CST-IV"),
    format = "list"
  )
  
  # Basic validation
  expect_is(example_data, "list")
  expect_true(all(c("counts", "taxonomy", "metadata") %in% names(example_data)))
  
  # Check counts
  expect_is(example_data$counts, "matrix")
  expect_equal(dim(example_data$counts), c(15, 5))
  
  # Check taxonomy
  expect_is(example_data$taxonomy, "data.frame")
  expect_equal(nrow(example_data$taxonomy), 15)
  
  # Check metadata
  expect_is(example_data$metadata, "data.frame")
  expect_equal(nrow(example_data$metadata), 5)
  
  # Check if metadata contains the group variable and clinical parameters
  expect_true("condition" %in% colnames(example_data$metadata))
  expect_true("pH" %in% colnames(example_data$metadata))
  
  # Check consistency
  expect_equal(rownames(example_data$counts), rownames(example_data$taxonomy))
  expect_equal(colnames(example_data$counts), rownames(example_data$metadata))
})