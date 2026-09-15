test_that("calculate_diversity correctly calculates Shannon diversity", {
  # Create a simple matrix
  counts <- matrix(c(
    10, 20, 30, 0, 0,   # Sample 1: more even distribution
    90, 5, 3, 1, 1,     # Sample 2: dominated by one taxon
    20, 20, 20, 20, 20  # Sample 3: perfectly even distribution
  ), nrow = 3, byrow = TRUE)
  colnames(counts) <- paste0("Feature", 1:5)
  rownames(counts) <- paste0("Sample", 1:3)
  
  # Calculate Shannon diversity
  shannon <- calculate_diversity(counts, method = "shannon")
  
  # Check results
  expect_equal(length(shannon), 3)
  expect_equal(names(shannon), paste0("Sample", 1:3))
  
  # Shannon diversity should be highest for sample 3 (even distribution)
  # and lowest for sample 2 (dominated by one taxon)
  expect_true(shannon[3] > shannon[1])
  expect_true(shannon[1] > shannon[2])
  
  # Calculate without normalization
  shannon_unnorm <- calculate_diversity(counts, method = "shannon", normalize = FALSE)
  
  # Results should be different from normalized version
  expect_false(all(shannon == shannon_unnorm))
  
  # Test on SummarizedExperiment
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(counts)))
    shannon_se <- calculate_diversity(se, method = "shannon")
    expect_equal(shannon, shannon_se)
  }
})

test_that("calculate_diversity correctly calculates Simpson diversity", {
  # Create a simple matrix
  counts <- matrix(c(
    10, 20, 30, 0, 0,   # Sample 1: more even distribution
    90, 5, 3, 1, 1,     # Sample 2: dominated by one taxon
    20, 20, 20, 20, 20  # Sample 3: perfectly even distribution
  ), nrow = 3, byrow = TRUE)
  colnames(counts) <- paste0("Feature", 1:5)
  rownames(counts) <- paste0("Sample", 1:3)
  
  # Calculate Simpson diversity
  simpson <- calculate_diversity(counts, method = "simpson")
  
  # Check results
  expect_equal(length(simpson), 3)
  expect_equal(names(simpson), paste0("Sample", 1:3))
  
  # Simpson diversity should be highest for sample 3 (even distribution)
  # and lowest for sample 2 (dominated by one taxon)
  expect_true(simpson[3] > simpson[1])
  expect_true(simpson[1] > simpson[2])
  
  # Check value for perfectly even distribution (should be close to 1-1/5 = 0.8)
  expect_equal(simpson[3], 0.8, tolerance = 0.01)
  
  # Check value for highly dominated sample (should be close to 0)
  expect_true(simpson[2] < 0.2)
})

test_that("calculate_diversity correctly calculates richness and evenness", {
  # Create a simple matrix with some zeros
  counts <- matrix(c(
    10, 20, 30, 0, 0,   # Sample 1: 3 features
    90, 5, 3, 1, 1,     # Sample 2: 5 features
    20, 20, 20, 20, 0   # Sample 3: 4 features
  ), nrow = 3, byrow = TRUE)
  colnames(counts) <- paste0("Feature", 1:5)
  rownames(counts) <- paste0("Sample", 1:3)
  
  # Calculate richness
  richness <- calculate_diversity(counts, method = "richness")
  
  # Check results
  expect_equal(richness, c(Sample1 = 3, Sample2 = 5, Sample3 = 4))
  
  # Calculate evenness
  evenness <- calculate_diversity(counts, method = "evenness")
  
  # Check results
  expect_equal(length(evenness), 3)
  
  # Evenness should be highest for sample 3 (most even non-zero distribution)
  # and lowest for sample 2 (dominated by one taxon)
  expect_true(evenness[3] > evenness[1])
  expect_true(evenness[1] > evenness[2])
  
  # Perfect evenness would be 1
  expect_true(all(evenness <= 1))
  expect_true(all(evenness > 0))
})

test_that("calculate_beta_diversity correctly calculates Bray-Curtis distances", {
  # Skip test if required packages aren't available
  skip_if_not_installed("SummarizedExperiment")
  
  # Create a simple matrix
  counts <- matrix(c(
    10, 20, 30, 0, 0,   # Sample 1: similar to sample 3
    90, 5, 3, 1, 1,     # Sample 2: very different
    15, 25, 35, 5, 0    # Sample 3: similar to sample 1
  ), nrow = 3, byrow = TRUE)
  colnames(counts) <- paste0("Feature", 1:5)
  rownames(counts) <- paste0("Sample", 1:3)
  
  # Calculate Bray-Curtis distances
  bray <- calculate_beta_diversity(counts, method = "bray")
  
  # Check results
  expect_equal(dim(bray), c(3, 3))
  expect_equal(rownames(bray), paste0("Sample", 1:3))
  expect_equal(colnames(bray), paste0("Sample", 1:3))
  
  # Distances should be symmetric
  expect_equal(bray[1, 2], bray[2, 1])
  expect_equal(bray[1, 3], bray[3, 1])
  expect_equal(bray[2, 3], bray[3, 2])
  
  # Diagonal should be zeros
  expect_equal(diag(bray), c(0, 0, 0))
  
  # Sample 1 and 3 should be more similar than either is to sample 2
  expect_true(bray[1, 3] < bray[1, 2])
  expect_true(bray[1, 3] < bray[2, 3])
  
  # Test on SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(counts)))
  bray_se <- calculate_beta_diversity(se, method = "bray")
  expect_equal(bray, bray_se)
})

test_that("calculate_beta_diversity correctly calculates Jaccard distances", {
  # Create a simple presence/absence matrix
  counts <- matrix(c(
    1, 1, 1, 0, 0,  # Sample 1: 3 features, 2 shared with sample 3
    1, 0, 0, 1, 1,  # Sample 2: 3 features, 1 shared with sample 1
    1, 1, 0, 1, 0   # Sample 3: 3 features, 2 shared with sample 1, 1 with sample 2
  ), nrow = 3, byrow = TRUE)
  colnames(counts) <- paste0("Feature", 1:5)
  rownames(counts) <- paste0("Sample", 1:3)
  
  # Calculate Jaccard distances
  jaccard <- calculate_beta_diversity(counts, method = "jaccard")
  
  # Check results
  expect_equal(dim(jaccard), c(3, 3))
  
  # For binary data, Jaccard distance = 1 - (A ∩ B) / (A ∪ B)
  # Sample 1 and 2: shared = 1, union = 5, distance = 1 - 1/5 = 0.8
  expect_equal(jaccard[1, 2], 0.8)
  
  # Sample 1 and 3: shared = 2, union = 4, distance = 1 - 2/4 = 0.5
  expect_equal(jaccard[1, 3], 0.5)
  
  # Sample 2 and 3: shared = 1, union = 5, distance = 1 - 1/5 = 0.8
  expect_equal(jaccard[2, 3], 0.8)
})

test_that("helper distance functions work correctly", {
  # Create a simple matrix
  counts <- matrix(c(
    0.5, 0.3, 0.2, 0.0, 0.0,  # Sample 1
    0.9, 0.1, 0.0, 0.0, 0.0,  # Sample 2
    0.3, 0.3, 0.3, 0.1, 0.0   # Sample 3
  ), nrow = 3, byrow = TRUE)
  colnames(counts) <- paste0("Feature", 1:5)
  rownames(counts) <- paste0("Sample", 1:3)
  
  # Calculate Bray-Curtis distances manually
  bc_dist <- bray_curtis_distance(counts)
  
  # Check dimensions and symmetry
  expect_equal(dim(bc_dist), c(3, 3))
  expect_equal(bc_dist[1, 2], bc_dist[2, 1])
  
  # Calculate Jaccard distances manually
  jc_dist <- jaccard_distance(counts)
  
  # Check dimensions and symmetry
  expect_equal(dim(jc_dist), c(3, 3))
  expect_equal(jc_dist[1, 2], jc_dist[2, 1])
  
  # Presence/absence should affect Jaccard but not necessarily Bray-Curtis the same way
  # Add a small test case to verify this
  binary_counts <- matrix(c(
    1, 1, 0,   # Sample A: 2 features
    1, 0, 1    # Sample B: 2 features, 1 shared
  ), nrow = 2, byrow = TRUE)
  
  uneven_counts <- matrix(c(
    10, 90, 0,    # Sample A: 2 features
    90, 0, 10     # Sample B: 2 features, 1 shared
  ), nrow = 2, byrow = TRUE)
  
  # Jaccard should be the same for both (1 - 1/3 = 0.667)
  expect_equal(
    jaccard_distance(binary_counts)[1, 2],
    jaccard_distance(uneven_counts)[1, 2]
  )
  
  # Bray-Curtis should be different
  expect_false(
    bray_curtis_distance(binary_counts)[1, 2] == 
    bray_curtis_distance(uneven_counts)[1, 2]
  )
})