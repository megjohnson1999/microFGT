# Tests for diversity functions with FGTExperiment class

library(testthat)
library(microFGT)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)

# Tests for calculate_diversity (alpha diversity)
test_that("calculate_diversity computes alpha diversity metrics correctly", {
  # Create test data with known diversity properties
  # Equal abundance - should have maximum diversity (Shannon = log(rows))
  equal_fgt <- create_test_fgt(rows = 10, cols = 5)
  counts <- assay(equal_fgt, "counts")
  counts[] <- 1  # Set all abundances equal
  assay(equal_fgt, "counts") <- counts
  
  # Very uneven abundance - should have low diversity
  uneven_fgt <- create_test_fgt(rows = 10, cols = 5)
  counts <- assay(uneven_fgt, "counts")
  counts[] <- 0
  for (j in 1:ncol(counts)) {
    # Make one dominant taxon per sample (90%)
    counts[j, j] <- 90
    # Distribute remaining 10% among other taxa
    others <- setdiff(1:nrow(counts), j)
    if (length(others) > 0) {
      counts[others, j] <- 10 / length(others)
    }
  }
  assay(uneven_fgt, "counts") <- counts
  
  # 1. Test Shannon diversity
  shannon_equal <- calculate_diversity(equal_fgt, method = "shannon")
  shannon_uneven <- calculate_diversity(uneven_fgt, method = "shannon")
  
  # Equal should have higher Shannon diversity
  expect_true(all(shannon_equal > shannon_uneven))
  
  # For equal abundance, Shannon should be log(rows)
  # But with small numerical differences due to normalization
  expect_equal(shannon_equal, rep(log(10), 5), tolerance = 1e-5)
  
  # 2. Test Simpson diversity
  simpson_equal <- calculate_diversity(equal_fgt, method = "simpson")
  simpson_uneven <- calculate_diversity(uneven_fgt, method = "simpson")
  
  # Equal should have higher Simpson diversity (closer to 1)
  expect_true(all(simpson_equal > simpson_uneven))
  
  # For equal abundance, Simpson should be 1 - 1/rows
  expect_equal(simpson_equal, rep(1 - 1/10, 5), tolerance = 1e-5)
  
  # 3. Test Richness
  # Create zero-inflated data
  sparse_fgt <- create_test_fgt(rows = 10, cols = 5)
  counts <- assay(sparse_fgt, "counts")
  # First sample has all taxa
  # Second sample is missing 1 taxon
  counts[1, 2] <- 0
  # Third sample is missing 2 taxa
  counts[1:2, 3] <- 0
  # Fourth sample is missing 5 taxa
  counts[1:5, 4] <- 0
  # Fifth sample has only 1 taxon
  counts[2:10, 5] <- 0
  assay(sparse_fgt, "counts") <- counts
  
  richness <- calculate_diversity(sparse_fgt, method = "richness")
  
  # Check richness values
  expect_equal(richness[1], 10)  # All 10 taxa
  expect_equal(richness[2], 9)   # Missing 1 taxon
  expect_equal(richness[3], 8)   # Missing 2 taxa
  expect_equal(richness[4], 5)   # Missing 5 taxa
  expect_equal(richness[5], 1)   # Only 1 taxon
  
  # 4. Test evenness (Pielou's evenness)
  evenness_equal <- calculate_diversity(equal_fgt, method = "evenness")
  evenness_uneven <- calculate_diversity(uneven_fgt, method = "evenness")
  
  # Equal should have maximum evenness (1.0)
  expect_true(all(evenness_equal > evenness_uneven))
  expect_equal(evenness_equal, rep(1.0, 5), tolerance = 1e-5)
  
  # 5. Test assay_name parameter
  # Create a transformed assay
  rel_fgt <- transformAbundance(equal_fgt, type = "relative")
  expect_error(calculate_diversity(rel_fgt, assay_name = "nonexistent"),
               "not found in the object")
  
  shannon_rel <- calculate_diversity(rel_fgt, method = "shannon", assay_name = "relative_counts")
  expect_equal(shannon_rel, shannon_equal)
  
  # 6. Test normalize parameter
  # Should get same result with normalize=FALSE on already normalized data
  shannon_no_norm <- calculate_diversity(rel_fgt, method = "shannon", 
                                        assay_name = "relative_counts",
                                        normalize = FALSE)
  expect_equal(shannon_no_norm, shannon_rel)
  
  # 7. Test error handling
  expect_error(calculate_diversity("not_a_matrix"),
               "must be a matrix or SummarizedExperiment object")
})

# Tests for calculate_beta_diversity
test_that("calculate_beta_diversity computes beta diversity metrics correctly", {
  # Create test data with known beta diversity properties
  # Create FGTExperiment with samples having different degrees of similarity
  fgt <- create_test_fgt(rows = 10, cols = 6)
  counts <- matrix(0, nrow = 10, ncol = 6)
  rownames(counts) <- rownames(assay(fgt))
  colnames(counts) <- colnames(assay(fgt))
  
  # Sample 1 and 2 are identical
  counts[1:5, 1] <- c(10, 20, 30, 40, 50)
  counts[1:5, 2] <- c(10, 20, 30, 40, 50)
  
  # Sample 3 is similar to 1 and 2
  counts[1:5, 3] <- c(15, 25, 35, 45, 55)
  
  # Sample 4 has different taxa from 1, 2, 3
  counts[6:10, 4] <- c(10, 20, 30, 40, 50)
  
  # Sample 5 shares half its taxa with 1/2/3 and half with 4
  counts[3:7, 5] <- c(30, 40, 50, 20, 10)
  
  # Sample 6 is completely different (all zeros)
  
  assay(fgt, "counts") <- counts
  
  # 1. Test Bray-Curtis dissimilarity
  bray <- calculate_beta_diversity(fgt, method = "bray")
  
  # Check matrix properties
  expect_true(is.matrix(bray))
  expect_equal(dim(bray), c(6, 6))
  expect_equal(rownames(bray), colnames(fgt))
  expect_equal(colnames(bray), colnames(fgt))
  
  # Check symmetry and diagonal
  expect_true(isSymmetric(bray))
  expect_equal(diag(bray), rep(0, 6))
  
  # Check specific distance relationships
  # Identical samples should have distance 0
  expect_equal(bray[1, 2], 0)
  
  # Similar samples should have small distance
  expect_lt(bray[1, 3], 0.2)
  
  # Different samples should have large distance
  expect_gt(bray[1, 4], 0.8)
  
  # Sample 5 should be intermediate between groups
  expect_lt(bray[5, 1], bray[4, 1])
  expect_lt(bray[5, 4], bray[1, 4])
  
  # Sample 6 (all zeros) should have maximal distance
  expect_equal(bray[6, 1:5], rep(1, 5))
  
  # 2. Test Jaccard distance
  jaccard <- calculate_beta_diversity(fgt, method = "jaccard")
  
  # Check matrix properties
  expect_true(is.matrix(jaccard))
  expect_equal(dim(jaccard), c(6, 6))
  
  # Check specific distance relationships similar to Bray-Curtis
  # Identical samples should have distance 0
  expect_equal(jaccard[1, 2], 0)
  
  # 3. Test UniFrac distance (requires tree)
  skip_if_not_installed("phyloseq")
  
  # Create FGTExperiment with tree
  tree_fgt <- create_test_fgt(rows = 10, cols = 6, include_tree = TRUE)
  assay(tree_fgt, "counts") <- counts
  
  unifrac <- calculate_beta_diversity(tree_fgt, method = "unifrac")
  
  # Check matrix properties
  expect_true(is.matrix(unifrac))
  expect_equal(dim(unifrac), c(6, 6))
  
  # Check symmetry and diagonal
  expect_true(isSymmetric(unifrac))
  expect_equal(diag(unifrac), rep(0, 6))
  
  # Identical samples should have distance 0
  expect_equal(unifrac[1, 2], 0)
  
  # 4. Test error handling
  # Missing assay
  expect_error(calculate_beta_diversity(fgt, assay_name = "nonexistent"),
               "not found in the object")
  
  # UniFrac without tree
  expect_error(calculate_beta_diversity(fgt, method = "unifrac"),
               "requires a phylogenetic tree")
  
  # Non-SummarizedExperiment object
  expect_error(calculate_beta_diversity("not_a_matrix"),
               "must be a matrix or SummarizedExperiment object")
})

# Mock ordination test for use with beta diversity metrics
test_that("beta diversity can be used for ordination", {
  skip_if_not_installed("vegan")
  
  # Create test data
  fgt <- create_test_fgt(rows = 20, cols = 10)
  
  # Calculate beta diversity matrix
  dist_matrix <- calculate_beta_diversity(fgt, method = "bray")
  
  # Perform PCoA/MDS ordination on the distance matrix
  pcoa <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  
  # Check ordination results
  expect_true(is.list(pcoa))
  expect_equal(dim(pcoa$points), c(10, 2))
  expect_true(length(pcoa$eig) >= 2)
  
  # Test with Jaccard distance
  dist_matrix <- calculate_beta_diversity(fgt, method = "jaccard")
  pcoa <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  expect_equal(dim(pcoa$points), c(10, 2))
})

# Test using diversity functions with transformed counts
test_that("diversity functions work correctly with transformed data", {
  # Create test data
  fgt <- create_test_fgt(rows = 15, cols = 8)
  
  # Create different transformations
  fgt_rel <- transformAbundance(fgt, type = "relative")
  fgt_log <- transformAbundance(fgt, type = "log")
  fgt_clr <- transformAbundance(fgt, type = "clr")
  
  # Calculate alpha diversity on transformed data
  alpha_orig <- calculate_diversity(fgt, method = "shannon")
  alpha_rel <- calculate_diversity(fgt_rel, assay_name = "relative_counts", normalize = FALSE)
  alpha_log <- calculate_diversity(fgt_log, assay_name = "log_counts", normalize = TRUE)
  
  # Relative abundance transformation with normalize=FALSE should match original with normalize=TRUE
  expect_equal(alpha_orig, alpha_rel)
  
  # Log-transformed data should give different results
  expect_false(isTRUE(all.equal(alpha_orig, alpha_log)))
  
  # Calculate beta diversity on transformed data
  beta_orig <- calculate_beta_diversity(fgt, method = "bray")
  beta_rel <- calculate_beta_diversity(fgt_rel, assay_name = "relative_counts", normalize = FALSE)
  beta_clr <- calculate_beta_diversity(fgt_clr, assay_name = "clr_counts", normalize = FALSE)
  
  # Relative abundance with normalize=FALSE should match original with normalize=TRUE
  expect_equal(beta_orig, beta_rel)
  
  # CLR-transformed data should give different results
  expect_false(isTRUE(all.equal(beta_orig, beta_clr)))
})

# Test diversity with sparse and edge-case data
test_that("diversity functions handle sparse and edge-case data correctly", {
  # Create FGTExperiment with sparse data
  sparse_fgt <- create_test_fgt(rows = 10, cols = 5)
  counts <- assay(sparse_fgt, "counts")
  counts[] <- 0
  
  # Sample 1: Just one non-zero taxon
  counts[1, 1] <- 100
  
  # Sample 2: 2 non-zero taxa
  counts[1:2, 2] <- c(50, 50)
  
  # Sample 3: 5 non-zero taxa with equal values
  counts[1:5, 3] <- 20
  
  # Sample 4: 5 non-zero taxa with different values
  counts[1:5, 4] <- c(10, 20, 30, 40, 50)
  
  # Sample 5: All zeros
  assay(sparse_fgt, "counts") <- counts
  
  # Alpha diversity
  richness <- calculate_diversity(sparse_fgt, method = "richness")
  expect_equal(richness, c(1, 2, 5, 5, 0))
  
  shannon <- calculate_diversity(sparse_fgt, method = "shannon")
  # Shannon should be 0 for sample 1, higher for others, 0 for all zeros
  expect_equal(shannon[1], 0)
  expect_equal(shannon[5], 0)
  expect_true(shannon[3] > shannon[2])  # More taxa = higher diversity
  expect_true(shannon[3] > shannon[4])  # Equal abundance > uneven for same richness
  
  simpson <- calculate_diversity(sparse_fgt, method = "simpson")
  # Simpson should be 0 for sample 1, higher for others, 0 for all zeros
  expect_equal(simpson[1], 0)
  expect_equal(simpson[5], 0)
  
  # Beta diversity
  bray <- calculate_beta_diversity(sparse_fgt, method = "bray")
  # Sample 5 (all zeros) vs any other sample should have distance 1
  expect_equal(bray[5, 1:4], rep(1, 4))
  
  # Sample 1 and 2 share one taxon
  expect_lt(bray[1, 2], 1)
  
  # Samples 3 and 4 have same taxa but different abundances
  expect_lt(bray[3, 4], bray[1, 3])
  
  # Jaccard (presence/absence)
  jaccard <- calculate_beta_diversity(sparse_fgt, method = "jaccard")
  # Samples 3 and 4 have identical presence/absence pattern
  expect_equal(jaccard[3, 4], 0)
})