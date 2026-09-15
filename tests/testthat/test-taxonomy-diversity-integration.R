# Tests for integration between taxonomy and diversity functions

library(testthat)
library(microFGT)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)

test_that("taxonomic aggregation affects diversity metrics in expected ways", {
  # Create FGTExperiment with taxonomy and controlled abundance patterns
  fgt <- create_test_fgt(rows = 20, cols = 5, include_taxonomy = TRUE)
  
  # Manual control of taxonomy to create nested structure
  # We'll create 4 phyla with 5 genera each (all unique)
  phyla <- rep(c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria"), each = 5)
  genera <- paste0("Genus", 1:20)
  
  # Manually set abundances to have equal values within each phylum
  counts <- matrix(0, nrow = 20, ncol = 5)
  rownames(counts) <- paste0("Feature", 1:20)
  colnames(counts) <- paste0("Sample", 1:5)
  
  # Set counts so that each phylum has different total abundance
  # Firmicutes: 40% of abundance
  counts[1:5,] <- 8
  # Bacteroidetes: 30% of abundance
  counts[6:10,] <- 6
  # Proteobacteria: 20% of abundance
  counts[11:15,] <- 4
  # Actinobacteria: 10% of abundance
  counts[16:20,] <- 2
  
  # Update FGTExperiment object
  rowData(fgt)$Phylum <- phyla
  rowData(fgt)$Genus <- genera
  assay(fgt, "counts") <- counts
  
  # Calculate alpha diversity on original data
  shannon_orig <- calculate_diversity(fgt, method = "shannon")
  richness_orig <- calculate_diversity(fgt, method = "richness")
  
  # Aggregate at phylum level
  phylum_fgt <- aggregate_taxa(fgt, rank = "Phylum")
  
  # Calculate diversity on aggregated data
  shannon_phylum <- calculate_diversity(phylum_fgt, method = "shannon")
  richness_phylum <- calculate_diversity(phylum_fgt, method = "richness")
  
  # Richness should be lower after aggregation
  expect_true(all(richness_phylum < richness_orig))
  expect_equal(richness_phylum, rep(4, 5))  # 4 phyla
  
  # Shannon diversity should be lower after aggregation
  expect_true(all(shannon_phylum < shannon_orig))
  
  # Calculate beta diversity on original and aggregated data
  beta_orig <- calculate_beta_diversity(fgt, method = "bray")
  beta_phylum <- calculate_beta_diversity(phylum_fgt, method = "bray")
  
  # Beta diversity should be affected by aggregation
  # Since we've made all samples identical, both beta matrices should be zeros
  expect_equal(beta_orig, matrix(0, nrow = 5, ncol = 5, dimnames = list(colnames(fgt), colnames(fgt))))
  expect_equal(beta_phylum, matrix(0, nrow = 5, ncol = 5, dimnames = list(colnames(phylum_fgt), colnames(phylum_fgt))))
})

test_that("taxonomy normalization preserves diversity metrics", {
  # Create FGTExperiment with messy taxonomy
  fgt <- create_test_fgt(rows = 10, cols = 5, include_taxonomy = TRUE)
  
  # Add prefixes and confidence values to some taxonomy levels
  rowData(fgt)$Phylum <- paste0("p__", rowData(fgt)$Phylum)
  rowData(fgt)$Genus <- paste0("g__", rowData(fgt)$Genus, " (90%)")
  
  # Calculate diversity before normalization
  shannon_before <- calculate_diversity(fgt, method = "shannon")
  beta_before <- calculate_beta_diversity(fgt, method = "bray")
  
  # Normalize taxonomy
  norm_fgt <- normalize_taxonomy(fgt)
  
  # Calculate diversity after normalization
  shannon_after <- calculate_diversity(norm_fgt, method = "shannon")
  beta_after <- calculate_beta_diversity(norm_fgt, method = "bray")
  
  # Normalization shouldn't affect abundance values, so diversity should be identical
  expect_equal(shannon_before, shannon_after)
  expect_equal(beta_before, beta_after)
  
  # Test aggregation with normalized taxonomy
  phylum_orig <- aggregate_taxa(fgt, rank = "Phylum")
  phylum_norm <- aggregate_taxa(norm_fgt, rank = "Phylum")
  
  # The results should have different row names (due to prefix removal)
  expect_false(identical(rownames(phylum_orig), rownames(phylum_norm)))
  
  # But abundance totals should be preserved
  expect_equal(colSums(assay(phylum_orig)), colSums(assay(phylum_norm)))
})

test_that("complete workflow from taxonomy to diversity analysis works", {
  # Create FGTExperiment with taxonomy
  fgt <- create_test_fgt(rows = 20, cols = 8, include_taxonomy = TRUE)
  
  # 1. Normalize taxonomy
  norm_fgt <- normalize_taxonomy(fgt)
  
  # 2. Create taxonomic strings
  lineage_fgt <- create_tax_strings(norm_fgt)
  
  # 3. Aggregate at phylum level
  phylum_fgt <- aggregate_taxa(lineage_fgt, rank = "Phylum")
  
  # 4. Transform to relative abundance
  rel_fgt <- transformAbundance(phylum_fgt, type = "relative")
  
  # 5. Calculate alpha diversity
  shannon <- calculate_diversity(rel_fgt, method = "shannon", 
                               assay_name = "relative_counts",
                               normalize = FALSE)
  
  # 6. Calculate beta diversity
  bray <- calculate_beta_diversity(rel_fgt, method = "bray",
                                 assay_name = "relative_counts",
                                 normalize = FALSE)
  
  # Verify results
  expect_true(is.numeric(shannon))
  expect_length(shannon, 8)
  expect_true(is.matrix(bray))
  expect_equal(dim(bray), c(8, 8))
  
  # Test workflow with CLR transformation
  clr_fgt <- transformAbundance(phylum_fgt, type = "clr")
  
  # Calculate diversity on CLR-transformed data
  shannon_clr <- calculate_diversity(clr_fgt, method = "shannon", 
                                   assay_name = "clr_counts",
                                   normalize = TRUE)
  
  # CLR-transformed data should give different shannon values
  expect_false(identical(shannon, shannon_clr))
})

test_that("diversity comparisons work between taxonomic groups", {
  # Create FGTExperiment with two distinct sample groups
  fgt <- create_test_fgt(rows = 20, cols = 10, include_taxonomy = TRUE)
  
  # Create two different sample types
  colData(fgt)$group <- rep(c("A", "B"), each = 5)
  
  # Modify counts to create systematic differences between groups
  # Group A: Firmicutes dominant
  # Group B: Bacteroidetes dominant
  counts <- assay(fgt, "counts")
  
  # Get indices for different phyla
  firmicutes_idx <- which(rowData(fgt)$Phylum == "Firmicutes")
  bacteroidetes_idx <- which(rowData(fgt)$Phylum == "Bacteroidetes")
  
  if (length(firmicutes_idx) > 0 && length(bacteroidetes_idx) > 0) {
    # Increase Firmicutes in group A
    counts[firmicutes_idx, colData(fgt)$group == "A"] <- 
      counts[firmicutes_idx, colData(fgt)$group == "A"] * 5
    
    # Increase Bacteroidetes in group B
    counts[bacteroidetes_idx, colData(fgt)$group == "B"] <- 
      counts[bacteroidetes_idx, colData(fgt)$group == "B"] * 5
    
    assay(fgt, "counts") <- counts
    
    # Aggregate at phylum level
    phylum_fgt <- aggregate_taxa(fgt, rank = "Phylum")
    
    # Calculate diversity by group
    shannon_values <- calculate_diversity(phylum_fgt, method = "shannon")
    
    # Group samples by group
    group_a_shannon <- shannon_values[colData(fgt)$group == "A"]
    group_b_shannon <- shannon_values[colData(fgt)$group == "B"]
    
    # Groups should have different diversity patterns
    # But we can't predict exactly how they'll differ due to random data generation
    # So just check that we can calculate diversity by group
    expect_length(group_a_shannon, 5)
    expect_length(group_b_shannon, 5)
    
    # Calculate beta diversity
    beta_div <- calculate_beta_diversity(phylum_fgt, method = "bray")
    
    # Within-group distances should be smaller than between-group distances
    # Get within-group and between-group distances
    within_a <- beta_div[colData(fgt)$group == "A", colData(fgt)$group == "A"]
    within_b <- beta_div[colData(fgt)$group == "B", colData(fgt)$group == "B"]
    between <- beta_div[colData(fgt)$group == "A", colData(fgt)$group == "B"]
    
    # Calculate average distances
    avg_within_a <- mean(within_a[lower.tri(within_a)])
    avg_within_b <- mean(within_b[lower.tri(within_b)])
    avg_between <- mean(between)
    
    # Between-group distances should be larger
    expect_true(avg_between > avg_within_a)
    expect_true(avg_between > avg_within_b)
  } else {
    skip("Not enough samples of different phyla for comparison test")
  }
})