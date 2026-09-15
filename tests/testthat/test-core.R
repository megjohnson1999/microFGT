test_that("FGTExperiment constructor works", {
  # Create simple test data
  counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("OTU", 1:10)
  colnames(counts) <- paste0("Sample", 1:10)
  
  # Create minimal FGTExperiment
  fgt_exp <- FGTExperiment(assays = list(counts = counts))
  
  # Tests
  expect_s4_class(fgt_exp, "FGTExperiment")
  expect_s4_class(fgt_exp, "TreeSummarizedExperiment")
  expect_equal(dim(fgt_exp), c(10, 10))
  expect_equal(experimentType(fgt_exp), "amplicon")
})

test_that("filter_taxa works correctly", {
  # Create test data
  counts <- matrix(c(
    # High prevalence, high abundance
    rep(10, 10),
    # Low prevalence, high abundance 
    c(100, rep(0, 9)),
    # High prevalence, low abundance
    rep(1, 10)
  ), nrow = 3, byrow = TRUE)
  
  rownames(counts) <- c("high_both", "low_prev", "low_abund")
  colnames(counts) <- paste0("S", 1:10)
  
  fgt_exp <- FGTExperiment(assays = list(counts = counts))
  
  # Filter with different thresholds
  res1 <- filter_taxa(fgt_exp, min_prevalence = 0.5, min_abundance = 0.01)
  res2 <- filter_taxa(fgt_exp, min_prevalence = 0.5, min_abundance = 0.2)
  res3 <- filter_taxa(fgt_exp, min_prevalence = 0.2, min_abundance = 0.01)
  
  # Tests
  expect_equal(dim(res1), c(2, 10))  # Both high_both and low_abund pass
  expect_equal(dim(res2), c(1, 10))  # Only high_both passes
  expect_equal(dim(res3), c(3, 10))  # All taxa pass
})

test_that("transform_abundance works correctly", {
  # Create test data
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- paste0("OTU", 1:3)
  colnames(counts) <- paste0("S", 1:4)
  
  fgt_exp <- FGTExperiment(assays = list(counts = counts))
  
  # Apply transformations
  rel_result <- transform_abundance(fgt_exp, type = "relative")
  log_result <- transform_abundance(fgt_exp, type = "log")
  clr_result <- transform_abundance(fgt_exp, type = "clr")
  presence_result <- transform_abundance(fgt_exp, type = "presence")
  
  # Tests
  expect_true("relative" %in% names(SummarizedExperiment::assays(rel_result)))
  expect_true("log" %in% names(SummarizedExperiment::assays(log_result)))
  expect_true("clr" %in% names(SummarizedExperiment::assays(clr_result)))
  expect_true("presence" %in% names(SummarizedExperiment::assays(presence_result)))
  
  # Check relative abundance sums to 1 per sample
  rel_abundances <- SummarizedExperiment::assays(rel_result)[["relative"]]
  expect_equal(colSums(rel_abundances), rep(1, ncol(rel_abundances)))
  
  # Check presence values are logical
  presence_values <- SummarizedExperiment::assays(presence_result)[["presence"]]
  expect_type(presence_values, "logical")
  expect_equal(presence_values, counts > 0)
})

test_that("conversion between formats works", {
  skip_if_not_installed("phyloseq")
  
  # Create test data
  counts <- matrix(sample(0:100, 60, replace = TRUE), nrow = 10, ncol = 6)
  rownames(counts) <- paste0("OTU", 1:10)
  colnames(counts) <- paste0("Sample", 1:6)
  
  taxa <- data.frame(
    Kingdom = rep("Bacteria", 10),
    Phylum = sample(c("Firmicutes", "Bacteroidetes"), 10, replace = TRUE),
    row.names = rownames(counts)
  )
  
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 3),
    row.names = colnames(counts)
  )
  
  # Create from dada2-like objects
  dada2_fgt <- import_dada2(counts, taxa, metadata)
  
  # Tests for dada2 import
  expect_s4_class(dada2_fgt, "FGTExperiment")
  expect_equal(dim(dada2_fgt), c(10, 6))
  expect_equal(colnames(SummarizedExperiment::colData(dada2_fgt)), "group")
  
  # Create phyloseq object
  if (requireNamespace("phyloseq", quietly = TRUE)) {
    otu_table <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
    tax_table <- phyloseq::tax_table(as.matrix(taxa))
    sample_data <- phyloseq::sample_data(metadata)
    
    physeq <- phyloseq::phyloseq(otu_table, tax_table, sample_data)
    
    # Convert to FGTExperiment
    physeq_fgt <- phyloseq_to_fgt(physeq)
    
    # Tests for phyloseq conversion
    expect_s4_class(physeq_fgt, "FGTExperiment")
    expect_equal(dim(physeq_fgt), c(10, 6))
    expect_equal(colnames(SummarizedExperiment::colData(physeq_fgt)), "group")
  }
})