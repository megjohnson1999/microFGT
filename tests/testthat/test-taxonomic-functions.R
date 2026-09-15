test_that("get_taxonomic_ranks works with standard ranks", {
  # Create a simple SummarizedExperiment with taxonomic ranks
  counts <- matrix(sample(0:100, 60, replace = TRUE), nrow = 10, ncol = 6)
  rownames(counts) <- paste0("Feature", 1:10)
  colnames(counts) <- paste0("Sample", 1:6)
  
  taxonomy <- data.frame(
    phylum = c("Firmicutes", "Bacteroidetes", "Firmicutes", "Proteobacteria", 
               "Bacteroidetes", "Firmicutes", "Actinobacteria", "Firmicutes", 
               "Proteobacteria", "Bacteroidetes"),
    genus = c("Lactobacillus", "Prevotella", "Streptococcus", "Escherichia", 
              "Bacteroides", "Staphylococcus", "Bifidobacterium", "Clostridium", 
              "Pseudomonas", "Porphyromonas"),
    row.names = rownames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(taxonomy)
  )
  
  # Test with standard_only = TRUE (default)
  ranks <- get_taxonomic_ranks(se)
  expect_equal(ranks, c("phylum", "genus"))
  
  # Test with standard_only = FALSE
  ranks_all <- get_taxonomic_ranks(se, standard_only = FALSE)
  expect_equal(ranks_all, c("phylum", "genus"))
  
  # Add non-standard column
  SummarizedExperiment::rowData(se)$other_info <- 1:10
  ranks <- get_taxonomic_ranks(se)
  expect_equal(ranks, c("phylum", "genus"))
  
  ranks_all <- get_taxonomic_ranks(se, standard_only = FALSE)
  expect_equal(ranks_all, c("phylum", "genus", "other_info"))
})

test_that("aggregate_taxa correctly aggregates at a taxonomic rank", {
  # Create a simple SummarizedExperiment with taxonomic ranks
  counts <- matrix(1:60, nrow = 10, ncol = 6)
  rownames(counts) <- paste0("Feature", 1:10)
  colnames(counts) <- paste0("Sample", 1:6)
  
  # Create taxonomy with some duplicates at genus level
  taxonomy <- data.frame(
    phylum = c(rep("Firmicutes", 4), rep("Bacteroidetes", 3), rep("Proteobacteria", 3)),
    genus = c(rep("Lactobacillus", 2), "Streptococcus", "Staphylococcus", 
              rep("Bacteroides", 2), "Prevotella", rep("Escherichia", 2), "Pseudomonas"),
    row.names = rownames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(taxonomy)
  )
  
  # Aggregate at genus level
  genus_se <- aggregate_taxa(se, rank = "genus")
  
  # Check dimensions
  expect_equal(nrow(genus_se), 6)  # 6 unique genera
  expect_equal(ncol(genus_se), 6)  # Same number of samples
  
  # Check that Lactobacillus (which had 2 features) has been summed correctly
  lacto_idx <- which(SummarizedExperiment::rowData(genus_se)$genus == "Lactobacillus")
  expect_equal(length(lacto_idx), 1)
  
  # The first two features were Lactobacillus, so we sum their counts
  expected_lacto_counts <- colSums(counts[1:2, ])
  expect_equal(SummarizedExperiment::assays(genus_se)$counts[lacto_idx, ], expected_lacto_counts)
  
  # Aggregate at phylum level
  phylum_se <- aggregate_taxa(se, rank = "phylum")
  
  # Check dimensions
  expect_equal(nrow(phylum_se), 3)  # 3 unique phyla
  expect_equal(ncol(phylum_se), 6)  # Same number of samples
  
  # Check that the metadata was added
  expect_true("aggregation_info" %in% names(SummarizedExperiment::metadata(phylum_se)))
  expect_equal(SummarizedExperiment::metadata(phylum_se)$aggregation_info$rank, "phylum")
  expect_equal(SummarizedExperiment::metadata(phylum_se)$aggregation_info$original_feature_count, 10)
  expect_equal(SummarizedExperiment::metadata(phylum_se)$aggregation_info$aggregated_feature_count, 3)
})

test_that("normalize_taxonomy correctly normalizes taxonomic names", {
  # Create a simple SummarizedExperiment with messy taxonomy
  counts <- matrix(1:30, nrow = 5, ncol = 6)
  rownames(counts) <- paste0("Feature", 1:5)
  colnames(counts) <- paste0("Sample", 1:6)
  
  taxonomy <- data.frame(
    phylum = c("p__Firmicutes", "p__Bacteroidetes", "Firmicutes", "p__Proteobacteria (99%)", NA),
    genus = c("g__Lactobacillus", "g__Bacteroides", "g__Streptococcus", "Escherichia", ""),
    row.names = rownames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(taxonomy)
  )
  
  # Normalize taxonomy
  norm_se <- normalize_taxonomy(se)
  
  # Check that prefixes were removed
  expect_equal(
    SummarizedExperiment::rowData(norm_se)$phylum,
    c("Firmicutes", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Unclassified_phylum")
  )
  
  # Check that confidence values were removed
  expect_false(any(grepl("\\(", SummarizedExperiment::rowData(norm_se)$phylum)))
  
  # Check that NA values were replaced
  expect_false(any(is.na(SummarizedExperiment::rowData(norm_se)$phylum)))
  expect_false(any(is.na(SummarizedExperiment::rowData(norm_se)$genus)))
  
  # Test with make_unique = FALSE (duplicate "Firmicutes" would remain)
  norm_se2 <- normalize_taxonomy(se, make_unique = FALSE)
  expect_equal(sum(SummarizedExperiment::rowData(norm_se2)$phylum == "Firmicutes"), 2)
})

test_that("create_tax_strings generates correct taxonomic strings", {
  # Create a simple SummarizedExperiment with taxonomic ranks
  counts <- matrix(1:24, nrow = 4, ncol = 6)
  rownames(counts) <- paste0("Feature", 1:4)
  colnames(counts) <- paste0("Sample", 1:6)
  
  taxonomy <- data.frame(
    kingdom = c("Bacteria", "Bacteria", "Bacteria", "Archaea"),
    phylum = c("Firmicutes", "Bacteroidetes", "Proteobacteria", NA),
    genus = c("Lactobacillus", "Bacteroides", NA, NA),
    species = c("L. crispatus", NA, NA, NA),
    row.names = rownames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(taxonomy)
  )
  
  # Create taxonomy strings with lineage format
  lineage_se <- create_tax_strings(se, format = "lineage")
  tax_strings <- SummarizedExperiment::rowData(lineage_se)$taxonomy
  
  # Check the first entry (complete lineage)
  expect_equal(
    tax_strings[1],
    "k__Bacteria;p__Firmicutes;g__Lactobacillus;s__L. crispatus"
  )
  
  # Second entry (missing species)
  expect_equal(
    tax_strings[2],
    "k__Bacteria;p__Bacteroidetes;g__Bacteroides"
  )
  
  # Create taxonomy strings with text format
  text_se <- create_tax_strings(se, format = "text")
  tax_text <- SummarizedExperiment::rowData(text_se)$taxonomy
  
  # Check the first entry (complete lineage in text format)
  expect_equal(
    tax_text[1],
    "Bacteria Firmicutes Lactobacillus L. crispatus"
  )
})

test_that("parse_tax_strings correctly splits taxonomic strings", {
  # Create a simple SummarizedExperiment with taxonomy strings
  counts <- matrix(1:18, nrow = 3, ncol = 6)
  rownames(counts) <- paste0("Feature", 1:3)
  colnames(counts) <- paste0("Sample", 1:6)
  
  taxonomy <- data.frame(
    taxonomy = c(
      "k__Bacteria;p__Firmicutes;g__Lactobacillus;s__L. crispatus",
      "k__Bacteria;p__Bacteroidetes;g__Bacteroides",
      "k__Bacteria;p__Proteobacteria"
    ),
    row.names = rownames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(taxonomy)
  )
  
  # Parse taxonomy strings
  parsed_se <- parse_tax_strings(se)
  
  # Check that the right columns were created
  rowdata <- SummarizedExperiment::rowData(parsed_se)
  expect_true(all(c("kingdom", "phylum", "genus", "species") %in% colnames(rowdata)))
  
  # Check values
  expect_equal(rowdata$kingdom, c("Bacteria", "Bacteria", "Bacteria"))
  expect_equal(rowdata$phylum, c("Firmicutes", "Bacteroidetes", "Proteobacteria"))
  expect_equal(rowdata$genus, c("Lactobacillus", "Bacteroides", NA))
  expect_equal(rowdata$species, c("L. crispatus", NA, NA))
  
  # Test with prefixes retained
  parsed_se2 <- parse_tax_strings(se, remove_prefixes = FALSE)
  rowdata2 <- SummarizedExperiment::rowData(parsed_se2)
  expect_equal(rowdata2$kingdom, c("k__Bacteria", "k__Bacteria", "k__Bacteria"))
  expect_equal(rowdata2$phylum, c("p__Firmicutes", "p__Bacteroidetes", "p__Proteobacteria"))
})