# Tests for taxonomic functions with FGTExperiment class

library(testthat)
library(microFGT)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)

# Helper function to create test taxonomy data
create_test_taxonomy <- function(rows = 10) {
  taxa <- data.frame(
    Kingdom = rep("Bacteria", rows),
    Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria"), 
                   rows, replace = TRUE),
    Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria"), 
                  rows, replace = TRUE),
    Order = sample(c("Lactobacillales", "Clostridiales", "Bacteroidales", "Bifidobacteriales"), 
                  rows, replace = TRUE),
    Family = sample(c("Lactobacillaceae", "Bifidobacteriaceae", "Prevotellaceae"), 
                   rows, replace = TRUE),
    Genus = sample(c("Lactobacillus", "Gardnerella", "Prevotella", "Sneathia"), 
                  rows, replace = TRUE),
    Species = paste0("species", seq_len(rows)),
    row.names = paste0("Feature", seq_len(rows)),
    stringsAsFactors = FALSE
  )
  return(taxa)
}

# Tests for get_taxonomic_ranks function
test_that("get_taxonomic_ranks identifies taxonomic ranks correctly", {
  # Create FGTExperiment with known taxonomy
  fgt <- create_test_fgt(rows = 10, cols = 5, include_taxonomy = TRUE)
  
  # Test with standard_only = TRUE (default)
  ranks <- get_taxonomic_ranks(fgt)
  expect_true(is.character(ranks))
  expect_true(all(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %in% ranks))
  
  # Test with standard_only = FALSE
  all_ranks <- get_taxonomic_ranks(fgt, standard_only = FALSE)
  expect_true(is.character(all_ranks))
  expect_true(all(ranks %in% all_ranks))
  
  # Test with object having no rowData
  no_taxa_fgt <- create_test_fgt(rows = 5, cols = 3, include_taxonomy = FALSE)
  empty_ranks <- get_taxonomic_ranks(no_taxa_fgt)
  expect_equal(length(empty_ranks), 0)
  
  # Test with non-SummarizedExperiment object
  expect_error(get_taxonomic_ranks(matrix(1:10)), "must be a SummarizedExperiment-like object")
})

# Tests for aggregate_taxa function
test_that("aggregate_taxa correctly aggregates at specified taxonomic ranks", {
  # Create FGTExperiment with known taxonomy
  fgt <- create_test_fgt(rows = 20, cols = 5, include_taxonomy = TRUE)
  
  # We need to ensure there are duplicates in the taxonomic ranks
  # Manually set some duplicate phyla
  phyla <- c("Firmicutes", "Bacteroidetes", "Actinobacteria")
  rowData(fgt)$Phylum <- rep(phyla, length.out = 20)
  
  # Aggregate at phylum level
  phylum_fgt <- aggregate_taxa(fgt, rank = "Phylum")
  
  # Check result is still an FGTExperiment
  expect_s4_class(phylum_fgt, "FGTExperiment")
  
  # Should have 3 rows (one per phylum)
  expect_equal(nrow(phylum_fgt), 3)
  
  # Check row names match phyla
  expect_setequal(rownames(phylum_fgt), phyla)
  
  # Check column count is unchanged
  expect_equal(ncol(phylum_fgt), 5)
  
  # Check counts aggregation: total counts should be preserved
  expect_equal(sum(assay(fgt)), sum(assay(phylum_fgt)))
  
  # Test with a different rank (Genus)
  genera <- unique(rowData(fgt)$Genus)
  genus_fgt <- aggregate_taxa(fgt, rank = "Genus")
  expect_s4_class(genus_fgt, "FGTExperiment")
  expect_equal(nrow(genus_fgt), length(genera))
  expect_setequal(rownames(genus_fgt), genera)
  
  # Test aggregation with empty values
  # Set some genus values to empty
  rowData(fgt)$Genus[c(1, 5, 10)] <- ""
  rowData(fgt)$Genus[c(2, 6)] <- NA
  genus_empty_fgt <- aggregate_taxa(fgt, rank = "Genus", empty_label = "Unknown")
  expect_true("Unknown" %in% rownames(genus_empty_fgt))
  
  # Test aggregation metadata
  expect_true(!is.null(metadata(genus_fgt)$aggregation_info))
  expect_equal(metadata(genus_fgt)$aggregation_info$rank, "Genus")
  expect_equal(metadata(genus_fgt)$aggregation_info$original_feature_count, 20)
  
  # Test error handling
  expect_error(aggregate_taxa(fgt, rank = "NonExistentRank"), 
               "not found in rowData")
  
  # Test with non-SummarizedExperiment object
  expect_error(aggregate_taxa(matrix(1:10)), 
               "must be a SummarizedExperiment-like object")
})

# Tests for normalize_taxonomy function
test_that("normalize_taxonomy correctly normalizes taxonomic names", {
  # Create FGTExperiment with messy taxonomy
  fgt <- create_test_fgt(rows = 10, cols = 5, include_taxonomy = TRUE)
  
  # Manually mess up the taxonomy with prefixes and confidence values
  rowData(fgt)$Phylum <- paste0("p__", rowData(fgt)$Phylum)
  rowData(fgt)$Genus <- paste0("g__", rowData(fgt)$Genus)
  rowData(fgt)$Family[1:5] <- paste0(rowData(fgt)$Family[1:5], " (97%)")
  rowData(fgt)$Genus[6:10] <- paste0(rowData(fgt)$Genus[6:10], " (88%)")
  
  # Add some duplicates
  rowData(fgt)$Genus[c(2, 4)] <- "g__Lactobacillus"
  
  # Run normalize_taxonomy
  normalized_fgt <- normalize_taxonomy(fgt)
  
  # Check result is still an FGTExperiment
  expect_s4_class(normalized_fgt, "FGTExperiment")
  
  # Check prefixes are removed
  phyla <- rowData(normalized_fgt)$Phylum
  expect_false(any(grepl("^p__", phyla)))
  
  genera <- rowData(normalized_fgt)$Genus
  expect_false(any(grepl("^g__", genera)))
  
  # Check confidence values are removed
  families <- rowData(normalized_fgt)$Family
  expect_false(any(grepl("\\([0-9]+%\\)", families)))
  
  genera <- rowData(normalized_fgt)$Genus
  expect_false(any(grepl("\\([0-9]+%\\)", genera)))
  
  # Test with remove_prefixes = FALSE
  keep_prefix_fgt <- normalize_taxonomy(fgt, remove_prefixes = FALSE)
  phyla <- rowData(keep_prefix_fgt)$Phylum
  expect_true(any(grepl("^p__", phyla)))
  
  # Test with remove_confidence = FALSE
  keep_conf_fgt <- normalize_taxonomy(fgt, remove_confidence = FALSE)
  families <- rowData(keep_conf_fgt)$Family
  expect_true(any(grepl("\\([0-9]+%\\)", families)))
  
  # Test with make_unique = TRUE (default)
  # Check that duplicated genera are now unique
  genera <- rowData(normalized_fgt)$Genus
  expect_equal(length(unique(genera)), length(genera))
  
  # Test with non-SummarizedExperiment object
  expect_error(normalize_taxonomy(matrix(1:10)), 
               "must be a SummarizedExperiment-like object")
})

# Tests for create_tax_strings function
test_that("create_tax_strings correctly creates taxonomic strings", {
  # Create FGTExperiment with taxonomy
  fgt <- create_test_fgt(rows = 10, cols = 5, include_taxonomy = TRUE)
  
  # Create lineage strings
  lineage_fgt <- create_tax_strings(fgt, format = "lineage")
  
  # Check result is still an FGTExperiment
  expect_s4_class(lineage_fgt, "FGTExperiment")
  
  # Check taxonomy column exists
  expect_true("taxonomy" %in% colnames(rowData(lineage_fgt)))
  
  # Check format of strings
  tax_strings <- rowData(lineage_fgt)$taxonomy
  expect_true(all(grepl("^k__Bacteria;p__", tax_strings)))
  expect_true(all(grepl(";g__[^;]+;s__species", tax_strings)))
  
  # Test with format = "text"
  text_fgt <- create_tax_strings(fgt, format = "text")
  text_strings <- rowData(text_fgt)$taxonomy
  expect_false(any(grepl(";", text_strings)))
  expect_true(all(grepl("Bacteria", text_strings)))
  
  # Test with add_prefixes = FALSE
  no_prefix_fgt <- create_tax_strings(fgt, add_prefixes = FALSE)
  no_prefix_strings <- rowData(no_prefix_fgt)$taxonomy
  expect_false(any(grepl("k__", no_prefix_strings)))
  expect_false(any(grepl("p__", no_prefix_strings)))
  
  # Test with custom new_column
  custom_col_fgt <- create_tax_strings(fgt, new_column = "full_taxonomy")
  expect_true("full_taxonomy" %in% colnames(rowData(custom_col_fgt)))
  
  # Test with specific ranks
  specific_ranks_fgt <- create_tax_strings(fgt, ranks = c("Phylum", "Genus"))
  specific_strings <- rowData(specific_ranks_fgt)$taxonomy
  
  # Should only contain phylum and genus
  expect_true(all(grepl("^p__", specific_strings)))
  expect_true(all(grepl(";g__", specific_strings)))
  expect_false(any(grepl("c__", specific_strings)))
  
  # Test error handling
  # Make a copy with missing ranks
  limited_fgt <- fgt
  rowData(limited_fgt) <- rowData(limited_fgt)[, c("Kingdom", "Phylum")]
  expect_error(create_tax_strings(limited_fgt, ranks = c("Kingdom", "Phylum", "Genus")), 
               "ranks are missing")
  
  # Test with non-SummarizedExperiment object
  expect_error(create_tax_strings(matrix(1:10)), 
               "must be a SummarizedExperiment-like object")
})

# Tests for parse_tax_strings function
test_that("parse_tax_strings correctly parses taxonomic strings", {
  # Create FGTExperiment with taxonomy
  fgt <- create_test_fgt(rows = 10, cols = 5, include_taxonomy = TRUE)
  
  # Create lineage strings first
  lineage_fgt <- create_tax_strings(fgt, format = "lineage")
  
  # Remove the individual rank columns to ensure they're recreated from the string
  tax_string_column <- rowData(lineage_fgt)$taxonomy
  rowData(lineage_fgt) <- S4Vectors::DataFrame(taxonomy = tax_string_column)
  
  # Parse the strings back to individual ranks
  parsed_fgt <- parse_tax_strings(lineage_fgt)
  
  # Check result is still an FGTExperiment
  expect_s4_class(parsed_fgt, "FGTExperiment")
  
  # Check taxonomic ranks were recreated
  rank_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  expect_true(all(rank_cols %in% tolower(colnames(rowData(parsed_fgt)))))
  
  # Check prefix removal (default is TRUE)
  phylum_vals <- rowData(parsed_fgt)$phylum
  expect_false(any(grepl("^p__", phylum_vals)))
  
  # Test with remove_prefixes = FALSE
  keep_prefix_fgt <- parse_tax_strings(lineage_fgt, remove_prefixes = FALSE)
  phylum_vals <- rowData(keep_prefix_fgt)$phylum
  expect_true(all(grepl("^p__", phylum_vals)))
  
  # Test with custom split character
  # Create custom format first
  custom_fgt <- lineage_fgt
  rowData(custom_fgt)$taxonomy <- gsub(";", "|", rowData(custom_fgt)$taxonomy)
  
  custom_parsed <- parse_tax_strings(custom_fgt, split = "|", remove_prefixes = FALSE)
  expect_true(all(rank_cols %in% tolower(colnames(rowData(custom_parsed)))))
  
  # Test error handling
  # No taxonomy column
  no_tax_fgt <- fgt
  rowData(no_tax_fgt) <- rowData(no_tax_fgt)[, -which(colnames(rowData(no_tax_fgt)) == "taxonomy")]
  expect_error(parse_tax_strings(no_tax_fgt), "not found in rowData")
  
  # Test with non-SummarizedExperiment object
  expect_error(parse_tax_strings(matrix(1:10)), 
               "must be a SummarizedExperiment-like object")
})

# Test combined operations
test_that("taxonomic functions work together in a workflow", {
  # Create FGTExperiment with taxonomy
  fgt <- create_test_fgt(rows = 20, cols = 5, include_taxonomy = TRUE)
  
  # 1. Add some prefixes and confidence values
  rowData(fgt)$Phylum <- paste0("p__", rowData(fgt)$Phylum)
  rowData(fgt)$Genus[1:10] <- paste0(rowData(fgt)$Genus[1:10], " (90%)")
  
  # 2. Normalize taxonomy
  normalized_fgt <- normalize_taxonomy(fgt)
  
  # 3. Create taxonomic strings
  lineage_fgt <- create_tax_strings(normalized_fgt)
  
  # 4. Aggregate at genus level
  genus_fgt <- aggregate_taxa(lineage_fgt, rank = "Genus")
  
  # Verify results
  expect_s4_class(genus_fgt, "FGTExperiment")
  expect_lt(nrow(genus_fgt), nrow(fgt))  # Should have fewer rows after aggregation
  expect_true("taxonomy" %in% colnames(rowData(lineage_fgt)))  # Should have taxonomy strings
  
  # 5. Now try parsing the strings in the aggregated object
  parsed_fgt <- parse_tax_strings(genus_fgt)
  
  # The genus level object should have taxonomy strings and now parsed columns
  expect_true("kingdom" %in% tolower(colnames(rowData(parsed_fgt))))
})