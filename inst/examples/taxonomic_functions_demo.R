library(microFGT)
library(SummarizedExperiment)

# Load example data (assuming small amplicon data is available)
counts_file <- system.file("extdata", "microFGT_example_small_amplicon_counts.rds", package = "microFGT")
taxonomy_file <- system.file("extdata", "microFGT_example_small_amplicon_taxonomy.rds", package = "microFGT")
metadata_file <- system.file("extdata", "microFGT_example_small_amplicon_metadata.rds", package = "microFGT")

counts <- readRDS(counts_file)
taxonomy <- readRDS(taxonomy_file)
metadata <- readRDS(metadata_file)

# Create a SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = counts),
  rowData = taxonomy,
  colData = metadata
)

# Print available taxonomic ranks
cat("Available taxonomic ranks:\n")
ranks <- get_taxonomic_ranks(se)
print(ranks)

# Normalize taxonomy (remove prefixes, clean up names)
se_norm <- normalize_taxonomy(se)
cat("\nNormalized taxonomy (first 5 entries):\n")
print(head(rowData(se_norm)[, ranks], 5))

# Create taxonomic strings
se_strings <- create_tax_strings(se_norm, format = "lineage")
cat("\nTaxonomic strings (first 5 entries):\n")
print(head(rowData(se_strings)$taxonomy, 5))

# Aggregate at genus level
genus_se <- aggregate_taxa(se_norm, rank = "Genus")
cat("\nDimensions after aggregation at genus level:\n")
cat("Original features:", nrow(se), "\n")
cat("Features after genus aggregation:", nrow(genus_se), "\n")

# Plot a simple feature count at phylum level
phylum_se <- aggregate_taxa(se_norm, rank = "Phylum")
phylum_counts <- colSums(assays(phylum_se)$counts)

cat("\nTotal counts per sample at phylum level:\n")
print(phylum_counts)

# Show aggregation info
cat("\nAggregation info:\n")
print(metadata(phylum_se)$aggregation_info)

# Demonstrate parsing taxonomic strings (create then parse back)
# First create strings if they don't exist
if (!"taxonomy" %in% colnames(rowData(se_norm))) {
  se_strings <- create_tax_strings(se_norm, format = "lineage")
} else {
  se_strings <- se_norm
}

# Then parse them
se_parsed <- parse_tax_strings(se_strings)
cat("\nParsed taxonomy (first 5 entries):\n")
print(head(rowData(se_parsed)[, ranks], 5))

cat("\nTaxonomic functions demo completed successfully!\n")