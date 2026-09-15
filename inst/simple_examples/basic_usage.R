#!/usr/bin/env Rscript

# Basic Usage Examples for microFGT
# This script demonstrates minimal working examples

# Load required packages
library(microFGT)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)

# =========================
# Example 1: Create a microFGT object
# =========================

# Create a simple count matrix
counts <- matrix(rpois(50, 20), nrow=10, ncol=5)
rownames(counts) <- paste0("Feature", 1:10)
colnames(counts) <- paste0("Sample", 1:5)

# Create simple taxonomy data
taxa <- data.frame(
  Kingdom = rep("Bacteria", 10),
  Phylum = sample(c("Firmicutes", "Bacteroidetes", "Proteobacteria"), 10, replace=TRUE),
  row.names = rownames(counts)
)

# Create sample metadata
metadata <- data.frame(
  Group = rep(c("A", "B"), c(2, 3)),
  pH = runif(5, min=3.8, max=4.5),
  row.names = colnames(counts)
)

# Create FGTExperiment using the workaround method
# Step 1: Create SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = counts),
  rowData = taxa,
  colData = metadata
)

# Step 2: Convert to TreeSummarizedExperiment
tse <- TreeSummarizedExperiment(se)

# Step 3: Create FGTExperiment
fgt <- methods::new("FGTExperiment",
                  tse,
                  experimentType = "amplicon",
                  fgtMetadata = S4Vectors::SimpleList())

# Step 4: Ensure assay is named "counts"
assayNames(fgt)[1] <- "counts"

# Display basic information about the object
print(fgt)

# =========================
# Example 2: Basic data manipulation
# =========================

# If transform_abundance is available, use it
if (exists("transform_abundance", mode="function", where="package:microFGT")) {
  # Transform to relative abundance
  # Note: we use tryCatch to handle potential errors
  tryCatch({
    fgt_rel <- transform_abundance(fgt, type = "relative", assay_name = "counts")
    print("Transformed to relative abundance:")
    print(head(assays(fgt_rel)$relative[, 1:3]))
  }, error = function(e) {
    cat("Error in transform_abundance:", conditionMessage(e), "\n")
  })
}

# If filter_taxa is available, use it
if (exists("filter_taxa", mode="function", where="package:microFGT")) {
  # Filter taxa
  tryCatch({
    fgt_filtered <- filter_taxa(fgt, min_prevalence = 0.5, min_abundance = 0.001, 
                               assay_name = "counts")
    cat("Filtered taxa:", nrow(fgt), "->", nrow(fgt_filtered), "\n")
  }, error = function(e) {
    cat("Error in filter_taxa:", conditionMessage(e), "\n")
  })
}

# If aggregate_taxa is available, use it
if (exists("aggregate_taxa", mode="function", where="package:microFGT")) {
  # Aggregate by phylum
  tryCatch({
    fgt_agg <- aggregate_taxa(fgt, rank = "Phylum", assay_name = "counts")
    cat("Aggregated to phylum level:", nrow(fgt), "->", nrow(fgt_agg), "\n")
    print(rownames(fgt_agg))
  }, error = function(e) {
    cat("Error in aggregate_taxa:", conditionMessage(e), "\n")
  })
}

# =========================
# Example 3: Basic visualization (if available)
# =========================

# If plot_taxa_composition is available, use it
if (exists("plot_taxa_composition", mode="function", where="package:microFGT")) {
  # First transform to relative abundance if needed
  if (!exists("fgt_rel")) {
    # Try manual relative abundance calculation
    counts_mat <- assays(fgt)$counts
    rel_counts <- t(t(counts_mat) / colSums(counts_mat))
    assays(fgt)[["relative"]] <- rel_counts
    fgt_rel <- fgt
  }
  
  # Plot taxa composition
  tryCatch({
    p <- plot_taxa_composition(fgt_rel, rank = "Phylum", top_n = 5)
    print("Taxa composition plot created successfully")
  }, error = function(e) {
    cat("Error in plot_taxa_composition:", conditionMessage(e), "\n")
  })
}

# If plot_alpha_diversity is available, use it
if (exists("plot_alpha_diversity", mode="function", where="package:microFGT")) {
  tryCatch({
    p <- plot_alpha_diversity(fgt, metrics = c("shannon"), group_var = "Group", 
                             assay_name = "counts")
    print("Alpha diversity plot created successfully")
  }, error = function(e) {
    cat("Error in plot_alpha_diversity:", conditionMessage(e), "\n")
  })
}

# =========================
# Example 4: Accessors
# =========================

# Extract assay data
counts_data <- assays(fgt)$counts
cat("Extracted count data dimensions:", dim(counts_data), "\n")

# Extract row data (taxonomy)
tax_data <- rowData(fgt)
cat("Taxonomy data columns:", colnames(tax_data), "\n")

# Extract column data (metadata)
meta_data <- colData(fgt)
cat("Metadata columns:", colnames(meta_data), "\n")

# Get experiment type
if (exists("experimentType", mode="function", where="package:microFGT")) {
  exp_type <- experimentType(fgt)
  cat("Experiment type:", exp_type, "\n")
} else {
  cat("Experiment type:", fgt@experimentType, "\n")
}

cat("\nExamples completed!\n")