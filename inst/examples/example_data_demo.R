#!/usr/bin/env Rscript

# Demonstration of microFGT example data generation
library(microFGT)

# 1. Generate a simple example dataset
cat("Generating a small example dataset...\n")
small_data <- generate_fgt_example_data(
  n_samples = 8,
  n_features = 30,
  sample_groups = c("Healthy", "BV"),
  group_proportions = c(0.75, 0.25),
  community_types = c("CST-I", "CST-III", "CST-IV"),
  format = "list"
)

# Examine the data
cat("\nCounts matrix dimensions:", dim(small_data$counts), "\n")
cat("First few counts:\n")
print(small_data$counts[1:5, 1:3])

cat("\nTaxonomic classifications (first 5 features):\n")
print(small_data$taxonomy[1:5, ])

cat("\nSample metadata (first 3 samples):\n")
print(small_data$metadata[1:3, ])

# 2. Generate an FGTExperiment object
cat("\nGenerating an FGTExperiment object...\n")
fgt_experiment <- generate_fgt_example_data(
  n_samples = 8,
  n_features = 30,
  sample_groups = c("Healthy", "BV"),
  group_proportions = c(0.75, 0.25),
  community_types = c("CST-I", "CST-III", "CST-IV"),
  format = "FGTExperiment"
)

# Print basic information about the object
cat("\nFGTExperiment object:\n")
print(fgt_experiment)

# 3. Generate a complex dataset with specific parameters
cat("\nGenerating a more complex dataset...\n")
complex_data <- generate_fgt_example_data(
  n_samples = 20,
  n_features = 100,
  sample_groups = c("Healthy", "BV", "Intermediate"),
  group_proportions = c(0.5, 0.3, 0.2),
  community_types = c("CST-I", "CST-III", "CST-IV", "CST-V"),
  sequencing_depth = c(10000, 50000),
  include_tree = TRUE,
  format = "list"
)

# Show clinical parameter distributions by group
cat("\nClinical parameter distributions by group:\n")

# pH distribution
by_group <- split(complex_data$metadata$pH, complex_data$metadata$condition)
cat("\npH summary by group:\n")
print(lapply(by_group, summary))

# Nugent score distribution
by_group <- split(complex_data$metadata$Nugent_Score, complex_data$metadata$condition)
cat("\nNugent score summary by group:\n")
print(lapply(by_group, summary))

# 4. Save example data to files
output_dir <- tempdir()
cat("\nSaving example data to directory:", output_dir, "\n")

files <- save_example_data(
  seq_tab = small_data$counts,
  taxonomy = small_data$taxonomy,
  metadata = small_data$metadata,
  output_dir = output_dir,
  base_name = "microFGT_example_demo"
)

cat("\nSaved files:\n")
print(files)

# 5. Check taxonomic composition
if (require(SummarizedExperiment) && require(TreeSummarizedExperiment)) {
  # Convert to relative abundance to see taxonomic proportions
  counts_rel <- t(t(complex_data$counts) / colSums(complex_data$counts))
  
  # Get mean abundance by phylum and condition
  phylum_abundances <- list()
  
  for (i in 1:ncol(counts_rel)) {
    sample_id <- colnames(counts_rel)[i]
    condition <- complex_data$metadata[sample_id, "condition"]
    
    for (phylum in unique(complex_data$taxonomy$Phylum)) {
      phylum_idx <- which(complex_data$taxonomy$Phylum == phylum)
      abundance <- sum(counts_rel[phylum_idx, i])
      
      key <- paste(condition, phylum, sep = "_")
      if (!key %in% names(phylum_abundances)) {
        phylum_abundances[[key]] <- c()
      }
      phylum_abundances[[key]] <- c(phylum_abundances[[key]], abundance)
    }
  }
  
  # Calculate mean abundance for each phylum by condition
  phylum_means <- lapply(phylum_abundances, mean)
  
  cat("\nMean phylum abundances by condition:\n")
  for (condition in unique(complex_data$metadata$condition)) {
    cat("\n", condition, ":\n", sep = "")
    for (phylum in unique(complex_data$taxonomy$Phylum)) {
      key <- paste(condition, phylum, sep = "_")
      if (key %in% names(phylum_means)) {
        cat("  ", phylum, ": ", round(phylum_means[[key]] * 100, 1), "%\n", sep = "")
      }
    }
  }
}

cat("\nExample data demonstration completed.\n")