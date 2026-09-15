# Script to generate example datasets for microFGT package
# The datasets created here will be stored in inst/extdata and will be accessible via
# the load_example_data() function when the package is installed.

# Make sure we're in the right directory
if (!file.exists("DESCRIPTION") || !grepl("microFGT", readLines("DESCRIPTION", n = 1))) {
  stop("Run this script from the package root directory.")
}

# Create directory if it doesn't exist
if (!dir.exists("inst/extdata")) {
  dir.create("inst/extdata", recursive = TRUE)
}

# Source necessary functions
devtools::load_all()

# Generate small amplicon dataset
small_amplicon <- generate_fgt_example_data(
  n_samples = 10,
  n_features = 50,
  community_types = c("CST-I", "CST-III", "CST-IV"),
  include_tree = TRUE,
  random_seed = 123,
  format = "list"
)

# Save components
saveRDS(small_amplicon$counts, "inst/extdata/microFGT_example_small_amplicon_counts.rds")
saveRDS(small_amplicon$taxonomy, "inst/extdata/microFGT_example_small_amplicon_taxonomy.rds")
saveRDS(small_amplicon$metadata, "inst/extdata/microFGT_example_small_amplicon_metadata.rds")
saveRDS(small_amplicon$tree, "inst/extdata/microFGT_example_small_amplicon_tree.rds")

# Generate medium amplicon dataset
medium_amplicon <- generate_fgt_example_data(
  n_samples = 30,
  n_features = 100,
  sample_groups = c("Healthy", "BV", "Intermediate"),
  group_proportions = c(0.5, 0.3, 0.2),
  community_types = c("CST-I", "CST-III", "CST-IV", "CST-V"),
  include_tree = TRUE,
  random_seed = 456,
  format = "list"
)

# Save components
saveRDS(medium_amplicon$counts, "inst/extdata/microFGT_example_medium_amplicon_counts.rds")
saveRDS(medium_amplicon$taxonomy, "inst/extdata/microFGT_example_medium_amplicon_taxonomy.rds")
saveRDS(medium_amplicon$metadata, "inst/extdata/microFGT_example_medium_amplicon_metadata.rds")
saveRDS(medium_amplicon$tree, "inst/extdata/microFGT_example_medium_amplicon_tree.rds")

# Create a very simple metagenomic example
small_metagenomic <- generate_fgt_example_data(
  n_samples = 5,
  n_features = 30,
  include_tree = FALSE,
  random_seed = 789,
  format = "list"
)

# Save components
saveRDS(small_metagenomic$counts, "inst/extdata/microFGT_example_small_metagenomic_counts.rds")
saveRDS(small_metagenomic$taxonomy, "inst/extdata/microFGT_example_small_metagenomic_taxonomy.rds")
saveRDS(small_metagenomic$metadata, "inst/extdata/microFGT_example_small_metagenomic_metadata.rds")

# Create README file for extdata
cat('# Example Data Files for microFGT

This directory contains example datasets for testing and demonstration.

Available datasets:
- microFGT_example_small_amplicon: Small amplicon dataset (10 samples, 50 features)
- microFGT_example_medium_amplicon: Medium amplicon dataset (30 samples, 100 features)
- microFGT_example_small_metagenomic: Small metagenomic dataset (5 samples, 30 features)

These datasets can be loaded using the `load_example_data()` function:

```r
# Load small amplicon dataset
data <- load_example_data("small", "amplicon")

# Load medium dataset but return as list instead of FGTExperiment
data_list <- load_example_data("medium", "amplicon", as_fgt_experiment = FALSE)
```

The datasets were generated using the `generate_fgt_example_data()` function.
', file = "inst/extdata/README.md")

message("Example datasets have been created in the inst/extdata directory.")