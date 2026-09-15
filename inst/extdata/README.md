# Example Data Files for microFGT

This directory contains example datasets for testing and demonstration.

## Available Datasets

- **microFGT_example_small_amplicon**: Small amplicon dataset (10 samples, 50 features)
- **microFGT_example_medium_amplicon**: Medium amplicon dataset (30 samples, 100 features)
- **microFGT_example_small_metagenomic**: Small metagenomic dataset (5 samples, 30 features)

## Usage

These datasets can be loaded using the `load_example_data()` function:

```r
# Load small amplicon dataset
data <- load_example_data("small", "amplicon")

# Load medium dataset but return as list instead of FGTExperiment
data_list <- load_example_data("medium", "amplicon", as_fgt_experiment = FALSE)
```

## Dataset Generation

The datasets were generated using the `generate_fgt_example_data()` function. The generation script is located in `data-raw/prepare_example_data.R`.

## Dataset Content

Each dataset consists of multiple components:

- **counts**: A count matrix of taxa abundance
- **taxonomy**: Taxonomic classifications for each feature
- **metadata**: Sample metadata including clinical parameters
- **tree**: Phylogenetic tree (if available)