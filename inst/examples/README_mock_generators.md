# Mock Data Generators for microFGT

This directory contains functions to generate mock data for various FGT microbiome analysis tools. These generators help with testing and development by creating synthetic data that mimics the format and content of real tool outputs.

## Current Mock Data Generators

### SpeciateIT

SpeciateIT is a tool for taxonomic classification of vaginal microbiome sequences. The mock data generator creates data that mimics the output format of SpeciateIT's `MC_order7_results.txt` file.

#### Functions

1. `generate_mock_speciateit()`: Generates a data frame with mock SpeciateIT classifications
2. `write_mock_speciateit()`: Writes mock SpeciateIT data to a file
3. `convert_speciateit_to_fgt()`: Converts mock SpeciateIT data to a format compatible with `add_speciateit()`

#### Usage

To use these functions in your own scripts:

```r
# Load the mock data generators
source(system.file("examples", "load_mock_generators.R", package = "microFGT"))

# Generate mock data with default parameters
mock_data <- generate_mock_speciateit(n_sequences = 100)

# Generate data with custom parameters
custom_data <- generate_mock_speciateit(
  n_sequences = 50,
  species_distribution = c(
    "Lactobacillus crispatus" = 0.7,
    "Gardnerella vaginalis" = 0.3
  ),
  include_unclassified = 0.1,
  seed = 123
)

# Write data to a file
output_file <- write_mock_speciateit("speciateit_mock.txt", n_sequences = 20)

# Convert to FGTMicrobiome-compatible format
fgt_format <- convert_speciateit_to_fgt(mock_data)

# Use with FGTMicrobiome
# fgt <- FGTMicrobiome(...)
# fgt_with_speciate <- add_speciateit(fgt, fgt_format)
```

For a complete example, see `complete_speciateIT_example.R` in this directory.

## Parameters

### `generate_mock_speciateit()`

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `n_sequences` | Number of sequences to generate | 100 |
| `seq_id_prefix` | Prefix for sequence IDs | "Seq" |
| `species_distribution` | Named vector of probabilities for different species | Realistic vaginal microbiome distribution |
| `confidence_range` | Range of confidence values | c(0.7, 1.0) |
| `decisions_range` | Range of decision counts | c(5, 20) |
| `include_unclassified` | Proportion of sequences that should be unclassified | 0.05 |
| `seed` | Random seed for reproducibility | NULL |

## Future Development

In the future, we plan to add mock data generators for:

1. VIRGO - Vaginal reference gene catalog for functional profiling
2. VALENCIA - Vaginal community state type classifier

## Contributing

Contributions to improve these mock data generators or add new ones are welcome. Please see the main package's CONTRIBUTING guidelines.