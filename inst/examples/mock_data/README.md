# Mock Data Generators for microFGT

This directory contains functions to generate realistic mock data for various FGT microbiome analysis tools. These generators help with testing, development, and demonstration by creating synthetic data that mimics the format and content of real tool outputs.

## Available Mock Data Generators

The package provides mock data generators for three key FGT microbiome analysis tools:

1. **SpeciateIT** - Taxonomic classification of vaginal microbiome sequences
2. **VIRGO** - Vaginal non-redundant gene catalog and functional profiling
3. **VALENCIA** - Vaginal community state type (CST) classification

## Getting Started

### Loading the Mock Data Generators

You can use the mock data generators in your R script:

```r
# Load the microFGT package
library(microFGT)

# Source the mock data generators directly from the package
source(system.file("examples", "mock_data_generators.R", package = "microFGT"))
```

Alternatively, during development:

```r
# When working directly within the package directory
source("R/mock_data_generators.R")
```

### Basic Usage

#### SpeciateIT Mock Data

```r
# Generate mock SpeciateIT data
speciateit_data <- generate_mock_speciateit(n_sequences = 100)

# Write to a file
output_file <- write_mock_speciateit("mock_speciateIT.txt", n_sequences = 100)

# Convert to FGTMicrobiome-compatible format
fgt_format <- convert_speciateit_to_fgt(speciateit_data)
```

#### VIRGO Mock Data

```r
# Generate mock VIRGO data
virgo_data <- generate_mock_virgo(n_genes = 1000, n_samples = 10)

# Write to files
output_dir <- write_mock_virgo("virgo_output", "mock_virgo", n_genes = 1000, n_samples = 10)

# Convert to FGTMicrobiome-compatible format
fgt_format <- convert_virgo_to_fgt(virgo_data)
```

#### VALENCIA Mock Data

```r
# Generate mock VALENCIA data
valencia_data <- generate_mock_valencia(
  n_samples = 20, 
  include_scores = TRUE,
  add_taxonomic_data = TRUE
)

# Write to files
output_dir <- write_mock_valencia("valencia_output", "mock_valencia", n_samples = 20)

# Convert to FGTMicrobiome-compatible format
fgt_format <- convert_valencia_to_fgt(valencia_data)
```

#### Coordinated Dataset

```r
# Generate a coordinated dataset for all three tools
mock_dataset <- generate_mock_fgt_dataset(
  n_samples = 10,
  n_sequences = 1000,
  n_genes = 1000,
  create_files = TRUE,
  output_dir = "mock_data_output"
)
```

## Function Reference

### SpeciateIT Mock Data

#### `generate_mock_speciateit()`

Generates a dataframe with mock SpeciateIT classifications.

**Parameters:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `n_sequences` | Number of sequences to generate | 100 |
| `seq_id_prefix` | Prefix for sequence IDs | "Seq" |
| `species_distribution` | Named vector of probabilities for different species | Realistic vaginal microbiome distribution |
| `confidence_range` | Range of confidence values | c(0.7, 1.0) |
| `decisions_range` | Range of decision counts | c(5, 20) |
| `include_unclassified` | Proportion of sequences that should be unclassified | 0.05 |
| `seed` | Random seed for reproducibility | NULL |

**Returns:** A data frame with the columns: "Sequence ID", "Classification", "posterior probability", and "number of Decisions"

#### `write_mock_speciateit()`

Writes mock SpeciateIT data to a file.

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `file_path` | Path to the output file |
| `...` | Additional parameters passed to `generate_mock_speciateit()` |

**Returns:** Path to the created file (invisibly)

#### `convert_speciateit_to_fgt()`

Converts mock SpeciateIT data to FGTMicrobiome format.

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `speciateit_results` | Data frame from `generate_mock_speciateit()` |

**Returns:** A list with an `assignments` component suitable for `add_speciateit()`

### VIRGO Mock Data

#### `generate_mock_virgo()`

Generates mock VIRGO gene abundance data.

**Parameters:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `n_genes` | Number of genes to include | 1000 |
| `n_samples` | Number of samples to generate | 10 |
| `sample_id_prefix` | Prefix for sample IDs | "Sample" |
| `include_zeros` | Proportion of zero counts (sparsity) | 0.7 |
| `abundance_distribution` | Distribution for non-zero counts | "lognormal" |
| `abundance_params` | Parameters for the abundance distribution | c(1, 2) |
| `simulation_type` | Type of simulation | "realistic" |
| `cst_distribution` | Named vector of CST probabilities | Realistic CST distribution |
| `seed` | Random seed for reproducibility | NULL |

**Returns:** A list with components: `counts`, `genes`, `metadata`, and `test_out`

#### `write_mock_virgo()`

Writes mock VIRGO data to files.

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `dir_path` | Directory path to write output files |
| `base_name` | Base name for output files |
| `...` | Additional parameters passed to `generate_mock_virgo()` |

**Returns:** Path to the created directory (invisibly)

#### `convert_virgo_to_fgt()`

Converts mock VIRGO data to FGTMicrobiome format.

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `virgo_results` | List from `generate_mock_virgo()` |

**Returns:** A list with components suitable for `add_virgo()`

### VALENCIA Mock Data

#### `generate_mock_valencia()`

Generates mock VALENCIA community state type classifications.

**Parameters:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `n_samples` | Number of samples to generate | 20 |
| `sample_id_prefix` | Prefix for sample IDs | "Sample" |
| `cst_distribution` | Named vector of CST probabilities | Realistic CST distribution |
| `include_scores` | Whether to include CST probability scores | TRUE |
| `add_taxonomic_data` | Whether to add taxonomic abundance data | TRUE |
| `seed` | Random seed for reproducibility | NULL |

**Returns:** A list with components: `cst`, `scores` (optional), and `abundance` (optional)

#### `write_mock_valencia()`

Writes mock VALENCIA data to files.

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `dir_path` | Directory path to write output files |
| `base_name` | Base name for output files |
| `...` | Additional parameters passed to `generate_mock_valencia()` |

**Returns:** Path to the created directory (invisibly)

#### `convert_valencia_to_fgt()`

Converts mock VALENCIA data to FGTMicrobiome format.

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `valencia_results` | List from `generate_mock_valencia()` |

**Returns:** A list with components suitable for `add_valencia()`

### Coordinated Dataset

#### `generate_mock_fgt_dataset()`

Generates a coordinated set of mock data for all three tools with consistent samples and CSTs.

**Parameters:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `n_samples` | Number of samples to generate | 10 |
| `n_sequences` | Number of sequences for SpeciateIT data | 1000 |
| `n_genes` | Number of genes for VIRGO data | 1000 |
| `sample_id_prefix` | Prefix for sample IDs | "Sample" |
| `cst_distribution` | Named vector of CST probabilities | Realistic CST distribution |
| `seed` | Random seed for reproducibility | NULL |
| `create_files` | Whether to write files to disk | FALSE |
| `output_dir` | Directory to write files to | NULL |

**Returns:** A list with components for each tool's data, metadata, and optional file paths

## Example Workflows

### Testing FGTMicrobiome Functions

```r
# Generate coordinated mock data
mock_data <- generate_mock_fgt_dataset(n_samples = 10, n_sequences = 200)

# Create a basic FGTMicrobiome object
counts <- matrix(sample(0:100, 200*10, replace=TRUE), nrow = 200, ncol = 10)
rownames(counts) <- mock_data$speciateit$raw$`Sequence ID`[1:200]
colnames(counts) <- mock_data$metadata$Sample

taxonomy <- data.frame(
  Feature = rownames(counts),
  Kingdom = rep("Bacteria", 200),
  Phylum = rep("Firmicutes", 200),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  Sample = mock_data$metadata$Sample,
  Subject = paste0("Subject_", rep(1:5, each=2)),
  row.names = mock_data$metadata$Sample,
  stringsAsFactors = FALSE
)

# Create FGTMicrobiome object
fgt <- FGTMicrobiome(counts, taxonomy, metadata)

# Add mock tool results
fgt_with_speciate <- add_speciateit(fgt, mock_data$speciateit$fgt_compatible)
fgt_with_virgo <- add_virgo(fgt_with_speciate, mock_data$virgo$fgt_compatible)
fgt_with_valencia <- add_valencia(fgt_with_virgo, mock_data$valencia$fgt_compatible)

# The object now contains all three types of data
fgt_complete <- fgt_with_valencia
```

### Generating Test Data with Specific Characteristics

```r
# Create L. crispatus dominated test data (CST I)
cst_i_dist <- c("I" = 1.0)  # 100% CST I

l_crispatus_data <- generate_mock_fgt_dataset(
  n_samples = 5,
  cst_distribution = cst_i_dist,
  seed = 123
)

# Create diverse test data with custom CST distribution
diverse_dist <- c(
  "I" = 0.2,     # L. crispatus
  "III" = 0.2,   # L. iners
  "IV-A" = 0.3,  # G. vaginalis
  "IV-B" = 0.3   # Diverse anaerobes
)

diverse_data <- generate_mock_fgt_dataset(
  n_samples = 20,
  cst_distribution = diverse_dist,
  seed = 456
)
```

## Complete Demo

For a comprehensive demonstration of all mock data generators, see the [`mock_data_generators_demo.R`](mock_data_generators_demo.R) script in this directory.