# microFGT Testing Strategy

This document outlines a comprehensive testing strategy for the microFGT package, covering both automated tests and manual verification steps.

## Testing Goals

1. **Ensure code correctness**: Verify that all functions behave as expected
2. **Prevent regressions**: Ensure new changes don't break existing functionality
3. **Validate scientific accuracy**: Confirm biological interpretations are correct
4. **Test performance**: Evaluate runtime and memory usage with realistic datasets
5. **Verify usability**: Ensure the package can be effectively used by researchers

## Automated Testing Framework

### Unit Tests (via testthat)

- **Core Class Tests**
  - FGTExperiment class initialization and validation
  - Accessor and setter methods
  - Coercion methods between classes
  - Validation of slot constraints

- **Data Processing Functions**
  - Filtering functions (prevalence, abundance)
  - Transformation functions (relative abundance, log, CLR)
  - Taxonomic aggregation

- **Analysis Functions**
  - Diversity calculations
  - Differential abundance tests
  - Community state typing

- **Visualization Functions**
  - Plot generation correctness
  - Consistent styling and labeling

### Integration Tests

- **End-to-end workflows**
  - Data import → processing → analysis → visualization
  - Processing pipelines with different parameter combinations
  - Format conversion between ecosystem tools

### Performance Tests

- **Benchmarking key functions**
  - Runtime with increasing dataset sizes
  - Memory usage profiling
  - Comparison of algorithm implementations

## Implementation Plan

### Phase 1: Expand Unit Test Coverage

1. **Create test fixtures**
   - Mock datasets of varying sizes (small, medium, large)
   - Sample data representing different FGT microbiome profiles
   - Create helper functions for test setup

2. **Implement comprehensive unit tests**
   - Target at least 90% code coverage
   - Test both expected behavior and edge cases
   - Include tests for error conditions

3. **Set up test organization**
   - Group tests by functionality area
   - Implement shared test fixtures

### Phase 2: Integration and Workflow Tests

1. **Create end-to-end test workflows**
   - Standard data processing pipeline tests
   - Format conversion tests
   - Test interoperability with other Bioconductor packages

2. **Implement snapshot testing for visualizations**
   - Test plot generation for consistency
   - Validate plot elements and styling

### Phase 3: Performance Testing

1. **Create benchmarking suite**
   - Test with datasets of increasing size
   - Profile memory usage of key functions
   - Compare performance of alternative implementations

2. **Implement continuous performance testing**
   - Track performance metrics over time
   - Alert on significant regressions

## Manual Testing Plan

For the package maintainer and users to verify functionality:

### Installation Testing

```r
# From GitHub
remotes::install_github("shandley/microFGT")

# Local installation from source
install.packages("path/to/microFGT", repos = NULL, type = "source")

# Verify installation
library(microFGT)
packageVersion("microFGT")
```

### Basic Functionality Testing

Use the example data to verify core functionality:

```r
# Load example data
data <- load_example_data()

# Create FGTExperiment
fgt_exp <- create_fgt_experiment(data$counts, data$taxonomy, data$metadata)

# Verify structure
print(fgt_exp)
experimentType(fgt_exp)
head(assay(fgt_exp))
```

### Workflow Testing

Test the standard workflow steps:

```r
# Preprocess data
fgt_filtered <- filter_taxa(fgt_exp, min_prevalence = 0.1, min_abundance = 0.001)
fgt_rel <- transform_abundance(fgt_filtered, type = "relative")

# Analyze data
diversity_results <- alpha_diversity(fgt_rel)
ordination_result <- beta_diversity(fgt_rel, method = "bray")

# Visualize results
plot_diversity(diversity_results, group_var = "condition")
plot_ordination(ordination_result, color_by = "condition")
```

### Interactive Testing Checklist

For each release, manually verify:

- [ ] Package loads without errors
- [ ] Example data can be accessed
- [ ] Core functions run without errors
- [ ] Visualizations render correctly
- [ ] Documentation is accessible via `?function_name`
- [ ] Vignettes run without errors

## Continuous Integration

Automated tests will run on GitHub Actions for:

- [ ] R CMD check
- [ ] Unit tests via testthat
- [ ] Code coverage reporting
- [ ] Documentation building
- [ ] R package linting

## Next Steps

1. Implement unit tests for all core functions
2. Set up test fixtures for different data types
3. Expand code coverage to at least 90%
4. Implement integration tests for key workflows
5. Create performance benchmarking suite