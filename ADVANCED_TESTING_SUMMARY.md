# Advanced Testing for microFGT: Taxonomic and Diversity Analysis Functions

This document summarizes the implementation of advanced testing for the microFGT package, focusing on taxonomic analysis and diversity calculation functions for the composition-based FGTExperiment class.

## Overview

We have implemented comprehensive tests for the taxonomic and diversity functions in the microFGT package. These tests validate core functionality, integration with FGTExperiment class, edge cases, and integrated workflows. Tests are organized in three main categories:

1. **Taxonomic Function Tests**: Testing functions to manipulate, aggregate, and normalize taxonomic information
2. **Diversity Function Tests**: Testing calculation of alpha and beta diversity metrics
3. **Integration Tests**: Testing how taxonomic and diversity functions work together in workflows

## Taxonomic Function Tests

File: `tests/testthat/test-taxonomic-functions-FGTExperiment.R`

This test file covers the following functions:

### `get_taxonomic_ranks()`
- Tests identification of standard taxonomic ranks
- Tests handling of non-standard taxonomic data
- Tests error conditions and input validation

### `aggregate_taxa()`
- Tests aggregation at different taxonomic levels (phylum, genus, etc.)
- Tests handling of empty or missing taxonomic values
- Tests aggregation metadata and proper count summation
- Tests integration with FGTExperiment class structure

### `normalize_taxonomy()`
- Tests removal of taxonomic prefixes (k__, p__, etc.)
- Tests removal of confidence values
- Tests handling of duplicate taxa
- Tests handling of various taxonomy formats

### `create_tax_strings()`
- Tests creation of lineage strings from individual rank columns
- Tests different output formats (lineage vs. text)
- Tests custom column naming
- Tests handling of missing data

### `parse_tax_strings()`
- Tests parsing of taxonomy strings into individual rank columns
- Tests handling of different separators
- Tests prefix preservation/removal
- Tests reconstruction of taxonomic information

## Diversity Function Tests

File: `tests/testthat/test-diversity-functions-FGTExperiment.R`

This test file covers the following functions:

### `calculate_diversity()` (Alpha Diversity)
- Tests Shannon diversity calculation
- Tests Simpson diversity calculation
- Tests Richness calculation
- Tests Evenness calculation
- Tests with different data distributions (even, uneven, sparse)
- Tests with transformed abundance data

### `calculate_beta_diversity()`
- Tests Bray-Curtis dissimilarity
- Tests Jaccard distance
- Tests UniFrac distance (when phylogenetic tree is available)
- Tests with different data similarity patterns
- Tests with sparse and zero-inflated data

## Integration Tests

File: `tests/testthat/test-taxonomy-diversity-integration.R`

These tests validate how taxonomic and diversity functions work together:

### Taxonomic Aggregation and Diversity
- Tests how taxonomic aggregation affects diversity metrics
- Verifies that aggregating taxa reduces both richness and Shannon diversity
- Tests effects on beta diversity matrices

### Taxonomy Normalization and Diversity Preservation
- Tests that normalizing taxonomy doesn't alter diversity metrics
- Validates that abundance totals remain consistent after normalization

### End-to-End Workflows
- Tests complete workflows from taxonomy normalization to diversity analysis
- Tests combination of transformations and diversity calculations
- Tests workflows with both alpha and beta diversity metrics

### Group Comparisons
- Tests diversity comparisons between taxonomic groups
- Validates that systematic differences between sample groups are detectable
- Tests between-group vs. within-group beta diversity patterns

## Test Runner

File: `run_advanced_tests.R`

A script that runs all advanced tests and reports results:
- Runs each test file independently
- Provides detailed output for each test
- Summarizes test results
- Returns appropriate exit code based on test success

## Next Steps

These tests establish a strong foundation for validating taxonomic and diversity functionality in the microFGT package. Future testing expansions could include:

1. **Performance Testing**: Testing with larger datasets (100+ samples, 1000+ features)
2. **Visualization Testing**: Testing functions that visualize diversity metrics
3. **Machine Learning Integration**: Testing functions that use diversity metrics for prediction
4. **Reference Datasets**: Testing against published datasets with known diversity values
5. **Import/Export Testing**: Testing round-trip workflows with different file formats

## Running the Tests

To run the advanced tests, use the provided script:

```bash
./run_advanced_tests.R
```

This will execute all test files and provide a summary of results.