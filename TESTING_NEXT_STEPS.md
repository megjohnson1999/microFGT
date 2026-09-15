# microFGT Testing Next Steps

## Overview

This document outlines the next steps for expanding the testing suite for the microFGT package, building upon the foundation of composition-based FGTExperiment class tests.

## Priority 1: Advanced Functionality Tests

### 1.1 Taxonomic Analysis Tests

- Create tests for taxonomic aggregation functions
- Test taxonomy parsing and manipulation
- Test hierarchical taxonomic operations
- Verify behavior with missing or malformed taxonomy data

```r
test_that("aggregate_taxa functions correctly", {
  # Test aggregation at different taxonomic levels
  # Test handling of missing taxa
  # Test proper name generation
  # Test preservation of row and column metadata
})
```

### 1.2 Diversity Analysis Tests

- Test alpha diversity calculations (Shannon, Simpson, etc.)
- Test beta diversity calculations (Bray-Curtis, UniFrac, etc.)
- Test ordination methods (PCoA, NMDS, etc.)
- Test visualization of diversity metrics

```r
test_that("diversity calculations are correct", {
  # Create test data with known diversity values
  # Verify calculated values match expectations
  # Test behavior with edge cases (all zeros, single taxon)
})
```

### 1.3 Import/Export Function Tests

- Test SpeciateIT data import/export
- Test VALENCIA data import/export
- Test VIRGO data import/export
- Test general-purpose formats (biom, phyloseq, etc.)

```r
test_that("import_speciateit correctly parses data", {
  # Test with both minimal and complete data
  # Verify proper taxonomy and count handling
  # Test error handling for malformed inputs
})
```

## Priority 2: Comprehensive Integration Tests

### 2.1 End-to-End Workflow Tests

- Test complete analysis workflows from import to visualization
- Test multi-platform integrations (SpeciateIT + VALENCIA + VIRGO)
- Test real-world data scenarios

```r
test_that("complete workflow functions correctly", {
  # Generate mock data
  # Import into FGTExperiment
  # Transform and analyze
  # Verify outputs at each step
})
```

### 2.2 Multi-Component Integration Tests

- Test interactions between transformation and analysis functions
- Test visualization functions with various inputs
- Test compatibility of different data types

```r
test_that("transformation and diversity work together", {
  # Test how diversity metrics change with transformations
  # Verify transformation effects on ordination
  # Test visualization of transformed data
})
```

## Priority 3: Performance and Robustness Tests

### 3.1 Performance Tests

- Test with large datasets (1000+ features, 100+ samples)
- Test memory usage with sparse matrices
- Test parallel processing capabilities

```r
test_that("large dataset performance is acceptable", {
  # Generate large mock dataset
  # Measure time for key operations
  # Verify memory doesn't exceed reasonable limits
})
```

### 3.2 Edge Case Tests

- Test with extreme data (very sparse, all zeros, etc.)
- Test with malformed inputs
- Test recovery from errors

```r
test_that("edge cases are handled properly", {
  # Test with all-zero matrices
  # Test with single-feature or single-sample data
  # Test with extreme values (very large, NaN, Inf)
})
```

### 3.3 Validation Tests

- Test data validation functions
- Test quality control processes
- Test handling of missing or corrupted data

```r
test_that("data validation catches problems", {
  # Test with corrupted taxonomy
  # Test with mismatched dimensions
  # Test with incorrect data types
})
```

## Priority 4: User-Facing Documentation Tests

### 4.1 Example Tests

- Test all examples in function documentation
- Test vignette code
- Test README examples

```r
test_that("documentation examples work", {
  # Extract examples from documentation
  # Run in clean environment
  # Verify outputs match expectations
})
```

### 4.2 Interface Tests

- Test function parameter interfaces
- Test deprecation warnings
- Test error messages for clarity

```r
test_that("function interfaces are consistent", {
  # Check parameter naming across functions
  # Verify parameter defaults
  # Test coherent error messages
})
```

## Implementation Plan

### Phase 1: Core Functionality (2-3 weeks)

- Implement taxonomic analysis tests
- Implement diversity analysis tests
- Implement import/export tests
- Extend transformation tests

### Phase 2: Integration (2-3 weeks)

- Implement end-to-end workflow tests
- Implement multi-component integration tests
- Create tests for real-world scenarios

### Phase 3: Performance and Edge Cases (1-2 weeks)

- Implement performance tests
- Implement edge case tests
- Implement validation tests

### Phase 4: Documentation and Cleanup (1 week)

- Implement documentation example tests
- Implement interface tests
- Consolidate and refactor test suite

## Success Criteria

1. **Coverage**: >80% test coverage for all package functionality
2. **Specificity**: Tests catch actual bugs and regressions
3. **Performance**: Tests run efficiently and complete in reasonable time
4. **Clarity**: Test errors provide clear indication of what failed
5. **Maintainability**: Test suite is well-organized and documented

This plan builds on the existing foundation of FGTExperiment class tests to create a comprehensive test suite covering all aspects of the microFGT package.