# microFGT Comprehensive Testing Strategy

This document outlines the testing strategy and infrastructure for the microFGT package, designed to ensure code quality, functionality, and performance.

## Testing Philosophy

The microFGT testing suite follows these core principles:

1. **Comprehensive Coverage**: Tests should cover all key functionality, edge cases, and error handling.
2. **Performance Monitoring**: Tests should evaluate and track performance with datasets of varying sizes.
3. **Integration Testing**: Tests should validate that components work correctly together.
4. **Documentation Verification**: Tests should ensure examples in documentation work as expected.
5. **Regression Prevention**: Tests should detect any regressions in functionality.

## Test Categories

The microFGT test suite is organized into the following categories:

### 1. Core Functionality Tests

Tests fundamental package functionality including:

- **Class Definition Tests**: Validate S4 class structures, inheritance, and validation methods.
- **Constructor Tests**: Ensure objects can be created with various input combinations.
- **Method Tests**: Verify methods work correctly for different object types.
- **Transformation Tests**: Test data transformation functions like relative abundance, log, CLR.

**Key files:**
- `test-FGTExperiment-class-definition.R`
- `test-FGTExperiment-transformations.R`
- `test-taxonomic-functions-FGTExperiment.R`
- `test-diversity-functions-FGTExperiment.R`
- `test-taxonomy-diversity-integration.R`

### 2. Integration Tests

Tests how components work together and end-to-end workflows:

- **Platform Integration**: Tests workflows for SpeciateIT, VIRGO, and VALENCIA data.
- **Cross-Component Tests**: Ensures different package components work together.
- **Import/Export Tests**: Validates data can be properly imported and exported.

**Key files:**
- `test-end-to-end-workflows.R`
- `test-multi-component-integration.R`

### 3. Performance and Robustness Tests

Tests package performance with various data sizes and edge cases:

- **Performance Benchmarks**: Tests execution time and memory usage with different data sizes.
- **Edge Cases**: Tests behavior with extreme inputs like empty matrices or large sparse datasets.
- **Data Validation**: Tests input validation and error handling.

**Key files:**
- `test-performance.R`
- `test-performance-benchmarks.R`
- `test-edge-cases.R`
- `test-validation.R`

### 4. Documentation and Interface Tests

Tests documentation examples and interface consistency:

- **Documentation Examples**: Validates examples in function documentation and vignettes.
- **Interface Consistency**: Ensures function interfaces follow consistent patterns.
- **Parameter Consistency**: Checks for consistent parameter names and defaults.

**Key files:**
- `test-documentation-examples.R`
- `test-interface-consistency.R`

## Test Infrastructure

The microFGT package uses several tools to support comprehensive testing:

### Test Runners

- **run_all_tests.R**: Master script to run all test suites with various options.
- **run_all_tests_enhanced.R**: Enhanced version with coverage, profiling, and reporting.
- **run_composition_tests.R**, **run_advanced_tests.R**, etc.: Category-specific test runners.
- **generate_test_report.R**: Generates HTML reports with test results, coverage, and performance metrics.

### Test Helpers

- **helper-microFGT.R**: Core test helpers including fixture generation and validation functions.
- **helper-test-reporting.R**: Advanced test helpers for performance testing and report generation.

### Test Fixtures

The package provides several types of test data fixtures:

- **Standard Test Data**: Small, predictable datasets for basic tests.
- **Sparse Test Data**: Datasets with varying sparsity to test efficiency.
- **Edge Case Data**: Empty, single-feature, or other edge cases.
- **CST-specific Data**: Data with community state types for testing taxonomic functions.
- **Time Series Data**: Longitudinal data for testing temporal analysis.
- **Large-scale Data**: Bigger datasets for performance testing.

## Running Tests

### Basic Test Execution

Run all tests (limited mode):
```r
Rscript run_all_tests.R
```

Run specific test categories:
```r
Rscript run_all_tests.R --core-only    # Run only core tests
Rscript run_all_tests.R --integration-only
Rscript run_all_tests.R --performance-only
Rscript run_all_tests.R --docs-only
```

Run full test suite (including large tests):
```r
Rscript run_all_tests.R --full
```

### Enhanced Test Execution

Run tests with enhanced reporting:
```r
Rscript run_all_tests_enhanced.R --coverage --profile --report
```

Generate test reports without running tests:
```r
Rscript generate_test_report.R --no-tests
```

## Continuous Integration

The microFGT package uses GitHub Actions for continuous integration, which:

1. Runs the test suite on multiple R versions
2. Calculates and tracks test coverage
3. Generates test reports
4. Performs package checks

## Test Maintenance

To maintain the test suite:

1. Add tests for any new functionality
2. Create regression tests for fixed bugs
3. Update test fixtures as needed for new data types
4. Review performance tests periodically to ensure they remain relevant

## Coverage Goals

The microFGT package has the following test coverage goals:

- Core Functionality: ≥ 95% coverage
- Exported Functions: 100% coverage
- Overall Package: ≥ 90% coverage

## Performance Benchmarks

Performance tests track:

1. Execution time scaling with dataset size
2. Memory usage with different operations
3. Function call overhead
4. Algorithmic efficiency

## Test Report Interpretation

The test reports provide several key metrics:

- **Pass Rate**: Percentage of tests that pass
- **Coverage**: Percentage of code covered by tests
- **Performance Metrics**: Execution time and memory usage
- **Health Score**: Composite score of test quality

## Contributing Tests

When contributing to the microFGT package:

1. Add tests for any new functionality
2. Run the test suite to ensure no regressions
3. Include performance tests for computationally intensive operations
4. Verify documentation examples work correctly