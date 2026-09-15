# microFGT Testing Infrastructure

## Overview

The microFGT package includes a comprehensive automated testing framework designed to ensure reliability, performance, and cross-platform compatibility.

## Test Categories

1. **Unit Tests**: Test individual functions and methods
2. **Edge Case Tests**: Test handling of unusual or extreme inputs
3. **Integration Tests**: Test interaction between components
4. **Performance Tests**: Benchmark performance and memory usage
5. **Platform Tests**: Test platform-specific functionality (SpeciateIT, VALENCIA, VIRGO)

## Running Tests

### Local Testing

```r
# Run all tests
testthat::test_local()

# Run specific test file
testthat::test_file("tests/testthat/test-edge-cases-comprehensive.R")

# Run automated test suite
source("tests/run_automated_tests.R")
run_automated_tests()

# Run specific categories
run_automated_tests(test_categories = c("edge_cases", "performance"))
```

### Command Line

```bash
# Run all tests
Rscript tests/run_automated_tests.R

# Run specific categories
Rscript tests/run_automated_tests.R edge_cases integration

# Run with coverage
R CMD check --as-cran
```

### Continuous Integration

Tests run automatically on:
- Every push to main/develop branches
- Every pull request
- Daily at 1 AM UTC
- Manual workflow dispatch

## Test Coverage

The package aims for >95% code coverage. Coverage reports are generated automatically and uploaded to Codecov.

To check coverage locally:
```r
covr::package_coverage()
covr::report()
```

## Performance Benchmarking

Performance tests run on multiple data sizes:
- Small: 100 × 20
- Medium: 1,000 × 50
- Large: 5,000 × 100
- Extra Large: 10,000 × 200

Benchmarks measure:
- Object creation time
- Transformation performance
- Memory usage
- Scalability

## Test Fixtures

Test fixtures are generated using helper functions:

```r
# Generate empty data
empty_data <- generate_test_fixture("empty")

# Generate sparse data
sparse_data <- generate_test_fixture("sparse", rows = 1000, cols = 100, sparsity = 0.95)

# Generate CST-based data
cst_data <- generate_test_fixture("CST", cols = 20)

# Generate time series data
ts_data <- generate_test_fixture("time_series", subjects = 5, timepoints = 4)
```

## Testing Best Practices

1. **Write tests first**: Follow TDD when adding new features
2. **Test edge cases**: Always test boundary conditions
3. **Use fixtures**: Use consistent test data generators
4. **Mock external dependencies**: Don't rely on external files/services
5. **Keep tests fast**: Individual tests should complete in <1 second
6. **Make tests independent**: Tests shouldn't depend on execution order
7. **Clear assertions**: Use descriptive expect_* messages

## Adding New Tests

1. Create test file: `tests/testthat/test-{feature-name}.R`
2. Include standard headers and helpers
3. Group related tests using `context()`
4. Write descriptive test names
5. Include edge cases and error conditions
6. Run locally before committing

Example test structure:
```r
context("Feature Name")

test_that("function handles normal input correctly", {
  # Arrange
  test_data <- generate_test_fixture("standard")
  
  # Act
  result <- my_function(test_data)
  
  # Assert
  expect_equal(dim(result), c(10, 5))
  expect_true(all(result >= 0))
})

test_that("function handles edge cases", {
  # Test empty input
  expect_error(my_function(NULL), "Input cannot be NULL")
  
  # Test invalid input
  expect_warning(my_function(-1), "Negative values")
})
```

## Debugging Failed Tests

1. Run test in isolation:
   ```r
   testthat::test_file("tests/testthat/test-failing.R")
   ```

2. Use browser() for interactive debugging:
   ```r
   test_that("debugging example", {
     browser()  # Execution stops here
     result <- my_function(input)
     expect_equal(result, expected)
   })
   ```

3. Check test reports in `test-results/` directory

4. Review CI logs for platform-specific issues

## Maintenance

- Review test coverage weekly
- Update performance baselines quarterly
- Review and refactor tests during major updates
- Monitor CI build times and optimize as needed

## Contact

For questions about testing, please file an issue on GitHub or contact the package maintainers.