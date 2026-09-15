#!/usr/bin/env Rscript

#' Automated Test Runner for microFGT
#' 
#' This script provides a comprehensive automated testing framework with
#' reporting, coverage analysis, and performance benchmarking.

# Load required packages
library(testthat)
library(covr)
library(bench)

# Source helper functions
source("tests/testthat/helper-test-reporting.R")

#' Run automated test suite
#' 
#' @param test_categories Character vector of test categories to run
#' @param generate_report Logical, whether to generate HTML report
#' @param check_coverage Logical, whether to check code coverage
#' @param min_coverage Numeric, minimum required coverage percentage
#' @param output_dir Character, directory for output files
#' @param parallel Logical, whether to run tests in parallel
#' 
#' @return List containing test results, coverage, and timing information
run_automated_tests <- function(test_categories = c("all"),
                               generate_report = TRUE,
                               check_coverage = TRUE,
                               min_coverage = 90,
                               output_dir = "test-results",
                               parallel = FALSE) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Start timing
  start_time <- Sys.time()
  
  # Initialize results
  results <- list(
    start_time = start_time,
    test_results = list(),
    coverage_results = NULL,
    performance_results = list(),
    summary = list()
  )
  
  # Define test categories
  all_categories <- list(
    unit = "test-[^-]+\\.R$",
    edge_cases = "test-edge-cases.*\\.R$",
    integration = "test-.*integration.*\\.R$",
    performance = "test-performance.*\\.R$",
    platform = "test-platform.*\\.R$",
    all = "test-.*\\.R$"
  )
  
  # Determine which tests to run
  if ("all" %in% test_categories) {
    test_patterns <- all_categories$all
  } else {
    test_patterns <- unlist(all_categories[test_categories])
  }
  
  # Set up parallel processing if requested
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    cores <- min(parallel::detectCores() - 1, 4)
    message(sprintf("Running tests in parallel using %d cores", cores))
  } else {
    cores <- 1
  }
  
  # Run tests
  message("Running automated tests...")
  test_results <- test_dir(
    "tests/testthat",
    filter = paste(test_patterns, collapse = "|"),
    reporter = "summary",
    stop_on_failure = FALSE
  )
  
  results$test_results <- test_results
  
  # Calculate test summary
  results$summary$total_tests <- length(test_results)
  results$summary$passed <- sum(sapply(test_results, function(x) x$passed))
  results$summary$failed <- sum(sapply(test_results, function(x) x$failed))
  results$summary$skipped <- sum(sapply(test_results, function(x) x$skipped))
  results$summary$warnings <- sum(sapply(test_results, function(x) x$warnings))
  
  # Check code coverage if requested
  if (check_coverage) {
    message("Calculating code coverage...")
    
    coverage <- package_coverage(
      type = "tests",
      combine_types = FALSE,
      line_exclusions = list(
        "R/zzz.R",
        "R/utils-pipe.R"
      )
    )
    
    coverage_pct <- percent_coverage(coverage)
    results$coverage_results <- list(
      coverage = coverage,
      percent = coverage_pct,
      meets_minimum = coverage_pct >= min_coverage
    )
    
    # Save coverage report
    coverage_report <- file.path(output_dir, "coverage.html")
    report(coverage, file = coverage_report)
    message(sprintf("Coverage report saved to: %s", coverage_report))
    
    # Check if coverage meets minimum
    if (coverage_pct < min_coverage) {
      warning(sprintf("Code coverage (%.1f%%) is below minimum threshold (%.1f%%)",
                     coverage_pct, min_coverage))
    }
  }
  
  # Run performance benchmarks for specific test categories
  if (any(c("performance", "all") %in% test_categories)) {
    message("Running performance benchmarks...")
    
    perf_results <- run_performance_benchmarks()
    results$performance_results <- perf_results
    
    # Save performance results
    perf_file <- file.path(output_dir, "performance_results.rds")
    saveRDS(perf_results, perf_file)
  }
  
  # End timing
  end_time <- Sys.time()
  results$end_time <- end_time
  results$duration <- difftime(end_time, start_time, units = "secs")
  
  # Generate report if requested
  if (generate_report) {
    message("Generating test report...")
    report_file <- generate_test_report(results, output_dir)
    message(sprintf("Test report saved to: %s", report_file))
  }
  
  # Print summary
  cat("\n===== TEST SUMMARY =====\n")
  cat(sprintf("Total tests: %d\n", results$summary$total_tests))
  cat(sprintf("Passed: %d\n", results$summary$passed))
  cat(sprintf("Failed: %d\n", results$summary$failed))
  cat(sprintf("Skipped: %d\n", results$summary$skipped))
  
  if (check_coverage) {
    cat(sprintf("Coverage: %.1f%%\n", results$coverage_results$percent))
  }
  
  cat(sprintf("Duration: %.2f seconds\n", as.numeric(results$duration)))
  cat("========================\n\n")
  
  # Return results
  invisible(results)
}

#' Run performance benchmarks
#' 
#' @return List of benchmark results
run_performance_benchmarks <- function() {
  # Define test sizes
  test_sizes <- list(
    small = c(rows = 100, cols = 20),
    medium = c(rows = 1000, cols = 50),
    large = c(rows = 5000, cols = 100),
    xlarge = c(rows = 10000, cols = 200)
  )
  
  results <- list()
  
  for (size_name in names(test_sizes)) {
    size <- test_sizes[[size_name]]
    
    # Generate test data
    test_data <- generate_test_fixture("large", 
                                     rows = size["rows"], 
                                     cols = size["cols"])
    
    # Benchmark FGTExperiment creation
    creation_bench <- mark(
      FGTExperiment = {
        FGTExperiment(
          assays = list(counts = test_data$counts),
          rowData = DataFrame(test_data$taxonomy),
          colData = DataFrame(test_data$metadata),
          experimentType = "amplicon"
        )
      },
      iterations = 3,
      check = FALSE
    )
    
    # Create object for transformation benchmarks
    fgt <- FGTExperiment(
      assays = list(counts = test_data$counts),
      rowData = DataFrame(test_data$taxonomy),
      colData = DataFrame(test_data$metadata),
      experimentType = "amplicon"
    )
    
    # Benchmark transformations
    transform_bench <- mark(
      relative = transformAbundance(fgt, type = "relative"),
      log = transformAbundance(fgt, type = "log"),
      clr = transformAbundance(fgt, type = "clr"),
      iterations = 3,
      check = FALSE
    )
    
    results[[size_name]] <- list(
      dimensions = size,
      creation = creation_bench,
      transformations = transform_bench
    )
  }
  
  return(results)
}

#' Generate HTML test report
#' 
#' @param results Test results object
#' @param output_dir Output directory
#' 
#' @return Path to generated report
generate_test_report <- function(results, output_dir) {
  report_file <- file.path(output_dir, "test_report.html")
  
  # Create HTML content
  html_content <- sprintf('
<!DOCTYPE html>
<html>
<head>
  <title>microFGT Test Report</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    .header { background-color: #f0f0f0; padding: 10px; }
    .summary { margin: 20px 0; }
    .passed { color: green; }
    .failed { color: red; }
    .warning { color: orange; }
    table { border-collapse: collapse; margin: 20px 0; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #f0f0f0; }
  </style>
</head>
<body>
  <div class="header">
    <h1>microFGT Automated Test Report</h1>
    <p>Generated: %s</p>
    <p>Duration: %.2f seconds</p>
  </div>
  
  <div class="summary">
    <h2>Test Summary</h2>
    <p>Total Tests: %d</p>
    <p class="passed">Passed: %d</p>
    <p class="failed">Failed: %d</p>
    <p>Skipped: %d</p>
    %s
  </div>
  
  %s
  
  %s
</body>
</html>',
    format(results$start_time, "%Y-%m-%d %H:%M:%S"),
    as.numeric(results$duration),
    results$summary$total_tests,
    results$summary$passed,
    results$summary$failed,
    results$summary$skipped,
    if (!is.null(results$coverage_results)) {
      sprintf('<p>Code Coverage: %.1f%% %s</p>',
              results$coverage_results$percent,
              ifelse(results$coverage_results$meets_minimum,
                     '<span class="passed">✓</span>',
                     '<span class="failed">✗</span>'))
    } else "",
    generate_performance_section(results$performance_results),
    generate_detailed_results_section(results$test_results)
  )
  
  writeLines(html_content, report_file)
  return(report_file)
}

#' Generate performance section for report
generate_performance_section <- function(perf_results) {
  if (length(perf_results) == 0) return("")
  
  section <- '<div class="performance">
    <h2>Performance Benchmarks</h2>
    <table>
      <tr>
        <th>Size</th>
        <th>Dimensions</th>
        <th>Creation Time</th>
        <th>Relative Transform</th>
        <th>Log Transform</th>
        <th>CLR Transform</th>
      </tr>'
  
  for (size_name in names(perf_results)) {
    result <- perf_results[[size_name]]
    section <- paste0(section, sprintf('
      <tr>
        <td>%s</td>
        <td>%d × %d</td>
        <td>%.3f s</td>
        <td>%.3f s</td>
        <td>%.3f s</td>
        <td>%.3f s</td>
      </tr>',
      toupper(size_name),
      result$dimensions["rows"],
      result$dimensions["cols"],
      median(result$creation$time) / 1e9,
      median(result$transformations$time[[1]]) / 1e9,
      median(result$transformations$time[[2]]) / 1e9,
      median(result$transformations$time[[3]]) / 1e9
    ))
  }
  
  section <- paste0(section, '</table></div>')
  return(section)
}

#' Generate detailed test results section
generate_detailed_results_section <- function(test_results) {
  # This would generate detailed test-by-test results
  # For now, return a placeholder
  return('<div class="details"><h2>Detailed Results</h2>
          <p>See individual test files for detailed results.</p></div>')
}

# Run tests if script is executed directly
if (sys.nframe() == 0) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default options
  categories <- if (length(args) > 0) args else "all"
  
  # Run tests
  results <- run_automated_tests(
    test_categories = categories,
    generate_report = TRUE,
    check_coverage = TRUE,
    parallel = TRUE
  )
  
  # Exit with appropriate code
  quit(status = ifelse(results$summary$failed > 0, 1, 0))
}