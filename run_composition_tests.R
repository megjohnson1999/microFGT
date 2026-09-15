#!/usr/bin/env Rscript

# Script to run composition-based FGTExperiment tests

# Load required packages
suppressPackageStartupMessages({
  library(microFGT)
  library(testthat)
  library(TreeSummarizedExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(methods)
})

# Set options
options(testthat.summary.max_reports = 10)

# Report start time
cat("Running comprehensive FGTExperiment tests at", format(Sys.time()), "\n\n")

# Define test files to run
test_files <- c(
  "test-FGTExperiment-class-definition.R",
  "test-FGTExperiment-transformations.R",
  "test-FGTExperiment-migration.R",
  "test-mock-data-integration.R"
)

# Create a reporter
reporter <- testthat::ProgressReporter$new()

# Run tests and collect results
results <- lapply(test_files, function(test_file) {
  cat("Running tests from", test_file, "...\n")
  
  # Construct file path
  file_path <- file.path("tests", "testthat", test_file)
  
  if (!file.exists(file_path)) {
    cat("  File not found:", file_path, "\n")
    return(NULL)
  }
  
  # Run tests in the file
  results <- tryCatch({
    testthat::test_file(file_path, reporter = reporter)
  }, error = function(e) {
    cat("  Error running tests:", e$message, "\n")
    NULL
  })
  
  cat("  Done with", test_file, "\n")
  return(results)
})

# Print summary
cat("\nTest Summary:\n")
cat("-------------\n")

# Count tests and failures
total_tests <- reporter$n_ok + reporter$n_skip + reporter$n_fail + reporter$n_warn
cat("Total tests run:", total_tests, "\n")
cat("Passed:", reporter$n_ok, "\n")
cat("Skipped:", reporter$n_skip, "\n")
cat("Warnings:", reporter$n_warn, "\n")
cat("Failures:", reporter$n_fail, "\n")

# Report end time
cat("\nTests completed at", format(Sys.time()), "\n")