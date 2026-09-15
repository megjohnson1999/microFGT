#!/usr/bin/env Rscript

# Script to run the advanced taxonomic and diversity function tests for the FGTExperiment class

# Set up environment
library(testthat)
library(microFGT)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(methods)

# Define test files to run
test_files <- c(
  "tests/testthat/test-taxonomic-functions-FGTExperiment.R",
  "tests/testthat/test-diversity-functions-FGTExperiment.R",
  "tests/testthat/test-taxonomy-diversity-integration.R"
)

# Run each test file
results <- list()
for (file in test_files) {
  cat("\n\n======================================================================\n")
  cat("Running tests from:", file, "\n")
  cat("======================================================================\n\n")
  
  # Run the tests
  result <- try({
    test_file(file, reporter = SummaryReporter$new())
  })
  
  # Store result
  results[[file]] <- !inherits(result, "try-error")
  
  if (inherits(result, "try-error")) {
    cat("\nERROR in", file, ":", attr(result, "condition")$message, "\n")
  }
}

# Print summary
cat("\n\n======================================================================\n")
cat("TEST SUMMARY\n")
cat("======================================================================\n\n")

all_passed <- TRUE
for (file in names(results)) {
  status <- if (results[[file]]) "PASSED" else "FAILED"
  cat(sprintf("%-50s: %s\n", file, status))
  all_passed <- all_passed && results[[file]]
}

cat("\n")
if (all_passed) {
  cat("All tests passed successfully!\n")
} else {
  cat("Some tests failed. Please check the output above for details.\n")
}

# Exit with appropriate status code
if (!all_passed) {
  quit(status = 1)
}