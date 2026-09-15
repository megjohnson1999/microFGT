#!/usr/bin/env Rscript

# Enhanced Master Script to run all test suites for the microFGT package
# Provides comprehensive test reporting, coverage analysis, and performance profiling

# Setup environment
suppressPackageStartupMessages({
  library(testthat)
  library(microFGT)
  library(TreeSummarizedExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(methods)
})

# Try to load optional packages for enhanced testing
has_covr <- requireNamespace("covr", quietly = TRUE)
has_profvis <- requireNamespace("profvis", quietly = TRUE)
has_crayon <- requireNamespace("crayon", quietly = TRUE)
has_DT <- requireNamespace("DT", quietly = TRUE)
has_htmlwidgets <- requireNamespace("htmlwidgets", quietly = TRUE)

# Setup colored output if available
if (has_crayon) {
  header <- crayon::bold
  success <- crayon::green
  warning <- crayon::yellow
  error <- crayon::red
  info <- crayon::blue
} else {
  header <- warning <- error <- success <- info <- function(x) x
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
run_core <- TRUE
run_integration <- TRUE
run_performance <- TRUE
run_documentation <- TRUE
run_coverage <- FALSE
run_profile <- FALSE
output_report <- FALSE
only_failed_tests <- FALSE

if (length(args) > 0) {
  if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: run_all_tests_enhanced.R [options]\n")
    cat("\nOptions:\n")
    cat("  --core-only          Run only core functionality tests\n")
    cat("  --integration-only   Run only integration tests\n")
    cat("  --performance-only   Run only performance and robustness tests\n")
    cat("  --docs-only          Run only documentation and interface tests\n")
    cat("  --full               Run all tests in full mode (may take a long time)\n")
    cat("  --coverage           Calculate test coverage (requires 'covr' package)\n")
    cat("  --profile            Profile test performance (requires 'profvis' package)\n")
    cat("  --report             Generate HTML test report (requires 'DT' and 'htmlwidgets')\n")
    cat("  --failed-only        After initial run, rerun only failed tests\n")
    cat("  --help, -h           Show this help message\n")
    quit(status = 0)
  }
  
  if ("--core-only" %in% args) {
    run_integration <- FALSE
    run_performance <- FALSE
    run_documentation <- FALSE
  } else if ("--integration-only" %in% args) {
    run_core <- FALSE
    run_performance <- FALSE
    run_documentation <- FALSE
  } else if ("--performance-only" %in% args) {
    run_core <- FALSE
    run_integration <- FALSE
    run_documentation <- FALSE
  } else if ("--docs-only" %in% args) {
    run_core <- FALSE
    run_integration <- FALSE
    run_performance <- FALSE
  }
  
  if ("--coverage" %in% args) {
    if (!has_covr) {
      cat(error("Error: Package 'covr' is required for coverage analysis. Please install it.\n"))
      quit(status = 1)
    }
    run_coverage <- TRUE
  }
  
  if ("--profile" %in% args) {
    if (!has_profvis) {
      cat(error("Error: Package 'profvis' is required for profiling. Please install it.\n"))
      quit(status = 1)
    }
    run_profile <- TRUE
  }
  
  if ("--report" %in% args) {
    if (!has_DT || !has_htmlwidgets) {
      cat(error("Error: Packages 'DT' and 'htmlwidgets' are required for HTML reports. Please install them.\n"))
      quit(status = 1)
    }
    output_report <- TRUE
  }
  
  if ("--failed-only" %in% args) {
    only_failed_tests <- TRUE
  }
  
  if ("--full" %in% args) {
    Sys.setenv(MICROBIOME_FULL_TEST = "true")
    Sys.setenv(MICROFGT_FULL_TEST = "true")
  }
}

# Print test plan
cat("\n", header("======================================================================"), "\n")
cat(header("microFGT COMPREHENSIVE TEST SUITE (ENHANCED)"), "\n")
cat(header("======================================================================"), "\n\n")

cat(info("Test plan:"), "\n")
if (run_core) cat("- Core functionality tests\n")
if (run_integration) cat("- Integration tests\n")
if (run_performance) cat("- Performance and robustness tests\n")
if (run_documentation) cat("- Documentation and interface tests\n")
if (run_coverage) cat("- " + info("Test coverage analysis"), "\n")
if (run_profile) cat("- " + info("Test performance profiling"), "\n")
if (output_report) cat("- " + info("Generating HTML test report"), "\n")

if (identical(Sys.getenv("MICROBIOME_FULL_TEST"), "true") && 
    identical(Sys.getenv("MICROFGT_FULL_TEST"), "true")) {
  cat("\n", warning("Running in FULL TEST mode. This may take a long time."), "\n")
} else {
  cat("\n", info("Running in LIMITED TEST mode. Use --full for comprehensive testing."), "\n")
}

cat("\n")

# Define test categories with metadata for better organization
test_categories <- list(
  core = list(
    name = "Core Functionality",
    description = "Basic functionality tests for FGTExperiment class and methods",
    files = c(
      "tests/testthat/test-FGTExperiment-class-definition.R",
      "tests/testthat/test-FGTExperiment-transformations.R",
      "tests/testthat/test-taxonomic-functions-FGTExperiment.R",
      "tests/testthat/test-diversity-functions-FGTExperiment.R",
      "tests/testthat/test-taxonomy-diversity-integration.R"
    )
  ),
  integration = list(
    name = "Integration",
    description = "Tests for component integration and end-to-end workflows",
    files = c(
      "tests/testthat/test-end-to-end-workflows.R",
      "tests/testthat/test-multi-component-integration.R"
    )
  ),
  performance = list(
    name = "Performance and Robustness",
    description = "Tests for performance, edge cases, and validation",
    files = c(
      "tests/testthat/test-performance.R",
      "tests/testthat/test-edge-cases.R",
      "tests/testthat/test-validation.R"
    )
  ),
  documentation = list(
    name = "Documentation and Interface",
    description = "Tests for documentation examples and interface consistency",
    files = c(
      "tests/testthat/test-documentation-examples.R",
      "tests/testthat/test-interface-consistency.R"
    )
  )
)

# Storage for test results
test_results <- list(
  passed = list(),
  failed = list(),
  skipped = list(),
  warnings = list(),
  timings = list(),
  errors = list()
)

# Custom reporter to capture more detailed results
DetailedReporter <- R6::R6Class("DetailedReporter",
  inherit = testthat::Reporter,
  public = list(
    test_results = NULL,
    test_file = NULL,
    
    initialize = function() {
      super$initialize()
      self$test_results <- list(
        passed = list(),
        failed = list(),
        skipped = list(),
        warnings = list()
      )
    },
    
    start_file = function(file) {
      self$test_file <- file
    },
    
    add_result = function(context, test, result) {
      if (expectation_skip(result)) {
        self$test_results$skipped[[length(self$test_results$skipped) + 1]] <- list(
          file = self$test_file,
          context = context,
          test = test,
          message = result$message
        )
      } else if (expectation_warning(result)) {
        self$test_results$warnings[[length(self$test_results$warnings) + 1]] <- list(
          file = self$test_file,
          context = context,
          test = test,
          message = result$message
        )
      } else if (expectation_success(result)) {
        self$test_results$passed[[length(self$test_results$passed) + 1]] <- list(
          file = self$test_file,
          context = context,
          test = test
        )
      } else {
        self$test_results$failed[[length(self$test_results$failed) + 1]] <- list(
          file = self$test_file,
          context = context,
          test = test,
          message = result$message
        )
      }
    }
  )
)

# Storage for all detailed test results
all_results <- list()

# Helper function to run test files with enhanced reporting
run_test_category <- function(category_info, rerun_failed = FALSE, failed_tests = NULL) {
  cat("\n\n", header("======================================================================"), "\n")
  cat(header(paste0("RUNNING ", toupper(category_info$name), " TESTS")), "\n")
  cat(header("======================================================================"), "\n")
  cat(info(category_info$description), "\n\n")
  
  category_results <- list(
    files = list(),
    timings = list(),
    passed = 0,
    failed = 0,
    skipped = 0,
    warnings = 0,
    detailed = list()
  )
  
  files_to_run <- category_info$files
  
  # If rerunning failed tests, filter to only those
  if (rerun_failed && !is.null(failed_tests)) {
    files_to_run <- intersect(files_to_run, failed_tests)
    if (length(files_to_run) == 0) {
      cat(info("No failed tests to rerun in this category.\n"))
      return(NULL)
    }
  }
  
  for (file in files_to_run) {
    # Skip files that don't exist
    if (!file.exists(file)) {
      cat(warning(paste("Skipping non-existent file:", file)), "\n")
      next
    }
    
    cat("\n", info(paste("Running tests from:", file)), "\n")
    cat("----------------------------------------------------------------------\n")
    
    # Create detailed reporter
    detailed_reporter <- DetailedReporter$new()
    
    # Run the tests and time them
    start_time <- Sys.time()
    result <- try({
      if (run_profile) {
        profile <- profvis::profvis({
          test_file(file, reporter = MultiReporter$new(list(
            SummaryReporter$new(),
            detailed_reporter
          )))
        }, interval = 0.01)
        
        # Save profile for this file
        profvis_file <- paste0("test_profile_", basename(file), ".html")
        htmlwidgets::saveWidget(profile, profvis_file)
        cat(info(paste("Performance profile saved to:", profvis_file)), "\n")
      } else {
        test_file(file, reporter = MultiReporter$new(list(
          SummaryReporter$new(),
          detailed_reporter
        )))
      }
    })
    end_time <- Sys.time()
    
    # Calculate timing
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
    category_results$timings[[file]] <- elapsed
    
    # Store result
    if (inherits(result, "try-error")) {
      category_results$failed <- category_results$failed + 1
      category_results$files[[file]] <- FALSE
      category_results$detailed[[file]] <- list(
        error = attr(result, "condition")$message,
        detailed_results = detailed_reporter$test_results
      )
      cat("\n", error(paste("ERROR in", file, ":", attr(result, "condition")$message)), "\n")
    } else {
      # Count results
      n_passed <- length(detailed_reporter$test_results$passed)
      n_failed <- length(detailed_reporter$test_results$failed)
      n_skipped <- length(detailed_reporter$test_results$skipped)
      n_warnings <- length(detailed_reporter$test_results$warnings)
      
      # Store results
      category_results$passed <- category_results$passed + n_passed
      category_results$failed <- category_results$failed + n_failed
      category_results$skipped <- category_results$skipped + n_skipped
      category_results$warnings <- category_results$warnings + n_warnings
      category_results$files[[file]] <- n_failed == 0
      category_results$detailed[[file]] <- detailed_reporter$test_results
      
      # Print summary for this file
      cat("\n", info("Results:"), " ")
      cat(success(paste0(n_passed, " passed")), ", ")
      if (n_failed > 0) cat(error(paste0(n_failed, " failed")), ", ")
      if (n_skipped > 0) cat(warning(paste0(n_skipped, " skipped")), ", ")
      if (n_warnings > 0) cat(warning(paste0(n_warnings, " warnings")), ", ")
      cat(info(paste0("(", format(elapsed, digits = 2), " seconds)")), "\n")
      
      # If there were failures, print details
      if (n_failed > 0) {
        cat("\n", error("Failed tests:"), "\n")
        for (failure in detailed_reporter$test_results$failed) {
          cat("  - ", error(failure$test), ": ", failure$message, "\n")
        }
      }
    }
  }
  
  return(category_results)
}

# Run initial test suite
all_results <- list()
failed_files <- c()

# Function to run all categories
run_all_categories <- function(rerun_failed = FALSE, failed_files = NULL) {
  results <- list()
  
  # Run each category of tests
  if (run_core) {
    results$core <- run_test_category(test_categories$core, rerun_failed, failed_files)
  }
  
  if (run_integration) {
    results$integration <- run_test_category(test_categories$integration, rerun_failed, failed_files)
  }
  
  if (run_performance) {
    results$performance <- run_test_category(test_categories$performance, rerun_failed, failed_files)
  }
  
  if (run_documentation) {
    results$documentation <- run_test_category(test_categories$documentation, rerun_failed, failed_files)
  }
  
  return(results)
}

# Run the initial test suite
all_results <- run_all_categories(FALSE)

# Collect failed files for rerun
if (only_failed_tests) {
  failed_files <- c()
  for (category in names(all_results)) {
    if (!is.null(all_results[[category]])) {
      for (file in names(all_results[[category]]$files)) {
        if (!all_results[[category]]$files[[file]]) {
          failed_files <- c(failed_files, file)
        }
      }
    }
  }
  
  if (length(failed_files) > 0) {
    cat("\n\n", header("======================================================================"), "\n")
    cat(header("RERUNNING FAILED TESTS"), "\n")
    cat(header("======================================================================"), "\n\n")
    
    # Rerun just the failed tests
    rerun_results <- run_all_categories(TRUE, failed_files)
    
    # Update results for files that were rerun
    for (category in names(rerun_results)) {
      if (!is.null(rerun_results[[category]])) {
        for (file in names(rerun_results[[category]]$files)) {
          all_results[[category]]$files[[file]] <- rerun_results[[category]]$files[[file]]
          all_results[[category]]$detailed[[file]] <- rerun_results[[category]]$detailed[[file]]
        }
      }
    }
  }
}

# Run test coverage if requested
if (run_coverage) {
  cat("\n\n", header("======================================================================"), "\n")
  cat(header("CALCULATING TEST COVERAGE"), "\n")
  cat(header("======================================================================"), "\n\n")
  
  # Determine which test files to include based on what was run
  coverage_files <- c()
  if (run_core) coverage_files <- c(coverage_files, test_categories$core$files)
  if (run_integration) coverage_files <- c(coverage_files, test_categories$integration$files)
  if (run_performance) coverage_files <- c(coverage_files, test_categories$performance$files)
  if (run_documentation) coverage_files <- c(coverage_files, test_categories$documentation$files)
  
  cat(info("Calculating coverage for:"), "\n")
  for (file in coverage_files) {
    cat("  - ", file, "\n")
  }
  
  # Calculate coverage
  coverage_result <- try({
    package_coverage <- covr::package_coverage(
      ".",
      type = "tests",
      test_files = coverage_files
    )
    
    # Print coverage summary
    cat("\n", info("Coverage Summary:"), "\n")
    print(covr::tally_coverage(package_coverage))
    
    # Save detailed HTML report
    coverage_report <- "coverage_report.html"
    covr::report(package_coverage, file = coverage_report)
    cat("\n", info(paste("Detailed coverage report saved to:", coverage_report)), "\n")
    
    # Return coverage data
    list(
      coverage = package_coverage,
      summary = covr::tally_coverage(package_coverage)
    )
  })
  
  if (inherits(coverage_result, "try-error")) {
    cat("\n", error("Error calculating coverage:"), "\n")
    cat(error(attr(coverage_result, "condition")$message), "\n")
    all_results$coverage <- NULL
  } else {
    all_results$coverage <- coverage_result
  }
}

# Generate HTML report if requested
if (output_report) {
  cat("\n\n", header("======================================================================"), "\n")
  cat(header("GENERATING HTML TEST REPORT"), "\n")
  cat(header("======================================================================"), "\n\n")
  
  # Create report data
  report_data <- data.frame(
    Category = character(),
    File = character(),
    Status = character(),
    Passed = numeric(),
    Failed = numeric(),
    Skipped = numeric(),
    Warnings = numeric(),
    Time_Seconds = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Collect data from all categories
  for (category_name in names(all_results)) {
    category <- all_results[[category_name]]
    if (is.null(category) || category_name == "coverage") next
    
    category_display_name <- test_categories[[category_name]]$name
    
    for (file in names(category$files)) {
      # Get counts for this file
      detailed <- category$detailed[[file]]
      n_passed <- if(is.list(detailed)) length(detailed$passed) else 0
      n_failed <- if(is.list(detailed)) length(detailed$failed) else 0
      n_skipped <- if(is.list(detailed)) length(detailed$skipped) else 0
      n_warnings <- if(is.list(detailed)) length(detailed$warnings) else 0
      
      # Get status
      if (is.logical(category$files[[file]])) {
        status <- if(category$files[[file]]) "Passed" else "Failed"
      } else {
        status <- "Error"
      }
      
      # Add to report data
      report_data <- rbind(report_data, data.frame(
        Category = category_display_name,
        File = file,
        Status = status,
        Passed = n_passed,
        Failed = n_failed,
        Skipped = n_skipped,
        Warnings = n_warnings,
        Time_Seconds = category$timings[[file]],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Create HTML widget
  report_widget <- DT::datatable(
    report_data,
    caption = paste("microFGT Test Results -", Sys.time()),
    rownames = FALSE,
    options = list(
      pageLength = 50,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    )
  ) %>%
    DT::formatStyle(
      'Status',
      backgroundColor = DT::styleEqual(
        c('Passed', 'Failed', 'Error'),
        c('#dff0d8', '#f2dede', '#fcf8e3')
      )
    ) %>%
    DT::formatRound('Time_Seconds', 2)
  
  # Save as HTML
  report_file <- "test_report.html"
  htmlwidgets::saveWidget(report_widget, report_file)
  cat(info(paste("Test report saved to:", report_file)), "\n")
}

# Print overall summary
cat("\n\n", header("======================================================================"), "\n")
cat(header("OVERALL TEST SUMMARY"), "\n")
cat(header("======================================================================"), "\n\n")

all_passed <- TRUE
total_passed <- 0
total_failed <- 0
total_skipped <- 0
total_warnings <- 0
total_time <- 0

for (category_name in names(all_results)) {
  category <- all_results[[category_name]]
  if (is.null(category) || category_name == "coverage") next
  
  category_display_name <- test_categories[[category_name]]$name
  
  # Calculate category time
  category_time <- sum(unlist(category$timings), na.rm = TRUE)
  total_time <- total_time + category_time
  
  # Calculate category status
  files_status <- unlist(category$files)
  category_passed <- length(files_status) > 0 && all(files_status)
  
  # Print category summary
  status_str <- if(category_passed) success("PASSED") else error("FAILED")
  time_str <- sprintf("(%.2f seconds)", category_time)
  cat(sprintf("%-25s: %s %s\n", 
              category_display_name, 
              status_str,
              info(time_str)))
  
  # Print category stats
  cat(sprintf("  Tests: %s passed, %s failed, %s skipped, %s warnings\n",
              success(category$passed),
              if(category$failed > 0) error(category$failed) else category$failed,
              if(category$skipped > 0) warning(category$skipped) else category$skipped,
              if(category$warnings > 0) warning(category$warnings) else category$warnings))
  
  # Add to totals
  total_passed <- total_passed + category$passed
  total_failed <- total_failed + category$failed
  total_skipped <- total_skipped + category$skipped
  total_warnings <- total_warnings + category$warnings
  
  # If not passed, list failed files
  if (!category_passed) {
    cat("  Failed files:\n")
    for (file in names(category$files)) {
      if (!is.logical(category$files[[file]]) || !category$files[[file]]) {
        cat("  - ", error(file), "\n")
        # Add first few failure messages if available
        if (is.list(category$detailed[[file]]) && 
            is.list(category$detailed[[file]]$failed) && 
            length(category$detailed[[file]]$failed) > 0) {
          for (i in 1:min(3, length(category$detailed[[file]]$failed))) {
            failure <- category$detailed[[file]]$failed[[i]]
            cat("    * ", error(failure$test), ": ", failure$message, "\n")
          }
          if (length(category$detailed[[file]]$failed) > 3) {
            cat("    * ", warning("...and ", length(category$detailed[[file]]$failed) - 3, " more failures"), "\n")
          }
        }
      }
    }
  }
  
  all_passed <- all_passed && category_passed
}

# Print totals
cat("\n", header("TOTALS:"), "\n")
cat(sprintf("  %s passed, %s failed, %s skipped, %s warnings\n",
            success(total_passed),
            if(total_failed > 0) error(total_failed) else total_failed,
            if(total_skipped > 0) warning(total_skipped) else total_skipped,
            if(total_warnings > 0) warning(total_warnings) else total_warnings))
cat(sprintf("  Total execution time: %s seconds\n", info(sprintf("%.2f", total_time))))

cat("\n")
if (all_passed) {
  cat(success("All test suites passed successfully!"), "\n")
} else {
  cat(error("Some test suites failed. See details above."), "\n")
}

# Additional resources
if (output_report) {
  cat("\n", info("Additional resources:"), "\n")
  cat("  - HTML Test Report: ", info("test_report.html"), "\n")
}

if (run_coverage) {
  cat("  - Coverage Report: ", info("coverage_report.html"), "\n")
}

if (run_profile) {
  cat("  - Performance profiles are saved as: ", info("test_profile_*.html"), "\n")
}

# Exit with appropriate status code
if (!all_passed) {
  quit(status = 1)
}