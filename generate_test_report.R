#!/usr/bin/env Rscript

# Comprehensive test report generator for microFGT package
# Generates HTML reports, coverage analysis, and performance visualizations

# Check for required packages
required_packages <- c("testthat", "covr", "DT", "htmlwidgets", "ggplot2", "dplyr", "knitr", "rmarkdown")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Missing packages required for report generation:\n")
  cat(paste(" -", missing_packages, collapse = "\n"), "\n")
  cat("\nPlease install these packages using:\n")
  cat(paste0("install.packages(c('", paste(missing_packages, collapse = "', '"), "')"), "\n")
  quit(status = 1)
}

# Load required packages
suppressPackageStartupMessages({
  library(testthat)
  library(covr)
  library(DT)
  library(htmlwidgets)
  library(ggplot2)
  library(dplyr)
  library(knitr)
  library(rmarkdown)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
output_dir <- "test_reports"
run_tests <- TRUE
run_coverage <- TRUE
include_performance <- TRUE

if (length(args) > 0) {
  if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: generate_test_report.R [options]\n")
    cat("\nOptions:\n")
    cat("  --output-dir=PATH     Output directory for reports (default: test_reports)\n")
    cat("  --no-tests            Skip running tests and use existing results\n")
    cat("  --no-coverage         Skip coverage analysis\n")
    cat("  --no-performance      Skip performance tests\n")
    cat("  --help, -h            Show this help message\n")
    quit(status = 0)
  }
  
  # Parse output directory
  out_dir_arg <- grep("^--output-dir=", args, value = TRUE)
  if (length(out_dir_arg) > 0) {
    output_dir <- sub("^--output-dir=", "", out_dir_arg)
  }
  
  # Check for other flags
  if ("--no-tests" %in% args) run_tests <- FALSE
  if ("--no-coverage" %in% args) run_coverage <- FALSE
  if ("--no-performance" %in% args) include_performance <- FALSE
}

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to get package information
get_package_info <- function() {
  pkg_name <- "microFGT"
  pkg_version <- as.character(packageVersion(pkg_name))
  pkg_description <- packageDescription(pkg_name)
  
  list(
    name = pkg_name,
    version = pkg_version,
    description = pkg_description$Description,
    date = Sys.Date(),
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = Sys.info()["sysname"]
  )
}

# Run tests and collect results
run_test_suite <- function() {
  cat("Running test suite...\n")
  
  # Create a custom reporter to capture detailed test results
  DetailedReporter <- R6::R6Class("DetailedReporter",
    inherit = testthat::Reporter,
    public = list(
      results = NULL,
      
      initialize = function() {
        super$initialize()
        self$results <- list(
          passed = list(),
          failed = list(),
          skipped = list(),
          warnings = list(),
          files = list(),
          contexts = list(),
          total_time = 0
        )
      },
      
      start_file = function(file) {
        self$results$files[[file]] <- list(
          tests = 0,
          passed = 0,
          failed = 0,
          skipped = 0,
          warnings = 0,
          contexts = list(),
          start_time = Sys.time()
        )
      },
      
      end_file = function() {
        current_file <- names(self$results$files)[length(self$results$files)]
        self$results$files[[current_file]]$end_time <- Sys.time()
        self$results$files[[current_file]]$duration <- 
          difftime(self$results$files[[current_file]]$end_time, 
                 self$results$files[[current_file]]$start_time, 
                 units = "secs")
      },
      
      start_context = function(context) {
        current_file <- names(self$results$files)[length(self$results$files)]
        self$results$contexts[[context]] <- list(
          file = current_file,
          tests = list()
        )
        self$results$files[[current_file]]$contexts[[context]] <- list(
          tests = 0,
          passed = 0,
          failed = 0,
          skipped = 0,
          warnings = 0
        )
      },
      
      add_result = function(context, test, result) {
        current_file <- names(self$results$files)[length(self$results$files)]
        
        # Increment file counters
        self$results$files[[current_file]]$tests <- 
          self$results$files[[current_file]]$tests + 1
        
        # Increment context counters
        if (!context %in% names(self$results$files[[current_file]]$contexts)) {
          self$results$files[[current_file]]$contexts[[context]] <- list(
            tests = 0,
            passed = 0,
            failed = 0,
            skipped = 0,
            warnings = 0
          )
        }
        
        self$results$files[[current_file]]$contexts[[context]]$tests <- 
          self$results$files[[current_file]]$contexts[[context]]$tests + 1
        
        # Process result
        if (expectation_skip(result)) {
          self$results$skipped[[length(self$results$skipped) + 1]] <- list(
            file = current_file,
            context = context,
            test = test,
            message = result$message
          )
          self$results$files[[current_file]]$skipped <- 
            self$results$files[[current_file]]$skipped + 1
          self$results$files[[current_file]]$contexts[[context]]$skipped <- 
            self$results$files[[current_file]]$contexts[[context]]$skipped + 1
        } else if (expectation_warning(result)) {
          self$results$warnings[[length(self$results$warnings) + 1]] <- list(
            file = current_file,
            context = context,
            test = test,
            message = result$message
          )
          self$results$files[[current_file]]$warnings <- 
            self$results$files[[current_file]]$warnings + 1
          self$results$files[[current_file]]$contexts[[context]]$warnings <- 
            self$results$files[[current_file]]$contexts[[context]]$warnings + 1
          self$results$files[[current_file]]$passed <- 
            self$results$files[[current_file]]$passed + 1
          self$results$files[[current_file]]$contexts[[context]]$passed <- 
            self$results$files[[current_file]]$contexts[[context]]$passed + 1
          self$results$passed[[length(self$results$passed) + 1]] <- list(
            file = current_file,
            context = context,
            test = test
          )
        } else if (expectation_success(result)) {
          self$results$passed[[length(self$results$passed) + 1]] <- list(
            file = current_file,
            context = context,
            test = test
          )
          self$results$files[[current_file]]$passed <- 
            self$results$files[[current_file]]$passed + 1
          self$results$files[[current_file]]$contexts[[context]]$passed <- 
            self$results$files[[current_file]]$contexts[[context]]$passed + 1
        } else {
          self$results$failed[[length(self$results$failed) + 1]] <- list(
            file = current_file,
            context = context,
            test = test,
            message = result$message
          )
          self$results$files[[current_file]]$failed <- 
            self$results$files[[current_file]]$failed + 1
          self$results$files[[current_file]]$contexts[[context]]$failed <- 
            self$results$files[[current_file]]$contexts[[context]]$failed + 1
        }
        
        # Store full test info
        self$results$contexts[[context]]$tests[[test]] <- list(
          result = if (expectation_success(result)) "passed" else if (expectation_skip(result)) "skipped" else "failed",
          message = if (expectation_success(result)) NULL else result$message
        )
      },
      
      end_reporter = function() {
        self$results$total_time <- sum(sapply(self$results$files, function(f) as.numeric(f$duration)))
      }
    )
  )
  
  # Create reporter
  detailed_reporter <- DetailedReporter$new()
  
  # Run tests
  start_time <- Sys.time()
  test_dir(
    "tests/testthat",
    reporter = MultiReporter$new(list(
      SummaryReporter$new(),
      detailed_reporter
    )),
    stop_on_failure = FALSE
  )
  end_time <- Sys.time()
  
  # Add overall time
  detailed_reporter$results$start_time <- start_time
  detailed_reporter$results$end_time <- end_time
  detailed_reporter$results$total_duration <- difftime(end_time, start_time, units = "secs")
  
  # Add summary
  detailed_reporter$results$summary <- list(
    n_tests = length(detailed_reporter$results$passed) + 
             length(detailed_reporter$results$failed) + 
             length(detailed_reporter$results$skipped),
    n_passed = length(detailed_reporter$results$passed),
    n_failed = length(detailed_reporter$results$failed),
    n_skipped = length(detailed_reporter$results$skipped),
    n_warnings = length(detailed_reporter$results$warnings),
    n_files = length(detailed_reporter$results$files),
    pass_rate = length(detailed_reporter$results$passed) / 
                (length(detailed_reporter$results$passed) + 
                 length(detailed_reporter$results$failed))
  )
  
  return(detailed_reporter$results)
}

# Run code coverage analysis
calculate_coverage <- function() {
  cat("Calculating code coverage...\n")
  
  # Calculate coverage
  coverage_result <- package_coverage(".")
  
  # Generate summary
  coverage_summary <- tally_coverage(coverage_result)
  
  # Generate per-file coverage
  file_coverage <- file_coverage(coverage_result)
  
  # Generate function coverage
  function_cov <- function_coverage(coverage_result)
  
  return(list(
    coverage = coverage_result,
    summary = coverage_summary,
    file_coverage = file_coverage,
    function_coverage = function_cov
  ))
}

# Generate test summary tables
generate_test_tables <- function(test_results) {
  # File summary table
  file_summary <- data.frame(
    File = names(test_results$files),
    Tests = sapply(test_results$files, `[[`, "tests"),
    Passed = sapply(test_results$files, `[[`, "passed"),
    Failed = sapply(test_results$files, `[[`, "failed"),
    Skipped = sapply(test_results$files, `[[`, "skipped"),
    Warnings = sapply(test_results$files, `[[`, "warnings"),
    Duration = sapply(test_results$files, function(f) round(as.numeric(f$duration), 2)),
    stringsAsFactors = FALSE
  )
  
  # Context summary table
  context_summary <- data.frame(
    Context = names(test_results$contexts),
    File = sapply(test_results$contexts, `[[`, "file"),
    Tests = sapply(test_results$contexts, function(c) length(c$tests)),
    stringsAsFactors = FALSE
  )
  
  # Add pass/fail info to context summary
  context_summary$Passed <- 0
  context_summary$Failed <- 0
  context_summary$Skipped <- 0
  
  for (i in seq_len(nrow(context_summary))) {
    context <- context_summary$Context[i]
    file <- context_summary$File[i]
    
    if (context %in% names(test_results$files[[file]]$contexts)) {
      context_data <- test_results$files[[file]]$contexts[[context]]
      context_summary$Passed[i] <- context_data$passed
      context_summary$Failed[i] <- context_data$failed
      context_summary$Skipped[i] <- context_data$skipped
    }
  }
  
  # Failed tests table
  failed_tests <- data.frame(
    File = sapply(test_results$failed, `[[`, "file"),
    Context = sapply(test_results$failed, `[[`, "context"),
    Test = sapply(test_results$failed, `[[`, "test"),
    Message = sapply(test_results$failed, `[[`, "message"),
    stringsAsFactors = FALSE
  )
  
  # Skipped tests table
  skipped_tests <- data.frame(
    File = sapply(test_results$skipped, `[[`, "file"),
    Context = sapply(test_results$skipped, `[[`, "context"),
    Test = sapply(test_results$skipped, `[[`, "test"),
    Reason = sapply(test_results$skipped, `[[`, "message"),
    stringsAsFactors = FALSE
  )
  
  return(list(
    file_summary = file_summary,
    context_summary = context_summary,
    failed_tests = failed_tests,
    skipped_tests = skipped_tests
  ))
}

# Generate coverage visualizations
generate_coverage_charts <- function(coverage_results) {
  # Create coverage summary chart
  coverage_df <- data.frame(
    Category = c("Total", "R", "src"),
    Coverage = c(
      coverage_results$summary$value,
      coverage_results$summary$value_by_type[["R"]],
      coverage_results$summary$value_by_type[["src"]]
    ),
    stringsAsFactors = FALSE
  )
  
  coverage_chart <- ggplot(coverage_df, aes(x = Category, y = Coverage * 100, fill = Category)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f%%", Coverage * 100)), vjust = -0.5) +
    ylim(0, 100) +
    theme_minimal() +
    labs(
      title = "Code Coverage by Category",
      x = "",
      y = "Coverage (%)"
    ) +
    theme(legend.position = "none")
  
  # File coverage chart (top 20 most important files)
  if (nrow(coverage_results$file_coverage) > 0) {
    # Filter to R files and sort by coverage
    file_cov_df <- as.data.frame(coverage_results$file_coverage)
    file_cov_df$filename <- rownames(file_cov_df)
    file_cov_df <- file_cov_df[grep("\\.R$", file_cov_df$filename), ]
    
    # Sort by lines and get top 20
    file_cov_df <- file_cov_df[order(-file_cov_df$lines), ]
    if (nrow(file_cov_df) > 20) file_cov_df <- file_cov_df[1:20, ]
    
    # Calculate coverage percentage
    file_cov_df$coverage_pct <- file_cov_df$value * 100
    
    # Create file chart if we have data
    if (nrow(file_cov_df) > 0) {
      file_cov_df$filename <- basename(file_cov_df$filename)
      file_cov_df$filename <- factor(file_cov_df$filename, 
                                   levels = file_cov_df$filename[order(file_cov_df$coverage_pct)])
      
      file_chart <- ggplot(file_cov_df, aes(x = filename, y = coverage_pct, fill = coverage_pct)) +
        geom_col() +
        geom_text(aes(label = sprintf("%.1f%%", coverage_pct)), hjust = -0.1) +
        scale_fill_gradient(low = "lightcoral", high = "forestgreen") +
        coord_flip() +
        theme_minimal() +
        labs(
          title = "Coverage by File (Top 20 by Lines)",
          x = "",
          y = "Coverage (%)"
        ) +
        theme(legend.position = "none")
    } else {
      file_chart <- NULL
    }
  } else {
    file_chart <- NULL
  }
  
  # Function coverage chart (top 20 by complexity)
  if (length(coverage_results$function_coverage) > 0) {
    # Convert to data frame
    fn_cov_df <- data.frame(
      function_name = names(coverage_results$function_coverage),
      coverage = sapply(coverage_results$function_coverage, `[[`, "value"),
      code = sapply(coverage_results$function_coverage, `[[`, "code")
    )
    
    # Get code length as approximation of complexity
    fn_cov_df$complexity <- nchar(fn_cov_df$code)
    
    # Sort by complexity and get top 20
    fn_cov_df <- fn_cov_df[order(-fn_cov_df$complexity), ]
    if (nrow(fn_cov_df) > 20) fn_cov_df <- fn_cov_df[1:20, ]
    
    # Calculate coverage percentage
    fn_cov_df$coverage_pct <- fn_cov_df$coverage * 100
    
    # Create function chart
    if (nrow(fn_cov_df) > 0) {
      fn_cov_df$function_name <- factor(fn_cov_df$function_name, 
                                     levels = fn_cov_df$function_name[order(fn_cov_df$coverage_pct)])
      
      function_chart <- ggplot(fn_cov_df, aes(x = function_name, y = coverage_pct, fill = coverage_pct)) +
        geom_col() +
        geom_text(aes(label = sprintf("%.1f%%", coverage_pct)), hjust = -0.1) +
        scale_fill_gradient(low = "lightcoral", high = "forestgreen") +
        coord_flip() +
        theme_minimal() +
        labs(
          title = "Coverage by Function (Top 20 by Complexity)",
          x = "",
          y = "Coverage (%)"
        ) +
        theme(legend.position = "none")
    } else {
      function_chart <- NULL
    }
  } else {
    function_chart <- NULL
  }
  
  return(list(
    coverage_chart = coverage_chart,
    file_chart = file_chart,
    function_chart = function_chart
  ))
}

# Create overall package health metrics
calculate_health_metrics <- function(test_results, coverage_results = NULL) {
  # Test metrics
  test_health <- list(
    test_count = test_results$summary$n_tests,
    pass_rate = test_results$summary$pass_rate * 100,
    skip_rate = test_results$summary$n_skipped / test_results$summary$n_tests * 100,
    test_file_count = test_results$summary$n_files
  )
  
  # Coverage metrics if available
  if (!is.null(coverage_results)) {
    coverage_health <- list(
      overall_coverage = coverage_results$summary$value * 100,
      r_coverage = coverage_results$summary$value_by_type[["R"]] * 100,
      src_coverage = if ("src" %in% names(coverage_results$summary$value_by_type)) 
                     coverage_results$summary$value_by_type[["src"]] * 100 else NA
    )
  } else {
    coverage_health <- list(
      overall_coverage = NA,
      r_coverage = NA,
      src_coverage = NA
    )
  }
  
  # Calculate performance tests metrics if available
  performance_metrics <- list(
    performance_tests = 0,
    edge_case_tests = 0
  )
  
  if ("test-performance.R" %in% names(test_results$files)) {
    performance_metrics$performance_tests <- test_results$files[["test-performance.R"]]$tests
  }
  
  if ("test-performance-benchmarks.R" %in% names(test_results$files)) {
    performance_metrics$performance_tests <- performance_metrics$performance_tests + 
      test_results$files[["test-performance-benchmarks.R"]]$tests
  }
  
  if ("test-edge-cases.R" %in% names(test_results$files)) {
    performance_metrics$edge_case_tests <- test_results$files[["test-edge-cases.R"]]$tests
  }
  
  # Calculate overall health score (simplified version)
  test_score <- min(100, test_health$pass_rate) * 0.4
  coverage_score <- if (!is.null(coverage_results)) min(100, coverage_health$overall_coverage) * 0.4 else 0
  performance_score <- min(10, performance_metrics$performance_tests + performance_metrics$edge_case_tests) * 2
  
  overall_score <- test_score + coverage_score + performance_score
  if (is.null(coverage_results)) overall_score <- overall_score * 10/6  # Scale if no coverage
  
  health_metrics <- list(
    test = test_health,
    coverage = coverage_health,
    performance = performance_metrics,
    overall_score = overall_score
  )
  
  return(health_metrics)
}

# Create a health dashboard visualization
create_health_dashboard <- function(health_metrics) {
  # Create score gauges
  score_data <- data.frame(
    metric = c("Overall Health", "Test Pass Rate", "Code Coverage"),
    score = c(
      health_metrics$overall_score,
      health_metrics$test$pass_rate,
      health_metrics$coverage$overall_coverage
    ),
    stringsAsFactors = FALSE
  )
  
  # Check for NA values
  score_data$score[is.na(score_data$score)] <- 0
  
  # Define score colors
  score_data$color <- sapply(score_data$score, function(score) {
    if (score >= 90) return("forestgreen")
    if (score >= 75) return("olivedrab")
    if (score >= 60) return("gold")
    if (score >= 40) return("orange")
    return("firebrick")
  })
  
  # Create gauge plot
  gauge_plot <- ggplot(score_data, aes(x = metric, y = score, fill = color)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f", score)), position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
    scale_fill_identity() +
    scale_y_continuous(limits = c(0, 100)) +
    theme_minimal() +
    labs(
      title = "Package Health Metrics",
      x = "",
      y = "Score (0-100)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create test metrics summary
  test_metrics <- data.frame(
    Metric = c(
      "Total Tests",
      "Passing Tests",
      "Failed Tests",
      "Skipped Tests",
      "Test Files"
    ),
    Value = c(
      health_metrics$test$test_count,
      round(health_metrics$test$test_count * health_metrics$test$pass_rate / 100),
      health_metrics$test$test_count - round(health_metrics$test$test_count * health_metrics$test$pass_rate / 100) - 
        round(health_metrics$test$test_count * health_metrics$test$skip_rate / 100),
      round(health_metrics$test$test_count * health_metrics$test$skip_rate / 100),
      health_metrics$test$test_file_count
    ),
    stringsAsFactors = FALSE
  )
  
  return(list(
    gauge_plot = gauge_plot,
    test_metrics = test_metrics
  ))
}

# Generate the HTML report
generate_html_report <- function(test_results, coverage_results = NULL, health_metrics, test_tables, 
                              coverage_charts = NULL, health_dashboard, package_info) {
  # Create the report file
  report_file <- file.path(output_dir, "test_report.html")
  
  # Create report content
  report_content <- c(
    "---",
    "title: \"microFGT Package Test Report\"",
    "date: \"`r format(Sys.time(), '%B %d, %Y %H:%M')`\"",
    "output:",
    "  html_document:",
    "    theme: cosmo",
    "    toc: true",
    "    toc_float: true",
    "    code_folding: hide",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)",
    "library(DT)",
    "library(ggplot2)",
    "```",
    "",
    "## Package Information",
    "",
    paste0("- **Package**: ", package_info$name),
    paste0("- **Version**: ", package_info$version),
    paste0("- **R Version**: ", package_info$R_version),
    paste0("- **Platform**: ", package_info$platform),
    paste0("- **Report Date**: ", format(Sys.Date(), "%B %d, %Y")),
    "",
    "### Description",
    "",
    paste0(strwrap(package_info$description, width = 80), collapse = "\n"),
    "",
    "## Health Dashboard",
    "",
    "```{r health-gauge, echo=FALSE, fig.width=8, fig.height=4}",
    "health_dashboard$gauge_plot",
    "```",
    "",
    "### Key Test Metrics",
    "",
    "```{r test-metrics, echo=FALSE}",
    "knitr::kable(health_dashboard$test_metrics)",
    "```",
    "",
    "## Test Results Summary",
    "",
    paste0("- **Total Tests**: ", test_results$summary$n_tests),
    paste0("- **Passed**: ", test_results$summary$n_passed, " (", round(test_results$summary$pass_rate * 100, 1), "%)"),
    paste0("- **Failed**: ", test_results$summary$n_failed),
    paste0("- **Skipped**: ", test_results$summary$n_skipped),
    paste0("- **Test Files**: ", test_results$summary$n_files),
    paste0("- **Total Duration**: ", round(as.numeric(test_results$total_duration), 2), " seconds"),
    "",
    "### Test Files",
    "",
    "```{r file-summary, echo=FALSE}",
    "datatable(test_tables$file_summary, options = list(",
    "  pageLength = 10,",
    "  dom = 'Bfrtip',",
    "  buttons = c('copy', 'csv')",
    ")) %>%",
    "  formatStyle('Failed', ",
    "              backgroundColor = styleEqual(c(0), c('none', '#f8d7da')),",
    "              fontWeight = styleEqual(c(0), c('normal', 'bold')))",
    "```",
    "",
    "### Test Contexts",
    "",
    "```{r context-summary, echo=FALSE}",
    "datatable(test_tables$context_summary, options = list(",
    "  pageLength = 10,",
    "  dom = 'Bfrtip',",
    "  buttons = c('copy', 'csv')",
    ")) %>%",
    "  formatStyle('Failed', ",
    "              backgroundColor = styleEqual(c(0), c('none', '#f8d7da')),",
    "              fontWeight = styleEqual(c(0), c('normal', 'bold')))",
    "```"
  )
  
  # Add failed tests if any
  if (nrow(test_tables$failed_tests) > 0) {
    report_content <- c(
      report_content,
      "",
      "### Failed Tests",
      "",
      "```{r failed-tests, echo=FALSE}",
      "datatable(test_tables$failed_tests, options = list(",
      "  pageLength = 10,",
      "  dom = 'Bfrtip',",
      "  buttons = c('copy', 'csv')",
      "))",
      "```"
    )
  }
  
  # Add skipped tests if any
  if (nrow(test_tables$skipped_tests) > 0) {
    report_content <- c(
      report_content,
      "",
      "### Skipped Tests",
      "",
      "```{r skipped-tests, echo=FALSE}",
      "datatable(test_tables$skipped_tests, options = list(",
      "  pageLength = 10,",
      "  dom = 'Bfrtip',",
      "  buttons = c('copy', 'csv')",
      "))",
      "```"
    )
  }
  
  # Add coverage results if available
  if (!is.null(coverage_results)) {
    report_content <- c(
      report_content,
      "",
      "## Coverage Analysis",
      "",
      paste0("- **Overall Coverage**: ", round(coverage_results$summary$value * 100, 1), "%"),
      paste0("- **R Code Coverage**: ", round(coverage_results$summary$value_by_type[["R"]] * 100, 1), "%"),
      if ("src" %in% names(coverage_results$summary$value_by_type)) 
        paste0("- **C/C++ Code Coverage**: ", round(coverage_results$summary$value_by_type[["src"]] * 100, 1), "%") else "",
      "",
      "### Coverage Summary",
      "",
      "```{r coverage-chart, echo=FALSE, fig.width=8, fig.height=4}",
      "coverage_charts$coverage_chart",
      "```"
    )
    
    # Add file coverage chart if available
    if (!is.null(coverage_charts$file_chart)) {
      report_content <- c(
        report_content,
        "",
        "### File Coverage",
        "",
        "```{r file-coverage-chart, echo=FALSE, fig.width=10, fig.height=8}",
        "coverage_charts$file_chart",
        "```"
      )
    }
    
    # Add function coverage chart if available
    if (!is.null(coverage_charts$function_chart)) {
      report_content <- c(
        report_content,
        "",
        "### Function Coverage",
        "",
        "```{r function-coverage-chart, echo=FALSE, fig.width=10, fig.height=8}",
        "coverage_charts$function_chart",
        "```"
      )
    }
  }
  
  # Add performance section if applicable
  if (include_performance) {
    report_content <- c(
      report_content,
      "",
      "## Performance Tests",
      "",
      paste0("- **Performance Tests**: ", health_metrics$performance$performance_tests),
      paste0("- **Edge Case Tests**: ", health_metrics$performance$edge_case_tests),
      "",
      "Performance tests are critical for ensuring the package maintains optimal speed and memory usage, especially with larger datasets."
    )
    
    # If we have performance test files, add more info
    perf_files <- names(test_results$files)[grep("performance", names(test_results$files))]
    if (length(perf_files) > 0) {
      # Extract performance test descriptions
      perf_tests <- data.frame(
        File = character(),
        Test = character(),
        Status = character(),
        stringsAsFactors = FALSE
      )
      
      for (file in perf_files) {
        for (context in names(test_results$files[[file]]$contexts)) {
          context_data <- test_results$contexts[[context]]
          for (test_name in names(context_data$tests)) {
            perf_tests <- rbind(perf_tests, data.frame(
              File = file,
              Test = test_name,
              Status = context_data$tests[[test_name]]$result,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      if (nrow(perf_tests) > 0) {
        report_content <- c(
          report_content,
          "",
          "### Performance Test Details",
          "",
          "```{r perf-tests, echo=FALSE}",
          "datatable(perf_tests, options = list(",
          "  pageLength = 10,",
          "  dom = 'Bfrtip',",
          "  buttons = c('copy', 'csv')",
          ")) %>%",
          "  formatStyle('Status', ",
          "              backgroundColor = styleEqual(c('passed', 'failed', 'skipped'), ",
          "                                           c('#d4edda', '#f8d7da', '#fff3cd')))",
          "```"
        )
      }
    }
  }
  
  # Write the report file
  writeLines(report_content, report_file)
  
  # Render the report
  rmarkdown::render(
    report_file,
    output_file = "test_report.html",
    output_dir = output_dir,
    quiet = TRUE
  )
  
  # Return the report file path
  file.path(output_dir, "test_report.html")
}

# Main function to run everything
main <- function() {
  package_info <- get_package_info()
  
  # Run tests
  if (run_tests) {
    test_results <- run_test_suite()
  } else {
    test_results <- list(
      summary = list(
        n_tests = 0,
        n_passed = 0,
        n_failed = 0,
        n_skipped = 0,
        n_warnings = 0,
        n_files = 0,
        pass_rate = 0
      ),
      files = list(),
      contexts = list(),
      failed = list(),
      skipped = list(),
      warnings = list(),
      passed = list(),
      total_duration = 0
    )
  }
  
  # Run coverage analysis
  if (run_coverage) {
    coverage_results <- calculate_coverage()
    # Save coverage report
    coverage_report_file <- file.path(output_dir, "coverage_report.html")
    covr::report(coverage_results$coverage, file = coverage_report_file)
  } else {
    coverage_results <- NULL
  }
  
  # Generate test tables
  test_tables <- generate_test_tables(test_results)
  
  # Generate coverage charts if coverage was run
  if (!is.null(coverage_results)) {
    coverage_charts <- generate_coverage_charts(coverage_results)
  } else {
    coverage_charts <- NULL
  }
  
  # Calculate health metrics
  health_metrics <- calculate_health_metrics(test_results, coverage_results)
  
  # Create health dashboard
  health_dashboard <- create_health_dashboard(health_metrics)
  
  # Generate HTML report
  report_file <- generate_html_report(
    test_results, 
    coverage_results, 
    health_metrics,
    test_tables,
    coverage_charts,
    health_dashboard,
    package_info
  )
  
  cat("\nTest report generated at:", report_file, "\n")
  
  if (run_coverage) {
    cat("Coverage report generated at:", file.path(output_dir, "coverage_report.html"), "\n")
  }
  
  return(invisible(NULL))
}

# Run the main function
main()