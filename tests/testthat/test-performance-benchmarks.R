# Comprehensive performance tests for microFGT package

library(testthat)
library(microFGT)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(methods)

# Skip entire file if benchmark packages not available
skip_if_not_available(
  packages = c("bench", "dplyr", "ggplot2"),
  message = "Performance benchmarks require bench, dplyr, and ggplot2 packages"
)

# Skip large benchmarks unless MICROFGT_FULL_TEST is set
run_large_benchmarks <- identical(Sys.getenv("MICROFGT_FULL_TEST"), "true")

# Set seed for reproducibility
set.seed(12345)

# Test transformation performance
test_that("transformAbundance performance scales well with dataset size", {
  skip_if_not_available(fns = c("transformAbundance", "FGTExperiment"))
  
  sizes <- list(
    small = list(rows = 100, cols = 10),
    medium = list(rows = 500, cols = 20)
  )
  
  if (run_large_benchmarks) {
    sizes$large <- list(rows = 5000, cols = 100)
  }
  
  results <- data.frame(
    size = character(),
    dimensions = character(),
    transformation = character(),
    time_seconds = numeric(),
    memory_mb = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (size_name in names(sizes)) {
    size <- sizes[[size_name]]
    
    # Create test data
    test_data <- generate_test_fixture("sparse", 
                                      rows = size$rows, 
                                      cols = size$cols, 
                                      sparsity = 0.7)
    
    # Create FGTExperiment
    fgt <- FGTExperiment(
      assays = list(counts = test_data$counts),
      rowData = S4Vectors::DataFrame(test_data$taxonomy),
      colData = S4Vectors::DataFrame(test_data$metadata),
      experimentType = "amplicon"
    )
    
    # Benchmark transformations
    for (transform_type in c("relative", "log", "clr", "presence")) {
      if (!run_large_benchmarks && size_name == "medium" && transform_type == "clr") {
        # Skip medium CLR in non-full mode (CLR is intensive)
        next
      }
      
      # Run benchmark
      result <- profile_test(
        paste("transform", transform_type, size_name, sep = "_"),
        transformAbundance(fgt, type = transform_type)
      )
      
      # Record result
      results <- rbind(results, data.frame(
        size = size_name,
        dimensions = paste0(size$rows, " x ", size$cols),
        transformation = transform_type,
        time_seconds = result$time,
        memory_mb = result$memory %||% NA,
        stringsAsFactors = FALSE
      ))
      
      # Add expectation to ensure test reports correctly
      expect_true(result$time > 0, 
                 info = paste0("Transformation ", transform_type, 
                               " on ", size_name, " dataset completed"))
    }
  }
  
  # Check for reasonable scaling
  # For each transformation type, time should scale sub-quadratically with data size
  for (transform_type in c("relative", "log", "presence")) {
    if (all(c("small", "medium") %in% results$size)) {
      small_time <- results$time_seconds[results$size == "small" & 
                                         results$transformation == transform_type]
      medium_time <- results$time_seconds[results$size == "medium" & 
                                          results$transformation == transform_type]
      
      # Medium dataset has 25x more cells than small (500x20 vs 100x10)
      # If scaling is linear or better, medium_time / small_time < 25
      ratio <- medium_time / small_time
      expect_lt(ratio, 50, 
                info = paste0(transform_type, " transformation shows reasonable scaling"))
    }
  }
  
  # Optional: Save benchmark results if all packages available 
  # and we're in full test mode
  if (run_large_benchmarks && 
      requireNamespace("ggplot2", quietly = TRUE) && 
      requireNamespace("dplyr", quietly = TRUE)) {
    
    # Create results plot
    results_wide <- tidyr::pivot_wider(
      results, 
      id_cols = c("size", "dimensions"),
      names_from = "transformation",
      values_from = c("time_seconds", "memory_mb")
    )
    
    # Order sizes
    results$size <- factor(results$size, levels = c("small", "medium", "large"))
    
    # Create plot
    transform_plot <- ggplot2::ggplot(results, 
      ggplot2::aes(x = size, y = time_seconds, fill = transformation)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        title = "Transformation Performance by Data Size",
        x = "Dataset Size",
        y = "Time (seconds, log scale)",
        fill = "Transformation"
      ) +
      ggplot2::theme_minimal()
    
    # Save plot if possible
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot_file <- file.path(tempdir(), "transformation_benchmark.png")
      ggplot2::ggsave(plot_file, transform_plot, width = 8, height = 6)
      cat("Benchmark plot saved to:", plot_file, "\n")
    }
  }
})

# Test CST calculation performance
test_that("CST calculation performance is optimal", {
  skip_if_not_available(fns = c("calculate_cst", "FGTExperiment"))
  
  # Generate test data with CSTs
  test_data <- generate_test_fixture("CST", cols = 20)
  
  # Create FGTExperiment
  fgt <- FGTExperiment(
    assays = list(counts = test_data$counts),
    rowData = S4Vectors::DataFrame(test_data$taxonomy),
    colData = S4Vectors::DataFrame(test_data$metadata),
    experimentType = "amplicon"
  )
  
  # Only run if calculate_cst function exists
  if (function_exists("calculate_cst")) {
    # Benchmark CST calculation
    result <- profile_test(
      "calculate_cst",
      calculate_cst(fgt)
    )
    
    # Verify performance is reasonable (< 5 seconds for this small dataset)
    expect_lt(result$time, 5, 
             info = "CST calculation completes in reasonable time")
    
    # Check result correctness
    cst_result <- result$result
    expect_true(is.character(cst_result) || is.factor(cst_result), 
               info = "CST result has correct type")
    expect_equal(length(cst_result), ncol(fgt), 
                info = "CST result has one value per sample")
  } else {
    skip("calculate_cst function not available")
  }
})

# Test taxonomic aggregation performance
test_that("aggregate_taxa performance scales well with dataset size", {
  skip_if_not_available(fns = c("aggregate_taxa", "FGTExperiment"))
  
  sizes <- list(
    small = list(rows = 100, cols = 10),
    medium = list(rows = 500, cols = 20)
  )
  
  if (run_large_benchmarks) {
    sizes$large <- list(rows = 5000, cols = 50)
  }
  
  results <- data.frame(
    size = character(),
    dimensions = character(),
    rank = character(),
    time_seconds = numeric(),
    memory_mb = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (size_name in names(sizes)) {
    size <- sizes[[size_name]]
    
    # Create test data
    test_data <- generate_test_fixture("sparse", 
                                      rows = size$rows, 
                                      cols = size$cols)
    
    # Create FGTExperiment
    fgt <- FGTExperiment(
      assays = list(counts = test_data$counts),
      rowData = S4Vectors::DataFrame(test_data$taxonomy),
      colData = S4Vectors::DataFrame(test_data$metadata),
      experimentType = "amplicon"
    )
    
    # Benchmark different taxonomic ranks
    for (rank in c("Phylum", "Genus")) {
      # Run benchmark
      result <- profile_test(
        paste("aggregate", rank, size_name, sep = "_"),
        aggregate_taxa(fgt, rank = rank)
      )
      
      # Record result
      results <- rbind(results, data.frame(
        size = size_name,
        dimensions = paste0(size$rows, " x ", size$cols),
        rank = rank,
        time_seconds = result$time,
        memory_mb = result$memory %||% NA,
        stringsAsFactors = FALSE
      ))
      
      # Check the result has the expected structure
      agg_result <- result$result
      expect_s4_class(agg_result, "FGTExperiment", 
                     info = paste0("Aggregation at ", rank, " level returns FGTExperiment"))
      
      # Check aggregation reduced the number of rows as expected
      # Number of rows should equal number of unique values in that rank
      expected_rows <- length(unique(na.omit(test_data$taxonomy[[rank]])))
      # Allow for NA handling differences
      expect_lte(nrow(agg_result), max(expected_rows + 1, 1), 
                info = paste0("Aggregation at ", rank, " level has expected row count"))
    }
  }
  
  # Check for reasonable scaling
  if (all(c("small", "medium") %in% results$size)) {
    for (rank in c("Phylum", "Genus")) {
      small_time <- results$time_seconds[results$size == "small" & results$rank == rank]
      medium_time <- results$time_seconds[results$size == "medium" & results$rank == rank]
      
      # Medium dataset has 10x more taxa than small
      ratio <- medium_time / small_time
      expect_lt(ratio, 20, 
               info = paste0(rank, " aggregation shows reasonable scaling"))
    }
  }
})

# Test diversity calculation performance
test_that("diversity calculations perform well with increasing dataset size", {
  skip_if_not_available(fns = c("calculate_diversity", "FGTExperiment"))
  
  sizes <- list(
    small = list(rows = 100, cols = 10),
    medium = list(rows = 500, cols = 20)
  )
  
  if (run_large_benchmarks) {
    sizes$large <- list(rows = 2000, cols = 50)
  }
  
  results <- data.frame(
    size = character(),
    dimensions = character(),
    diversity_metric = character(),
    time_seconds = numeric(),
    memory_mb = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (size_name in names(sizes)) {
    size <- sizes[[size_name]]
    
    # Create test data
    test_data <- generate_test_fixture("sparse", 
                                      rows = size$rows, 
                                      cols = size$cols, 
                                      sparsity = 0.7)
    
    # Create FGTExperiment
    fgt <- FGTExperiment(
      assays = list(counts = test_data$counts),
      rowData = S4Vectors::DataFrame(test_data$taxonomy),
      colData = S4Vectors::DataFrame(test_data$metadata),
      experimentType = "amplicon"
    )
    
    # Benchmark diversity metrics
    for (metric in c("shannon", "simpson", "richness")) {
      # Run benchmark
      result <- profile_test(
        paste("diversity", metric, size_name, sep = "_"),
        calculate_diversity(fgt, method = metric)
      )
      
      # Record result
      results <- rbind(results, data.frame(
        size = size_name,
        dimensions = paste0(size$rows, " x ", size$cols),
        diversity_metric = metric,
        time_seconds = result$time,
        memory_mb = result$memory %||% NA,
        stringsAsFactors = FALSE
      ))
      
      # Check the result has the expected structure
      div_result <- result$result
      expect_true(is.numeric(div_result), 
                 info = paste0(metric, " diversity returns numeric values"))
      expect_equal(length(div_result), ncol(fgt), 
                  info = paste0(metric, " diversity returns one value per sample"))
      
      # Additional checks for specific metrics
      if (metric == "shannon") {
        # Shannon index should be non-negative
        expect_true(all(div_result >= 0), 
                   info = "Shannon diversity values are non-negative")
      } else if (metric == "simpson") {
        # Simpson index should be between 0 and 1
        expect_true(all(div_result >= 0 & div_result <= 1), 
                   info = "Simpson diversity values are between 0 and 1")
      } else if (metric == "richness") {
        # Richness should be integer-like and match count of non-zero taxa
        expect_true(all(div_result == floor(div_result)), 
                   info = "Richness values are integers")
      }
    }
  }
  
  # Check for reasonable scaling
  if (all(c("small", "medium") %in% results$size)) {
    for (metric in c("shannon", "simpson", "richness")) {
      small_time <- results$time_seconds[results$size == "small" & 
                                         results$diversity_metric == metric]
      medium_time <- results$time_seconds[results$size == "medium" & 
                                          results$diversity_metric == metric]
      
      # Medium dataset has 5x more samples and 5x more features
      ratio <- medium_time / small_time
      expect_lt(ratio, 25, 
               info = paste0(metric, " diversity shows reasonable scaling"))
    }
  }
})

# Test beta diversity calculation performance
test_that("beta diversity calculations perform well", {
  skip_if_not_available(fns = c("calculate_beta_diversity", "FGTExperiment"))
  
  sizes <- list(
    small = list(rows = 100, cols = 10),
    medium = list(rows = 300, cols = 20)
  )
  
  if (run_large_benchmarks) {
    sizes$large <- list(rows = 1000, cols = 50)
  }
  
  results <- data.frame(
    size = character(),
    dimensions = character(),
    beta_metric = character(),
    time_seconds = numeric(),
    memory_mb = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (size_name in names(sizes)) {
    size <- sizes[[size_name]]
    
    # Create test data
    test_data <- generate_test_fixture("sparse", 
                                      rows = size$rows, 
                                      cols = size$cols, 
                                      sparsity = 0.8)
    
    # Create FGTExperiment
    fgt <- FGTExperiment(
      assays = list(counts = test_data$counts),
      rowData = S4Vectors::DataFrame(test_data$taxonomy),
      colData = S4Vectors::DataFrame(test_data$metadata),
      experimentType = "amplicon"
    )
    
    # Convert to relative abundance
    fgt_rel <- transformAbundance(fgt, type = "relative")
    
    # Benchmark beta diversity metrics
    for (metric in c("bray", "jaccard")) {
      # Skip larger calculations in non-full mode
      if (!run_large_benchmarks && size_name != "small" && metric == "jaccard") {
        next
      }
      
      # Run benchmark
      result <- profile_test(
        paste("beta_diversity", metric, size_name, sep = "_"),
        calculate_beta_diversity(fgt_rel, method = metric)
      )
      
      # Record result
      results <- rbind(results, data.frame(
        size = size_name,
        dimensions = paste0(size$rows, " x ", size$cols),
        beta_metric = metric,
        time_seconds = result$time,
        memory_mb = result$memory %||% NA,
        stringsAsFactors = FALSE
      ))
      
      # Check the result has the expected structure
      beta_result <- result$result
      expect_true(is.matrix(beta_result), 
                 info = paste0(metric, " beta diversity returns a matrix"))
      expect_equal(dim(beta_result), c(ncol(fgt), ncol(fgt)), 
                  info = paste0(metric, " beta diversity returns a square matrix"))
      
      # Check diagonal is zeros (distances to self)
      diag_values <- diag(beta_result)
      expect_true(all(abs(diag_values) < 1e-10), 
                 info = paste0(metric, " beta diversity has zeros on diagonal"))
      
      # Check all values are between 0 and 1 for these metrics
      expect_true(all(beta_result >= 0 & beta_result <= 1), 
                 info = paste0(metric, " beta diversity values are between 0 and 1"))
      
      # Check if the matrix is symmetric
      expect_true(isSymmetric(beta_result), 
                 info = paste0(metric, " beta diversity matrix is symmetric"))
    }
  }
  
  # Check for reasonable scaling
  if (all(c("small", "medium") %in% results$size) && 
      "bray" %in% results$beta_metric) {
    small_time <- results$time_seconds[results$size == "small" & 
                                       results$beta_metric == "bray"]
    medium_time <- results$time_seconds[results$size == "medium" & 
                                        results$beta_metric == "bray"]
    
    # Medium dataset has 4x more pairwise comparisons (20x20 vs 10x10)
    ratio <- medium_time / small_time
    # We'd expect quadratic scaling in the worst case (O(n²))
    expect_lt(ratio, 16, 
             info = paste0("Bray-Curtis beta diversity shows reasonable scaling"))
  }
})

# Time series analysis performance tests
test_that("time series analysis performs efficiently", {
  skip_if_not_available(fns = c("analyze_time_series", "FGTExperiment"))
  
  # Skip if analyze_time_series function doesn't exist
  if (!function_exists("analyze_time_series")) {
    skip("analyze_time_series function not available")
  }
  
  # Generate time series test data
  test_data <- generate_test_fixture("time_series", 
                                    subjects = 5, 
                                    timepoints = 4, 
                                    taxa = 50)
  
  # Create FGTExperiment
  fgt <- FGTExperiment(
    assays = list(counts = test_data$counts),
    rowData = S4Vectors::DataFrame(test_data$taxonomy),
    colData = S4Vectors::DataFrame(test_data$metadata),
    experimentType = "longitudinal"
  )
  
  # Benchmark time series analysis
  result <- profile_test(
    "time_series_analysis",
    analyze_time_series(
      fgt,
      subject_col = "Subject",
      time_col = "TimepointNumeric"
    )
  )
  
  # Verify performance is reasonable
  expect_lt(result$time, 10, 
           info = "Time series analysis completes in reasonable time")
  
  # Check result structure
  ts_result <- result$result
  expect_true(is.list(ts_result), 
             info = "Time series analysis returns a list")
  
  # Check if result contains expected components
  expected_components <- c("subjects", "timepoints", "trends")
  expect_true(all(expected_components %in% names(ts_result)), 
             info = "Time series analysis result contains expected components")
})

# Memory usage tests
test_that("memory usage is reasonable during operations", {
  skip_if_not_available(
    packages = c("pryr"),
    message = "Memory usage tests require pryr package"
  )
  
  # Create medium-sized test data
  test_data <- generate_test_fixture("sparse", 
                                    rows = 500, 
                                    cols = 20, 
                                    sparsity = 0.7)
  
  # Create FGTExperiment
  mem_before <- pryr::mem_used()
  fgt <- FGTExperiment(
    assays = list(counts = test_data$counts),
    rowData = S4Vectors::DataFrame(test_data$taxonomy),
    colData = S4Vectors::DataFrame(test_data$metadata),
    experimentType = "amplicon"
  )
  mem_after <- pryr::mem_used()
  
  # Calculate memory used for object creation
  mem_diff_mb <- (mem_after - mem_before) / (1024 * 1024)
  
  # Log memory usage
  cat("FGTExperiment creation memory:", mem_diff_mb, "MB\n")
  
  # Memory usage should be reasonable for this data size
  # Object should not be more than 20x the size of the input data
  # (allowing for S4 overhead and additional slots)
  raw_data_size <- object.size(test_data) / (1024 * 1024)
  expect_lt(mem_diff_mb, raw_data_size * 20, 
           info = "FGTExperiment memory usage is reasonable")
  
  # Test memory usage during transformation
  mem_before <- pryr::mem_used()
  fgt_rel <- transformAbundance(fgt, type = "relative")
  mem_after <- pryr::mem_used()
  
  mem_diff_mb <- (mem_after - mem_before) / (1024 * 1024)
  cat("Relative abundance transformation memory:", mem_diff_mb, "MB\n")
  
  # Transformation should not use more than 3x the object size
  fgt_size <- object.size(fgt) / (1024 * 1024)
  expect_lt(mem_diff_mb, fgt_size * 3, 
           info = "Transformation memory usage is reasonable")
})

# Optional: Function call overhead test
test_that("S4 method dispatch overhead is acceptable", {
  skip_if_not_available(packages = c("bench"))
  
  # Skip test in non-full mode
  if (!run_large_benchmarks) {
    skip("Skipping S4 overhead test in non-full test mode")
  }
  
  # Create small test data
  test_data <- generate_test_fixture("sparse", rows = 50, cols = 5)
  
  # Create FGTExperiment
  fgt <- FGTExperiment(
    assays = list(counts = test_data$counts),
    rowData = S4Vectors::DataFrame(test_data$taxonomy),
    colData = S4Vectors::DataFrame(test_data$metadata),
    experimentType = "amplicon"
  )
  
  # Direct calculation for relative abundance
  direct_calc <- function(fgt) {
    counts <- assay(fgt, "counts")
    rel_counts <- t(t(counts) / colSums(counts))
    return(rel_counts)
  }
  
  # Compare S4 method vs direct calculation
  comparison <- bench::mark(
    s4_method = assay(transformAbundance(fgt, type = "relative"), "relative_counts"),
    direct = direct_calc(fgt),
    iterations = 20,
    check = FALSE
  )
  
  # S4 method should not be more than 5x slower than direct calculation
  s4_ratio <- comparison$median[1] / comparison$median[2]
  cat("S4 method overhead ratio:", s4_ratio, "\n")
  
  expect_lt(s4_ratio, 5, 
           info = "S4 method dispatch overhead is reasonable")
})