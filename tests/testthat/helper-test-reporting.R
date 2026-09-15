# Advanced test utilities for microFGT tests
library(testthat)

#' Generate test fixtures with specific properties for targeted testing
#'
#' @param test_case The specific test case to generate data for
#' @param ... Additional parameters specific to the test case
#'
#' @return A list containing test fixtures
generate_test_fixture <- function(test_case, ...) {
  switch(test_case,
    "empty" = {
      # Generate empty test data
      counts <- matrix(0, nrow = 10, ncol = 5)
      rownames(counts) <- paste0("Feature", seq_len(10))
      colnames(counts) <- paste0("Sample", seq_len(5))
      
      list(
        counts = counts,
        taxonomy = data.frame(
          Kingdom = rep("Bacteria", 10),
          Phylum = rep(NA, 10),
          Class = rep(NA, 10),
          Order = rep(NA, 10),
          Family = rep(NA, 10),
          Genus = rep(NA, 10),
          Species = rep(NA, 10),
          row.names = rownames(counts)
        ),
        metadata = data.frame(
          Sample = colnames(counts),
          Group = rep("Unknown", 5),
          row.names = colnames(counts)
        )
      )
    },
    "sparse" = {
      # Generate sparse test data
      rows <- list(...)$rows %||% 100
      cols <- list(...)$cols %||% 10
      sparsity <- list(...)$sparsity %||% 0.9
      
      # Create mostly empty matrix
      counts <- matrix(0, nrow = rows, ncol = cols)
      
      # Fill with sparse data
      idx <- sample(rows * cols, size = ceiling(rows * cols * (1 - sparsity)))
      counts[idx] <- rpois(length(idx), lambda = 20)
      
      rownames(counts) <- paste0("Feature", seq_len(rows))
      colnames(counts) <- paste0("Sample", seq_len(cols))
      
      list(
        counts = counts,
        taxonomy = data.frame(
          Kingdom = rep("Bacteria", rows),
          Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Proteobacteria", "Fusobacteria"), 
                          rows, replace = TRUE),
          Class = rep(NA, rows),
          Order = rep(NA, rows),
          Family = rep(NA, rows),
          Genus = rep(NA, rows),
          Species = rep(NA, rows),
          row.names = rownames(counts)
        ),
        metadata = data.frame(
          Sample = colnames(counts),
          Group = sample(c("Control", "Treatment"), cols, replace = TRUE),
          row.names = colnames(counts)
        )
      )
    },
    "large" = {
      # Generate large test data
      rows <- list(...)$rows %||% 1000
      cols <- list(...)$cols %||% 50
      
      # Create large matrix
      counts <- matrix(rpois(rows * cols, lambda = 10), nrow = rows, ncol = cols)
      rownames(counts) <- paste0("Feature", seq_len(rows))
      colnames(counts) <- paste0("Sample", seq_len(cols))
      
      list(
        counts = counts,
        taxonomy = data.frame(
          Kingdom = rep("Bacteria", rows),
          Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Proteobacteria", "Fusobacteria"), 
                          rows, replace = TRUE),
          Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria", "Gammaproteobacteria"), 
                         rows, replace = TRUE),
          Order = sample(LETTERS[1:10], rows, replace = TRUE),
          Family = sample(paste0("Family", 1:20), rows, replace = TRUE),
          Genus = sample(paste0("Genus", 1:30), rows, replace = TRUE),
          Species = paste0("Species", seq_len(rows)),
          row.names = rownames(counts)
        ),
        metadata = data.frame(
          Sample = colnames(counts),
          Group = sample(c("Control", "Treatment"), cols, replace = TRUE),
          Subgroup = sample(LETTERS[1:5], cols, replace = TRUE),
          Age = sample(20:80, cols, replace = TRUE),
          Value1 = rnorm(cols),
          Value2 = rnorm(cols),
          row.names = colnames(counts)
        )
      )
    },
    "CST" = {
      # Generate data with community state types
      cols <- list(...)$cols %||% 20
      
      # Define CSTs and their characteristic taxa
      csts <- c("I", "II", "III", "IV-A", "IV-B", "V")
      cst_dist <- c(0.3, 0.1, 0.2, 0.2, 0.1, 0.1)  # Distribution of CSTs
      
      # Assign CSTs to samples
      sample_csts <- sample(csts, cols, replace = TRUE, prob = cst_dist)
      
      # Define key taxa
      taxa <- c(
        "Lactobacillus_crispatus",     # Dominant in CST I
        "Lactobacillus_gasseri",       # Dominant in CST II
        "Lactobacillus_iners",         # Dominant in CST III
        "Gardnerella_vaginalis",       # Abundant in CST IV
        "Atopobium_vaginae",           # Abundant in CST IV
        "Prevotella_bivia",            # Present in CST IV
        "Sneathia_amnii",              # Present in CST IV
        "Megasphaera_sp",              # Present in CST IV
        "Lactobacillus_jensenii",      # Dominant in CST V
        "Other_bacteria"               # Background taxa
      )
      
      # Create count matrix
      counts <- matrix(0, nrow = length(taxa), ncol = cols)
      rownames(counts) <- taxa
      colnames(counts) <- paste0("Sample", seq_len(cols))
      
      # Fill in counts based on CST patterns
      for (i in seq_len(cols)) {
        cst <- sample_csts[i]
        
        if (cst == "I") {
          # CST I: L. crispatus dominated
          counts["Lactobacillus_crispatus", i] <- rbeta(1, 8, 2) * 10000
          counts["Lactobacillus_iners", i] <- rbeta(1, 1, 8) * 1000
          counts["Other_bacteria", i] <- rbeta(1, 1, 5) * 500
        } else if (cst == "II") {
          # CST II: L. gasseri dominated
          counts["Lactobacillus_gasseri", i] <- rbeta(1, 8, 2) * 10000
          counts["Lactobacillus_crispatus", i] <- rbeta(1, 1, 8) * 1000
          counts["Other_bacteria", i] <- rbeta(1, 1, 5) * 500
        } else if (cst == "III") {
          # CST III: L. iners dominated
          counts["Lactobacillus_iners", i] <- rbeta(1, 8, 2) * 10000
          counts["Gardnerella_vaginalis", i] <- rbeta(1, 2, 8) * 2000
          counts["Lactobacillus_crispatus", i] <- rbeta(1, 1, 8) * 500
          counts["Other_bacteria", i] <- rbeta(1, 1, 5) * 500
        } else if (cst == "IV-A") {
          # CST IV-A: G. vaginalis dominated diverse community
          counts["Gardnerella_vaginalis", i] <- rbeta(1, 8, 2) * 8000
          counts["Atopobium_vaginae", i] <- rbeta(1, 5, 2) * 4000
          counts["Prevotella_bivia", i] <- rbeta(1, 3, 2) * 2000
          counts["Sneathia_amnii", i] <- rbeta(1, 2, 5) * 1000
          counts["Megasphaera_sp", i] <- rbeta(1, 2, 5) * 1000
          counts["Lactobacillus_iners", i] <- rbeta(1, 1, 5) * 500
          counts["Other_bacteria", i] <- rbeta(1, 3, 3) * 2000
        } else if (cst == "IV-B") {
          # CST IV-B: A. vaginae dominated diverse community
          counts["Atopobium_vaginae", i] <- rbeta(1, 8, 2) * 8000
          counts["Gardnerella_vaginalis", i] <- rbeta(1, 5, 2) * 4000
          counts["Prevotella_bivia", i] <- rbeta(1, 5, 2) * 3000
          counts["Sneathia_amnii", i] <- rbeta(1, 3, 3) * 2000
          counts["Megasphaera_sp", i] <- rbeta(1, 3, 3) * 2000
          counts["Lactobacillus_iners", i] <- rbeta(1, 1, 8) * 200
          counts["Other_bacteria", i] <- rbeta(1, 3, 3) * 2000
        } else if (cst == "V") {
          # CST V: L. jensenii dominated
          counts["Lactobacillus_jensenii", i] <- rbeta(1, 8, 2) * 10000
          counts["Lactobacillus_iners", i] <- rbeta(1, 1, 8) * 1000
          counts["Other_bacteria", i] <- rbeta(1, 1, 5) * 500
        }
      }
      
      # Round to integers
      counts <- round(counts)
      
      # Create taxonomy data
      taxonomy <- data.frame(
        Taxon = rownames(counts),
        Genus = sapply(strsplit(rownames(counts), "_"), `[`, 1),
        Species = sapply(strsplit(rownames(counts), "_"), function(x) {
          if (length(x) > 1) paste(x[-1], collapse = "_") else NA
        }),
        row.names = rownames(counts),
        stringsAsFactors = FALSE
      )
      
      # Define higher taxonomy
      genus_to_taxonomy <- list(
        "Lactobacillus" = list(Phylum = "Firmicutes", Class = "Bacilli", Order = "Lactobacillales", Family = "Lactobacillaceae"),
        "Gardnerella" = list(Phylum = "Actinobacteria", Class = "Actinobacteria", Order = "Bifidobacteriales", Family = "Bifidobacteriaceae"),
        "Prevotella" = list(Phylum = "Bacteroidetes", Class = "Bacteroidia", Order = "Bacteroidales", Family = "Prevotellaceae"),
        "Atopobium" = list(Phylum = "Actinobacteria", Class = "Coriobacteriia", Order = "Coriobacteriales", Family = "Atopobiaceae"),
        "Sneathia" = list(Phylum = "Fusobacteria", Class = "Fusobacteriia", Order = "Fusobacteriales", Family = "Leptotrichiaceae"),
        "Megasphaera" = list(Phylum = "Firmicutes", Class = "Negativicutes", Order = "Selenomonadales", Family = "Veillonellaceae"),
        "Other" = list(Phylum = "Mixed", Class = "Mixed", Order = "Mixed", Family = "Mixed")
      )
      
      # Fill in higher taxonomy
      for (genus in unique(taxonomy$Genus)) {
        if (genus %in% names(genus_to_taxonomy)) {
          genus_ranks <- genus_to_taxonomy[[genus]]
          indices <- which(taxonomy$Genus == genus)
          for (rank in names(genus_ranks)) {
            taxonomy[indices, rank] <- genus_ranks[[rank]]
          }
        }
      }
      
      # Add Kingdom
      taxonomy$Kingdom <- "Bacteria"
      
      # Ensure columns are in proper order
      taxonomy <- taxonomy[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Taxon")]
      
      # Create sample metadata
      metadata <- data.frame(
        Sample = colnames(counts),
        CST = sample_csts,
        Nugent = c(0:10)[findInterval(runif(cols), seq(0, 1, length.out = 12))],
        pH = runif(cols, 3.8, 5.5),
        row.names = colnames(counts),
        stringsAsFactors = FALSE
      )
      
      list(
        counts = counts,
        taxonomy = taxonomy,
        metadata = metadata,
        cst = sample_csts
      )
    },
    "time_series" = {
      # Generate a time series dataset
      subjects <- list(...)$subjects %||% 5
      timepoints <- list(...)$timepoints %||% 4
      taxa <- list(...)$taxa %||% 20
      
      # Create sample IDs
      sample_ids <- c()
      subject_ids <- c()
      timepoint_ids <- c()
      
      for (s in 1:subjects) {
        for (t in 1:timepoints) {
          sample_ids <- c(sample_ids, paste0("Subject", s, "_T", t))
          subject_ids <- c(subject_ids, paste0("Subject", s))
          timepoint_ids <- c(timepoint_ids, paste0("T", t))
        }
      }
      
      # Create taxon IDs
      taxon_ids <- paste0("Taxon", 1:taxa)
      
      # Create count matrix with temporal patterns
      counts <- matrix(0, nrow = taxa, ncol = length(sample_ids))
      rownames(counts) <- taxon_ids
      colnames(counts) <- sample_ids
      
      # Function to simulate transition between communities
      simulate_transition <- function(start_state, end_state, steps) {
        t <- seq(0, 1, length.out = steps)
        states <- matrix(NA, nrow = length(start_state), ncol = steps)
        
        for (i in 1:length(start_state)) {
          # Linear interpolation with some noise
          states[i, ] <- start_state[i] + t * (end_state[i] - start_state[i]) + rnorm(steps, 0, 0.1)
        }
        
        # Ensure non-negative values
        states[states < 0] <- 0
        
        # Normalize to maintain total count
        totals <- colSums(states)
        for (j in 1:steps) {
          states[, j] <- states[, j] / totals[j] * sum(start_state)
        }
        
        return(states)
      }
      
      # Generate time series for each subject
      for (s in 1:subjects) {
        # Generate random start and end states
        start_state <- rexp(taxa, rate = 1)
        start_state <- start_state / sum(start_state) * 10000
        
        end_state <- rexp(taxa, rate = 1)
        end_state <- end_state / sum(end_state) * 10000
        
        # Simulate transition
        states <- simulate_transition(start_state, end_state, timepoints)
        
        # Fill in count matrix
        for (t in 1:timepoints) {
          sample_idx <- which(colnames(counts) == paste0("Subject", s, "_T", t))
          counts[, sample_idx] <- round(states[, t])
        }
      }
      
      # Create metadata
      metadata <- data.frame(
        Sample = sample_ids,
        Subject = subject_ids,
        Timepoint = timepoint_ids,
        TimepointNumeric = as.numeric(gsub("T", "", timepoint_ids)),
        Treatment = rep(sample(c("Control", "Treatment"), subjects, replace = TRUE), each = timepoints),
        row.names = sample_ids,
        stringsAsFactors = FALSE
      )
      
      # Create taxonomy
      taxonomy <- data.frame(
        Taxon = taxon_ids,
        Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Proteobacteria"), taxa, replace = TRUE),
        row.names = taxon_ids,
        stringsAsFactors = FALSE
      )
      
      list(
        counts = counts,
        taxonomy = taxonomy,
        metadata = metadata
      )
    },
    # Default case
    {
      stop("Unknown test case: ", test_case)
    }
  )
}

#' Create a test profile for performance tracking
#'
#' @param test_name Name of the test being profiled
#' @param expr Expression to profile
#'
#' @return A list containing timing and memory information
profile_test <- function(test_name, expr) {
  if (!requireNamespace("bench", quietly = TRUE)) {
    warning("Package 'bench' not available. Using basic profiling.")
    
    # Perform basic timing
    start_time <- Sys.time()
    result <- expr
    end_time <- Sys.time()
    
    return(list(
      test = test_name,
      time = as.numeric(end_time - start_time),
      memory = NA,
      result = result
    ))
  }
  
  # Use bench for detailed profiling
  bench_result <- bench::mark(
    result = expr,
    iterations = 1,
    check = FALSE
  )
  
  return(list(
    test = test_name,
    time = as.numeric(bench_result$time) / 1e9,  # Convert to seconds
    memory = as.numeric(bench_result$mem_alloc) / 1024^2, # Convert to MB
    result = bench_result$result[[1]]
  ))
}

#' Compare performance of multiple implementations
#'
#' @param implementations List of functions to compare
#' @param test_cases List of test cases
#' @param ... Additional arguments to pass to implementations
#'
#' @return A data frame with performance comparison
benchmark_implementations <- function(implementations, test_cases, ...) {
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("Package 'bench' is required for benchmarking")
  }
  
  results <- list()
  
  for (test_name in names(test_cases)) {
    test_case <- test_cases[[test_name]]
    
    # Run benchmarks
    benchmark <- bench::mark(
      min_iterations = 5,
      check = FALSE,
      filter_gc = FALSE,
      !!!lapply(implementations, function(fn) call(fn, test_case, ...))
    )
    
    # Extract results
    for (i in seq_along(implementations)) {
      results[[paste(test_name, names(implementations)[i], sep = "_")]] <- list(
        test = test_name,
        implementation = names(implementations)[i],
        time_median = as.numeric(benchmark$median[i]) / 1e9, # Convert to seconds
        time_min = as.numeric(benchmark$min[i]) / 1e9,
        time_max = as.numeric(benchmark$max[i]) / 1e9,
        memory = as.numeric(benchmark$mem_alloc[i]) / 1024^2, # Convert to MB
        gc_count = sum(benchmark$gc[[i]])
      )
    }
  }
  
  # Convert to data frame
  do.call(rbind, lapply(results, function(x) {
    data.frame(
      test = x$test,
      implementation = x$implementation,
      time_median_s = x$time_median,
      time_min_s = x$time_min,
      time_max_s = x$time_max,
      memory_mb = x$memory,
      gc_count = x$gc_count,
      stringsAsFactors = FALSE
    )
  }))
}

#' Check if two objects are identical within numerical tolerance
#'
#' @param a First object
#' @param b Second object
#' @param tolerance Numerical tolerance
#'
#' @return TRUE if objects are identical within tolerance, FALSE otherwise
check_equivalent <- function(a, b, tolerance = 1e-10) {
  if (is.numeric(a) && is.numeric(b)) {
    # For numeric objects, check if they're close enough
    identical(dim(a), dim(b)) && 
      all(abs(a - b) < tolerance, na.rm = TRUE) && 
      identical(is.na(a), is.na(b))
  } else {
    # For non-numeric objects, require exact equality
    identical(a, b)
  }
}

#' Track test coverage and report missing test coverage
#'
#' @param package_name Name of the package
#' @param functions Vector of function names to check
#' @param test_files Vector of test files to run
#'
#' @return A list with coverage information
check_test_coverage <- function(package_name, functions, test_files) {
  if (!requireNamespace("covr", quietly = TRUE)) {
    stop("Package 'covr' is required for coverage analysis")
  }
  
  # Get coverage for test files
  cov <- covr::package_coverage(package_name, type = "tests", test_files = test_files)
  
  # Extract coverage by function
  fn_coverage <- covr::function_coverage(cov, functions)
  
  # Format results
  results <- list(
    functions = setNames(lapply(functions, function(fn) {
      if (fn %in% names(fn_coverage)) {
        list(
          covered = fn_coverage[[fn]]$coverage,
          coverage_pct = fn_coverage[[fn]]$value * 100,
          tested = fn_coverage[[fn]]$value > 0
        )
      } else {
        list(
          covered = 0,
          coverage_pct = 0,
          tested = FALSE
        )
      }
    }), functions),
    summary = covr::tally_coverage(cov),
    missing = functions[!functions %in% names(fn_coverage) | 
                          sapply(fn_coverage[functions[functions %in% names(fn_coverage)]], 
                                 function(x) x$value == 0)]
  )
  
  return(results)
}

#' Check if a function with given name exists in a package
#'
#' @param fn_name Function name
#' @param package Package name (default: loaded namespace)
#'
#' @return TRUE if function exists, FALSE otherwise
function_exists <- function(fn_name, package = NULL) {
  if (is.null(package)) {
    exists(fn_name, mode = "function")
  } else {
    tryCatch({
      !is.null(getExportedValue(package, fn_name))
    }, error = function(e) FALSE)
  }
}

#' Safely skip a test if dependencies are not available
#'
#' @param packages Vector of package names
#' @param fns Vector of function names
#' @param message Custom skip message
#'
#' @return Invisible TRUE if test can proceed, otherwise skips test
skip_if_not_available <- function(packages = NULL, fns = NULL, message = NULL) {
  # Check packages
  if (!is.null(packages)) {
    missing_pkgs <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      skip(message %||% paste0("Required packages not available: ", 
                             paste(missing_pkgs, collapse = ", ")))
      return(invisible(FALSE))
    }
  }
  
  # Check functions
  if (!is.null(fns)) {
    missing_fns <- fns[!sapply(fns, function_exists)]
    if (length(missing_fns) > 0) {
      skip(message %||% paste0("Required functions not available: ", 
                             paste(missing_fns, collapse = ", ")))
      return(invisible(FALSE))
    }
  }
  
  invisible(TRUE)
}

# Utility function for handling NULL defaults (like %||% in purrr)
"%||%" <- function(x, y) {
  if (is.null(x)) y else x
}