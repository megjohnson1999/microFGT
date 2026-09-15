#' Validate NAMESPACE Exports
#'
#' Checks that all required exports are present in the NAMESPACE file,
#' including S4 class exports, method exports, and function exports.
#'
#' @param namespace_file Path to the NAMESPACE file (default: "NAMESPACE")
#' @param verbose Logical indicating whether to print detailed information
#'
#' @return Logical indicating if all checks passed
#' @export
validate_namespace <- function(namespace_file = "NAMESPACE", verbose = TRUE) {
  # Check if file exists
  if (!file.exists(namespace_file)) {
    if (verbose) cat("✗ NAMESPACE file not found:", namespace_file, "\n")
    return(FALSE)
  }
  
  # Read NAMESPACE content
  ns_content <- readLines(namespace_file)
  if (verbose) cat("Reading NAMESPACE file:", namespace_file, "\n")
  
  result <- TRUE
  
  # Required S4 class exports
  s4_classes <- c("FGTExperiment")
  
  if (verbose) cat("\nChecking S4 class exports:\n")
  for (class_name in s4_classes) {
    pattern <- paste0("^exportClasses\\(", class_name, "\\)")
    if (!any(grepl(pattern, ns_content))) {
      if (verbose) cat("✗ Missing exportClasses(", class_name, ") in NAMESPACE\n", sep = "")
      result <- FALSE
    } else {
      if (verbose) cat("✓ Found exportClasses(", class_name, ") in NAMESPACE\n", sep = "")
    }
  }
  
  # Required constructor exports
  constructors <- c("FGTExperiment")
  
  if (verbose) cat("\nChecking constructor exports:\n")
  for (constructor in constructors) {
    pattern <- paste0("^export\\(", constructor, "\\)")
    if (!any(grepl(pattern, ns_content))) {
      if (verbose) cat("✗ Missing export(", constructor, ") in NAMESPACE\n", sep = "")
      result <- FALSE
    } else {
      if (verbose) cat("✓ Found export(", constructor, ") in NAMESPACE\n", sep = "")
    }
  }
  
  # Required method exports - add more as needed
  methods <- list(
    c("experimentType", FALSE),      # Regular name
    c("experimentType<-", TRUE),     # Needs backticks
    c("fgtMetadata", FALSE),         # Regular name
    c("fgtMetadata<-", TRUE)         # Needs backticks
  )
  
  if (verbose) cat("\nChecking method exports:\n")
  for (method_info in methods) {
    method_name <- method_info[1]
    needs_backticks <- as.logical(method_info[2])
    
    if (needs_backticks) {
      pattern <- paste0("^export\\(`", method_name, "`\\)")
    } else {
      pattern <- paste0("^export\\(", method_name, "\\)")
    }
    
    if (!any(grepl(pattern, ns_content))) {
      if (verbose) cat("✗ Missing export for method", method_name, "in NAMESPACE\n")
      result <- FALSE
    } else {
      if (verbose) cat("✓ Found export for method", method_name, "in NAMESPACE\n")
    }
  }
  
  # Required function exports - add core functions here
  functions <- c(
    "transform_abundance",
    "filter_taxa",
    "classify_with_speciateit",
    "default_speciateit_db",
    "is_speciateit_available",
    "run_speciateit"
  )
  
  if (verbose) cat("\nChecking function exports:\n")
  for (function_name in functions) {
    pattern <- paste0("^export\\(", function_name, "\\)")
    if (!any(grepl(pattern, ns_content))) {
      if (verbose) cat("✗ Missing export(", function_name, ") in NAMESPACE\n", sep = "")
      result <- FALSE
    } else {
      if (verbose) cat("✓ Found export(", function_name, ") in NAMESPACE\n", sep = "")
    }
  }
  
  # Required imports for S4 functionality
  imports <- c(
    "^import\\(methods\\)",
    "^importClassesFrom\\(TreeSummarizedExperiment,TreeSummarizedExperiment\\)",
    "^importFrom\\(S4Vectors,SimpleList\\)"
  )
  
  if (verbose) cat("\nChecking essential imports:\n")
  for (import in imports) {
    if (!any(grepl(import, ns_content))) {
      if (verbose) cat("✗ Missing required import:", import, "in NAMESPACE\n")
      result <- FALSE
    } else {
      if (verbose) cat("✓ Found required import:", import, "in NAMESPACE\n")
    }
  }
  
  # Return overall result
  if (verbose) {
    if (result) {
      cat("\n✓ All NAMESPACE checks passed\n")
    } else {
      cat("\n✗ Some NAMESPACE checks failed\n")
    }
  }
  
  return(result)
}

# If script is executed directly, run the validation
if (identical(environment(), globalenv())) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Set defaults
  namespace_file <- "NAMESPACE"
  verbose <- TRUE
  
  # Process arguments
  if (length(args) > 0) {
    namespace_file <- args[1]
  }
  
  # Run validation
  result <- validate_namespace(namespace_file, verbose)
  
  # Set exit code based on result
  if (!result) {
    quit(status = 1)
  }
}