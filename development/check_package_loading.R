#' Check Package Loading and Dependencies
#'
#' Tests if the package loads correctly and if all dependencies are available.
#' This helps catch issues early in the CI pipeline.
#'
#' @param package_name Name of the package to check
#' @param verbose Logical indicating whether to print detailed information
#'
#' @return Logical indicating if all checks passed
#' @export
check_package_loading <- function(package_name = "microFGT", verbose = TRUE) {
  result <- TRUE
  
  # Check if package can be loaded
  if (verbose) cat("Checking if package", package_name, "can be loaded...\n")
  pkg_loaded <- FALSE
  
  tryCatch({
    # Try to load the package
    library(package_name, character.only = TRUE)
    pkg_loaded <- TRUE
    if (verbose) cat("✓ Package", package_name, "loaded successfully\n")
  }, error = function(e) {
    if (verbose) cat("✗ Failed to load package:", conditionMessage(e), "\n")
    result <- FALSE
  })
  
  if (!pkg_loaded) {
    return(result)
  }
  
  # Get package description
  pkg_desc <- utils::packageDescription(package_name)
  if (is.null(pkg_desc)) {
    if (verbose) cat("✗ Could not get package description\n")
    return(FALSE)
  }
  
  # Check package version
  if (verbose) cat("Package version:", pkg_desc$Version, "\n")
  
  # Extract dependencies
  depends <- pkg_desc$Depends
  imports <- pkg_desc$Imports
  suggests <- pkg_desc$Suggests
  
  # Parse dependencies
  parse_deps <- function(deps) {
    if (is.null(deps)) return(character(0))
    
    # Split by comma
    deps_list <- strsplit(deps, ",")[[1]]
    
    # Clean up
    deps_list <- gsub("\\s*\\(.*\\)", "", deps_list)  # Remove version requirements
    deps_list <- gsub("^\\s+|\\s+$", "", deps_list)   # Trim whitespace
    
    # Remove R dependency
    deps_list <- deps_list[!grepl("^R$", deps_list)]
    
    return(deps_list)
  }
  
  # Get all dependencies
  all_depends <- parse_deps(depends)
  all_imports <- parse_deps(imports)
  all_suggests <- parse_deps(suggests)
  
  if (verbose) {
    cat("\nDependencies:\n")
    cat("- Depends:", paste(all_depends, collapse = ", "), "\n")
    cat("- Imports:", paste(all_imports, collapse = ", "), "\n")
    cat("- Suggests:", paste(all_suggests, collapse = ", "), "\n")
  }
  
  # Check required dependencies (Depends + Imports)
  required_deps <- c(all_depends, all_imports)
  
  if (verbose) cat("\nChecking required dependencies:\n")
  for (dep in required_deps) {
    if (dep == "") next
    
    # Check if package is installed
    is_installed <- requireNamespace(dep, quietly = TRUE)
    
    if (!is_installed) {
      if (verbose) cat("✗ Required dependency", dep, "is not installed\n")
      result <- FALSE
    } else {
      if (verbose) cat("✓ Required dependency", dep, "is installed\n")
      
      # Try to load the package to ensure it works
      tryCatch({
        # Load namespace without attaching
        loadNamespace(dep)
        if (verbose) cat("  ✓ Dependency", dep, "namespace loaded successfully\n")
      }, error = function(e) {
        if (verbose) cat("  ✗ Error loading dependency", dep, "namespace:", conditionMessage(e), "\n")
        result <- FALSE
      })
    }
  }
  
  # Check suggested dependencies
  if (verbose) cat("\nChecking suggested dependencies:\n")
  for (dep in all_suggests) {
    if (dep == "") next
    
    # Check if package is installed
    is_installed <- requireNamespace(dep, quietly = TRUE)
    
    if (!is_installed) {
      if (verbose) cat("! Suggested dependency", dep, "is not installed\n")
    } else {
      if (verbose) cat("✓ Suggested dependency", dep, "is installed\n")
    }
  }
  
  # Check if imported functions are available
  # This is a basic check that tries to get a few functions from each imported package
  if (verbose) cat("\nChecking critical imports:\n")
  
  # List of critical imports to check
  critical_imports <- list(
    methods = c("is", "isClass", "new"),
    S4Vectors = c("SimpleList"),
    SummarizedExperiment = c("SummarizedExperiment"),
    TreeSummarizedExperiment = c("TreeSummarizedExperiment")
  )
  
  for (pkg in names(critical_imports)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (verbose) cat("✗ Critical package", pkg, "is not available\n")
      result <- FALSE
      next
    }
    
    for (fun in critical_imports[[pkg]]) {
      fun_exists <- exists(fun, envir = asNamespace(pkg), mode = "function") ||
                   exists(fun, envir = asNamespace(pkg))
      
      if (!fun_exists) {
        if (verbose) cat("✗ Function", fun, "from package", pkg, "not found\n")
        result <- FALSE
      } else {
        if (verbose) cat("✓ Function", fun, "from package", pkg, "is available\n")
      }
    }
  }
  
  # Return overall result
  if (verbose) {
    if (result) {
      cat("\n✓ All package loading checks passed\n")
    } else {
      cat("\n✗ Some package loading checks failed\n")
    }
  }
  
  return(result)
}

# If script is executed directly, run the check
if (identical(environment(), globalenv())) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Set defaults
  package_name <- "microFGT"
  verbose <- TRUE
  
  # Process arguments
  if (length(args) > 0) {
    package_name <- args[1]
  }
  
  # Run check
  result <- check_package_loading(package_name, verbose)
  
  # Set exit code based on result
  if (!result) {
    quit(status = 1)
  }
}