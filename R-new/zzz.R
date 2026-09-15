#' Package initialization function
#'
#' This function is called when the package is loaded with library() or require()
#'
#' @param libname Library name
#' @param pkgname Package name
#'
#' @return None
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Add any package-level initialization here
  
  # Register S4 class validity methods
  # Already done in AllClasses.R, but this ensures it happens at load time
  setValidity("FGTExperiment", .validate_FGTExperiment)
  
  invisible()
}

#' Package attachment function
#'
#' This function is called when the package is attached with library() or require()
#'
#' @param libname Library name
#' @param pkgname Package name
#'
#' @return None
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  # Ensure required packages are available
  required_packages <- c(
    "methods",
    "S4Vectors",
    "SummarizedExperiment",
    "TreeSummarizedExperiment"
  )
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    packageStartupMessage("Warning: The following required packages are missing: ", 
                          paste(missing_packages, collapse = ", "), 
                          ". Some functionality may not work.")
  }
  
  # Add startup message
  packageStartupMessage("microFGT: Comprehensive Analysis of Female Genital Tract Microbiome Data")
}

#' Namespace unload function
#'
#' This function is called when the package namespace is unloaded
#'
#' @param libpath Library path
#'
#' @return None
#' @keywords internal
.onUnload <- function(libpath) {
  # Clean up any resources here
  invisible()
}