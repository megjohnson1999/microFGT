#' @import methods
#' @importFrom SummarizedExperiment SummarizedExperiment assay assays assayNames
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
NULL

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
  # Register S4 class validity methods if needed
  # This ensures validation happens at load time
  if (exists(".validate_FGTExperiment", envir = asNamespace(pkgname), inherits = FALSE)) {
    setValidity("FGTExperiment", get(".validate_FGTExperiment", envir = asNamespace(pkgname)))
  }
  
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
  # Check for required packages
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
  
  # Display package information
  version <- utils::packageVersion("microFGT")
  msg <- paste0(
    "\n",
    "╔══════════════════════════════════════════════════════════════╗\n",
    "║ microFGT ", version, " - Female Genital Tract Microbiome Analysis  ║\n",
    "╠══════════════════════════════════════════════════════════════╣\n",
    "║ Comprehensive tools for amplicon and metagenomic FGT analysis ║\n",
    "║ Type 'browseVignettes(\"microFGT\")' for documentation          ║\n",
    "╚══════════════════════════════════════════════════════════════╝\n"
  )
  
  packageStartupMessage(msg)
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
  # Clean up any resources if needed
  invisible()
}