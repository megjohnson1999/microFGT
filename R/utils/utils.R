#' Download and install speciateIT
#'
#' Downloads and installs the speciateIT taxonomic classifier tool.
#'
#' @param force Force reinstallation if already installed
#' @param version Version to install
#' @param dest_dir Directory to install speciateIT to
#'
#' @return Logical indicating success
#' @export
install_speciateit <- function(force = FALSE, version = "latest", dest_dir = NULL) {
  # Check if already installed and force is FALSE
  if (!force && is_speciateit_available()) {
    message("speciateIT is already installed. Use force=TRUE to reinstall.")
    return(TRUE)
  }
  
  # This is a placeholder function that would need to be implemented
  # based on how speciateIT is actually distributed and installed
  message("Installation functionality not implemented yet.")
  message("Please manually install speciateIT following the instructions at:")
  message("https://github.com/username/speciateit")
  
  return(FALSE)
}

#' Download and configure speciateIT reference database
#'
#' Downloads and configures the reference database for speciateIT.
#'
#' @param destination Directory to store the database
#' @param force Force download even if database already exists
#'
#' @return Path to the installed database
#' @export
setup_speciateit_database <- function(destination = file.path(system.file(package = "microFGT"), 
                                                         "extdata", "reference_dbs", 
                                                         "speciateit"),
                                    force = FALSE) {
  # Check if destination directory exists
  if (!dir.exists(dirname(destination))) {
    dir.create(dirname(destination), recursive = TRUE)
  }
  
  # Check if database already exists
  if (!force && dir.exists(destination) && length(list.files(destination)) > 0) {
    message("speciateIT database already exists at: ", destination)
    message("Use force=TRUE to re-download.")
    return(destination)
  }
  
  # This is a placeholder function that would need to be implemented
  # based on how the speciateIT database is actually distributed
  message("Database download functionality not implemented yet.")
  message("Please manually download the database following the instructions at:")
  message("https://github.com/username/speciateit")
  
  # Set as default database if download succeeded
  options(microFGT.speciateit_db = destination)
  
  return(destination)
}

#' Normalize file path for cross-platform compatibility
#'
#' @param path File path to normalize
#'
#' @return Normalized path
#' @keywords internal
normalize_path <- function(path) {
  if (is.null(path)) return(NULL)
  
  # Convert to absolute path if not already
  if (!file.exists(path) && !dir.exists(path)) {
    # Return as is if file doesn't exist
    return(path)
  }
  
  # Get absolute path
  abs_path <- normalizePath(path, mustWork = FALSE)
  
  # Normalize path separators for cross-platform compatibility
  norm_path <- gsub("\\\\", "/", abs_path)
  
  return(norm_path)
}

#' Check if a package is installed
#'
#' @param package Name of the package to check
#' @param min_version Minimum required version (optional)
#'
#' @return Logical indicating if the package is installed (and meets version requirement)
#' @keywords internal
check_package <- function(package, min_version = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return(FALSE)
  }
  
  if (!is.null(min_version)) {
    pkg_version <- as.character(utils::packageVersion(package))
    if (utils::compareVersion(pkg_version, min_version) < 0) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Require packages with helpful error messages
#'
#' @param ... Package names to check
#' @param quietly Suppress messages
#'
#' @return Invisible logical vector indicating which packages are available
#' @keywords internal
require_packages <- function(..., quietly = FALSE) {
  pkgs <- c(...)
  status <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status) && !quietly) {
    missing_pkgs <- pkgs[!status]
    msg <- paste0(
      "The following packages are required but not installed: ",
      paste(missing_pkgs, collapse = ", "), ".\n",
      "Please install with:\n",
      "install.packages(c('", paste(missing_pkgs, collapse = "', '"), "'))"
    )
    
    # Check if missing packages are from Bioconductor
    bioc_pkgs <- c("TreeSummarizedExperiment", "SummarizedExperiment", 
                 "MultiAssayExperiment", "S4Vectors", "Biostrings", "phyloseq")
    
    missing_bioc <- intersect(missing_pkgs, bioc_pkgs)
    if (length(missing_bioc) > 0) {
      msg <- paste0(
        msg, "\n\n",
        "For Bioconductor packages (", paste(missing_bioc, collapse = ", "), "), use:\n",
        "if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n",
        "    install.packages(\"BiocManager\")\n",
        "BiocManager::install(c('", paste(missing_bioc, collapse = "', '"), "'))"
      )
    }
    
    message(msg)
  }
  
  invisible(status)
}