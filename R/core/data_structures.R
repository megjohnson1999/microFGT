#' Data structure utility functions
#'
#' This file contains utility functions for working with FGTExperiment and related data structures.
#' The FGTExperiment class itself is defined in FGTExperiment-class.R.
#'
#' @name data_structures
#' @keywords internal
NULL

#' Check if an object is an FGTExperiment
#'
#' @param x Object to check
#'
#' @return Logical indicating if x is an FGTExperiment
#' @export
is_fgt_experiment <- function(x) {
  inherits(x, "FGTExperiment")
}

#' Create a new FGTExperiment with modified data
#'
#' This function creates a new FGTExperiment by modifying an existing one.
#' It's useful for internal operations that need to return modified objects.
#'
#' @param fgt_exp Original FGTExperiment
#' @param new_assays List of new assays (optional)
#' @param new_rowData New rowData (optional)
#' @param new_colData New colData (optional)
#' @param new_rowTree New rowTree (optional)
#' @param new_experimentType New experimentType (optional)
#' @param new_fgtMetadata New fgtMetadata (optional)
#'
#' @return A new FGTExperiment with updated data
#' @keywords internal
update_fgt_experiment <- function(fgt_exp, 
                                new_assays = NULL,
                                new_rowData = NULL,
                                new_colData = NULL,
                                new_rowTree = NULL,
                                new_experimentType = NULL,
                                new_fgtMetadata = NULL) {
  
  # Check input
  if (!is_fgt_experiment(fgt_exp)) {
    stop(format_error("Object is not an FGTExperiment", "update_fgt_experiment"))
  }
  
  # Create a copy of the original
  result <- fgt_exp
  
  # Update assays if provided
  if (!is.null(new_assays)) {
    # Implementation depends on specifics of how to replace assays
    # This is just a placeholder
    for (name in names(new_assays)) {
      assay(result, name) <- new_assays[[name]]
    }
  }
  
  # Update rowData if provided
  if (!is.null(new_rowData)) {
    rowData(result) <- new_rowData
  }
  
  # Update colData if provided
  if (!is.null(new_colData)) {
    colData(result) <- new_colData
  }
  
  # Update rowTree if provided
  if (!is.null(new_rowTree)) {
    rowTree(result) <- new_rowTree
  }
  
  # Update experimentType if provided
  if (!is.null(new_experimentType)) {
    experimentType(result) <- new_experimentType
  }
  
  # Update fgtMetadata if provided
  if (!is.null(new_fgtMetadata)) {
    fgtMetadata(result) <- new_fgtMetadata
  }
  
  return(result)
}