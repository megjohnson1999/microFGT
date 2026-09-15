#' Legacy constructor for FGTExperiment objects
#'
#' @details
#' This file is maintained for backward compatibility. It provides an alias to the
#' main FGTExperiment constructor. Please use FGTExperiment directly instead of 
#' this function in new code.
#'
#' @keywords internal
NULL

#' Fixed Constructor for FGTExperiment objects (Deprecated)
#'
#' @description
#' This function is maintained for backward compatibility. Please use the 
#' standard `FGTExperiment()` constructor in all new code.
#'
#' @param assays List of matrices or similar objects containing count or abundance data
#' @param rowData DataFrame of feature metadata (e.g., taxonomic classifications)
#' @param colData DataFrame of sample metadata
#' @param rowTree phylo object or NULL for feature tree (e.g., phylogenetic tree)
#' @param experimentType Type of experiment ("amplicon", "metagenomic", "integrated")
#' @param ... Additional arguments passed to TreeSummarizedExperiment constructor
#'
#' @return An FGTExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a simple FGTExperiment (deprecated method)
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' 
#' # Create simple FGTExperiment
#' fgt_exp <- create_FGTExperiment(
#'   assays = list(counts = counts)
#' )
#' }
create_FGTExperiment <- function(assays, rowData = NULL, colData = NULL, rowTree = NULL,
                              experimentType = "amplicon", ...) {
  
  # Display deprecation warning
  warning("create_FGTExperiment() is deprecated. Please use FGTExperiment() instead.")
  
  # Forward to the standard constructor
  FGTExperiment(assays = assays, 
               rowData = rowData, 
               colData = colData, 
               rowTree = rowTree,
               experimentType = experimentType, 
               ...)