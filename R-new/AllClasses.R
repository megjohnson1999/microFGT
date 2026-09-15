#' @import methods
NULL

#' FGTExperiment Class
#'
#' A S4 class extending TreeSummarizedExperiment to provide specialized
#' functionality for female genital tract microbiome data.
#'
#' @slot experimentType Character string indicating the experiment type
#' @slot fgtMetadata SimpleList containing additional metadata
#'
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
#' @exportClass FGTExperiment
setClass("FGTExperiment",
         contains = "TreeSummarizedExperiment",
         slots = list(
           experimentType = "character",
           fgtMetadata = "SimpleList"
         ),
         prototype = list(
           experimentType = "amplicon",
           fgtMetadata = S4Vectors::SimpleList()
         ))

#' Validity method for FGTExperiment class
#'
#' @param object FGTExperiment object to validate
#' @return TRUE if object is valid, otherwise a character vector of error messages
#' @keywords internal
.validate_FGTExperiment <- function(object) {
  msg <- NULL
  
  # Check experimentType
  validTypes <- c("amplicon", "metagenomic", "integrated") 
  if (!(object@experimentType %in% validTypes)) {
    msg <- c(msg, paste("experimentType must be one of:", 
                       paste(validTypes, collapse=", ")))
  }
  
  # Check fgtMetadata
  if (!is(object@fgtMetadata, "SimpleList")) {
    msg <- c(msg, "fgtMetadata must be a SimpleList object")
  }
  
  # Return messages or TRUE
  if (is.null(msg)) TRUE else msg
}

# Register the validity method
setValidity("FGTExperiment", .validate_FGTExperiment)