#' FGTExperiment Class
#'
#' The FGTExperiment class provides specialized functionality for female genital tract (FGT) 
#' microbiome data analysis. It contains a TreeSummarizedExperiment object and adds slots 
#' for experiment type and custom metadata.
#'
#' @slot experimentData Internal TreeSummarizedExperiment object
#' @slot experimentType Character string indicating the type of experiment
#'   ("amplicon", "metagenomic", or "integrated").
#' @slot fgtMetadata SimpleList containing experiment-specific metadata.
#'
#' @name FGTExperiment-class
#' @rdname FGTExperiment-class
#' @exportClass FGTExperiment
#'
#' @examples
#' \dontrun{
#' # Create a simple FGTExperiment
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#'
#' # Create a basic FGTExperiment
#' fgt_exp <- FGTExperiment(assays = list(counts = counts))
#' }
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom methods setClass setMethod setGeneric setValidity show new is
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames assayNames<-
setClass("FGTExperiment",
         slots = list(
           experimentData = "TreeSummarizedExperiment",
           experimentType = "character",
           fgtMetadata = "SimpleList"
         ),
         prototype = list(
           experimentType = "amplicon",
           fgtMetadata = S4Vectors::SimpleList()
         ))

#' Validity method for FGTExperiment objects
#'
#' @param object FGTExperiment object to validate
#'
#' @return TRUE if the object is valid, otherwise an error message
setValidity("FGTExperiment", function(object) {
  msg <- NULL
  
  # Validate experimentType
  valid_types <- c("amplicon", "metagenomic", "integrated")
  if (length(object@experimentType) != 1 || 
      !object@experimentType %in% valid_types) {
    msg <- c(msg, paste("experimentType must be one of:",
                        paste(valid_types, collapse = ", ")))
  }
  
  # Validate fgtMetadata
  if (!is(object@fgtMetadata, "SimpleList")) {
    msg <- c(msg, "fgtMetadata must be a SimpleList object")
  }
  
  # Return TRUE if no errors, otherwise return error messages
  if (length(msg) == 0) TRUE else msg
})

#' Show method for FGTExperiment objects
#'
#' @param object FGTExperiment object to display
#'
#' @return None, prints to console
#'
#' @importFrom methods show
#' @export
setMethod("show", "FGTExperiment", function(object) {
  cat("FGTExperiment object with", nrow(object), "features and", ncol(object), "samples\n")
  cat("experimentType:", object@experimentType, "\n")
  
  fgt_meta <- object@fgtMetadata
  if (length(fgt_meta) > 0) {
    cat("fgtMetadata: SimpleList with", length(fgt_meta), "elements\n")
    for (i in seq_along(fgt_meta)) {
      cat("  $", names(fgt_meta)[i], ": ", 
          paste(class(fgt_meta[[i]])[1], collapse=", "), 
          if (length(fgt_meta[[i]]) == 1) {
            paste0(" [", as.character(fgt_meta[[i]]), "]")
          } else {
            paste0(" [length: ", length(fgt_meta[[i]]), "]")
          },
          "\n", sep="")
    }
  } else {
    cat("fgtMetadata: Empty SimpleList\n")
  }
  
  # Show the underlying TreeSummarizedExperiment
  cat("\nUnderlying TreeSummarizedExperiment:\n")
  show(object@experimentData)
})