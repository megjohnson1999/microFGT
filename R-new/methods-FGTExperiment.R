#' @import methods
NULL

#' Methods for the FGTExperiment class
#'
#' Implementation of methods for accessing and manipulating FGTExperiment objects.
#'
#' @name FGTExperiment-methods
#' @rdname FGTExperiment-methods
NULL

#' @describeIn FGTExperiment-methods Get experiment type
#' @export
setMethod("experimentType", "FGTExperiment", function(object) {
  object@experimentType
})

#' @describeIn FGTExperiment-methods Set experiment type
#' @export
setMethod("experimentType<-", "FGTExperiment", function(object, value) {
  validTypes <- c("amplicon", "metagenomic", "integrated")
  if (!value %in% validTypes) {
    stop("experimentType must be one of: ", paste(validTypes, collapse=", "))
  }
  
  object@experimentType <- value
  validObject(object)
  return(object)
})

#' @describeIn FGTExperiment-methods Get FGT metadata
#' @export
setMethod("fgtMetadata", "FGTExperiment", function(object) {
  object@fgtMetadata
})

#' @describeIn FGTExperiment-methods Set FGT metadata
#' @export
setMethod("fgtMetadata<-", "FGTExperiment", function(object, value) {
  if (!is(value, "SimpleList")) {
    stop("fgtMetadata must be a SimpleList")
  }
  
  object@fgtMetadata <- value
  validObject(object)
  return(object)
})

#' @describeIn FGTExperiment-methods Show method for FGTExperiment objects
#' @importFrom methods show
#' @param object The FGTExperiment object to display
#' @export
setMethod("show", "FGTExperiment", function(object) {
  callNextMethod()
  cat("experimentType:", object@experimentType, "\n")
  cat("fgtMetadata names(", length(object@fgtMetadata), "): ", 
      paste(names(object@fgtMetadata), collapse = ", "), "\n", sep = "")
})