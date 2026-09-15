#' @describeIn FGTExperiment Get experiment type
#' @export
setMethod("experimentType", "FGTExperiment", function(object) {
  object@experimentType
})

#' @describeIn FGTExperiment Set experiment type
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

#' @describeIn FGTExperiment Get FGT metadata
#' @export
setMethod("fgtMetadata", "FGTExperiment", function(object) {
  object@fgtMetadata
})

#' @describeIn FGTExperiment Set FGT metadata
#' @export
setMethod("fgtMetadata<-", "FGTExperiment", function(object, value) {
  if (!methods::is(value, "SimpleList")) {
    stop("fgtMetadata must be a SimpleList")
  }
  
  object@fgtMetadata <- value
  validObject(object)
  return(object)
})

# Proxy methods to pass through to the underlying TreeSummarizedExperiment
#' @export
setMethod("dim", "FGTExperiment", function(x) {
  dim(x@experimentData)
})

#' @export
setMethod("dimnames", "FGTExperiment", function(x) {
  dimnames(x@experimentData)
})

#' @export
setMethod("rownames", "FGTExperiment", function(x) {
  rownames(x@experimentData)
})

#' @export
setMethod("colnames", "FGTExperiment", function(x) {
  colnames(x@experimentData)
})

#' @export
setMethod("assays", "FGTExperiment", function(x, ...) {
  assays(x@experimentData, ...)
})

#' Get an assay
#' @param x FGTExperiment object
#' @param i Assay name or index
#' @param withDimnames Whether to include dimension names
#' @param ... Additional arguments
#' @export
setMethod("assay", signature(x = "FGTExperiment", i = "ANY"), function(x, i = 1, withDimnames = TRUE, ...) {
  SummarizedExperiment::assay(x@experimentData, i = i, withDimnames = withDimnames, ...)
})

#' @export
setMethod("assayNames", "FGTExperiment", function(x, ...) {
  SummarizedExperiment::assayNames(x@experimentData)
})

#' @export
setMethod("rowData", "FGTExperiment", function(x, ...) {
  rowData(x@experimentData, ...)
})

#' @export
setMethod("colData", "FGTExperiment", function(x, ...) {
  colData(x@experimentData, ...)
})

#' @export
setMethod("rowTree", "FGTExperiment", function(x, ...) {
  rowTree(x@experimentData, ...)
})

# Define coercion methods
#' @describeIn FGTExperiment Coerce from FGTExperiment to TreeSummarizedExperiment
#' @export
setAs("FGTExperiment", "TreeSummarizedExperiment", function(from) {
  from@experimentData
})

#' @describeIn FGTExperiment Coerce from FGTExperiment to SummarizedExperiment
#' @export
setAs("FGTExperiment", "SummarizedExperiment", function(from) {
  as(from@experimentData, "SummarizedExperiment")
})

#' @describeIn FGTExperiment Coerce from TreeSummarizedExperiment to FGTExperiment
#' @export
setAs("TreeSummarizedExperiment", "FGTExperiment", function(from) {
  FGTExperiment(from)
})

#' @describeIn FGTExperiment Coerce from SummarizedExperiment to FGTExperiment
#' @export
setAs("SummarizedExperiment", "FGTExperiment", function(from) {
  FGTExperiment(from)
})