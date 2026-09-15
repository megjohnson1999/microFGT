#' @import methods
NULL

#' Get the experiment type
#'
#' @param object FGTExperiment object
#'
#' @return Character string indicating the experiment type
#' @export
setGeneric("experimentType", function(object) standardGeneric("experimentType"))

#' Set the experiment type
#'
#' @param object FGTExperiment object
#' @param value New experiment type value
#'
#' @return Updated FGTExperiment object
#' @export
setGeneric("experimentType<-", function(object, value) standardGeneric("experimentType<-"))

#' Get the FGT metadata
#'
#' @param object FGTExperiment object
#'
#' @return SimpleList containing the metadata
#' @export
setGeneric("fgtMetadata", function(object) standardGeneric("fgtMetadata"))

#' Set the FGT metadata
#'
#' @param object FGTExperiment object
#' @param value SimpleList containing the new metadata
#'
#' @return Updated FGTExperiment object
#' @export
setGeneric("fgtMetadata<-", function(object, value) standardGeneric("fgtMetadata<-"))

#' Transform abundance values
#'
#' @param object FGTExperiment or SummarizedExperiment object
#' @param type Transformation type: "relative", "log", "clr", "presence"
#' @param assay_name Name of assay to transform
#' @param pseudocount Value to add before log transformations
#'
#' @return FGTExperiment or SummarizedExperiment with added assay containing the transformed values
#' @export
setGeneric("transformAbundance", function(object, type = "relative", assay_name = "counts", pseudocount = 1) 
  standardGeneric("transformAbundance"))

# Import generics from SummarizedExperiment and TreeSummarizedExperiment
#' @importFrom SummarizedExperiment assay assays rowData colData assayNames
#' @importFrom TreeSummarizedExperiment rowTree

# Define TreeSummarizedExperiment generics with proper signatures
#' @export
setGeneric("rowTree", function(x) standardGeneric("rowTree"))

# Create generics for base functions if they don't exist
if (!isGeneric("rownames")) {
  setGeneric("rownames", function(x) base::rownames(x))
}

if (!isGeneric("colnames")) {
  setGeneric("colnames", function(x) base::colnames(x))
}

if (!isGeneric("dim")) {
  setGeneric("dim", function(x) base::dim(x))
}

if (!isGeneric("dimnames")) {
  setGeneric("dimnames", function(x) base::dimnames(x))
}
