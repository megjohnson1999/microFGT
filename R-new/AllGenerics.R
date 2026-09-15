#' @import methods
NULL

#' Get the experiment type from a FGTExperiment object
#'
#' @param object FGTExperiment object
#'
#' @return Character string indicating the experiment type
#' @export
setGeneric("experimentType", function(object) standardGeneric("experimentType"))

#' Set the experiment type for a FGTExperiment object
#'
#' @param object FGTExperiment object
#' @param value New experiment type value (character)
#'
#' @return Updated FGTExperiment object
#' @export
setGeneric("experimentType<-", function(object, value) standardGeneric("experimentType<-"))

#' Get the FGT metadata from a FGTExperiment object
#'
#' @param object FGTExperiment object
#'
#' @return SimpleList containing the metadata
#' @export
setGeneric("fgtMetadata", function(object) standardGeneric("fgtMetadata"))

#' Set the FGT metadata for a FGTExperiment object
#'
#' @param object FGTExperiment object
#' @param value SimpleList containing the new metadata
#'
#' @return Updated FGTExperiment object
#' @export
setGeneric("fgtMetadata<-", function(object, value) standardGeneric("fgtMetadata<-"))

#' Transform abundance values in a FGTExperiment object
#'
#' @param object FGTExperiment or SummarizedExperiment object
#' @param type Transformation type: "relative", "log", "clr", "presence"
#' @param assay_name Name of assay to transform
#' @param pseudocount Value to add before log transformations
#'
#' @return A new object with transformed assay data
#' @export
setGeneric("transformAbundance", function(object, type = "relative", assay_name = "counts", pseudocount = 1) standardGeneric("transformAbundance"))