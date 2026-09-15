#' Legacy wrapper functions
#'
#' These functions are kept for backward compatibility
#' @name legacy_exports
#' @import methods
NULL

#' Test transform function
#'
#' @param data A matrix of data to transform
#' @param method Transformation method
#' @return Transformed matrix
#' @export
test_transform <- function(data, method = c("log", "relative", "sqrt")) {
  method <- match.arg(method)
  
  switch(method,
         "log" = log1p(data),
         "relative" = {
           t(apply(data, 2, function(x) x / sum(x)))
         },
         "sqrt" = sqrt(data))
}

#' Test if an object is an FGTExperiment
#'
#' @param x Object to test
#' @return Logical
#' @export
is_fgt <- function(x) {
  methods::is(x, "FGTExperiment")
}

#' Get experiment type
#'
#' @param x FGTExperiment object
#' @return Character string
#' @export
get_experiment_type <- function(x) {
  if (is_fgt(x)) {
    return(experimentType(x))
  } else {
    return(NA_character_)
  }
}

#' Get FGT metadata
#'
#' @param x FGTExperiment object
#' @return SimpleList
#' @export
get_fgt_metadata <- function(x) {
  if (is_fgt(x)) {
    return(fgtMetadata(x))
  } else {
    return(S4Vectors::SimpleList())
  }
}

#' Set experiment type
#'
#' @param x FGTExperiment object
#' @param value New value
#' @return Updated object
#' @export
set_experiment_type <- function(x, value) {
  if (!is_fgt(x)) {
    stop("Object must be an FGTExperiment")
  }
  
  experimentType(x) <- value
  return(x)
}

#' Set FGT metadata
#'
#' @param x FGTExperiment object
#' @param value New value
#' @return Updated object
#' @export
set_fgt_metadata <- function(x, value) {
  if (!is_fgt(x)) {
    stop("Object must be an FGTExperiment")
  }
  
  fgtMetadata(x) <- value
  return(x)
}

#' Create FGTExperiment object
#'
#' @param counts Counts matrix
#' @param taxonomy Taxonomy data
#' @param sample_data Sample metadata
#' @param experiment_type Experiment type
#' @param ... Additional arguments
#' @return FGTExperiment object
#' @export
create_fgt <- function(counts, taxonomy = NULL, sample_data = NULL, 
                      experiment_type = "amplicon", ...) {
  
  # Check inputs
  if (!is.matrix(counts)) {
    stop("counts must be a matrix")
  }
  
  # Prepare rowData if provided
  if (!is.null(taxonomy)) {
    rowData <- taxonomy
  } else {
    rowData <- data.frame(row.names = rownames(counts))
  }
  
  # Prepare colData if provided
  if (!is.null(sample_data)) {
    colData <- sample_data
  } else {
    colData <- data.frame(row.names = colnames(counts))
  }
  
  # Call the main constructor
  FGTExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData,
    experimentType = experiment_type,
    ...
  )
}