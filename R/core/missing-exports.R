#' Test Transform Function
#'
#' @param data A matrix of data to transform
#' @param method The transformation method to use
#'
#' @return A transformed matrix
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

#' Check if an object is an FGTExperiment
#'
#' @param x Object to test
#'
#' @return Logical indicating if x is an FGTExperiment
#' @export
is_fgt <- function(x) {
  methods::is(x, "FGTExperiment")
}

#' Get the experiment type from an object
#'
#' @param x An FGTExperiment or compatible object
#'
#' @return Character string with the experiment type
#' @export
get_experiment_type <- function(x) {
  if (is_fgt(x)) {
    return(experimentType(x))
  } else {
    return(NA_character_)
  }
}

#' Get the FGT metadata from an object
#'
#' @param x An FGTExperiment or compatible object
#'
#' @return SimpleList containing the FGT metadata
#' @export
get_fgt_metadata <- function(x) {
  if (is_fgt(x)) {
    return(fgtMetadata(x))
  } else {
    return(S4Vectors::SimpleList())
  }
}

#' Set the experiment type
#'
#' @param x An FGTExperiment object
#' @param value New experiment type value
#'
#' @return Updated FGTExperiment object
#' @export
set_experiment_type <- function(x, value) {
  if (!is_fgt(x)) {
    stop("Object must be an FGTExperiment")
  }
  
  experimentType(x) <- value
  return(x)
}

#' Set the FGT metadata
#'
#' @param x An FGTExperiment object
#' @param value SimpleList of new metadata
#'
#' @return Updated FGTExperiment object
#' @export
set_fgt_metadata <- function(x, value) {
  if (!is_fgt(x)) {
    stop("Object must be an FGTExperiment")
  }
  
  fgtMetadata(x) <- value
  return(x)
}

#' Create a new FGTExperiment object
#'
#' This is a simplified constructor for creating FGTExperiment objects.
#' 
#' @param counts A matrix of counts data
#' @param taxonomy Optional data frame with taxonomic classifications
#' @param sample_data Optional data frame with sample metadata
#' @param experiment_type Type of experiment (amplicon, metagenomic, integrated)
#' @param ... Additional arguments passed to TreeSummarizedExperiment constructor
#'
#' @return A new FGTExperiment object
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