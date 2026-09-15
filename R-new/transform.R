#' @import methods
NULL

#' Transformation methods for FGTExperiment objects
#'
#' Implementations of methods to transform abundance data in FGTExperiment objects.
#'
#' @name transformAbundance-methods
#' @rdname transformAbundance-methods
NULL

#' @describeIn transformAbundance-methods Transform abundance for FGTExperiment objects
#' @param object FGTExperiment object to transform
#' @param type Transformation type: "relative", "log", "clr", "presence"
#' @param assay_name Name of assay to transform
#' @param pseudocount Value to add before log transformations
#' @return A new FGTExperiment object with transformed data
#' @importFrom SummarizedExperiment assay assays assayNames
#' @export
setMethod("transformAbundance", "FGTExperiment", function(object, type = "relative", assay_name = "counts", pseudocount = 1) {
  # Validate inputs
  validTypes <- c("relative", "log", "clr", "presence")
  if (!type %in% validTypes) {
    stop("type must be one of: ", paste(validTypes, collapse=", "))
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
    stop("assay '", assay_name, "' not found")
  }
  
  # Get assay data
  counts <- SummarizedExperiment::assay(object, assay_name)
  
  # Perform transformation
  result <- switch(
    type,
    "relative" = {
      # Relative abundance (proportions)
      t(apply(counts, 2, function(x) {
        if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) {
          return(rep(NA, length(x)))
        }
        return(x / sum(x, na.rm = TRUE))
      }))
    },
    "log" = {
      # Log transformation (log(x + pseudocount))
      log1p(counts)
    },
    "clr" = {
      # Centered log-ratio transformation
      t(apply(counts + pseudocount, 2, function(x) {
        log_x <- log(x)
        log_x - mean(log_x)
      }))
    },
    "presence" = {
      # Presence/absence (logical)
      counts > 0
    }
  )
  
  # Ensure dimensions and dimnames are preserved
  dimnames(result) <- dimnames(counts)
  
  # Create new assay name
  new_assay_name <- paste0(type, "_", assay_name)
  
  # Create a copy with the new assay
  assay_list <- as.list(SummarizedExperiment::assays(object))
  assay_list[[new_assay_name]] <- result
  
  # Build a new object with the updated assays
  # This is safer than trying to modify the existing object
  new_object <- FGTExperiment(
    assays = assay_list,
    rowData = SummarizedExperiment::rowData(object),
    colData = SummarizedExperiment::colData(object),
    experimentType = experimentType(object),
    fgtMetadata = fgtMetadata(object)
  )
  
  return(new_object)
})

#' @describeIn transformAbundance-methods Transform abundance for SummarizedExperiment objects
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
setMethod("transformAbundance", "SummarizedExperiment", function(object, type = "relative", assay_name = "counts", pseudocount = 1) {
  # Validate inputs
  validTypes <- c("relative", "log", "clr", "presence")
  if (!type %in% validTypes) {
    stop("type must be one of: ", paste(validTypes, collapse=", "))
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(object)) {
    stop("assay '", assay_name, "' not found")
  }
  
  # Get assay data
  counts <- SummarizedExperiment::assay(object, assay_name)
  
  # Perform transformation
  result <- switch(
    type,
    "relative" = {
      # Relative abundance (proportions)
      t(apply(counts, 2, function(x) {
        if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) {
          return(rep(NA, length(x)))
        }
        return(x / sum(x, na.rm = TRUE))
      }))
    },
    "log" = {
      # Log transformation (log(x + pseudocount))
      log1p(counts)
    },
    "clr" = {
      # Centered log-ratio transformation
      t(apply(counts + pseudocount, 2, function(x) {
        log_x <- log(x)
        log_x - mean(log_x)
      }))
    },
    "presence" = {
      # Presence/absence (logical)
      counts > 0
    }
  )
  
  # Ensure dimensions and dimnames are preserved
  dimnames(result) <- dimnames(counts)
  
  # Create new assay name
  new_assay_name <- paste0(type, "_", assay_name)
  
  # Create a copy with the new assay
  assay_list <- as.list(SummarizedExperiment::assays(object))
  assay_list[[new_assay_name]] <- result
  
  # Build a new object with the updated assays
  new_object <- SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    rowData = SummarizedExperiment::rowData(object),
    colData = SummarizedExperiment::colData(object)
  )
  
  return(new_object)
})