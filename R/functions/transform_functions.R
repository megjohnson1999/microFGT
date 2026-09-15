#' Transform count data using various methods
#'
#' This function transforms count or abundance data in SummarizedExperiment objects
#' using various methods like relative abundance, log, CLR, etc.
#'
#' @param x A SummarizedExperiment or FGTExperiment object
#' @param method Transformation method: "relative", "log", "clr", "presence"
#' @param assay_name Name of the assay to transform
#' @param pseudocount Value to add before log transformations
#'
#' @return The input object with a new assay containing transformed values
#' @export
transform_counts <- function(x, method = c("relative", "log", "clr", "presence"),
                            assay_name = "counts", pseudocount = 1) {
  
  # Match the method argument
  method <- match.arg(method)
  
  # Check if the object has assays
  if (!methods::is(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment or FGTExperiment object")
  }
  
  # Check if the requested assay exists
  if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
    stop("Assay '", assay_name, "' not found in the object")
  }
  
  # Get the counts/abundance matrix
  counts <- SummarizedExperiment::assays(x)[[assay_name]]
  
  # Perform the transformation
  transformed <- switch(
    method,
    "relative" = {
      # Relative abundance (proportions)
      t(apply(counts, 2, function(col) {
        if (all(is.na(col)) || sum(col, na.rm = TRUE) == 0) {
          return(rep(NA, length(col)))
        }
        return(col / sum(col, na.rm = TRUE))
      }))
    },
    "log" = {
      # Log transformation (log(x + pseudocount))
      log1p(counts)
    },
    "clr" = {
      # Centered log-ratio transformation
      t(apply(counts + pseudocount, 2, function(col) {
        log_col <- log(col)
        return(log_col - mean(log_col))
      }))
    },
    "presence" = {
      # Presence/absence (logical)
      counts > 0
    }
  )
  
  # Add the transformed data as a new assay
  SummarizedExperiment::assays(x)[[method]] <- transformed
  
  return(x)
}

#' Filter features based on prevalence and abundance
#'
#' @param x A SummarizedExperiment or FGTExperiment object
#' @param min_prevalence Minimum prevalence (proportion of samples)
#' @param min_abundance Minimum abundance threshold 
#' @param abundance_type Type of abundance filter: "absolute" or "relative"
#' @param assay_name Name of the assay to use
#' @param detection_threshold Minimum value to consider a feature present
#'
#' @return Filtered SummarizedExperiment with fewer features
#' @export
filter_features <- function(x, min_prevalence = 0.1, min_abundance = 0.001,
                           abundance_type = c("relative", "absolute"),
                           assay_name = "counts", detection_threshold = 0) {
  
  # Match the abundance_type argument
  abundance_type <- match.arg(abundance_type)
  
  # Check if the object has assays
  if (!methods::is(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment or FGTExperiment object")
  }
  
  # Check if the requested assay exists
  if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
    stop("Assay '", assay_name, "' not found in the object")
  }
  
  # Check prevalence values
  if (min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be between 0 and 1")
  }
  
  # Get the counts/abundance matrix
  counts <- SummarizedExperiment::assays(x)[[assay_name]]
  
  # Calculate prevalence (proportion of samples where feature is present)
  prevalence <- rowSums(counts > detection_threshold, na.rm = TRUE) / ncol(counts)
  keep_prevalence <- prevalence >= min_prevalence
  
  # Filter by abundance if required
  if (min_abundance > 0) {
    if (abundance_type == "relative") {
      # Convert to relative abundance for filtering
      rel_counts <- t(apply(counts, 2, function(col) {
        if (all(is.na(col)) || sum(col, na.rm = TRUE) == 0) {
          return(rep(NA, length(col)))
        }
        return(col / sum(col, na.rm = TRUE))
      }))
      
      # Keep features with relative abundance >= threshold in at least one sample
      keep_abundance <- apply(rel_counts, 1, function(row) {
        any(row >= min_abundance, na.rm = TRUE)
      })
    } else {
      # Keep features with absolute abundance >= threshold in at least one sample
      keep_abundance <- apply(counts, 1, function(row) {
        any(row >= min_abundance, na.rm = TRUE)
      })
    }
  } else {
    # No abundance filtering
    keep_abundance <- rep(TRUE, nrow(counts))
  }
  
  # Combine filters
  keep <- keep_prevalence & keep_abundance
  
  # Subset the object
  return(x[keep, ])
}