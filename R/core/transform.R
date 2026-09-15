#' @describeIn transformAbundance Transform abundance for FGTExperiment objects
#' @export
setMethod("transformAbundance", "FGTExperiment", function(object, type = "relative", assay_name = "counts", pseudocount = 1) {
  # Get the underlying TreeSummarizedExperiment
  tse <- object@experimentData
  
  # Transform the TreeSummarizedExperiment
  transformed_tse <- transformAbundance(tse, type, assay_name, pseudocount)
  
  # Update the experimentData slot with the transformed TSE
  object@experimentData <- transformed_tse
  
  return(object)
})

#' @describeIn transformAbundance Transform abundance for TreeSummarizedExperiment objects
#' @export
setMethod("transformAbundance", "TreeSummarizedExperiment", function(object, type = "relative", assay_name = "counts", pseudocount = 1) {
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
      rel_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
      for (j in seq_len(ncol(counts))) {
        col_sum <- sum(counts[, j], na.rm = TRUE)
        if (col_sum > 0) {
          rel_counts[, j] <- counts[, j] / col_sum
        } else {
          rel_counts[, j] <- NA
        }
      }
      rel_counts
    },
    "log" = {
      # Log transformation (log(x + pseudocount))
      log1p(counts)
    },
    "clr" = {
      # Centered log-ratio transformation
      clr_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
      for (j in seq_len(ncol(counts))) {
        x <- counts[, j] + pseudocount
        log_x <- log(x)
        clr_counts[, j] <- log_x - mean(log_x)
      }
      clr_counts
    },
    "presence" = {
      # Presence/absence (logical)
      counts > 0
    }
  )
  
  # Ensure dimensions and names are preserved
  if (!identical(dim(result), dim(counts))) {
    stop("Transformation resulted in inconsistent dimensions")
  }
  
  rownames(result) <- rownames(counts)
  colnames(result) <- colnames(counts)
  
  # Create a copy with the new assay
  new_object <- object
  SummarizedExperiment::assays(new_object)[[paste0(type, "_", assay_name)]] <- result
  
  return(new_object)
})

#' @describeIn transformAbundance Transform abundance for SummarizedExperiment objects
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
      rel_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
      for (j in seq_len(ncol(counts))) {
        col_sum <- sum(counts[, j], na.rm = TRUE)
        if (col_sum > 0) {
          rel_counts[, j] <- counts[, j] / col_sum
        } else {
          rel_counts[, j] <- NA
        }
      }
      rel_counts
    },
    "log" = {
      # Log transformation (log(x + pseudocount))
      log1p(counts)
    },
    "clr" = {
      # Centered log-ratio transformation
      clr_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
      for (j in seq_len(ncol(counts))) {
        x <- counts[, j] + pseudocount
        log_x <- log(x)
        clr_counts[, j] <- log_x - mean(log_x)
      }
      clr_counts
    },
    "presence" = {
      # Presence/absence (logical)
      counts > 0
    }
  )
  
  # Ensure dimensions and names are preserved
  if (!identical(dim(result), dim(counts))) {
    stop("Transformation resulted in inconsistent dimensions")
  }
  
  rownames(result) <- rownames(counts)
  colnames(result) <- colnames(counts)
  
  # Create a copy with the new assay
  new_object <- object
  SummarizedExperiment::assays(new_object)[[paste0(type, "_", assay_name)]] <- result
  
  return(new_object)
})