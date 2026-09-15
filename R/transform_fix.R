#' Transform abundance values
#' 
#' @param x Input data - must be an FGTExperiment or SummarizedExperiment
#' @param type Transformation type: "relative", "log", "clr", "presence"
#' @param assay_name Name of assay to transform
#' @param pseudocount Value to add before log transformations
#' @param ... Additional arguments (not used)
#' @return Transformed data
#' @export
transform_abundance <- function(x, type = "relative", assay_name = "counts", pseudocount = 1, ...) {
  # Check input
  if (!methods::is(x, "FGTExperiment") && !methods::is(x, "SummarizedExperiment")) {
    stop("x must be a FGTExperiment or SummarizedExperiment object")
  }
  
  # Check if the requested assay exists
  if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
    stop(paste0("assay '", assay_name, "' not found"))
  }
  
  # Get the assay matrix
  counts <- SummarizedExperiment::assay(x, assay_name)
  
  # Perform the requested transformation
  transformed <- switch(
    type,
    "relative" = {
      # Relative abundance (proportions)
      rel_counts <- matrix(0, nrow = nrow(counts), ncol = ncol(counts))
      for (i in 1:ncol(counts)) {
        col_sum <- sum(counts[, i], na.rm = TRUE)
        if (col_sum > 0) {
          rel_counts[, i] <- counts[, i] / col_sum
        } else {
          rel_counts[, i] <- rep(NA, nrow(counts))
        }
      }
      dimnames(rel_counts) <- dimnames(counts)
      rel_counts
    },
    "log" = {
      # Log transformation (log(x + pseudocount))
      log_counts <- log1p(counts)
      dimnames(log_counts) <- dimnames(counts)
      log_counts
    },
    "clr" = {
      # Centered log-ratio transformation
      clr_counts <- matrix(0, nrow = nrow(counts), ncol = ncol(counts))
      for (i in 1:ncol(counts)) {
        col_log <- log(counts[, i] + pseudocount)
        clr_counts[, i] <- col_log - mean(col_log)
      }
      dimnames(clr_counts) <- dimnames(counts)
      clr_counts
    },
    "presence" = {
      # Presence/absence (logical)
      pres_counts <- counts > 0
      dimnames(pres_counts) <- dimnames(counts)
      pres_counts
    },
    stop(paste0("Unknown transformation type: ", type))
  )
  
  # Create a new object with the same structure
  if (methods::is(x, "FGTExperiment")) {
    # For FGTExperiment we need to handle the special slots
    assays_list <- as.list(SummarizedExperiment::assays(x))
    assays_list[[type]] <- transformed
    
    # Create a new FGTExperiment with all the same data but updated assays
    result <- FGTExperiment(
      assays = assays_list,
      rowData = SummarizedExperiment::rowData(x),
      colData = SummarizedExperiment::colData(x),
      experimentType = experimentType(x),
      fgtMetadata = fgtMetadata(x)
    )
  } else {
    # For SummarizedExperiment we can just update the assays
    assays_list <- as.list(SummarizedExperiment::assays(x))
    assays_list[[type]] <- transformed
    
    # Create a new SummarizedExperiment
    result <- SummarizedExperiment::SummarizedExperiment(
      assays = assays_list,
      rowData = SummarizedExperiment::rowData(x),
      colData = SummarizedExperiment::colData(x)
    )
  }
  
  return(result)
}