#' Direct implementation of FGTExperiment constructor
#'
#' This is a direct implementation of the FGTExperiment constructor that is compatible
#' with the new compositional class design.
#'
#' @param assays List of matrices or similar objects containing count data,
#'        or a SummarizedExperiment or TreeSummarizedExperiment object
#' @param rowData DataFrame of feature metadata
#' @param colData DataFrame of sample metadata
#' @param rowTree Phylogenetic tree of type phylo
#' @param experimentType Type of experiment ("amplicon", "metagenomic", "integrated")
#' @param fgtMetadata Additional metadata as a SimpleList
#' @param ... Additional arguments passed to TreeSummarizedExperiment
#'
#' @return A new FGTExperiment object
#' @export
FGTExperiment_direct <- function(assays = list(), rowData = NULL, colData = NULL, 
                              rowTree = NULL, experimentType = "amplicon", 
                              fgtMetadata = S4Vectors::SimpleList(), ...) {
  
  # Create the TreeSummarizedExperiment
  if (methods::is(assays, "TreeSummarizedExperiment")) {
    tse <- assays
  } else if (methods::is(assays, "SummarizedExperiment")) {
    tse <- tryCatch({
      TreeSummarizedExperiment::TreeSummarizedExperiment(
        SummarizedExperiment::assays(assays),
        rowData = SummarizedExperiment::rowData(assays),
        colData = SummarizedExperiment::colData(assays),
        rowTree = rowTree
      )
    }, error = function(e) {
      warning("Error converting SummarizedExperiment to TreeSummarizedExperiment: ", conditionMessage(e))
      # Create an empty TSE with the same dimensions
      tse_empty <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(counts = matrix(0, nrow=nrow(assays), ncol=ncol(assays))),
        rowData = SummarizedExperiment::rowData(assays),
        colData = SummarizedExperiment::colData(assays)
      )
      # Copy the assays over
      SummarizedExperiment::assays(tse_empty) <- SummarizedExperiment::assays(assays)
      tse_empty
    })
  } else {
    # Input validation for creating from scratch
    if (is.list(assays) && length(assays) > 0) {
      # Verify the first assay 
      first_assay <- assays[[1]]
      if (!is.matrix(first_assay) && !methods::is(first_assay, "Matrix")) {
        stop("The first element of 'assays' must be a matrix or Matrix object")
      }
      
      # Ensure row and column names are present
      if (is.null(rownames(first_assay))) {
        rownames(first_assay) <- paste0("Feature", seq_len(nrow(first_assay)))
        assays[[1]] <- first_assay
      }
      
      if (is.null(colnames(first_assay))) {
        colnames(first_assay) <- paste0("Sample", seq_len(ncol(first_assay)))
        assays[[1]] <- first_assay
      }
      
      # Prepare rowData if needed
      if (is.null(rowData)) {
        rowData <- S4Vectors::DataFrame(row.names = rownames(first_assay))
      } else if (!is.data.frame(rowData) && !methods::is(rowData, "DataFrame")) {
        rowData <- S4Vectors::DataFrame(rowData, row.names = rownames(first_assay))
      }
      
      # Prepare colData if needed
      if (is.null(colData)) {
        colData <- S4Vectors::DataFrame(row.names = colnames(first_assay))
      } else if (!is.data.frame(colData) && !methods::is(colData, "DataFrame")) {
        colData <- S4Vectors::DataFrame(colData, row.names = colnames(first_assay))
      }
      
      # Create TreeSummarizedExperiment
      tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = assays,
        rowData = rowData,
        colData = colData,
        rowTree = rowTree,
        ...
      )
    } else {
      # For empty constructor, create a simple example TreeSummarizedExperiment
      counts <- matrix(1:12, nrow = 3, ncol = 4)
      rownames(counts) <- paste0("Gene", 1:3)
      colnames(counts) <- paste0("Sample", 1:4)
      tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(counts = counts)
      )
    }
  }
  
  # Convert fgtMetadata to SimpleList if necessary
  if (!methods::is(fgtMetadata, "SimpleList")) {
    fgtMetadata <- S4Vectors::SimpleList(fgtMetadata)
  }
  
  # Validate experimentType
  validTypes <- c("amplicon", "metagenomic", "integrated")
  if (!experimentType %in% validTypes) {
    stop("experimentType must be one of: ", paste(validTypes, collapse=", "))
  }
  
  # Create FGTExperiment object using new() with the experimentData slot
  fgtObj <- methods::new("FGTExperiment",
                       experimentData = tse,
                       experimentType = experimentType,
                       fgtMetadata = fgtMetadata)
  
  # Validate and return the object
  methods::validObject(fgtObj)
  return(fgtObj)
}

#' Direct implementation of transformAbundance for FGTExperiment
#'
#' @param object FGTExperiment object to transform
#' @param type Transformation type ("relative", "log", "clr", "presence")
#' @param assay_name Name of assay to transform
#' @param pseudocount Value to add for log transformations
#'
#' @return Transformed FGTExperiment object
#' @export
transformAbundance_direct <- function(object, type = "relative", assay_name = "counts", pseudocount = 1) {
  if (!methods::is(object, "FGTExperiment")) {
    stop("Object must be an FGTExperiment")
  }

  # Validate inputs
  valid_types <- c("relative", "log", "clr", "presence")
  if (!type %in% valid_types) {
    stop("type must be one of: ", paste(valid_types, collapse = ", "))
  }

  # Extract the TSE
  tse <- object@experimentData

  if (!assay_name %in% SummarizedExperiment::assayNames(tse)) {
    stop("assay_name '", assay_name, "' not found in assays")
  }

  # Get the count matrix
  counts <- SummarizedExperiment::assay(tse, assay_name)

  # Apply transformation
  if (type == "relative") {
    # Relative abundance (proportions)
    result <- t(t(counts) / colSums(counts))
    new_name <- paste0("relative_", assay_name)
  } else if (type == "log") {
    # Log transformation
    result <- log(counts + pseudocount)
    new_name <- paste0("log_", assay_name)
  } else if (type == "clr") {
    # Centered log-ratio transformation
    counts_adj <- counts + pseudocount
    # Calculate geometric mean for each sample
    geo_means <- apply(counts_adj, 2, function(x) exp(mean(log(x))))
    # Apply CLR transformation
    result <- log(t(t(counts_adj) / geo_means))
    new_name <- paste0("clr_", assay_name)
  } else if (type == "presence") {
    # Presence/absence (binary)
    result <- counts > 0
    mode(result) <- "numeric"  # Convert logical to numeric
    new_name <- paste0("presence_", assay_name)
  }

  # Add the transformed matrix as a new assay
  assays_list <- SummarizedExperiment::assays(tse)
  assays_list[[new_name]] <- result
  SummarizedExperiment::assays(tse) <- assays_list

  # Update the FGTExperiment
  object@experimentData <- tse

  return(object)
}