#' Create a new FGTExperiment object
#'
#' @param assays List of matrices or similar objects containing count data
#' @param rowData DataFrame of feature metadata
#' @param colData DataFrame of sample metadata
#' @param rowTree Phylogenetic tree of type phylo
#' @param experimentType Type of experiment ("amplicon", "metagenomic", "integrated")
#' @param fgtMetadata Additional metadata as a SimpleList
#' @param ... Additional arguments passed to TreeSummarizedExperiment
#'
#' @return A new FGTExperiment object
#' @export
FGTExperiment <- function(assays = list(), rowData = NULL, colData = NULL, rowTree = NULL,
                          experimentType = "amplicon", fgtMetadata = S4Vectors::SimpleList(), ...) {
  
  # Create the TreeSummarizedExperiment
  if (methods::is(assays, "TreeSummarizedExperiment")) {
    tse <- assays
  } else if (methods::is(assays, "SummarizedExperiment")) {
    tse <- tryCatch({
      TreeSummarizedExperiment::TreeSummarizedExperiment(assays, rowTree = rowTree)
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
      stop("'assays' must be a non-empty list, SummarizedExperiment, or TreeSummarizedExperiment")
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
  
  # Create FGTExperiment object using new()
  fgtObj <- methods::new("FGTExperiment",
                         experimentData = tse,
                         experimentType = experimentType,
                         fgtMetadata = fgtMetadata)
  
  # Validate and return the object
  methods::validObject(fgtObj)
  return(fgtObj)
}