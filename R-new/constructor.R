#' @import methods
NULL

#' Create a new FGTExperiment object
#'
#' This constructor creates an FGTExperiment object from count data and metadata.
#'
#' @param assays List of matrices or similar objects containing count data
#' @param rowData DataFrame of feature metadata. If NULL, a minimal DataFrame is created.
#' @param colData DataFrame of sample metadata. If NULL, a minimal DataFrame is created.
#' @param rowTree Phylogenetic tree of type phylo. Optional.
#' @param experimentType Type of experiment: "amplicon", "metagenomic", or "integrated"
#' @param fgtMetadata Additional metadata as a SimpleList
#' @param ... Additional arguments passed to TreeSummarizedExperiment constructor
#'
#' @return A new FGTExperiment object
#' @export
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @examples
#' # Create sample data
#' counts <- matrix(rpois(30, lambda = 10), nrow = 6, ncol = 5)
#' rownames(counts) <- paste0("Feature", 1:6)
#' colnames(counts) <- paste0("Sample", 1:5)
#'
#' # Create a simple FGTExperiment object
#' fgt <- FGTExperiment(assays = list(counts = counts))
#' fgt
FGTExperiment <- function(assays, rowData = NULL, colData = NULL, rowTree = NULL,
                         experimentType = "amplicon", fgtMetadata = S4Vectors::SimpleList(), ...) {
  # Input validation
  if (!is.list(assays) || length(assays) == 0) {
    stop("'assays' must be a non-empty list")
  }

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

  # Create SummarizedExperiment 
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = rowData,
    colData = colData,
    ...
  )

  # Convert to TreeSummarizedExperiment with error handling
  tse <- tryCatch({
    TreeSummarizedExperiment::TreeSummarizedExperiment(se, rowTree = rowTree)
  }, error = function(e) {
    warning("Error creating TreeSummarizedExperiment: ", conditionMessage(e),
            ". Using SummarizedExperiment instead.")
    se
  })

  # Create FGTExperiment object
  obj <- methods::new("FGTExperiment",
                     tse,
                     experimentType = experimentType,
                     fgtMetadata = fgtMetadata)

  # Validate the object
  methods::validObject(obj)
  return(obj)
}