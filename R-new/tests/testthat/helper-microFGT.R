#' Create a test FGTExperiment object
#'
#' @param rows Number of rows
#' @param cols Number of columns
#' @param withTree Include a phylogenetic tree
#' @param experimentType Type of experiment
#'
#' @return FGTExperiment object
create_test_object <- function(rows=10, cols=5, withTree=FALSE, experimentType="amplicon") {
  # Create count matrix
  counts <- matrix(rpois(rows*cols, lambda=10), nrow=rows, ncol=cols)
  rownames(counts) <- paste0("Feature", seq_len(rows))
  colnames(counts) <- paste0("Sample", seq_len(cols))
  
  # Create row data
  rowData <- S4Vectors::DataFrame(
    FeatureID = rownames(counts),
    row.names = rownames(counts)
  )
  
  # Create column data
  colData <- S4Vectors::DataFrame(
    SampleID = colnames(counts),
    Condition = rep(c("A", "B"), length.out=cols),
    row.names = colnames(counts)
  )
  
  # Create tree if requested
  rowTree <- NULL
  if (withTree) {
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Package 'ape' needed to create test tree")
    }
    rowTree <- ape::rtree(rows)
    rowTree$tip.label <- rownames(counts)
  }
  
  # Create FGTExperiment
  FGTExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData,
    rowTree = rowTree,
    experimentType = experimentType,
    fgtMetadata = S4Vectors::SimpleList(
      CreatedBy = "test_helper",
      Date = Sys.Date()
    )
  )
}