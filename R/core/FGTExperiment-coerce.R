#' Coercion methods for FGTExperiment class
#'
#' This file contains methods for converting between FGTExperiment and related classes.
#'
#' @name FGTExperiment-coerce
#' @rdname FGTExperiment-coerce
#'
#' @importFrom methods setAs as
#' @importFrom S4Vectors SimpleList
NULL

#' Convert TreeSummarizedExperiment to FGTExperiment
#'
#' @param from A TreeSummarizedExperiment object
#' @param to Class to convert to
#'
#' @return An FGTExperiment object
#' @export
#' @examples
#' \dontrun{
#' # Create a simple TreeSummarizedExperiment
#' library(TreeSummarizedExperiment)
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' tse <- TreeSummarizedExperiment(assays = list(counts = counts))
#' 
#' # Convert to FGTExperiment
#' fgt_exp <- as(tse, "FGTExperiment")
#' }
setAs("TreeSummarizedExperiment", "FGTExperiment", function(from) {
  obj <- new("FGTExperiment", 
            from,
            experimentType = "amplicon",
            fgtMetadata = S4Vectors::SimpleList())
  
  # Ensure assay names are set
  if (length(assayNames(obj)) > 0 && assayNames(obj)[1] == "") {
    assayNames(obj)[1] <- "counts"
  }
  
  return(obj)
})

#' Convert SummarizedExperiment to FGTExperiment
#'
#' @param from A SummarizedExperiment object
#' @param to Class to convert to
#'
#' @return An FGTExperiment object
#' @export
#' @examples
#' \dontrun{
#' # Create a simple SummarizedExperiment
#' library(SummarizedExperiment)
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' 
#' # Convert to FGTExperiment
#' fgt_exp <- as(se, "FGTExperiment")
#' }
setAs("SummarizedExperiment", "FGTExperiment", function(from) {
  # First convert to TreeSummarizedExperiment
  tse <- tryCatch({
    TreeSummarizedExperiment::TreeSummarizedExperiment(from)
  }, error = function(e) {
    message("Error creating TreeSummarizedExperiment: ", conditionMessage(e))
    # Return the SummarizedExperiment if TSE creation fails
    from
  })
  
  # Then convert to FGTExperiment
  obj <- new("FGTExperiment", 
            tse,
            experimentType = "amplicon",
            fgtMetadata = S4Vectors::SimpleList())
  
  # Ensure assay names are set
  if (length(assayNames(obj)) > 0 && assayNames(obj)[1] == "") {
    assayNames(obj)[1] <- "counts"
  }
  
  return(obj)
})

#' Convert FGTExperiment to TreeSummarizedExperiment
#'
#' @param from An FGTExperiment object
#' @param to Class to convert to
#'
#' @return A TreeSummarizedExperiment object
#' @export
#' @examples
#' \dontrun{
#' # Create a simple FGTExperiment
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' fgt_exp <- FGTExperiment(assays = list(counts = counts))
#' 
#' # Convert to TreeSummarizedExperiment
#' tse <- as(fgt_exp, "TreeSummarizedExperiment")
#' }
setAs("FGTExperiment", "TreeSummarizedExperiment", function(from) {
  # Create a new TreeSummarizedExperiment with the same data but without extra slots
  tse <- new("TreeSummarizedExperiment",
             from)
  
  return(tse)
})

#' Convert FGTExperiment to SummarizedExperiment
#'
#' @param from An FGTExperiment object
#' @param to Class to convert to 
#'
#' @return A SummarizedExperiment object
#' @export
#' @examples
#' \dontrun{
#' # Create a simple FGTExperiment
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' fgt_exp <- FGTExperiment(assays = list(counts = counts))
#' 
#' # Convert to SummarizedExperiment
#' se <- as(fgt_exp, "SummarizedExperiment")
#' }
setAs("FGTExperiment", "SummarizedExperiment", function(from) {
  # First convert to TreeSummarizedExperiment
  tse <- as(from, "TreeSummarizedExperiment")
  
  # Then to SummarizedExperiment
  se <- new("SummarizedExperiment",
           tse)
  
  return(se)
})

#' Convert phyloseq to FGTExperiment
#'
#' This function is only available if the phyloseq package is installed.
#'
#' @param from A phyloseq object
#' @param to Class to convert to
#'
#' @return An FGTExperiment object
#' @export
#' @examples
#' \dontrun{
#' # Only if phyloseq is installed
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   # Create a simple phyloseq object
#'   library(phyloseq)
#'   counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#'   rownames(counts) <- paste0("Feature", 1:10)
#'   colnames(counts) <- paste0("Sample", 1:6)
#'   otu <- otu_table(counts, taxa_are_rows = TRUE)
#'   ps <- phyloseq(otu)
#'   
#'   # Convert to FGTExperiment
#'   fgt_exp <- as(ps, "FGTExperiment")
#' }
#' }
#' @importFrom methods setClass setAs existsMethod
setAs("phyloseq", "FGTExperiment", function(from) {
  # Check if phyloseq is available
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required for this coercion but is not installed.")
  }
  
  # Extract components from phyloseq object
  # OTU table (required)
  otu_table <- phyloseq::otu_table(from)
  
  # Check if taxa are rows, which is the expected orientation
  taxa_are_rows <- phyloseq::taxa_are_rows(otu_table)
  
  # Get counts matrix
  counts <- as(otu_table, "matrix")
  
  # Transpose if taxa are not rows
  if (!taxa_are_rows) {
    counts <- t(counts)
  }
  
  # Sample data (optional)
  if (phyloseq::nsamples(from) > 0 && !is.null(phyloseq::sample_data(from, FALSE))) {
    col_data <- data.frame(phyloseq::sample_data(from))
  } else {
    col_data <- data.frame(row.names = colnames(counts))
  }
  
  # Taxa data (optional)
  if (phyloseq::ntaxa(from) > 0 && !is.null(phyloseq::tax_table(from, FALSE))) {
    row_data <- data.frame(phyloseq::tax_table(from))
  } else {
    row_data <- data.frame(row.names = rownames(counts))
  }
  
  # Phylogenetic tree (optional)
  if (!is.null(phyloseq::phy_tree(from, FALSE))) {
    row_tree <- phyloseq::phy_tree(from)
  } else {
    row_tree <- NULL
  }
  
  # Create FGTExperiment object
  fgt_exp <- FGTExperiment(
    assays = list(counts = counts),
    rowData = row_data,
    colData = col_data,
    rowTree = row_tree,
    experimentType = "amplicon"
  )
  
  return(fgt_exp)
})

#' Convert FGTExperiment to phyloseq
#'
#' This function is only available if the phyloseq package is installed.
#'
#' @param from An FGTExperiment object
#' @param to Class to convert to
#'
#' @return A phyloseq object
#' @export
#' @examples
#' \dontrun{
#' # Only if phyloseq is installed
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   # Create a simple FGTExperiment
#'   counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#'   rownames(counts) <- paste0("Feature", 1:10)
#'   colnames(counts) <- paste0("Sample", 1:6)
#'   fgt_exp <- FGTExperiment(assays = list(counts = counts))
#'   
#'   # Convert to phyloseq
#'   ps <- as(fgt_exp, "phyloseq")
#' }
#' }
setAs("FGTExperiment", "phyloseq", function(from) {
  # Check if phyloseq is available
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required for this coercion but is not installed.")
  }
  
  # Extract components from FGTExperiment object
  # Default to the first assay (usually counts)
  assay_name <- assayNames(from)[1]
  counts <- assay(from, assay_name)
  
  # Create components for phyloseq object
  otu <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  
  # Create empty structure for other components
  components <- list(otu)
  names(components) <- c("otu_table")
  
  # Add sample data if available
  if (ncol(colData(from)) > 0) {
    sample_data <- phyloseq::sample_data(data.frame(colData(from)))
    components$sam_data <- sample_data
  }
  
  # Add taxonomy data if available
  if (ncol(rowData(from)) > 0) {
    tax_mat <- as.matrix(rowData(from))
    if (ncol(tax_mat) > 0) {
      tax_tab <- phyloseq::tax_table(tax_mat)
      components$tax_table <- tax_tab
    }
  }
  
  # Add tree if available
  row_tree <- rowTree(from)
  if (!is.null(row_tree)) {
    components$phy_tree <- row_tree
  }
  
  # Create phyloseq object
  ps <- do.call(phyloseq::phyloseq, components)
  
  return(ps)
})