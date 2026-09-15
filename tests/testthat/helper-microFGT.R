# Standardized test helpers for microFGT tests
library(microFGT)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(methods)

#' Create a test FGTExperiment object with standard parameters
#'
#' @param rows Number of features (rows)
#' @param cols Number of samples (columns)
#' @param experiment_type Type of experiment
#' @param include_taxonomy Whether to include taxonomy data
#' @param include_tree Whether to include a phylogenetic tree
#'
#' @return A standardized FGTExperiment object for testing
create_test_fgt <- function(rows = 10, 
                           cols = 5, 
                           experiment_type = "amplicon", 
                           include_taxonomy = TRUE,
                           include_tree = FALSE) {
  # Create count matrix
  counts <- matrix(rpois(rows * cols, lambda = 20), nrow = rows, ncol = cols)
  rownames(counts) <- paste0("Feature", seq_len(rows))
  colnames(counts) <- paste0("Sample", seq_len(cols))
  
  # Create taxonomy data if requested
  row_data <- NULL
  if (include_taxonomy) {
    taxa <- data.frame(
      Kingdom = rep("Bacteria", rows),
      Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria"), 
                     rows, replace = TRUE),
      Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria"), 
                    rows, replace = TRUE),
      Order = sample(c("Lactobacillales", "Clostridiales", "Bacteroidales", "Bifidobacteriales"), 
                    rows, replace = TRUE),
      Family = sample(c("Lactobacillaceae", "Bifidobacteriaceae", "Prevotellaceae"), 
                     rows, replace = TRUE),
      Genus = sample(c("Lactobacillus", "Gardnerella", "Prevotella", "Sneathia", "Atopobium"), 
                    rows, replace = TRUE),
      Species = paste0("Species", seq_len(rows)),
      row.names = rownames(counts)
    )
    row_data <- S4Vectors::DataFrame(taxa)
  }
  
  # Create sample metadata
  col_data <- S4Vectors::DataFrame(
    group = rep(c("A", "B"), length.out = cols),
    condition = sample(c("Healthy", "BV"), cols, replace = TRUE),
    pH = runif(cols, 3.8, 7.0),
    Nugent_Score = sample(0:10, cols, replace = TRUE),
    row.names = colnames(counts)
  )
  
  # Create phylogenetic tree if requested
  row_tree <- NULL
  if (include_tree) {
    # Create a simple tree structure
    edge <- matrix(c(
      rows + 1, 1,
      rows + 1, rows + 2,
      rows + 2, 2,
      rows + 2, 3,
      rows + 3, 4,
      rows + 3, 5,
      rows + 4, 6,
      rows + 4, 7,
      rows + 5, 8,
      rows + 5, 9,
      rows + 6, 10,
      rows + 6, rows + 3,
      rows + 7, rows + 4,
      rows + 7, rows + 5,
      rows + 8, rows + 6, 
      rows + 8, rows + 7 
    ), ncol = 2, byrow = TRUE)
    
    # Simple edge lengths
    edge.length <- rep(1, nrow(edge))
    
    # Create tip labels
    tip.label <- rownames(counts)
    
    # Create the tree
    row_tree <- structure(
      list(
        edge = edge,
        edge.length = edge.length,
        tip.label = tip.label,
        Nnode = rows - 1
      ),
      class = "phylo"
    )
  }
  
  # Create the FGTExperiment object
  fgt_exp <- FGTExperiment(
    assays = list(counts = counts),
    rowData = row_data,
    colData = col_data,
    rowTree = row_tree,
    experimentType = experiment_type,
    fgtMetadata = S4Vectors::SimpleList(
      creation_date = Sys.Date(),
      test_object = TRUE
    )
  )
  
  return(fgt_exp)
}

#' Check if transformation was correctly applied to FGTExperiment
#'
#' @param original Original FGTExperiment object before transformation
#' @param transformed Transformed FGTExperiment object
#' @param type Type of transformation applied
#' @param assay_name Name of the assay that was transformed
#' @param expected_assay_name Expected name of transformed assay
#'
#' @return TRUE if the transformation appears correct, errors otherwise
verify_transformation <- function(original, 
                                transformed, 
                                type = "relative", 
                                assay_name = "counts",
                                expected_assay_name = NULL) {
  # Determine expected assay name if not provided
  if (is.null(expected_assay_name)) {
    expected_assay_name <- paste0(type, "_", assay_name)
  }
  
  # Check that transformed assay exists
  if (!expected_assay_name %in% assayNames(transformed)) {
    stop("Expected assay '", expected_assay_name, "' not found in transformed object")
  }
  
  # Get original and transformed assays
  orig_assay <- assay(original, assay_name)
  trans_assay <- assay(transformed, expected_assay_name)
  
  # Check dimensions match
  if (!identical(dim(orig_assay), dim(trans_assay))) {
    stop("Transformed assay dimensions don't match original")
  }
  
  # Check rownames and colnames match
  if (!identical(rownames(orig_assay), rownames(trans_assay)) || 
      !identical(colnames(orig_assay), colnames(trans_assay))) {
    stop("Transformed assay names don't match original")
  }
  
  # Check transformation-specific properties
  if (type == "relative") {
    # Column sums should be 1 (or very close)
    col_sums <- colSums(trans_assay)
    if (!all(abs(col_sums - 1) < 1e-10)) {
      stop("Relative abundance transformation doesn't sum to 1 for each sample")
    }
  } else if (type == "log") {
    # Log values should be log(original + 1)
    expected <- log1p(orig_assay)
    if (!all(abs(trans_assay - expected) < 1e-10)) {
      stop("Log transformation doesn't match expected values")
    }
  } else if (type == "clr") {
    # Check that column means are close to zero
    col_means <- colMeans(trans_assay)
    if (!all(abs(col_means) < 1e-10)) {
      stop("CLR transformation doesn't have column means near zero")
    }
  } else if (type == "presence") {
    # Should be 1 where counts > 0 and 0 otherwise
    expected <- as.numeric(orig_assay > 0)
    dim(expected) <- dim(orig_assay)
    rownames(expected) <- rownames(orig_assay)
    colnames(expected) <- colnames(orig_assay)
    if (!identical(trans_assay, expected)) {
      stop("Presence/absence transformation doesn't match expected values")
    }
  }
  
  # If we made it here, the transformation seems correct
  return(TRUE)
}

#' Create a list of test assay matrices with different properties
#'
#' @param rows Number of features (rows)
#' @param cols Number of samples (columns)
#'
#' @return A list of named matrices for testing
create_test_assays <- function(rows = 10, cols = 5) {
  # Feature and sample names
  feature_names <- paste0("Feature", seq_len(rows))
  sample_names <- paste0("Sample", seq_len(cols))
  
  # Create standard count matrix
  counts <- matrix(rpois(rows * cols, lambda = 20), nrow = rows, ncol = cols)
  rownames(counts) <- feature_names
  colnames(counts) <- sample_names
  
  # Create a matrix with zeros
  counts_with_zeros <- counts
  counts_with_zeros[sample(seq_len(rows * cols), size = floor(rows * cols * 0.3))] <- 0
  rownames(counts_with_zeros) <- feature_names
  colnames(counts_with_zeros) <- sample_names
  
  # Create a sparse count matrix
  sparse_counts <- counts_with_zeros
  sparse_counts[sample(seq_len(rows * cols), size = floor(rows * cols * 0.6))] <- 0
  rownames(sparse_counts) <- feature_names
  colnames(sparse_counts) <- sample_names
  
  # Create a relative abundance matrix
  rel_abundance <- t(t(counts) / colSums(counts))
  rownames(rel_abundance) <- feature_names
  colnames(rel_abundance) <- sample_names
  
  # Return list of test matrices
  return(list(
    counts = counts,
    counts_with_zeros = counts_with_zeros,
    sparse_counts = sparse_counts,
    rel_abundance = rel_abundance
  ))
}

#' Helper function to check if all necessary packages are available
#'
#' @param packages Character vector of package names to check
#'
#' @return TRUE if all packages are available, skips test otherwise
check_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste0("Package '", pkg, "' not available"))
      return(FALSE)
    }
  }
  return(TRUE)
}