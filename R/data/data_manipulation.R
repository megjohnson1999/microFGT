#' Filter taxa based on prevalence and abundance
#'
#' Filters taxa in an FGTExperiment object based on prevalence and abundance thresholds.
#'
#' @param fgt_exp An FGTExperiment object
#' @param min_prevalence Minimum prevalence (proportion of samples)
#' @param min_abundance Minimum relative abundance
#' @param assay_name Name of the count assay to use
#'
#' @return Filtered FGTExperiment object
#' @export
#'
#' @examples
#' # Create a simple FGTExperiment
#' counts <- matrix(sample(0:100, 60, replace = TRUE), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' 
#' fgt_exp <- FGTExperiment(assays = list(counts = counts))
#' 
#' # Filter to keep only features present in at least 50% of samples
#' # and with at least 1% relative abundance
#' filtered_exp <- filter_taxa(fgt_exp, min_prevalence = 0.5, min_abundance = 0.01)
filter_taxa <- function(fgt_exp, min_prevalence = 0.1, min_abundance = 0.001, 
                        assay_name = "counts") {
  # Input validation
  if (!is(fgt_exp, "FGTExperiment")) {
    stop("Input must be an FGTExperiment object")
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(fgt_exp)) {
    stop(paste0("Assay '", assay_name, "' not found in the FGTExperiment object"))
  }
  
  # Extract count matrix
  counts <- SummarizedExperiment::assays(fgt_exp)[[assay_name]]
  
  # Calculate prevalence (proportion of samples where taxa are present)
  prevalence <- rowSums(counts > 0) / ncol(counts)
  
  # Calculate relative abundance (per sample, then take mean)
  rel_abundances <- t(t(counts) / colSums(counts))
  rel_abundances[is.na(rel_abundances)] <- 0
  mean_rel_abundance <- rowMeans(rel_abundances)
  
  # Apply filters
  keep_taxa <- prevalence >= min_prevalence & mean_rel_abundance >= min_abundance
  
  # Subset the object
  fgt_exp_filtered <- fgt_exp[keep_taxa, ]
  
  # Add filtering information to metadata
  metadata(fgt_exp_filtered)$filter_info <- list(
    min_prevalence = min_prevalence,
    min_abundance = min_abundance,
    original_taxa_count = nrow(fgt_exp),
    filtered_taxa_count = sum(keep_taxa),
    assay_used = assay_name
  )
  
  return(fgt_exp_filtered)
}

#' Transform abundance data
#'
#' Transforms count data in an FGTExperiment object to different abundance metrics.
#'
#' @param fgt_exp An FGTExperiment object
#' @param type Type of transformation: "relative", "clr" (centered log-ratio), 
#'             "log" (natural log), or "presence" (binary presence/absence)
#' @param assay_name Name of the count assay to use
#' @param new_assay_name Name for the new assay (defaults to value of 'type')
#' @param pseudocount Pseudocount to add for log and clr transformations
#'
#' @return FGTExperiment with additional assay
#' @export
#'
#' @examples
#' # Create a simple FGTExperiment
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' 
#' fgt_exp <- FGTExperiment(assays = list(counts = counts))
#' 
#' # Transform to relative abundance
#' fgt_exp <- transform_abundance(fgt_exp, type = "relative")
#' 
#' # Access the transformed data
#' rel_abundance <- assays(fgt_exp)[["relative"]]
transform_abundance <- function(fgt_exp, type = "relative", assay_name = "counts",
                              new_assay_name = NULL, pseudocount = 1) {
  # Input validation
  if (!is(fgt_exp, "FGTExperiment")) {
    stop("Input must be an FGTExperiment object")
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(fgt_exp)) {
    stop(paste0("Assay '", assay_name, "' not found in the FGTExperiment object"))
  }
  
  if (!type %in% c("relative", "clr", "log", "presence")) {
    stop("Type must be one of: 'relative', 'clr', 'log', 'presence'")
  }
  
  # Extract count matrix
  counts <- SummarizedExperiment::assays(fgt_exp)[[assay_name]]
  
  # Set new assay name if not provided
  if (is.null(new_assay_name)) {
    new_assay_name <- type
  }
  
  # Apply transformation with R implementations
  transformed <- switch(
    type,
    "relative" = {
      # Relative abundance
      rel_abundance <- t(t(counts) / colSums(counts))
      rel_abundance[is.na(rel_abundance)] <- 0
      rel_abundance
    },
    "clr" = {
      # Centered log-ratio transformation
      log_counts <- log(counts + pseudocount)
      t(t(log_counts) - colMeans(log_counts))
    },
    "log" = {
      # Log transformation
      log(counts + pseudocount)
    },
    "presence" = {
      # Presence/absence
      matrix(as.numeric(counts > 0), nrow = nrow(counts), ncol = ncol(counts),
             dimnames = dimnames(counts))
    }
  )
  
  # Add new assay to object
  SummarizedExperiment::assays(fgt_exp)[[new_assay_name]] <- transformed
  
  # Add transformation information to metadata
  if (is.null(metadata(fgt_exp)$transformations)) {
    metadata(fgt_exp)$transformations <- list()
  }
  
  metadata(fgt_exp)$transformations[[new_assay_name]] <- list(
    type = type,
    source_assay = assay_name,
    pseudocount = if(type %in% c("clr", "log")) pseudocount else NULL,
    date = Sys.time()
  )
  
  return(fgt_exp)
}

#' Aggregate taxa at a specific taxonomic rank
#'
#' Aggregates feature abundances at a specified taxonomic rank (e.g., Genus, Family).
#'
#' @param fgt_exp An FGTExperiment object
#' @param rank Taxonomic rank to aggregate at (e.g., "Genus", "Family")
#' @param assay_name Name of the assay to use for aggregation
#' @param empty_label Label to use for features without classification at the specified rank
#'
#' @return An FGTExperiment object with aggregated features
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data with taxonomy
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' 
#' taxonomy <- data.frame(
#'   Phylum = sample(c("Firmicutes", "Bacteroidetes"), 10, replace = TRUE),
#'   Genus = sample(c("Lactobacillus", "Gardnerella"), 10, replace = TRUE),
#'   row.names = rownames(counts)
#' )
#' 
#' fgt_exp <- FGTExperiment(
#'   assays = list(counts = counts),
#'   rowData = taxonomy
#' )
#' 
#' # Aggregate counts at genus level
#' genus_level <- aggregate_taxa(fgt_exp, rank = "Genus")
#' }
aggregate_taxa <- function(fgt_exp, rank, assay_name = "counts", empty_label = "Unclassified") {
  # Input validation
  if (!is(fgt_exp, "FGTExperiment")) {
    stop("Input must be an FGTExperiment object")
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(fgt_exp)) {
    stop(paste0("Assay '", assay_name, "' not found in the FGTExperiment object"))
  }
  
  # Check if rank column exists in rowData
  rowdata <- SummarizedExperiment::rowData(fgt_exp)
  if (!rank %in% colnames(rowdata)) {
    stop(paste0("Taxonomic rank '", rank, "' not found in rowData"))
  }
  
  # Extract count matrix and taxonomic information
  counts <- SummarizedExperiment::assays(fgt_exp)[[assay_name]]
  taxa_info <- rowdata[[rank]]
  
  # Replace NA or empty values with Unclassified
  taxa_info[is.na(taxa_info) | taxa_info == ""] <- empty_label
  
  # Get unique taxa at the specified rank
  unique_taxa <- unique(taxa_info)
  
  # Create new count matrix for aggregated taxa
  agg_counts <- matrix(0, nrow = length(unique_taxa), ncol = ncol(counts))
  rownames(agg_counts) <- unique_taxa
  colnames(agg_counts) <- colnames(counts)
  
  # Aggregate counts by taxonomic rank
  for (i in seq_along(unique_taxa)) {
    taxon <- unique_taxa[i]
    rows_to_sum <- which(taxa_info == taxon)
    agg_counts[i, ] <- colSums(counts[rows_to_sum, , drop = FALSE])
  }
  
  # Create new rowData for aggregated taxa
  agg_rowdata <- data.frame(row.names = unique_taxa)
  agg_rowdata[[rank]] <- unique_taxa
  
  # Create new FGTExperiment with aggregated data
  agg_fgt_exp <- FGTExperiment(
    assays = list(counts = agg_counts),
    rowData = agg_rowdata,
    colData = SummarizedExperiment::colData(fgt_exp),
    experimentType = experimentType(fgt_exp)
  )
  
  # Add aggregation information to metadata
  metadata(agg_fgt_exp)$aggregation_info <- list(
    original_feature_count = nrow(fgt_exp),
    aggregated_feature_count = length(unique_taxa),
    rank = rank,
    assay_used = assay_name,
    date = Sys.time()
  )
  
  return(agg_fgt_exp)
}

#' Get taxonomic ranks available in an FGTExperiment
#'
#' @param fgt_exp An FGTExperiment object
#'
#' @return Character vector of available taxonomic ranks
#' @export
#'
#' @examples
#' \dontrun{
#' # Get available ranks
#' ranks <- get_taxonomic_ranks(fgt_exp)
#' }
get_taxonomic_ranks <- function(fgt_exp) {
  # Input validation
  if (!is(fgt_exp, "FGTExperiment")) {
    stop("Input must be an FGTExperiment object")
  }
  
  # Get rowData column names that are likely taxonomic ranks
  rowdata <- SummarizedExperiment::rowData(fgt_exp)
  
  if (ncol(rowdata) == 0) {
    return(character(0))
  }
  
  # Standard taxonomic ranks to look for
  standard_ranks <- c(
    "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    "kingdom", "phylum", "class", "order", "family", "genus", "species"
  )
  
  # Find which columns match standard ranks
  rank_cols <- colnames(rowdata)[colnames(rowdata) %in% standard_ranks]
  
  # If no standard ranks found, return all column names
  if (length(rank_cols) == 0) {
    return(colnames(rowdata))
  }
  
  return(rank_cols)
}