#!/usr/bin/env Rscript

# Direct Solution for microFGT Example Data
# This standalone script can be run directly without relying on package exports

#' Direct Example Data Loader and Transformer
#' 
#' This script provides direct functions for:
#' 1. Loading example data from the microFGT package
#' 2. Transforming abundance data with pure R implementation
#' 3. Creating simple visualizations
#'
#' Copy this entire file to your working directory and source it:
#' source("direct_solution.R")
#' 
#' Then use the functions directly:
#' fgt_data <- direct_load_example()
#' rel_data <- direct_transform_abundance(fgt_data, "relative")

# Load required packages
library(microFGT)
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("Package 'SummarizedExperiment' is required. Please install with BiocManager::install('SummarizedExperiment')")
}
if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
  stop("Package 'TreeSummarizedExperiment' is required. Please install with BiocManager::install('TreeSummarizedExperiment')")
}
if (!requireNamespace("S4Vectors", quietly = TRUE)) {
  stop("Package 'S4Vectors' is required. Please install with BiocManager::install('S4Vectors')")
}

#' Direct example data loader
#' 
#' @param as_fgt_experiment Convert to FGTExperiment (default TRUE)
#' @return FGTExperiment object or list of data components
direct_load_example <- function(as_fgt_experiment = TRUE) {
  # Find package data directory
  data_dir <- system.file("extdata", package = "microFGT")
  if (data_dir == "") {
    stop("Package 'microFGT' not properly installed or example data not found.")
  }
  
  # File paths
  count_file <- file.path(data_dir, "microFGT_example_small_amplicon_counts.rds")
  tax_file <- file.path(data_dir, "microFGT_example_small_amplicon_taxonomy.rds")
  meta_file <- file.path(data_dir, "microFGT_example_small_amplicon_metadata.rds")
  tree_file <- file.path(data_dir, "microFGT_example_small_amplicon_tree.rds")
  
  # Verify files exist
  if (!file.exists(count_file)) stop("Count file not found: ", count_file)
  if (!file.exists(tax_file)) stop("Taxonomy file not found: ", tax_file)
  if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)
  if (!file.exists(tree_file)) stop("Tree file not found: ", tree_file)
  
  # Load the data
  counts <- readRDS(count_file)
  taxonomy <- readRDS(tax_file)
  metadata <- readRDS(meta_file)
  tree <- readRDS(tree_file)
  
  message("Successfully loaded example data")
  
  # Return as list if not converting to FGTExperiment
  if (!as_fgt_experiment) {
    return(list(
      counts = counts,
      taxonomy = taxonomy, 
      metadata = metadata,
      tree = tree
    ))
  }
  
  # Create a SummarizedExperiment
  message("Creating SummarizedExperiment...")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = taxonomy,
    colData = metadata
  )
  
  # Convert to TreeSummarizedExperiment
  message("Converting to TreeSummarizedExperiment...")
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(se, rowTree = tree)
  
  # Try to create FGTExperiment
  message("Creating FGTExperiment...")
  
  fgt <- tryCatch({
    # Try to create an FGTExperiment object
    obj <- methods::new("FGTExperiment", 
                      tse,
                      experimentType = "amplicon",
                      fgtMetadata = S4Vectors::SimpleList())
    
    # Set assay names
    SummarizedExperiment::assayNames(obj)[1] <- "counts"
    obj
  }, error = function(e) {
    message("Error creating FGTExperiment: ", conditionMessage(e))
    message("Returning TreeSummarizedExperiment instead")
    tse
  })
  
  return(fgt)
}

#' Direct transform abundance
#' 
#' @param obj FGTExperiment or TreeSummarizedExperiment object
#' @param type Transformation type: "relative", "clr", "log", or "presence"
#' @param assay_name Name of assay to transform
#' @param pseudocount Value to add before log transformation
#' @return Object with transformed assay added
direct_transform_abundance <- function(obj, type = "relative", 
                                     assay_name = "counts",
                                     pseudocount = 1) {
  # Check input
  if (!any(c("SummarizedExperiment", "TreeSummarizedExperiment", "FGTExperiment") %in% 
           class(obj))) {
    stop("Input must be a SummarizedExperiment, TreeSummarizedExperiment, or FGTExperiment")
  }
  
  # Check assay exists
  if (!assay_name %in% SummarizedExperiment::assayNames(obj)) {
    stop("Assay '", assay_name, "' not found in object")
  }
  
  # Check type
  if (!type %in% c("relative", "clr", "log", "presence")) {
    stop("Type must be one of: 'relative', 'clr', 'log', 'presence'")
  }
  
  # Get count matrix
  counts <- SummarizedExperiment::assays(obj)[[assay_name]]
  
  # Apply transformation using pure R
  transformed <- switch(
    type,
    "relative" = {
      # Relative abundance
      rel_abundance <- t(t(counts) / colSums(counts))
      rel_abundance[is.na(rel_abundance)] <- 0
      rel_abundance
    },
    "clr" = {
      # Centered log-ratio
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
  
  # Add transformed data as new assay
  SummarizedExperiment::assays(obj)[[type]] <- transformed
  
  message("Added transformed data as assay '", type, "'")
  return(obj)
}

#' Direct plot taxa composition
#'
#' @param obj FGTExperiment or similar object
#' @param assay_name Name of assay to use (should be relative abundance)
#' @param rank Taxonomic rank to aggregate at
#' @param top_n Number of top taxa to show individually
#' @return ggplot object
direct_plot_taxa <- function(obj, assay_name = "relative", rank = "Phylum", top_n = 5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for data manipulation")
  }
  
  # Get abundance data
  counts <- SummarizedExperiment::assays(obj)[[assay_name]]
  
  # Get taxonomy
  taxa <- SummarizedExperiment::rowData(obj)
  
  # Check if rank exists
  if (!rank %in% colnames(taxa)) {
    stop("Rank '", rank, "' not found in taxonomy")
  }
  
  # Create data frame for plotting
  plot_data <- data.frame()
  
  for (i in 1:ncol(counts)) {
    sample_id <- colnames(counts)[i]
    sample_data <- data.frame(
      Sample = sample_id,
      Taxon = taxa[[rank]],
      Abundance = counts[, i]
    )
    plot_data <- rbind(plot_data, sample_data)
  }
  
  # Aggregate by taxon
  agg_data <- dplyr::group_by(plot_data, Sample, Taxon)
  agg_data <- dplyr::summarize(agg_data, Abundance = sum(Abundance))
  
  # Calculate total abundance per sample
  totals <- dplyr::group_by(agg_data, Sample)
  totals <- dplyr::summarize(totals, Total = sum(Abundance))
  
  # Join to get relative abundance
  agg_data <- dplyr::left_join(agg_data, totals, by = "Sample")
  agg_data$RelativeAbundance <- agg_data$Abundance / agg_data$Total
  
  # Calculate mean abundance per taxon
  taxon_means <- dplyr::group_by(agg_data, Taxon)
  taxon_means <- dplyr::summarize(taxon_means, MeanAbundance = mean(RelativeAbundance))
  taxon_means <- dplyr::arrange(taxon_means, desc(MeanAbundance))
  
  # Get top taxa
  top_taxa <- dplyr::slice(taxon_means, 1:top_n)$Taxon
  
  # Label other taxa
  agg_data$TaxonGroup <- ifelse(agg_data$Taxon %in% top_taxa, 
                              as.character(agg_data$Taxon), "Other")
  agg_data$TaxonGroup <- factor(agg_data$TaxonGroup, 
                               levels = c(top_taxa, "Other"))
  
  # Create plot
  p <- ggplot2::ggplot(agg_data, ggplot2::aes(x = Sample, y = RelativeAbundance, 
                                          fill = TaxonGroup)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(title = paste("Taxonomic Composition at", rank, "Level"),
                x = "Sample", y = "Relative Abundance", fill = rank) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Usage examples (commented out)
# 
# # Load example data
# fgt_data <- direct_load_example()
# 
# # Transform to relative abundance
# fgt_rel <- direct_transform_abundance(fgt_data, "relative")
# 
# # Plot taxa composition
# plot <- direct_plot_taxa(fgt_rel)
# print(plot)

message("Direct solution functions loaded successfully.")
message("Use direct_load_example() to load example data")
message("Use direct_transform_abundance() to transform data")
message("Use direct_plot_taxa() to create taxonomy plots")