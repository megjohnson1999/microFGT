#' Calculate alpha diversity metrics
#'
#' Calculates various alpha diversity metrics for microbiome data.
#'
#' @param x A matrix, SummarizedExperiment or FGTExperiment object
#' @param method Diversity method: "shannon", "simpson", "richness", "evenness"
#' @param assay_name Name of the assay to use (for SummarizedExperiment objects)
#' @param normalize Logical indicating whether to normalize counts before calculation
#'
#' @return A named vector of diversity values for each sample
#' @export
calculate_diversity <- function(x, method = c("shannon", "simpson", "richness", "evenness"),
                               assay_name = "counts", normalize = TRUE) {
  
  # Match the method argument
  method <- match.arg(method)
  
  # Extract the abundance matrix if x is a SummarizedExperiment
  if (methods::is(x, "SummarizedExperiment")) {
    if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
      stop("Assay '", assay_name, "' not found in the object")
    }
    abundance <- SummarizedExperiment::assays(x)[[assay_name]]
  } else if (is.matrix(x)) {
    abundance <- x
  } else {
    stop("Input must be a matrix or SummarizedExperiment object")
  }
  
  # Normalize if requested
  if (normalize) {
    # Convert to relative abundance
    abundance <- t(apply(abundance, 2, function(col) {
      if (sum(col, na.rm = TRUE) == 0) {
        return(rep(0, length(col)))
      }
      return(col / sum(col, na.rm = TRUE))
    }))
  }
  
  # Calculate the requested diversity metric
  diversity_values <- switch(
    method,
    "shannon" = {
      # Shannon diversity: -sum(p_i * log(p_i))
      -colSums(ifelse(abundance > 0, abundance * log(abundance), 0), na.rm = TRUE)
    },
    "simpson" = {
      # Simpson diversity: 1 - sum(p_i^2)
      1 - colSums(abundance^2, na.rm = TRUE)
    },
    "richness" = {
      # Richness: number of non-zero values
      colSums(abundance > 0, na.rm = TRUE)
    },
    "evenness" = {
      # Pielou's evenness: shannon / log(richness)
      shannon <- -colSums(ifelse(abundance > 0, abundance * log(abundance), 0), na.rm = TRUE)
      richness <- colSums(abundance > 0, na.rm = TRUE)
      shannon / log(richness)
    }
  )
  
  return(diversity_values)
}

#' Calculate beta diversity metrics
#'
#' Calculates pairwise beta diversity metrics between samples.
#'
#' @param x A matrix, SummarizedExperiment or FGTExperiment object
#' @param method Distance method: "bray", "jaccard", "unifrac"
#' @param assay_name Name of the assay to use (for SummarizedExperiment objects)
#' @param normalize Logical indicating whether to normalize counts before calculation
#' @param tree Phylogenetic tree (needed for UniFrac)
#'
#' @return A distance matrix of beta diversity values
#' @export
calculate_beta_diversity <- function(x, method = c("bray", "jaccard", "unifrac"),
                                    assay_name = "counts", normalize = TRUE, tree = NULL) {
  
  # Match the method argument
  method <- match.arg(method)
  
  # Extract the abundance matrix if x is a SummarizedExperiment
  if (methods::is(x, "SummarizedExperiment")) {
    if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
      stop("Assay '", assay_name, "' not found in the object")
    }
    abundance <- SummarizedExperiment::assays(x)[[assay_name]]
  } else if (is.matrix(x)) {
    abundance <- x
  } else {
    stop("Input must be a matrix or SummarizedExperiment object")
  }
  
  # Normalize if requested
  if (normalize) {
    # Convert to relative abundance
    abundance <- t(apply(abundance, 2, function(col) {
      if (sum(col, na.rm = TRUE) == 0) {
        return(rep(0, length(col)))
      }
      return(col / sum(col, na.rm = TRUE))
    }))
  }
  
  # Extract tree if available and required
  if (method == "unifrac" && is.null(tree)) {
    if (methods::is(x, "TreeSummarizedExperiment")) {
      tree <- TreeSummarizedExperiment::rowTree(x)
      if (is.null(tree)) {
        stop("UniFrac requires a phylogenetic tree, but none was provided and none found in the object")
      }
    } else {
      stop("UniFrac requires a phylogenetic tree")
    }
  }
  
  # Calculate the requested beta diversity metric
  dist_matrix <- switch(
    method,
    "bray" = {
      # Bray-Curtis dissimilarity
      bray_curtis_distance(abundance)
    },
    "jaccard" = {
      # Jaccard distance
      jaccard_distance(abundance)
    },
    "unifrac" = {
      # UniFrac distance - requires a tree
      if (requireNamespace("phyloseq", quietly = TRUE)) {
        unifrac_distance(abundance, tree)
      } else {
        stop("UniFrac calculation requires the phyloseq package")
      }
    }
  )
  
  return(dist_matrix)
}

# Helper function for Bray-Curtis dissimilarity
bray_curtis_distance <- function(abundance) {
  n_samples <- ncol(abundance)
  dist_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  colnames(dist_matrix) <- colnames(abundance)
  rownames(dist_matrix) <- colnames(abundance)
  
  for (i in 1:(n_samples-1)) {
    for (j in (i+1):n_samples) {
      # Bray-Curtis formula: sum(abs(x_i - y_i)) / sum(x_i + y_i)
      sample_i <- abundance[, i]
      sample_j <- abundance[, j]
      
      numerator <- sum(abs(sample_i - sample_j), na.rm = TRUE)
      denominator <- sum(sample_i + sample_j, na.rm = TRUE)
      
      if (denominator == 0) {
        dist_matrix[i, j] <- dist_matrix[j, i] <- 0
      } else {
        bray_dist <- numerator / denominator
        dist_matrix[i, j] <- dist_matrix[j, i] <- bray_dist
      }
    }
  }
  
  return(dist_matrix)
}

# Helper function for Jaccard distance
jaccard_distance <- function(abundance) {
  # Convert to presence/absence
  presence <- abundance > 0
  
  n_samples <- ncol(presence)
  dist_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  colnames(dist_matrix) <- colnames(presence)
  rownames(dist_matrix) <- colnames(presence)
  
  for (i in 1:(n_samples-1)) {
    for (j in (i+1):n_samples) {
      # Jaccard formula: 1 - (A ∩ B) / (A ∪ B)
      sample_i <- presence[, i]
      sample_j <- presence[, j]
      
      shared <- sum(sample_i & sample_j, na.rm = TRUE)
      total <- sum(sample_i | sample_j, na.rm = TRUE)
      
      if (total == 0) {
        dist_matrix[i, j] <- dist_matrix[j, i] <- 0
      } else {
        jaccard_dist <- 1 - (shared / total)
        dist_matrix[i, j] <- dist_matrix[j, i] <- jaccard_dist
      }
    }
  }
  
  return(dist_matrix)
}

# Helper function for UniFrac distance (requires phyloseq)
unifrac_distance <- function(abundance, tree) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("UniFrac calculation requires the phyloseq package")
  }
  
  # Create a phyloseq object
  otu_table <- phyloseq::otu_table(t(abundance), taxa_are_rows = FALSE)
  
  if (!methods::is(tree, "phylo")) {
    stop("Tree must be a phylo object")
  }
  
  phy_tree <- tree
  
  # Create phyloseq object
  ps <- phyloseq::phyloseq(otu_table, phy_tree)
  
  # Calculate UniFrac distance
  unifrac_dist <- phyloseq::distance(ps, method = "unifrac")
  
  # Convert to matrix
  dist_matrix <- as.matrix(unifrac_dist)
  
  return(dist_matrix)
}