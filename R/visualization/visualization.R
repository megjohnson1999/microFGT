#' Plot taxonomic composition
#'
#' Creates a stacked bar plot of taxonomic composition.
#'
#' @param fgt_exp An FGTExperiment object
#' @param rank Taxonomic rank to display (e.g., "Genus", "Phylum")
#' @param top_n Number of top taxa to display individually
#' @param group_var Optional grouping variable from colData
#' @param assay_name Name of the assay to use (preferably relative abundance)
#' @param color_palette Optional vector of colors for plotting
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data
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
#' metadata <- data.frame(
#'   group = rep(c("A", "B"), each = 3),
#'   row.names = colnames(counts)
#' )
#' 
#' fgt_exp <- FGTExperiment(
#'   assays = list(counts = counts),
#'   rowData = taxonomy,
#'   colData = metadata
#' )
#' 
#' # Transform to relative abundance
#' fgt_rel <- transform_abundance(fgt_exp, type = "relative")
#' 
#' # Plot taxonomic composition
#' plot_taxa_composition(fgt_rel, rank = "Phylum")
#' 
#' # Group by metadata variable
#' plot_taxa_composition(fgt_rel, rank = "Genus", group_var = "group")
#' }
plot_taxa_composition <- function(fgt_exp, rank, top_n = 10, group_var = NULL,
                                assay_name = "relative", color_palette = NULL) {
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  # Input validation
  if (!is(fgt_exp, "FGTExperiment")) {
    stop("Input must be an FGTExperiment object")
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(fgt_exp)) {
    stop(paste0("Assay '", assay_name, "' not found in the FGTExperiment object"))
  }
  
  # Get taxonomic rank
  rowdata <- SummarizedExperiment::rowData(fgt_exp)
  if (!rank %in% colnames(rowdata)) {
    stop(paste0("Taxonomic rank '", rank, "' not found in rowData"))
  }
  
  # Extract data for plotting
  abundances <- SummarizedExperiment::assays(fgt_exp)[[assay_name]]
  taxa_ranks <- rowdata[[rank]]
  
  # Replace NA or empty values
  taxa_ranks[is.na(taxa_ranks) | taxa_ranks == ""] <- "Unknown"
  
  # Check if data is in relative abundance
  if (assay_name != "relative") {
    col_sums <- colSums(abundances)
    if (any(col_sums > 1.1)) {  # Allow for small floating point differences
      warning("Data doesn't appear to be in relative abundance. Results may be misleading.")
    }
  }
  
  # Prepare data frame for plotting
  plot_data <- data.frame()
  for (sample_idx in 1:ncol(abundances)) {
    sample_id <- colnames(abundances)[sample_idx]
    for (taxon_idx in 1:nrow(abundances)) {
      taxon <- taxa_ranks[taxon_idx]
      abundance <- abundances[taxon_idx, sample_idx]
      
      # Skip zero abundances
      if (abundance > 0) {
        plot_data <- rbind(plot_data, data.frame(
          sample = sample_id,
          taxon = taxon,
          abundance = abundance,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Aggregate by taxon
  agg_data <- stats::aggregate(
    abundance ~ sample + taxon, 
    data = plot_data,
    FUN = sum
  )
  
  # Add group variable if specified
  if (!is.null(group_var)) {
    if (!group_var %in% colnames(SummarizedExperiment::colData(fgt_exp))) {
      stop(paste0("Group variable '", group_var, "' not found in colData"))
    }
    
    # Extract group information
    groups <- SummarizedExperiment::colData(fgt_exp)[[group_var]]
    group_df <- data.frame(
      sample = colnames(abundances),
      group = groups,
      stringsAsFactors = FALSE
    )
    
    # Merge with plot data
    agg_data <- merge(agg_data, group_df, by = "sample")
  }
  
  # Find top taxa
  taxon_totals <- stats::aggregate(
    abundance ~ taxon, 
    data = agg_data,
    FUN = sum
  )
  taxon_totals <- taxon_totals[order(taxon_totals$abundance, decreasing = TRUE), ]
  top_taxa <- taxon_totals$taxon[1:min(top_n, nrow(taxon_totals))]
  
  # Label other taxa
  agg_data$taxon_group <- ifelse(agg_data$taxon %in% top_taxa, 
                               as.character(agg_data$taxon), 
                               "Other")
  
  # Convert to factor with top taxa ordered by abundance
  agg_data$taxon_group <- factor(agg_data$taxon_group,
                               levels = c(top_taxa, "Other"))
  
  # Create plot
  if (is.null(group_var)) {
    p <- ggplot2::ggplot(agg_data, ggplot2::aes(x = sample, y = abundance, 
                                            fill = taxon_group)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(title = paste("Taxonomic composition at", rank, "level"),
                x = "Sample",
                y = "Relative Abundance",
                fill = rank) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    p <- ggplot2::ggplot(agg_data, ggplot2::aes(x = sample, y = abundance, 
                                            fill = taxon_group)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(title = paste("Taxonomic composition at", rank, "level"),
                x = "Sample",
                y = "Relative Abundance",
                fill = rank) +
      ggplot2::facet_grid(. ~ group, scales = "free_x", space = "free") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  # Apply custom color palette if provided
  if (!is.null(color_palette)) {
    if (length(color_palette) < length(levels(agg_data$taxon_group))) {
      warning("Not enough colors in palette. Using default colors.")
    } else {
      p <- p + ggplot2::scale_fill_manual(values = color_palette)
    }
  }
  
  return(p)
}

#' Plot alpha diversity
#'
#' Creates box or violin plots of alpha diversity metrics.
#'
#' @param fgt_exp An FGTExperiment object
#' @param metrics Vector of diversity metrics to calculate ("shannon", "simpson", "richness")
#' @param group_var Optional grouping variable from colData
#' @param assay_name Name of the assay to use (raw counts recommended)
#' @param plot_type Type of plot ("box" or "violin")
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' 
#' metadata <- data.frame(
#'   group = rep(c("A", "B"), each = 3),
#'   row.names = colnames(counts)
#' )
#' 
#' fgt_exp <- FGTExperiment(
#'   assays = list(counts = counts),
#'   colData = metadata
#' )
#' 
#' # Plot alpha diversity
#' plot_alpha_diversity(fgt_exp, metrics = c("shannon"), group_var = "group")
#' }
plot_alpha_diversity <- function(fgt_exp, metrics = c("shannon"), group_var = NULL,
                               assay_name = "counts", plot_type = "box") {
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  # Input validation
  if (!is(fgt_exp, "FGTExperiment")) {
    stop("Input must be an FGTExperiment object")
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(fgt_exp)) {
    stop(paste0("Assay '", assay_name, "' not found in the FGTExperiment object"))
  }
  
  valid_metrics <- c("shannon", "simpson", "richness")
  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop("Invalid metrics: ", paste(invalid_metrics, collapse = ", "), 
         ". Valid metrics are: ", paste(valid_metrics, collapse = ", "))
  }
  
  plot_type <- match.arg(plot_type, c("box", "violin"))
  
  # Extract count matrix
  counts <- SummarizedExperiment::assays(fgt_exp)[[assay_name]]
  
  # Calculate diversity metrics
  diversity_df <- data.frame(
    sample = colnames(counts),
    stringsAsFactors = FALSE
  )
  
  # Calculate Shannon diversity
  if ("shannon" %in% metrics) {
    # Convert to relative abundance if not already
    rel_counts <- t(t(counts) / colSums(counts))
    rel_counts[is.na(rel_counts)] <- 0
    
    # Shannon index: -sum(p_i * log(p_i))
    shannon <- apply(rel_counts, 2, function(x) {
      p <- x[x > 0]  # Remove zeros
      -sum(p * log(p))
    })
    
    diversity_df$shannon <- shannon
  }
  
  # Calculate Simpson diversity
  if ("simpson" %in% metrics) {
    # Convert to relative abundance if not already
    rel_counts <- t(t(counts) / colSums(counts))
    rel_counts[is.na(rel_counts)] <- 0
    
    # Simpson index: 1 - sum(p_i^2)
    simpson <- apply(rel_counts, 2, function(x) {
      1 - sum(x^2)
    })
    
    diversity_df$simpson <- simpson
  }
  
  # Calculate richness (number of observed features)
  if ("richness" %in% metrics) {
    richness <- colSums(counts > 0)
    diversity_df$richness <- richness
  }
  
  # Add group variable if specified
  if (!is.null(group_var)) {
    if (!group_var %in% colnames(SummarizedExperiment::colData(fgt_exp))) {
      stop(paste0("Group variable '", group_var, "' not found in colData"))
    }
    
    # Extract group information
    groups <- SummarizedExperiment::colData(fgt_exp)[[group_var]]
    diversity_df$group <- groups
  }
  
  # Convert to long format for plotting
  plot_data <- stats::reshape(
    diversity_df,
    direction = "long",
    varying = metrics,
    v.names = "value",
    timevar = "metric",
    times = metrics,
    idvar = "sample"
  )
  
  # Create plot
  if (is.null(group_var)) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = metric, y = value)) +
      ggplot2::labs(title = "Alpha Diversity",
                x = "Metric",
                y = "Value") +
      ggplot2::theme_minimal()
    
    if (plot_type == "box") {
      p <- p + ggplot2::geom_boxplot()
    } else {
      p <- p + ggplot2::geom_violin() + ggplot2::geom_jitter(width = 0.1, alpha = 0.5)
    }
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = value, fill = group)) +
      ggplot2::labs(title = "Alpha Diversity",
                x = group_var,
                y = "Value") +
      ggplot2::facet_wrap(~ metric, scales = "free_y") +
      ggplot2::theme_minimal()
    
    if (plot_type == "box") {
      p <- p + ggplot2::geom_boxplot()
    } else {
      p <- p + ggplot2::geom_violin() + ggplot2::geom_jitter(width = 0.1, alpha = 0.5)
    }
  }
  
  return(p)
}