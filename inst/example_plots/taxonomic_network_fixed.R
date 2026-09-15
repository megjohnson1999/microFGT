#' Taxonomic Network Visualization
#' 
#' This script demonstrates how to create and visualize co-occurrence
#' networks from microbiome data using the microFGT package.

# Create output directory if it doesn't exist
dir.create("inst/example_plots/images", recursive = TRUE, showWarnings = FALSE)

# Load required libraries
library(microFGT)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(tidyr)

# Helper function to create mock FGT experiment
create_mock_fgt_experiment <- function(counts, sample_metadata, taxonomy, experiment_name) {
  # Create a simple struct to mimic the FGTExperiment object
  fgt_exp <- list(
    counts = counts,
    metadata = sample_metadata,
    taxonomy = taxonomy,
    experiment_name = experiment_name
  )
  
  # Define a simple class
  class(fgt_exp) <- "mock_fgt_experiment"
  
  return(fgt_exp)
}

# Function for transforming abundance in mock FGT experiments
transformAbundance <- function(fgt_exp, type = "relative") {
  if (inherits(fgt_exp, "mock_fgt_experiment")) {
    counts_data <- fgt_exp$counts
    
    if (type == "relative") {
      # Calculate relative abundance
      counts_total <- colSums(counts_data)
      transformed <- t(t(counts_data) / counts_total)
    } else if (type == "log") {
      # Log transformation (log10(x+1))
      transformed <- log10(counts_data + 1)
    } else if (type == "clr") {
      # CLR transformation
      transformed <- counts_data
      for (i in 1:ncol(counts_data)) {
        # Replace zeros with small value
        sample_data <- counts_data[, i]
        sample_data[sample_data == 0] <- min(sample_data[sample_data > 0]) / 2
        
        # Calculate geometric mean
        geom_mean <- exp(mean(log(sample_data)))
        
        # Apply CLR transformation
        transformed[, i] <- log(sample_data / geom_mean)
      }
    } else {
      # Default - no transformation
      transformed <- counts_data
    }
    
    # Create new object with transformed counts
    result <- fgt_exp
    result$counts <- transformed
    return(result)
  } else {
    # Try to handle real FGTExperiment objects
    stop("Real FGTExperiment not implemented")
  }
}

# Generate mock data
generate_mock_data <- function() {
  set.seed(42)
  
  # Create sample metadata
  n_samples <- 100
  metadata <- data.frame(
    sample_id = paste0("S", 1:n_samples),
    group = sample(c("Group A", "Group B", "Group C"), n_samples, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Create taxonomy
  n_taxa <- 50
  genera <- c("Lactobacillus", "Gardnerella", "Atopobium", "Prevotella", "Sneathia", 
              "Megasphaera", "Streptococcus", "Bifidobacterium", "Ureaplasma", "Mobiluncus")
  
  taxonomy <- data.frame(
    taxon_id = paste0("Taxon_", 1:n_taxa),
    Kingdom = "Bacteria",
    Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Fusobacteria"), n_taxa, replace = TRUE),
    Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria", "Fusobacteriia"), n_taxa, replace = TRUE),
    Order = sample(c("Lactobacillales", "Clostridiales", "Bacteroidales", "Bifidobacteriales", "Fusobacteriales"), n_taxa, replace = TRUE),
    Family = sample(c("Lactobacillaceae", "Bifidobacteriaceae", "Prevotellaceae", "Streptococcaceae"), n_taxa, replace = TRUE),
    Genus = sample(genera, n_taxa, replace = TRUE),
    Species = paste0(sample(genera, n_taxa, replace = TRUE), " ", sample(letters, n_taxa, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  
  # Generate correlated counts
  generate_correlated_counts <- function(n_samples, n_taxa, n_blocks = 5) {
    # Create empty count matrix
    counts <- matrix(0, nrow = n_taxa, ncol = n_samples)
    
    # Create blocks of correlated taxa
    taxa_per_block <- floor(n_taxa / n_blocks)
    remaining <- n_taxa - (taxa_per_block * n_blocks)
    
    block_sizes <- rep(taxa_per_block, n_blocks)
    if (remaining > 0) {
      block_sizes[1:remaining] <- block_sizes[1:remaining] + 1
    }
    
    start_idx <- 1
    for (b in 1:n_blocks) {
      end_idx <- start_idx + block_sizes[b] - 1
      if (end_idx > n_taxa) end_idx <- n_taxa
      
      block_taxa <- start_idx:end_idx
      
      # Generate positively correlated base pattern
      base_pattern <- rexp(n_samples, rate = 1/1000)
      
      # Add variation for each taxon in the block
      for (t in block_taxa) {
        # Base pattern plus noise
        taxon_pattern <- base_pattern * (runif(1, 0.8, 1.2)) + rnorm(n_samples, 0, 100)
        counts[t, ] <- pmax(0, round(taxon_pattern))  # Ensure non-negative integers
      }
      
      start_idx <- end_idx + 1
      if (start_idx > n_taxa) break
    }
    
    # Add anti-correlated patterns
    if (n_blocks >= 2) {
      # Choose a few pairs of blocks to be anti-correlated
      for (i in 1:min(2, n_blocks-1)) {
        block1 <- sample(1:n_blocks, 1)
        remaining_blocks <- (1:n_blocks)[-block1]
        block2 <- sample(remaining_blocks, 1)
        
        # Get taxa indices for these blocks
        block1_start <- 1 + sum(block_sizes[1:(block1-1)])
        block1_end <- block1_start + block_sizes[block1] - 1
        block1_taxa <- block1_start:block1_end
        
        block2_start <- 1 + sum(block_sizes[1:(block2-1)])
        block2_end <- block2_start + block_sizes[block2] - 1
        block2_taxa <- block2_start:block2_end
        
        # Create anti-correlation by inverse scaling
        for (t1 in sample(block1_taxa, min(3, length(block1_taxa)))) {
          for (t2 in sample(block2_taxa, min(3, length(block2_taxa)))) {
            total <- rowSums(counts[c(t1, t2), ])
            anti_factor <- runif(1, 0.7, 0.9)
            counts[t1, ] <- round(total * anti_factor)
            counts[t2, ] <- round(total * (1 - anti_factor))
          }
        }
      }
    }
    
    colnames(counts) <- paste0("Sample_", 1:n_samples)
    rownames(counts) <- paste0("Taxon_", 1:n_taxa)
    
    return(counts)
  }
  
  # Generate counts with correlations
  counts <- generate_correlated_counts(n_samples, n_taxa)
  colnames(counts) <- metadata$sample_id
  
  return(create_mock_fgt_experiment(
    counts = counts,
    sample_metadata = metadata,
    taxonomy = taxonomy,
    experiment_name = "Mock Dataset with Correlated Taxa"
  ))
}

# Generate the data
fgt_exp <- generate_mock_data()

# Function to create taxa correlation network
create_taxa_network <- function(fgt_exp, correlation_threshold = 0.6, p_value_threshold = 0.05) {
  # Convert to relative abundance
  rel_abund <- transformAbundance(fgt_exp, type = "relative")
  
  # Get abundance data
  counts_data <- rel_abund$counts
  
  # Filter low-abundance taxa
  min_prevalence <- 0.1  # Present in at least 10% of samples
  min_rel_abund <- 0.001  # At least 0.1% relative abundance on average
  
  prevalence <- rowSums(counts_data > 0) / ncol(counts_data)
  mean_abundance <- rowMeans(counts_data)
  
  keep_taxa <- which(prevalence >= min_prevalence & mean_abundance >= min_rel_abund)
  
  if (length(keep_taxa) < 2) {
    stop("Not enough taxa meet the abundance and prevalence thresholds")
  }
  
  filtered_counts <- counts_data[keep_taxa, ]
  
  # Calculate correlation matrix
  cor_matrix <- cor(t(filtered_counts), method = "spearman")
  
  # Calculate p-values
  n <- ncol(filtered_counts)
  t_stat <- cor_matrix * sqrt((n - 2) / (1 - cor_matrix^2))
  p_values <- 2 * pt(-abs(t_stat), df = n - 2)
  
  # Create edge list for significant correlations
  edge_list <- data.frame()
  
  for (i in 1:(nrow(cor_matrix) - 1)) {
    for (j in (i + 1):nrow(cor_matrix)) {
      cor_val <- cor_matrix[i, j]
      p_val <- p_values[i, j]
      
      if (abs(cor_val) >= correlation_threshold && p_val <= p_value_threshold) {
        edge_list <- rbind(edge_list, data.frame(
          from = rownames(filtered_counts)[i],
          to = rownames(filtered_counts)[j],
          correlation = cor_val,
          p_value = p_val,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (nrow(edge_list) == 0) {
    stop("No significant correlations found with the specified thresholds")
  }
  
  # Get taxonomy information for nodes
  taxa_info <- fgt_exp$taxonomy
  node_list <- data.frame(
    taxon_id = rownames(filtered_counts),
    abundance = rowMeans(filtered_counts),
    stringsAsFactors = FALSE
  )
  
  # Add taxonomy information
  node_list <- merge(node_list, taxa_info, by = "taxon_id", all.x = TRUE)
  
  # Create network
  g <- graph_from_data_frame(edge_list, vertices = node_list, directed = FALSE)
  
  # Add edge attributes
  E(g)$weight <- abs(edge_list$correlation)
  E(g)$correlation <- edge_list$correlation
  E(g)$type <- ifelse(edge_list$correlation > 0, "positive", "negative")
  
  # Add node attributes
  V(g)$abundance <- node_list$abundance
  V(g)$genus <- node_list$Genus
  V(g)$phylum <- node_list$Phylum
  
  return(g)
}

# Function to plot taxa network
plot_taxa_network <- function(network, detect_communities = TRUE, node_size_by = "abundance") {
  # Check if network is empty
  if (ecount(network) == 0) {
    stop("Network has no edges")
  }
  
  # Community detection
  if (detect_communities) {
    # Try different algorithms and use the one with highest modularity
    comm_methods <- list(
      cluster_fast_greedy(network),
      cluster_walktrap(network),
      cluster_louvain(network)
    )
    
    # Choose method with highest modularity
    modularity_scores <- sapply(comm_methods, modularity)
    best_method <- which.max(modularity_scores)
    communities <- comm_methods[[best_method]]
    
    # Add community information to nodes
    V(network)$community <- membership(communities)
  } else {
    # Default community - use phylum or genus
    if (all(!is.na(V(network)$phylum))) {
      phylum_levels <- as.integer(as.factor(V(network)$phylum))
      V(network)$community <- phylum_levels
    } else {
      V(network)$community <- 1
    }
  }
  
  # Set node size based on attribute
  if (node_size_by == "abundance" && !is.null(V(network)$abundance)) {
    # Log scale for abundance
    V(network)$size <- 1 + 3 * log10(1 + V(network)$abundance * 1000)
  } else if (node_size_by == "degree") {
    V(network)$size <- 1 + degree(network) / 2
  } else {
    V(network)$size <- 3
  }
  
  # Set node labels
  if (all(!is.na(V(network)$genus))) {
    V(network)$label <- V(network)$genus
  } else {
    V(network)$label <- names(V(network))
  }
  
  # Create the plot
  set.seed(42)  # For reproducible layout
  p <- ggraph(network, layout = "fr") +
    geom_edge_link(aes(width = weight, colour = type, alpha = weight),
                  show.legend = TRUE) +
    scale_edge_colour_manual(values = c("positive" = "forestgreen", "negative" = "firebrick")) +
    scale_edge_width(range = c(0.1, 2)) +
    scale_edge_alpha(range = c(0.1, 1)) +
    geom_node_point(aes(size = V(network)$size, fill = as.factor(V(network)$community)), shape = 21, color = "white") +
    geom_node_text(aes(label = V(network)$label), repel = TRUE, size = 3, alpha = 0.7) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Taxonomic Co-occurrence Network",
      subtitle = "Edge color indicates positive (green) or negative (red) correlation",
      edge_colour = "Correlation",
      edge_width = "Strength",
      edge_alpha = "Strength",
      fill = "Community"
    ) +
    theme_graph() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  return(p)
}

# Create the co-occurrence network
taxa_network <- create_taxa_network(fgt_exp, 
                                   correlation_threshold = 0.5,
                                   p_value_threshold = 0.05)

# Plot the network and save
network_plot <- plot_taxa_network(taxa_network, 
                                detect_communities = TRUE,
                                node_size_by = "abundance")

# Save the plot
ggsave("inst/example_plots/images/taxa_network.png", network_plot, width = 12, height = 10, dpi = 300)

# Create additional visualization with specific layout
network_plot2 <- ggraph(taxa_network, layout = "circle") +
  geom_edge_link(aes(width = weight, colour = type),
                alpha = 0.7, show.legend = TRUE) +
  scale_edge_colour_manual(values = c("positive" = "forestgreen", "negative" = "firebrick")) +
  scale_edge_width(range = c(0.1, 2)) +
  geom_node_point(aes(size = V(taxa_network)$size, fill = as.factor(V(taxa_network)$community)), shape = 21, color = "white") +
  geom_node_text(aes(label = V(taxa_network)$label), repel = TRUE, size = 3, check_overlap = TRUE) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Circular Taxonomic Co-occurrence Network",
    subtitle = "Organized by community structure",
    edge_colour = "Correlation",
    edge_width = "Strength",
    fill = "Community"
  ) +
  theme_graph() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

# Save the second plot
ggsave("inst/example_plots/images/taxa_network_circular.png", network_plot2, width = 12, height = 10, dpi = 300)

# Print success message
cat("Taxonomic network plots created successfully!\n")
cat("Plots saved to inst/example_plots/images/\n")