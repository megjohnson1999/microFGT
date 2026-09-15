#' Microbiome Fingerprints Visualization
#'
#' This script demonstrates how to create unique visual "fingerprints" 
#' for microbiome samples that can be used for quick visual identification
#' and comparison.

# Create output directory if it doesn't exist
dir.create("inst/example_plots/images", recursive = TRUE, showWarnings = FALSE)

# Load required libraries
library(microFGT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggforce)
library(grid)
library(gridExtra)

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

# Generate mock data
generate_mock_data <- function() {
  set.seed(42)
  
  # Create sample metadata
  n_samples <- 20
  metadata <- data.frame(
    sample_id = paste0("S", 1:n_samples),
    group = sample(c("Healthy", "BV", "Mixed"), n_samples, replace = TRUE),
    subject_id = sample(paste0("Subject_", 1:10), n_samples, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Generate mock counts with distinct patterns by group
  n_taxa <- 50
  counts <- matrix(0, nrow = n_taxa, ncol = n_samples)
  rownames(counts) <- paste0("Taxon_", 1:n_taxa)
  colnames(counts) <- metadata$sample_id
  
  # Create taxonomy
  genera <- c("Lactobacillus", "Gardnerella", "Atopobium", "Prevotella", "Sneathia", 
              "Megasphaera", "Streptococcus", "Bifidobacterium", "Ureaplasma", "Mobiluncus")
  
  taxonomy <- data.frame(
    taxon_id = rownames(counts),
    Kingdom = "Bacteria",
    Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Fusobacteria"), n_taxa, replace = TRUE),
    Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria", "Fusobacteriia"), n_taxa, replace = TRUE),
    Order = sample(c("Lactobacillales", "Clostridiales", "Bacteroidales", "Bifidobacteriales", "Fusobacteriales"), n_taxa, replace = TRUE),
    Family = sample(c("Lactobacillaceae", "Bifidobacteriaceae", "Prevotellaceae", "Streptococcaceae"), n_taxa, replace = TRUE),
    Genus = sample(genera, n_taxa, replace = TRUE),
    Species = paste0(sample(genera, n_taxa, replace = TRUE), " ", sample(letters, n_taxa, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  
  # Generate counts based on group
  for (i in 1:n_samples) {
    group <- metadata$group[i]
    
    # Base pattern with low counts for all taxa
    counts[, i] <- sample(0:5, n_taxa, replace = TRUE)
    
    if (group == "Healthy") {
      # Lactobacillus dominated
      lacto_indices <- which(taxonomy$Genus == "Lactobacillus")
      if (length(lacto_indices) > 0) {
        main_lacto <- sample(lacto_indices, 1)
        counts[main_lacto, i] <- sample(5000:20000, 1)
        
        # A few other lactobacilli at lower abundance
        other_lacto <- setdiff(lacto_indices, main_lacto)
        if (length(other_lacto) > 0) {
          for (j in sample(other_lacto, min(2, length(other_lacto)))) {
            counts[j, i] <- sample(100:1000, 1)
          }
        }
      }
    } else if (group == "BV") {
      # Diverse anaerobes
      bv_genera <- c("Gardnerella", "Atopobium", "Prevotella", "Sneathia", "Megasphaera")
      
      for (genus in bv_genera) {
        genus_indices <- which(taxonomy$Genus == genus)
        if (length(genus_indices) > 0) {
          for (j in sample(genus_indices, min(2, length(genus_indices)))) {
            counts[j, i] <- sample(1000:5000, 1)
          }
        }
      }
    } else {  # Mixed
      # Mix of lactobacilli and anaerobes
      lacto_indices <- which(taxonomy$Genus == "Lactobacillus")
      if (length(lacto_indices) > 0) {
        main_lacto <- sample(lacto_indices, 1)
        counts[main_lacto, i] <- sample(2000:10000, 1)
      }
      
      # Some anaerobes
      anaerobe_genera <- c("Gardnerella", "Prevotella")
      for (genus in anaerobe_genera) {
        genus_indices <- which(taxonomy$Genus == genus)
        if (length(genus_indices) > 0) {
          for (j in sample(genus_indices, min(1, length(genus_indices)))) {
            counts[j, i] <- sample(1000:3000, 1)
          }
        }
      }
    }
  }
  
  return(create_mock_fgt_experiment(
    counts = counts,
    sample_metadata = metadata,
    taxonomy = taxonomy,
    experiment_name = "Mock Dataset for Fingerprinting"
  ))
}

# Generate mock data
fgt_exp <- generate_mock_data()

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

# Function to generate unique microbiome fingerprints
generate_microbiome_fingerprint <- function(fgt_exp, method = "circular", n_taxa = 20, 
                                           color_palette = "viridis", sample_subset = NULL) {
  # Convert to relative abundance
  rel_abund <- transformAbundance(fgt_exp, type = "relative")
  
  # Get abundance data
  counts_data <- rel_abund$counts
  
  # Get taxonomy
  taxa_info <- rel_abund$taxonomy
  
  # Get metadata
  metadata <- rel_abund$metadata
  
  # Select sample subset if specified
  if (!is.null(sample_subset)) {
    if (is.numeric(sample_subset)) {
      # Use as count
      sample_ids <- sample(colnames(counts_data), min(sample_subset, ncol(counts_data)))
    } else {
      # Use as vector of sample IDs
      sample_ids <- intersect(sample_subset, colnames(counts_data))
      if (length(sample_ids) == 0) {
        stop("None of the specified sample IDs were found in the dataset")
      }
    }
    counts_data <- counts_data[, sample_ids, drop = FALSE]
    metadata <- metadata[metadata$sample_id %in% sample_ids, , drop = FALSE]
  }
  
  # Select top taxa
  taxa_means <- rowMeans(counts_data)
  top_taxa_idx <- order(taxa_means, decreasing = TRUE)[1:min(n_taxa, nrow(counts_data))]
  top_taxa_data <- counts_data[top_taxa_idx, , drop = FALSE]
  
  # Generate fingerprints
  fingerprints <- list()
  
  if (method == "circular") {
    # Circular fingerprint with taxa as segments
    for (sample_id in colnames(top_taxa_data)) {
      # Extract data for this sample
      sample_data <- top_taxa_data[, sample_id]
      
      # Normalize to sum to 1
      if (sum(sample_data) > 0) {
        sample_data <- sample_data / sum(sample_data)
      }
      
      # Create data frame for plotting
      n <- length(sample_data)
      plot_data <- data.frame(
        taxon = factor(rownames(top_taxa_data), levels = rownames(top_taxa_data)),
        abundance = sample_data,
        genus = taxa_info$Genus[match(rownames(top_taxa_data), taxa_info$taxon_id)],
        start = seq(0, 2*pi*(1-1/n), length.out = n),
        end = seq(2*pi/n, 2*pi, length.out = n)
      )
      
      # Sort by phylum then abundance
      if ("Phylum" %in% colnames(taxa_info)) {
        plot_data$phylum <- taxa_info$Phylum[match(rownames(top_taxa_data), taxa_info$taxon_id)]
        plot_data <- plot_data[order(plot_data$phylum, -plot_data$abundance), ]
      } else {
        plot_data <- plot_data[order(-plot_data$abundance), ]
      }
      
      # Create color palette
      if (color_palette == "viridis") {
        colors <- viridis(n)
      } else if (color_palette == "plasma") {
        colors <- viridis(n, option = "plasma")
      } else if (color_palette == "genus") {
        # Color by genus
        unique_genera <- unique(plot_data$genus)
        colors <- setNames(rainbow(length(unique_genera)), unique_genera)
        colors <- colors[plot_data$genus]
      } else {
        colors <- rainbow(n)
      }
      
      # Create plot
      p <- ggplot() +
        geom_arc_bar(data = plot_data, 
                    aes(x0 = 0, y0 = 0, r0 = 0.5, r = 1 + abundance*2, 
                        start = start, end = end, fill = taxon),
                    color = "white", size = 0.1) +
        scale_fill_manual(values = colors) +
        coord_fixed() +
        theme_void() +
        theme(legend.position = "none",
              plot.background = element_rect(fill = "white", color = NA),
              plot.margin = margin(5, 5, 5, 5)) +
        labs(title = sample_id)
      
      # Store the plot
      fingerprints[[sample_id]] <- p
    }
  } else if (method == "heatmap") {
    # Heatmap fingerprint
    for (sample_id in colnames(top_taxa_data)) {
      # Extract data for this sample
      sample_data <- top_taxa_data[, sample_id]
      
      # Normalize to sum to 1
      if (sum(sample_data) > 0) {
        sample_data <- sample_data / sum(sample_data)
      }
      
      # Create a square matrix from the vector
      size <- ceiling(sqrt(length(sample_data)))
      matrix_data <- matrix(0, nrow = size, ncol = size)
      
      # Fill the matrix (reshape vector to matrix)
      for (i in 1:length(sample_data)) {
        row <- ceiling(i / size)
        col <- ((i - 1) %% size) + 1
        if (row <= size && col <= size) {
          matrix_data[row, col] <- sample_data[i]
        }
      }
      
      # Create data frame for plotting
      plot_data <- expand.grid(x = 1:size, y = 1:size)
      plot_data$value <- as.vector(matrix_data)
      
      # Create plot
      p <- ggplot(plot_data, aes(x = x, y = y, fill = value)) +
        geom_tile(color = "white", size = 0.1) +
        scale_fill_viridis(option = ifelse(color_palette == "viridis", "viridis", "plasma")) +
        coord_fixed() +
        theme_void() +
        theme(legend.position = "none",
              plot.background = element_rect(fill = "white", color = NA),
              plot.margin = margin(5, 5, 5, 5)) +
        labs(title = sample_id)
      
      # Store the plot
      fingerprints[[sample_id]] <- p
    }
  }
  
  return(fingerprints)
}

# Function to calculate visual similarity between fingerprints
calculate_fingerprint_similarity <- function(fingerprints_data) {
  # Get taxa and samples info
  taxa <- rownames(fingerprints_data)
  samples <- colnames(fingerprints_data)
  n_samples <- length(samples)
  
  # Calculate similarity matrix
  similarity_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  rownames(similarity_matrix) <- samples
  colnames(similarity_matrix) <- samples
  
  # Calculate similarity for each pair of samples
  for (i in 1:(n_samples-1)) {
    for (j in (i+1):n_samples) {
      # Get profiles
      profile1 <- fingerprints_data[, i]
      profile2 <- fingerprints_data[, j]
      
      # Calculate cosine similarity
      similarity <- sum(profile1 * profile2) / 
        (sqrt(sum(profile1^2)) * sqrt(sum(profile2^2)))
      
      # Store in matrix
      similarity_matrix[i, j] <- similarity
      similarity_matrix[j, i] <- similarity
    }
  }
  
  # Set diagonal to 1 (self-similarity)
  diag(similarity_matrix) <- 1
  
  return(similarity_matrix)
}

# Function to arrange fingerprints in a grid with labels
plot_fingerprint_grid <- function(fingerprints, ncol = 4, group_by = NULL, metadata = NULL) {
  # If grouping is specified and metadata is provided
  if (!is.null(group_by) && !is.null(metadata)) {
    # Ensure group_by column exists in metadata
    if (!group_by %in% colnames(metadata)) {
      warning("Group column not found in metadata. Ignoring grouping.")
      group_by <- NULL
    } else {
      # Get group for each sample
      sample_groups <- metadata[[group_by]]
      names(sample_groups) <- metadata$sample_id
      
      # Only keep samples that have fingerprints
      sample_groups <- sample_groups[names(sample_groups) %in% names(fingerprints)]
      
      # Rearrange fingerprints by group
      grouped_fingerprints <- list()
      for (group in unique(sample_groups)) {
        group_samples <- names(sample_groups)[sample_groups == group]
        group_fingerprints <- fingerprints[names(fingerprints) %in% group_samples]
        
        # Create group title
        group_title <- ggplot() + 
          annotate("text", x = 0, y = 0, label = group, size = 6) +
          theme_void()
        
        # Add to list
        grouped_fingerprints[[paste0("title_", group)]] <- group_title
        grouped_fingerprints <- c(grouped_fingerprints, group_fingerprints)
      }
      
      # Replace original fingerprints with grouped ones
      fingerprints <- grouped_fingerprints
    }
  }
  
  # Arrange plots in a grid
  grid_plot <- gridExtra::grid.arrange(
    grobs = fingerprints,
    ncol = ncol
  )
  
  return(grid_plot)
}

# Generate microbiome fingerprints
fingerprints_circular <- generate_microbiome_fingerprint(
  fgt_exp, 
  method = "circular", 
  n_taxa = 20, 
  color_palette = "viridis"
)

fingerprints_heatmap <- generate_microbiome_fingerprint(
  fgt_exp, 
  method = "heatmap", 
  n_taxa = 25, 
  color_palette = "plasma"
)

# Get sample metadata for grouping
metadata <- fgt_exp$metadata

# Save circular fingerprints
grid_plot_circular <- plot_fingerprint_grid(
  fingerprints_circular, 
  ncol = 4,
  group_by = "group",
  metadata = metadata
)

# Save the circular fingerprint grid
ggsave("inst/example_plots/images/microbiome_fingerprints_circular.png", 
       grid_plot_circular, width = 12, height = 10, dpi = 300)

# Save heatmap fingerprints
grid_plot_heatmap <- plot_fingerprint_grid(
  fingerprints_heatmap, 
  ncol = 4,
  group_by = "group",
  metadata = metadata
)

# Save the heatmap fingerprint grid
ggsave("inst/example_plots/images/microbiome_fingerprints_heatmap.png", 
       grid_plot_heatmap, width = 12, height = 10, dpi = 300)

# Calculate fingerprint similarity
rel_abund <- transformAbundance(fgt_exp, type = "relative")
counts_data <- rel_abund$counts
top_taxa_idx <- order(rowMeans(counts_data), decreasing = TRUE)[1:20]
top_taxa_data <- counts_data[top_taxa_idx, ]

similarity_matrix <- calculate_fingerprint_similarity(top_taxa_data)

# Create similarity heatmap
similarity_df <- as.data.frame(as.table(similarity_matrix))
names(similarity_df) <- c("Sample1", "Sample2", "Similarity")

# Plot similarity heatmap
similarity_plot <- ggplot(similarity_df, aes(x = Sample1, y = Sample2, fill = Similarity)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Microbiome Fingerprint Similarity",
    subtitle = "Cosine similarity between sample profiles",
    x = NULL,
    y = NULL
  ) +
  coord_fixed()

# Save similarity plot
ggsave("inst/example_plots/images/fingerprint_similarity.png", 
       similarity_plot, width = 10, height = 9, dpi = 300)

# Print success message
cat("Microbiome fingerprint plots created successfully!\n")
cat("Plots saved to inst/example_plots/images/\n")