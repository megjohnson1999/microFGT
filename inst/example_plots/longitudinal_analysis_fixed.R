#' Longitudinal Analysis Example
#' 
#' This script demonstrates how to perform and visualize longitudinal analysis
#' of microbiome stability using the microFGT package.

# Create output directory if it doesn't exist
dir.create("inst/example_plots/images", recursive = TRUE, showWarnings = FALSE)

# Load required libraries
library(microFGT)
library(ggplot2)
library(dplyr)
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

# Generate mock time series data
generate_mock_data <- function() {
  set.seed(42)
  n_samples <- 30
  n_timepoints <- 5
  n_subjects <- 6

  # Create sample metadata
  metadata <- data.frame(
    sample_id = paste0("S", 1:(n_samples)),
    subject_id = rep(paste0("Subject_", 1:n_subjects), each = n_timepoints),
    timepoint = rep(1:n_timepoints, n_subjects),
    visit_week = rep(c(0, 4, 12, 24, 48), n_subjects),
    age = sample(21:45, n_samples, replace = TRUE),
    treatment = rep(rep(c("A", "B"), each = n_timepoints), n_subjects/2)
  )

  # Generate counts
  n_taxa <- 50

  # Generate mock counts with time continuity
  # Create base subject profiles
  subject_profiles <- matrix(
    rexp(n_subjects * n_taxa, rate = 1),
    nrow = n_subjects,
    ncol = n_taxa
  )

  # Add time variation
  counts <- matrix(0, nrow = n_taxa, ncol = n_samples)
  for (s in 1:n_subjects) {
    for (t in 1:n_timepoints) {
      idx <- (s-1)*n_timepoints + t
      # Base profile plus time-dependent noise
      sample_profile <- subject_profiles[s,] +
        rnorm(n_taxa, mean = 0, sd = 0.1 * t)
      # Ensure non-negative
      sample_profile[sample_profile < 0] <- 0
      # Scale to counts
      counts[, idx] <- round(sample_profile * 1000)
    }
  }

  # Set dimnames
  taxa_names <- paste0("Taxon_", 1:n_taxa)
  rownames(counts) <- taxa_names
  colnames(counts) <- metadata$sample_id

  # Create taxonomy
  taxonomy <- data.frame(
    taxon_id = taxa_names,
    Kingdom = "Bacteria",
    Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Proteobacteria"),
                   n_taxa, replace = TRUE),
    Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria", "Gammaproteobacteria"),
                  n_taxa, replace = TRUE),
    Order = sample(c("Lactobacillales", "Clostridiales", "Bacteroidales", "Bifidobacteriales", "Enterobacteriales"),
                  n_taxa, replace = TRUE),
    Family = sample(c("Lactobacillaceae", "Clostridiaceae", "Bacteroidaceae", "Bifidobacteriaceae", "Enterobacteriaceae"),
                   n_taxa, replace = TRUE),
    Genus = sample(c("Lactobacillus", "Clostridium", "Bacteroides", "Bifidobacterium", "Escherichia"),
                  n_taxa, replace = TRUE),
    Species = paste0(sample(c("Lactobacillus", "Clostridium", "Bacteroides", "Bifidobacterium", "Escherichia"),
                           n_taxa, replace = TRUE), " ", sample(letters, n_taxa, replace = TRUE))
  )

  # Create a FGT experiment object
  return(create_mock_fgt_experiment(
    counts = counts,
    sample_metadata = metadata,
    taxonomy = taxonomy,
    experiment_name = "Mock Time Series Dataset"
  ))
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

# Generate the time series data
fgt_ts <- generate_mock_data()

# Function to plot taxa stability over time
plot_taxa_stability <- function(fgt_exp, top_taxa = 5, group_by = "subject_id") {
  # Handle different class types
  if (inherits(fgt_exp, "mock_fgt_experiment")) {
    # Simple handling for our mock class
    counts_data <- fgt_exp$counts
    metadata <- fgt_exp$metadata
    taxa_info <- fgt_exp$taxonomy

    # Convert to relative abundance
    counts_total <- colSums(counts_data)
    rel_counts <- t(t(counts_data) / counts_total)
  } else {
    # Try to handle FGTExperiment objects
    tryCatch({
      # Get relative abundance data
      rel_abund <- transformAbundance(fgt_exp, type = "relative")
      counts_data <- assays(experimentData(rel_abund))$counts
      metadata <- as.data.frame(colData(experimentData(rel_abund)))
      taxa_info <- as.data.frame(rowData(experimentData(rel_abund)))
      rel_counts <- counts_data
    }, error = function(e) {
      # Fallback for direct matrices
      counts_data <- fgt_exp
      metadata <- data.frame(
        sample_id = colnames(counts_data),
        subject_id = rep(paste0("Subject_", seq_len(ncol(counts_data) %/% 5)), each = 5),
        timepoint = rep(1:5, ncol(counts_data) %/% 5),
        visit_week = rep(c(0, 4, 12, 24, 48), ncol(counts_data) %/% 5)
      )
      taxa_info <- data.frame(
        taxon_id = rownames(counts_data),
        Genus = rep("Unknown", nrow(counts_data)),
        Species = rownames(counts_data)
      )

      # Convert to relative abundance
      counts_total <- colSums(counts_data)
      rel_counts <- t(t(counts_data) / counts_total)
    })
  }

  # Identify top taxa across all samples
  taxa_means <- rowMeans(rel_counts)
  top_taxa_idx <- order(taxa_means, decreasing = TRUE)[1:min(top_taxa, length(taxa_means))]

  # Create a data frame for plotting
  plot_data <- data.frame(rel_counts[top_taxa_idx,])
  colnames(plot_data) <- metadata$sample_id
  plot_data$taxon_id <- rownames(rel_counts)[top_taxa_idx]

  # Add taxonomy information
  plot_data$Genus <- taxa_info$Genus[match(plot_data$taxon_id, taxa_info$taxon_id)]
  plot_data$Species <- taxa_info$Species[match(plot_data$taxon_id, taxa_info$taxon_id)]
  plot_data$taxon_name <- ifelse(is.na(plot_data$Species),
                                plot_data$Genus,
                                plot_data$Species)

  # Reshape for plotting
  plot_data_long <- tidyr::pivot_longer(
    plot_data,
    cols = metadata$sample_id,
    names_to = "sample_id",
    values_to = "abundance"
  )

  # Add metadata
  plot_data_long <- merge(
    plot_data_long,
    metadata,
    by = "sample_id"
  )

  # Sort by timepoint
  plot_data_long <- plot_data_long[order(plot_data_long[[group_by]], plot_data_long$timepoint),]

  # Create the plot
  p <- ggplot(plot_data_long, aes(x = visit_week, y = abundance,
                                 color = taxon_name, group = taxon_name)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~get(group_by)) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Stability of Top Taxa Over Time",
      subtitle = paste("Top", top_taxa, "most abundant taxa"),
      x = "Visit Week",
      y = "Relative Abundance",
      color = "Taxon"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "lightblue", color = "gray"),
      strip.text = element_text(face = "bold")
    )

  return(p)
}

# Function to calculate microbiome volatility index
calculate_volatility_index <- function(fgt_exp, method = "jensen_shannon", by_subject = TRUE) {
  # Handle different class types
  if (inherits(fgt_exp, "mock_fgt_experiment")) {
    # Simple handling for our mock class
    counts_data <- fgt_exp$counts
    metadata <- fgt_exp$metadata

    # Convert to relative abundance
    counts_total <- colSums(counts_data)
    rel_counts <- t(t(counts_data) / counts_total)
  } else {
    # Try to handle FGTExperiment objects
    tryCatch({
      # Get relative abundance data
      rel_abund <- transformAbundance(fgt_exp, type = "relative")
      counts_data <- assays(experimentData(rel_abund))$counts
      metadata <- as.data.frame(colData(experimentData(rel_abund)))
      rel_counts <- counts_data
    }, error = function(e) {
      # Fallback for direct matrices
      counts_data <- fgt_exp
      metadata <- data.frame(
        sample_id = colnames(counts_data),
        subject_id = rep(paste0("Subject_", seq_len(ncol(counts_data) %/% 5)), each = 5),
        timepoint = rep(1:5, ncol(counts_data) %/% 5),
        visit_week = rep(c(0, 4, 12, 24, 48), ncol(counts_data) %/% 5)
      )

      # Convert to relative abundance
      counts_total <- colSums(counts_data)
      rel_counts <- t(t(counts_data) / counts_total)
    })
  }

  # Calculate distances between consecutive timepoints
  subjects <- unique(metadata$subject_id)
  result <- data.frame()

  for (subject in subjects) {
    subject_samples <- metadata$sample_id[metadata$subject_id == subject]
    subject_timepoints <- metadata$timepoint[metadata$subject_id == subject]

    # Sort by timepoint
    idx <- order(subject_timepoints)
    subject_samples <- subject_samples[idx]
    subject_timepoints <- subject_timepoints[idx]

    if (length(subject_samples) < 2) next

    # Calculate distances between consecutive samples
    for (i in 1:(length(subject_samples)-1)) {
      s1 <- subject_samples[i]
      s2 <- subject_samples[i+1]
      t1 <- subject_timepoints[i]
      t2 <- subject_timepoints[i+1]

      # Calculate distance based on method
      if (method == "jensen_shannon") {
        # Get profiles
        p1 <- rel_counts[, s1]
        p2 <- rel_counts[, s2]

        # Calculate JS distance (simplified)
        m <- (p1 + p2)/2
        p1_norm <- p1 / sum(p1)
        p2_norm <- p2 / sum(p2)
        m_norm <- m / sum(m)

        # Handle zeros
        p1_norm[p1_norm == 0] <- 1e-10
        p2_norm[p2_norm == 0] <- 1e-10
        m_norm[m_norm == 0] <- 1e-10

        d1 <- sum(p1_norm * log(p1_norm / m_norm), na.rm = TRUE)
        d2 <- sum(p2_norm * log(p2_norm / m_norm), na.rm = TRUE)
        dist <- sqrt((d1 + d2)/2)
      } else {
        # Euclidean distance
        p1 <- rel_counts[, s1]
        p2 <- rel_counts[, s2]
        dist <- sqrt(sum((p1 - p2)^2, na.rm = TRUE))
      }

      # Store result
      result <- rbind(result, data.frame(
        subject_id = subject,
        timepoint1 = t1,
        timepoint2 = t2,
        distance = dist
      ))
    }
  }

  # Calculate volatility index
  if (by_subject) {
    volatility <- aggregate(distance ~ subject_id, data = result, FUN = mean)
    colnames(volatility)[2] <- "volatility_index"
  } else {
    volatility <- data.frame(
      overall_volatility = mean(result$distance, na.rm = TRUE)
    )
  }

  return(volatility)
}

# Plot taxa stability
stability_plot <- plot_taxa_stability(fgt_ts, top_taxa = 5, group_by = "subject_id")

# Save the plot
ggsave("inst/example_plots/images/taxa_stability.png", stability_plot, width = 12, height = 8, dpi = 300)

# Calculate volatility index
volatility <- calculate_volatility_index(fgt_ts, method = "jensen_shannon", by_subject = TRUE)

# Plot volatility index
volatility_plot <- ggplot(volatility, aes(x = reorder(subject_id, volatility_index), y = volatility_index)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(volatility_index, 3)), vjust = -0.5) +
  labs(
    title = "Microbiome Volatility Index by Subject",
    subtitle = "Higher values indicate more temporal instability",
    x = "Subject ID",
    y = "Volatility Index (Jensen-Shannon Distance)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

# Save the plot
ggsave("inst/example_plots/images/volatility_index.png", volatility_plot, width = 10, height = 6, dpi = 300)

# Print success message
cat("Longitudinal analysis plots created successfully!\n")
cat("Plots saved to inst/example_plots/images/\n")