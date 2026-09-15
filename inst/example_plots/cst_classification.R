#' Community State Type (CST) Classification Example
#' 
#' This script demonstrates how to classify samples into community state types
#' and visualize CST distributions and transitions using the microFGT package.

# Create output directory if it doesn't exist
dir.create("inst/example_plots/images", recursive = TRUE, showWarnings = FALSE)

# Load required libraries
library(microFGT)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load example data or create mock data
# Try to load example data, otherwise create mock data
fgt_exp <- NULL
fgt_ts <- NULL

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

# Generate mock data directly
generate_mock_data <- function() {
  set.seed(42)
  
  # Create sample metadata
  n_samples <- 150
  metadata <- data.frame(
    sample_id = paste0("S", 1:n_samples),
    age_group = sample(c("18-25", "26-35", "36-45", "46+"), n_samples, replace = TRUE),
    race = sample(c("White", "Black", "Asian", "Hispanic", "Other"), n_samples, replace = TRUE),
    pH = round(rnorm(n_samples, 4.5, 0.8), 1),
    Nugent_score = sample(0:10, n_samples, replace = TRUE)
  )
  
  # Generate mock counts with realistic patterns
  n_taxa <- 50
  counts <- matrix(0, nrow = n_taxa, ncol = n_samples)
  rownames(counts) <- paste0("Taxon_", 1:n_taxa)
  colnames(counts) <- metadata$sample_id
  
  # Create taxonomy
  genera <- c("Lactobacillus", "Gardnerella", "Atopobium", "Prevotella", "Sneathia", 
              "Megasphaera", "Streptococcus", "Bifidobacterium", "Ureaplasma", "Mobiluncus")
  
  species <- c(
    "Lactobacillus_crispatus", "Lactobacillus_iners", "Lactobacillus_gasseri", "Lactobacillus_jensenii",
    "Gardnerella_vaginalis", "Atopobium_vaginae", "Prevotella_bivia", "Sneathia_amnii",
    "Megasphaera_type_1", "Streptococcus_agalactiae", "Bifidobacterium_breve",
    "Ureaplasma_urealyticum", "Mobiluncus_curtisii"
  )
  
  taxonomy <- data.frame(
    taxon_id = rownames(counts),
    Kingdom = "Bacteria",
    Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Fusobacteria"), n_taxa, replace = TRUE),
    Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria", "Fusobacteriia"), n_taxa, replace = TRUE),
    Order = sample(c("Lactobacillales", "Clostridiales", "Bacteroidales", "Bifidobacteriales", "Fusobacteriales"), n_taxa, replace = TRUE),
    Family = sample(c("Lactobacillaceae", "Bifidobacteriaceae", "Prevotellaceae", "Streptococcaceae"), n_taxa, replace = TRUE),
    Genus = sample(genera, n_taxa, replace = TRUE),
    Species = sample(species, n_taxa, replace = TRUE)
  )
  
  # Create patterns for 5 CSTs
  # CST I: L. crispatus dominated
  # CST II: L. gasseri dominated  
  # CST III: L. iners dominated
  # CST IV: Diverse anaerobes
  # CST V: L. jensenii dominated
  
  # Assign random CSTs to samples
  cst_assignments <- sample(c("CST_I", "CST_II", "CST_III", "CST_IV", "CST_V"), 
                           n_samples, replace = TRUE, 
                           prob = c(0.3, 0.1, 0.25, 0.25, 0.1))
  
  # Generate counts based on CST
  for (i in 1:n_samples) {
    cst <- cst_assignments[i]
    
    # Base profile with low counts for all taxa
    counts[, i] <- sample(0:10, n_taxa, replace = TRUE)
    
    # CST-specific dominant taxa
    if (cst == "CST_I") {
      # Find Lactobacillus crispatus indices
      l_crispatus_idx <- which(taxonomy$Species == "Lactobacillus_crispatus")
      if (length(l_crispatus_idx) > 0) {
        # Make L. crispatus dominant
        counts[l_crispatus_idx[1], i] <- sample(5000:20000, 1)
      }
    } else if (cst == "CST_II") {
      # Find Lactobacillus gasseri indices
      l_gasseri_idx <- which(taxonomy$Species == "Lactobacillus_gasseri")
      if (length(l_gasseri_idx) > 0) {
        # Make L. gasseri dominant
        counts[l_gasseri_idx[1], i] <- sample(5000:15000, 1)
      }
    } else if (cst == "CST_III") {
      # Find Lactobacillus iners indices  
      l_iners_idx <- which(taxonomy$Species == "Lactobacillus_iners")
      if (length(l_iners_idx) > 0) {
        # Make L. iners dominant
        counts[l_iners_idx[1], i] <- sample(5000:15000, 1)
      }
    } else if (cst == "CST_IV") {
      # Diverse anaerobes, no single dominant species
      anaerobe_species <- c("Gardnerella_vaginalis", "Atopobium_vaginae", 
                           "Prevotella_bivia", "Sneathia_amnii", "Megasphaera_type_1")
      
      for (species_name in anaerobe_species) {
        species_idx <- which(taxonomy$Species == species_name)
        if (length(species_idx) > 0) {
          counts[species_idx[1], i] <- sample(500:5000, 1)
        }
      }
    } else if (cst == "CST_V") {
      # Find Lactobacillus jensenii indices
      l_jensenii_idx <- which(taxonomy$Species == "Lactobacillus_jensenii")
      if (length(l_jensenii_idx) > 0) {
        # Make L. jensenii dominant
        counts[l_jensenii_idx[1], i] <- sample(5000:15000, 1)
      }
    }
  }
  
  # Create a time series dataset
  n_subjects <- 10
  n_timepoints <- 5
  n_ts_samples <- n_subjects * n_timepoints
  
  ts_metadata <- data.frame(
    sample_id = paste0("TS_", 1:n_ts_samples),
    subject_id = rep(paste0("Subject_", 1:n_subjects), each = n_timepoints),
    timepoint = rep(1:n_timepoints, n_subjects),
    visit_week = rep(c(0, 4, 12, 24, 48), n_subjects),
    age_group = sample(c("18-25", "26-35", "36-45", "46+"), n_ts_samples, replace = TRUE)
  )
  
  # Generate time series counts with CST transitions
  ts_counts <- matrix(0, nrow = n_taxa, ncol = n_ts_samples)
  rownames(ts_counts) <- paste0("Taxon_", 1:n_taxa)
  colnames(ts_counts) <- ts_metadata$sample_id
  
  # Assign initial CSTs to subjects
  subject_cst <- sample(c("CST_I", "CST_II", "CST_III", "CST_IV", "CST_V"),
                        n_subjects, replace = TRUE)
  
  # Generate time series with potential transitions
  for (s in 1:n_subjects) {
    for (t in 1:n_timepoints) {
      idx <- (s-1)*n_timepoints + t
      
      # Possible CST transition
      if (t > 1 && runif(1) < 0.3) {  # 30% chance of transition
        subject_cst[s] <- sample(c("CST_I", "CST_II", "CST_III", "CST_IV", "CST_V"), 1)
      }
      
      # Base profile with low counts for all taxa
      ts_counts[, idx] <- sample(0:10, n_taxa, replace = TRUE)
      
      # Apply CST-specific pattern
      cst <- subject_cst[s]
      
      if (cst == "CST_I") {
        l_crispatus_idx <- which(taxonomy$Species == "Lactobacillus_crispatus")
        if (length(l_crispatus_idx) > 0) {
          ts_counts[l_crispatus_idx[1], idx] <- sample(5000:20000, 1)
        }
      } else if (cst == "CST_II") {
        l_gasseri_idx <- which(taxonomy$Species == "Lactobacillus_gasseri")
        if (length(l_gasseri_idx) > 0) {
          ts_counts[l_gasseri_idx[1], idx] <- sample(5000:15000, 1)
        }
      } else if (cst == "CST_III") {
        l_iners_idx <- which(taxonomy$Species == "Lactobacillus_iners")
        if (length(l_iners_idx) > 0) {
          ts_counts[l_iners_idx[1], idx] <- sample(5000:15000, 1)
        }
      } else if (cst == "CST_IV") {
        anaerobe_species <- c("Gardnerella_vaginalis", "Atopobium_vaginae",
                             "Prevotella_bivia", "Sneathia_amnii", "Megasphaera_type_1")
        
        for (species_name in anaerobe_species) {
          species_idx <- which(taxonomy$Species == species_name)
          if (length(species_idx) > 0) {
            ts_counts[species_idx[1], idx] <- sample(500:5000, 1)
          }
        }
      } else if (cst == "CST_V") {
        l_jensenii_idx <- which(taxonomy$Species == "Lactobacillus_jensenii")
        if (length(l_jensenii_idx) > 0) {
          ts_counts[l_jensenii_idx[1], idx] <- sample(5000:15000, 1)
        }
      }
    }
  }
  
  return(list(
    fgt_exp = create_mock_fgt_experiment(
      counts = counts,
      sample_metadata = metadata,
      taxonomy = taxonomy,
      experiment_name = "Mock FGT Dataset with CSTs"
    ),
    fgt_ts = create_mock_fgt_experiment(
      counts = ts_counts,
      sample_metadata = ts_metadata,
      taxonomy = taxonomy,
      experiment_name = "Mock Time Series Dataset with CST Transitions"
    )
  ))
}

# Try to load example data, if fails use mock data
# Generate the data
mock_data <- generate_mock_data()
fgt_exp <- mock_data$fgt_exp
fgt_ts <- mock_data$fgt_ts
    race = sample(c("White", "Black", "Asian", "Hispanic", "Other"), n_samples, replace = TRUE),
    pH = round(rnorm(n_samples, 4.5, 0.8), 1),
    Nugent_score = sample(0:10, n_samples, replace = TRUE)
  )
  
  # Generate mock counts with realistic patterns
  n_taxa <- 50
  counts <- matrix(0, nrow = n_taxa, ncol = n_samples)
  rownames(counts) <- paste0("Taxon_", 1:n_taxa)
  colnames(counts) <- metadata$sample_id
  
  # Create taxonomy
  genera <- c("Lactobacillus", "Gardnerella", "Atopobium", "Prevotella", "Sneathia", 
              "Megasphaera", "Streptococcus", "Bifidobacterium", "Ureaplasma", "Mobiluncus")
  
  species <- c(
    "Lactobacillus_crispatus", "Lactobacillus_iners", "Lactobacillus_gasseri", "Lactobacillus_jensenii",
    "Gardnerella_vaginalis", "Atopobium_vaginae", "Prevotella_bivia", "Sneathia_amnii",
    "Megasphaera_type_1", "Streptococcus_agalactiae", "Bifidobacterium_breve",
    "Ureaplasma_urealyticum", "Mobiluncus_curtisii"
  )
  
  taxonomy <- data.frame(
    taxon_id = rownames(counts),
    Kingdom = "Bacteria",
    Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Fusobacteria"), n_taxa, replace = TRUE),
    Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria", "Fusobacteriia"), n_taxa, replace = TRUE),
    Order = sample(c("Lactobacillales", "Clostridiales", "Bacteroidales", "Bifidobacteriales", "Fusobacteriales"), n_taxa, replace = TRUE),
    Family = sample(c("Lactobacillaceae", "Bifidobacteriaceae", "Prevotellaceae", "Streptococcaceae"), n_taxa, replace = TRUE),
    Genus = sample(genera, n_taxa, replace = TRUE),
    Species = sample(species, n_taxa, replace = TRUE)
  )
  
  # Create patterns for 5 CSTs
  # CST I: L. crispatus dominated
  # CST II: L. gasseri dominated  
  # CST III: L. iners dominated
  # CST IV: Diverse anaerobes
  # CST V: L. jensenii dominated
  
  # Assign random CSTs to samples
  cst_assignments <- sample(c("CST_I", "CST_II", "CST_III", "CST_IV", "CST_V"), 
                           n_samples, replace = TRUE, 
                           prob = c(0.3, 0.1, 0.25, 0.25, 0.1))
  
  # Generate counts based on CST
  for (i in 1:n_samples) {
    cst <- cst_assignments[i]
    
    # Base profile with low counts for all taxa
    counts[, i] <- sample(0:10, n_taxa, replace = TRUE)
    
    # CST-specific dominant taxa
    if (cst == "CST_I") {
      # Find Lactobacillus crispatus indices
      l_crispatus_idx <- which(taxonomy$Species == "Lactobacillus_crispatus")
      if (length(l_crispatus_idx) > 0) {
        # Make L. crispatus dominant
        counts[l_crispatus_idx[1], i] <- sample(5000:20000, 1)
      }
    } else if (cst == "CST_II") {
      # Find Lactobacillus gasseri indices
      l_gasseri_idx <- which(taxonomy$Species == "Lactobacillus_gasseri")
      if (length(l_gasseri_idx) > 0) {
        # Make L. gasseri dominant
        counts[l_gasseri_idx[1], i] <- sample(5000:15000, 1)
      }
    } else if (cst == "CST_III") {
      # Find Lactobacillus iners indices  
      l_iners_idx <- which(taxonomy$Species == "Lactobacillus_iners")
      if (length(l_iners_idx) > 0) {
        # Make L. iners dominant
        counts[l_iners_idx[1], i] <- sample(5000:15000, 1)
      }
    } else if (cst == "CST_IV") {
      # Diverse anaerobes, no single dominant species
      anaerobe_species <- c("Gardnerella_vaginalis", "Atopobium_vaginae", 
                           "Prevotella_bivia", "Sneathia_amnii", "Megasphaera_type_1")
      
      for (species_name in anaerobe_species) {
        species_idx <- which(taxonomy$Species == species_name)
        if (length(species_idx) > 0) {
          counts[species_idx[1], i] <- sample(500:5000, 1)
        }
      }
    } else if (cst == "CST_V") {
      # Find Lactobacillus jensenii indices
      l_jensenii_idx <- which(taxonomy$Species == "Lactobacillus_jensenii")
      if (length(l_jensenii_idx) > 0) {
        # Make L. jensenii dominant
        counts[l_jensenii_idx[1], i] <- sample(5000:15000, 1)
      }
    }
  }
  
  # Create a simple struct to mimic the FGTExperiment object
  fgt_exp <- list(
    counts = counts,
    metadata = metadata,
    taxonomy = taxonomy,
    experiment_name = "Mock FGT Dataset with CSTs"
  )

  # Define a simple class
  class(fgt_exp) <- "mock_fgt_experiment"

  # Create time series data
  n_subjects <- 10
  n_timepoints <- 5
  n_ts_samples <- n_subjects * n_timepoints

  ts_metadata <- data.frame(
    sample_id = paste0("TS_", 1:n_ts_samples),
    subject_id = rep(paste0("Subject_", 1:n_subjects), each = n_timepoints),
    timepoint = rep(1:n_timepoints, n_subjects),
    visit_week = rep(c(0, 4, 12, 24, 48), n_subjects),
    age_group = sample(c("18-25", "26-35", "36-45", "46+"), n_ts_samples, replace = TRUE)
  )

  # Generate time series counts with CST transitions
  ts_counts <- matrix(0, nrow = n_taxa, ncol = n_ts_samples)
  rownames(ts_counts) <- paste0("Taxon_", 1:n_taxa)
  colnames(ts_counts) <- ts_metadata$sample_id

  # Assign initial CSTs to subjects
  subject_cst <- sample(c("CST_I", "CST_II", "CST_III", "CST_IV", "CST_V"),
                        n_subjects, replace = TRUE)

  # Generate time series with potential transitions
  for (s in 1:n_subjects) {
    for (t in 1:n_timepoints) {
      idx <- (s-1)*n_timepoints + t

      # Possible CST transition
      if (t > 1 && runif(1) < 0.3) {  # 30% chance of transition
        subject_cst[s] <- sample(c("CST_I", "CST_II", "CST_III", "CST_IV", "CST_V"), 1)
      }

      # Base profile with low counts for all taxa
      ts_counts[, idx] <- sample(0:10, n_taxa, replace = TRUE)

      # Apply CST-specific pattern
      cst <- subject_cst[s]

      if (cst == "CST_I") {
        l_crispatus_idx <- which(taxonomy$Species == "Lactobacillus_crispatus")
        if (length(l_crispatus_idx) > 0) {
          ts_counts[l_crispatus_idx[1], idx] <- sample(5000:20000, 1)
        }
      } else if (cst == "CST_II") {
        l_gasseri_idx <- which(taxonomy$Species == "Lactobacillus_gasseri")
        if (length(l_gasseri_idx) > 0) {
          ts_counts[l_gasseri_idx[1], idx] <- sample(5000:15000, 1)
        }
      } else if (cst == "CST_III") {
        l_iners_idx <- which(taxonomy$Species == "Lactobacillus_iners")
        if (length(l_iners_idx) > 0) {
          ts_counts[l_iners_idx[1], idx] <- sample(5000:15000, 1)
        }
      } else if (cst == "CST_IV") {
        anaerobe_species <- c("Gardnerella_vaginalis", "Atopobium_vaginae",
                             "Prevotella_bivia", "Sneathia_amnii", "Megasphaera_type_1")

        for (species_name in anaerobe_species) {
          species_idx <- which(taxonomy$Species == species_name)
          if (length(species_idx) > 0) {
            ts_counts[species_idx[1], idx] <- sample(500:5000, 1)
          }
        }
      } else if (cst == "CST_V") {
        l_jensenii_idx <- which(taxonomy$Species == "Lactobacillus_jensenii")
        if (length(l_jensenii_idx) > 0) {
          ts_counts[l_jensenii_idx[1], idx] <- sample(5000:15000, 1)
        }
      }
    }
  }

  # Create a simple struct to mimic the FGTExperiment object for time series
  fgt_ts <- list(
    counts = ts_counts,
    metadata = ts_metadata,
    taxonomy = taxonomy,
    experiment_name = "Mock Time Series Dataset with CST Transitions"
  )

  # Define a simple class
  class(fgt_ts) <- "mock_fgt_experiment"
})

# Function to classify CSTs
classify_cst <- function(fgt_exp) {
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
        stringsAsFactors = FALSE
      )
      taxa_info <- data.frame(
        taxon_id = rownames(counts_data),
        Species = rownames(counts_data),
        stringsAsFactors = FALSE
      )

      # Convert to relative abundance
      counts_total <- colSums(counts_data)
      rel_counts <- t(t(counts_data) / counts_total)
    })
  }

  # Check for key Lactobacillus species
  l_crispatus_idx <- grep("Lactobacillus.*crispatus", taxa_info$Species)
  l_gasseri_idx <- grep("Lactobacillus.*gasseri", taxa_info$Species)
  l_iners_idx <- grep("Lactobacillus.*iners", taxa_info$Species)
  l_jensenii_idx <- grep("Lactobacillus.*jensenii", taxa_info$Species)

  # CST classifications
  cst <- rep("CST_IV", ncol(rel_counts))  # Default to CST IV (diverse)

  # Apply CST classification rules
  for (i in 1:ncol(rel_counts)) {
    # Calculate abundance of each species
    if (length(l_crispatus_idx) > 0) {
      l_crispatus_abund <- sum(rel_counts[l_crispatus_idx, i])
    } else {
      l_crispatus_abund <- 0
    }

    if (length(l_gasseri_idx) > 0) {
      l_gasseri_abund <- sum(rel_counts[l_gasseri_idx, i])
    } else {
      l_gasseri_abund <- 0
    }

    if (length(l_iners_idx) > 0) {
      l_iners_abund <- sum(rel_counts[l_iners_idx, i])
    } else {
      l_iners_abund <- 0
    }

    if (length(l_jensenii_idx) > 0) {
      l_jensenii_abund <- sum(rel_counts[l_jensenii_idx, i])
    } else {
      l_jensenii_abund <- 0
    }

    # Apply classification rules
    # CST I: L. crispatus dominated (>50%)
    if (l_crispatus_abund > 0.5) {
      cst[i] <- "CST_I"
    }
    # CST II: L. gasseri dominated (>50%)
    else if (l_gasseri_abund > 0.5) {
      cst[i] <- "CST_II"
    }
    # CST III: L. iners dominated (>50%)
    else if (l_iners_abund > 0.5) {
      cst[i] <- "CST_III"
    }
    # CST V: L. jensenii dominated (>50%)
    else if (l_jensenii_abund > 0.5) {
      cst[i] <- "CST_V"
    }
    # Else remains CST IV (diverse)
  }

  # Add CST to metadata
  metadata$CST <- cst

  # Return results based on class
  if (inherits(fgt_exp, "mock_fgt_experiment")) {
    # Update the mock object
    fgt_exp$metadata <- metadata
    return(fgt_exp)
  } else {
    # Try to update FGTExperiment
    tryCatch({
      # Update the experiment object
      colData(experimentData(fgt_exp)) <- DataFrame(metadata)
      return(fgt_exp)
    }, error = function(e) {
      # Return the updated metadata
      return(list(data = rel_counts, metadata = metadata))
    })
  }
}

# Function to plot CST distribution
plot_cst_distribution <- function(fgt_exp, by = NULL) {
  # Get metadata with CST based on object type
  if (inherits(fgt_exp, "mock_fgt_experiment")) {
    metadata <- fgt_exp$metadata
  } else if (is.list(fgt_exp) && "metadata" %in% names(fgt_exp)) {
    metadata <- fgt_exp$metadata
  } else {
    # Try to get metadata from FGTExperiment
    tryCatch({
      metadata <- as.data.frame(colData(experimentData(fgt_exp)))
    }, error = function(e) {
      stop("Could not extract metadata from the provided object")
    })
  }

  if (!("CST" %in% colnames(metadata))) {
    stop("CST classification not found. Run classify_cst() first.")
  }

  # Create plot data
  if (!is.null(by) && by %in% colnames(metadata)) {
    # Group by specified variable
    plot_data <- as.data.frame(table(metadata[[by]], metadata$CST))
    colnames(plot_data) <- c("Group", "CST", "Count")

    # Calculate percentage within each group
    plot_data <- plot_data %>%
      group_by(Group) %>%
      mutate(Percentage = Count / sum(Count) * 100)

    # Create plot
    p <- ggplot(plot_data, aes(x = Group, y = Percentage, fill = CST)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_brewer(palette = "Set1") +
      labs(
        title = "Community State Type (CST) Distribution",
        subtitle = paste("Grouped by", by),
        x = by,
        y = "Percentage (%)",
        fill = "CST"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      )
  } else {
    # Overall distribution
    plot_data <- as.data.frame(table(metadata$CST))
    colnames(plot_data) <- c("CST", "Count")
    plot_data$Percentage <- plot_data$Count / sum(plot_data$Count) * 100

    # Create plot
    p <- ggplot(plot_data, aes(x = CST, y = Percentage, fill = CST)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.5) +
      scale_fill_brewer(palette = "Set1") +
      labs(
        title = "Community State Type (CST) Distribution",
        x = "Community State Type",
        y = "Percentage (%)"
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "none"
      )
  }

  return(p)
}

# Function to calculate CST transitions
calculate_cst_transitions <- function(fgt_ts_with_cst) {
  # Get metadata with CST based on object type
  if (inherits(fgt_ts_with_cst, "mock_fgt_experiment")) {
    metadata <- fgt_ts_with_cst$metadata
  } else if (is.list(fgt_ts_with_cst) && "metadata" %in% names(fgt_ts_with_cst)) {
    metadata <- fgt_ts_with_cst$metadata
  } else {
    # Try to get metadata from FGTExperiment
    tryCatch({
      metadata <- as.data.frame(colData(experimentData(fgt_ts_with_cst)))
    }, error = function(e) {
      stop("Could not extract metadata from the provided object")
    })
  }

  if (!("CST" %in% colnames(metadata))) {
    stop("CST classification not found. Run classify_cst() first.")
  }

  if (!("subject_id" %in% colnames(metadata)) || !("timepoint" %in% colnames(metadata))) {
    stop("Required columns 'subject_id' and 'timepoint' not found in metadata.")
  }

  # Calculate transitions
  transition_counts <- matrix(0, nrow = 5, ncol = 5)
  rownames(transition_counts) <- paste0("CST_", c("I", "II", "III", "IV", "V"))
  colnames(transition_counts) <- paste0("CST_", c("I", "II", "III", "IV", "V"))

  # Get unique subjects
  subjects <- unique(metadata$subject_id)

  for (subject in subjects) {
    # Get data for this subject
    subject_data <- metadata[metadata$subject_id == subject, ]

    # Sort by timepoint
    subject_data <- subject_data[order(subject_data$timepoint), ]

    # Calculate transitions
    if (nrow(subject_data) >= 2) {
      for (i in 1:(nrow(subject_data)-1)) {
        from_cst <- subject_data$CST[i]
        to_cst <- subject_data$CST[i+1]

        # Skip if either CST is NA
        if (is.na(from_cst) || is.na(to_cst)) next

        # Count transition
        transition_counts[from_cst, to_cst] <- transition_counts[from_cst, to_cst] + 1
      }
    }
  }

  return(transition_counts)
}

# Function to plot CST transitions
plot_cst_transitions <- function(transition_matrix) {
  # Convert matrix to data frame for plotting
  transition_df <- as.data.frame(as.table(transition_matrix))
  colnames(transition_df) <- c("From", "To", "Count")
  
  # Calculate percentage
  total_transitions <- sum(transition_df$Count)
  transition_df$Percentage <- transition_df$Count / total_transitions * 100
  
  # Add self-transition indicator
  transition_df$Type <- ifelse(as.character(transition_df$From) == as.character(transition_df$To), 
                               "Self", "Transition")
  
  # Create plot
  p <- ggplot(transition_df, aes(x = From, y = To, fill = Count)) +
    geom_tile() +
    geom_text(aes(label = ifelse(Count > 0, sprintf("%d\n(%.1f%%)", Count, Percentage), "")), 
              color = "black", size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      title = "Community State Type (CST) Transitions",
      subtitle = "Number and percentage of transitions between CSTs",
      x = "From CST",
      y = "To CST",
      fill = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}

# Run CST classification
fgt_with_cst <- classify_cst(fgt_exp)
fgt_ts_with_cst <- classify_cst(fgt_ts)

# Create and save CST distribution plot
cst_plot <- plot_cst_distribution(fgt_with_cst)
ggsave("inst/example_plots/images/cst_distribution.png", cst_plot, width = 8, height = 6, dpi = 300)

# Create and save CST distribution by group
cst_by_group_plot <- plot_cst_distribution(fgt_with_cst, by = "age_group")
ggsave("inst/example_plots/images/cst_by_age_group.png", cst_by_group_plot, width = 10, height = 6, dpi = 300)

# Calculate and plot CST transitions
transition_matrix <- calculate_cst_transitions(fgt_ts_with_cst)
transition_plot <- plot_cst_transitions(transition_matrix)
ggsave("inst/example_plots/images/cst_transitions.png", transition_plot, width = 8, height = 7, dpi = 300)

# Print success message
cat("CST classification plots created successfully!\n")
cat("Plots saved to inst/example_plots/images/\n")