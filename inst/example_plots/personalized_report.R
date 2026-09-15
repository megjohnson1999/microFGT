#' Personalized Microbiome Report Example
#'
#' This script demonstrates how to generate personalized microbiome "report cards"
#' for clinical reporting, with visual summaries of key microbiome metrics.

# Create output directory if it doesn't exist
dir.create("inst/example_plots/images", recursive = TRUE, showWarnings = FALSE)
dir.create("inst/example_plots/reports", recursive = TRUE, showWarnings = FALSE)

# Load required libraries
library(microFGT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(gridExtra)
library(grid)

# Load example data or create mock data
tryCatch({
  fgt_exp <- load_example_data()
}, error = function(e) {
  # Generate mock data with clinical metadata
  set.seed(42)
  
  # Create sample metadata
  n_samples <- 100
  n_patients <- 20
  
  # Create sample metadata with clinical variables
  patient_ids <- paste0("Patient_", 1:n_patients)
  metadata <- data.frame(
    sample_id = paste0("S", 1:n_samples),
    patient_id = sample(patient_ids, n_samples, replace = TRUE),
    visit_type = sample(c("Baseline", "Follow-up"), n_samples, replace = TRUE),
    visit_date = sample(seq(as.Date("2022-01-01"), as.Date("2023-01-01"), by = "day"), n_samples),
    pH = round(rnorm(n_samples, mean = 4.5, sd = 0.8), 1),
    Nugent_score = sample(0:10, n_samples, replace = TRUE),
    Amsel_criteria = sample(0:4, n_samples, replace = TRUE),
    symptoms = sample(c("None", "Mild", "Moderate", "Severe"), n_samples, replace = TRUE),
    age = sample(18:45, n_samples, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Add some logical clinical relationships
  metadata$BV_status <- ifelse(metadata$Nugent_score >= 7, "BV", 
                              ifelse(metadata$Nugent_score <= 3, "Normal", "Intermediate"))
  
  metadata$symptoms_score <- ifelse(metadata$symptoms == "None", 0,
                                  ifelse(metadata$symptoms == "Mild", 1,
                                        ifelse(metadata$symptoms == "Moderate", 2, 3)))
  
  # Create taxonomy
  n_taxa <- 50
  genera <- c("Lactobacillus", "Gardnerella", "Atopobium", "Prevotella", "Sneathia", 
              "Megasphaera", "Streptococcus", "Bifidobacterium", "Ureaplasma", "Mobiluncus")
  
  key_species <- c(
    "Lactobacillus_crispatus", "Lactobacillus_iners", "Lactobacillus_gasseri", "Lactobacillus_jensenii",
    "Gardnerella_vaginalis", "Atopobium_vaginae", "Prevotella_bivia", "Sneathia_amnii",
    "Megasphaera_type_1", "Streptococcus_agalactiae"
  )
  
  taxonomy <- data.frame(
    taxon_id = paste0("Taxon_", 1:n_taxa),
    Kingdom = "Bacteria",
    Phylum = sample(c("Firmicutes", "Bacteroidetes", "Actinobacteria", "Fusobacteria"), n_taxa, replace = TRUE),
    Class = sample(c("Bacilli", "Clostridia", "Bacteroidia", "Actinobacteria", "Fusobacteriia"), n_taxa, replace = TRUE),
    Order = sample(c("Lactobacillales", "Clostridiales", "Bacteroidales", "Bifidobacteriales", "Fusobacteriales"), n_taxa, replace = TRUE),
    Family = sample(c("Lactobacillaceae", "Bifidobacteriaceae", "Prevotellaceae", "Streptococcaceae"), n_taxa, replace = TRUE),
    Genus = sample(genera, n_taxa, replace = TRUE),
    Species = sample(c(key_species, paste0(sample(genera, n_taxa - length(key_species), replace = TRUE), 
                                          "_", sample(letters, n_taxa - length(key_species), replace = TRUE))),
                    n_taxa, replace = FALSE),
    stringsAsFactors = FALSE
  )
  
  # Create counts with patterns that align with clinical variables
  counts <- matrix(0, nrow = n_taxa, ncol = n_samples)
  rownames(counts) <- taxonomy$taxon_id
  colnames(counts) <- metadata$sample_id
  
  # Create realistic distributions based on clinical variables
  for (i in 1:n_samples) {
    # Get clinical status
    bv_status <- metadata$BV_status[i]
    
    # Base counts for all taxa
    counts[, i] <- sample(0:5, n_taxa, replace = TRUE)
    
    if (bv_status == "Normal") {
      # Lactobacillus dominated
      lactobacillus_idx <- which(taxonomy$Genus == "Lactobacillus")
      if (length(lactobacillus_idx) > 0) {
        # Choose one dominant Lactobacillus species
        dominant_idx <- sample(lactobacillus_idx, 1)
        counts[dominant_idx, i] <- sample(5000:20000, 1)
      }
    } else if (bv_status == "BV") {
      # BV associated bacteria
      bv_genera <- c("Gardnerella", "Atopobium", "Prevotella", "Sneathia", "Megasphaera")
      for (genus in bv_genera) {
        genus_idx <- which(taxonomy$Genus == genus)
        if (length(genus_idx) > 0) {
          for (j in sample(genus_idx, min(2, length(genus_idx)))) {
            counts[j, i] <- sample(1000:5000, 1)
          }
        }
      }
    } else {  # Intermediate
      # Mix of lactobacilli and BV bacteria
      lactobacillus_idx <- which(taxonomy$Genus == "Lactobacillus")
      if (length(lactobacillus_idx) > 0) {
        dominant_idx <- sample(lactobacillus_idx, 1)
        counts[dominant_idx, i] <- sample(2000:8000, 1)
      }
      
      # Some BV bacteria
      bv_genus <- sample(c("Gardnerella", "Prevotella"), 1)
      genus_idx <- which(taxonomy$Genus == bv_genus)
      if (length(genus_idx) > 0) {
        for (j in sample(genus_idx, min(1, length(genus_idx)))) {
          counts[j, i] <- sample(1000:3000, 1)
        }
      }
    }
  }
  
  # Create FGTExperiment object
  fgt_exp <- create_mock_fgt_experiment(
    counts = counts,
    sample_metadata = metadata,
    taxonomy = taxonomy,
    experiment_name = "Mock Dataset with Clinical Data"
  )
})

# Function to generate a microbiome report for a sample
generate_microbiome_report <- function(fgt_exp, sample_id, reference_population = NULL, 
                                      include_metrics = c("diversity", "key_taxa", "health_index")) {
  # Convert to relative abundance
  rel_abund <- transformAbundance(fgt_exp, type = "relative")
  
  # Get data
  counts_data <- assays(experimentData(rel_abund))$counts
  
  # Get metadata
  metadata <- as.data.frame(colData(experimentData(rel_abund)))
  
  # Get taxonomy
  taxa_info <- as.data.frame(rowData(experimentData(rel_abund)))
  
  # Check if sample exists
  if (!sample_id %in% colnames(counts_data)) {
    stop("Sample ID not found in dataset")
  }
  
  # Get sample data
  sample_data <- counts_data[, sample_id]
  sample_metadata <- metadata[metadata$sample_id == sample_id, ]
  
  # Define reference population if not provided
  if (is.null(reference_population)) {
    reference_data <- counts_data[, colnames(counts_data) != sample_id, drop = FALSE]
  } else if (is.character(reference_population)) {
    # Assume it's a column name in metadata
    if (!reference_population %in% colnames(metadata)) {
      stop("Reference population column not found in metadata")
    }
    
    reference_value <- sample_metadata[[reference_population]]
    reference_samples <- metadata$sample_id[metadata[[reference_population]] == reference_value & 
                                           metadata$sample_id != sample_id]
    
    if (length(reference_samples) == 0) {
      warning("No matching reference samples found. Using all other samples.")
      reference_samples <- setdiff(colnames(counts_data), sample_id)
    }
    
    reference_data <- counts_data[, reference_samples, drop = FALSE]
  } else if (is.numeric(reference_population) || is.logical(reference_population)) {
    # Assume it's a vector of sample indices or logical mask
    reference_samples <- colnames(counts_data)[reference_population]
    reference_samples <- setdiff(reference_samples, sample_id)
    reference_data <- counts_data[, reference_samples, drop = FALSE]
  } else {
    stop("Invalid reference_population format")
  }
  
  # Initialize report
  report <- list()
  
  # Basic information
  report$sample_id <- sample_id
  report$metadata <- sample_metadata
  report$timestamp <- Sys.time()
  
  # Calculate metrics
  report$metrics <- list()
  
  # Function to calculate alpha diversity
  calculate_alpha_diversity <- function(x, method = "shannon") {
    x <- x / sum(x)  # Normalize to relative abundance
    x <- x[x > 0]  # Remove zeros
    
    if (method == "shannon") {
      -sum(x * log(x))
    } else if (method == "simpson") {
      1 - sum(x^2)
    } else if (method == "invsimpson") {
      1 / sum(x^2)
    } else if (method == "richness") {
      sum(x > 0)
    } else {
      stop("Unknown diversity method")
    }
  }
  
  # Diversity metrics
  if ("diversity" %in% include_metrics) {
    # Calculate Shannon diversity
    report$metrics$shannon <- calculate_alpha_diversity(sample_data, "shannon")
    report$metrics$simpson <- calculate_alpha_diversity(sample_data, "simpson")
    report$metrics$richness <- calculate_alpha_diversity(sample_data, "richness")
    
    # Reference diversities
    ref_shannon <- apply(reference_data, 2, calculate_alpha_diversity, "shannon")
    ref_simpson <- apply(reference_data, 2, calculate_alpha_diversity, "simpson")
    ref_richness <- apply(reference_data, 2, calculate_alpha_diversity, "richness")
    
    # Percentiles
    report$metrics$shannon_percentile <- mean(report$metrics$shannon >= ref_shannon) * 100
    report$metrics$simpson_percentile <- mean(report$metrics$simpson >= ref_simpson) * 100
    report$metrics$richness_percentile <- mean(report$metrics$richness >= ref_richness) * 100
    
    # Reference stats
    report$metrics$ref_shannon_mean <- mean(ref_shannon)
    report$metrics$ref_shannon_sd <- sd(ref_shannon)
    report$metrics$ref_simpson_mean <- mean(ref_simpson)
    report$metrics$ref_simpson_sd <- sd(ref_simpson)
    report$metrics$ref_richness_mean <- mean(ref_richness)
    report$metrics$ref_richness_sd <- sd(ref_richness)
  }
  
  # Key taxa analysis
  if ("key_taxa" %in% include_metrics) {
    # Define key taxa of interest for FGT
    key_genera <- c("Lactobacillus", "Gardnerella", "Atopobium", "Prevotella", "Sneathia")
    
    # Extract abundance of key genera
    key_taxa_abundance <- list()
    
    for (genus in key_genera) {
      # Find index of the genus
      genus_idx <- which(taxa_info$Genus == genus)
      
      if (length(genus_idx) > 0) {
        # Sum abundance across all species of this genus
        abundance <- sum(sample_data[genus_idx])
        key_taxa_abundance[[genus]] <- abundance
        
        # Reference distribution
        ref_abundance <- colSums(reference_data[genus_idx, , drop = FALSE])
        key_taxa_abundance[[paste0(genus, "_percentile")]] <- mean(abundance >= ref_abundance) * 100
        key_taxa_abundance[[paste0(genus, "_ref_mean")]] <- mean(ref_abundance)
        key_taxa_abundance[[paste0(genus, "_ref_sd")]] <- sd(ref_abundance)
      } else {
        key_taxa_abundance[[genus]] <- 0
        key_taxa_abundance[[paste0(genus, "_percentile")]] <- NA
        key_taxa_abundance[[paste0(genus, "_ref_mean")]] <- NA
        key_taxa_abundance[[paste0(genus, "_ref_sd")]] <- NA
      }
    }
    
    report$metrics$key_taxa <- key_taxa_abundance
    
    # Calculate ratio of Lactobacillus to Gardnerella (L/G ratio)
    l_abundance <- key_taxa_abundance[["Lactobacillus"]]
    g_abundance <- key_taxa_abundance[["Gardnerella"]]
    
    if (g_abundance > 0) {
      report$metrics$lg_ratio <- l_abundance / g_abundance
    } else if (l_abundance > 0) {
      report$metrics$lg_ratio <- Inf
    } else {
      report$metrics$lg_ratio <- NA
    }
    
    # Calculate ratios for reference
    l_ref <- sapply(1:ncol(reference_data), function(i) {
      sum(reference_data[which(taxa_info$Genus == "Lactobacillus"), i])
    })
    
    g_ref <- sapply(1:ncol(reference_data), function(i) {
      sum(reference_data[which(taxa_info$Genus == "Gardnerella"), i])
    })
    
    lg_ratio_ref <- l_ref / g_ref
    lg_ratio_ref[is.infinite(lg_ratio_ref)] <- max(lg_ratio_ref[is.finite(lg_ratio_ref)]) * 2
    lg_ratio_ref[is.na(lg_ratio_ref)] <- 0
    
    report$metrics$lg_ratio_percentile <- mean(report$metrics$lg_ratio >= lg_ratio_ref, na.rm = TRUE) * 100
    report$metrics$lg_ratio_ref_mean <- mean(lg_ratio_ref, na.rm = TRUE)
    report$metrics$lg_ratio_ref_sd <- sd(lg_ratio_ref, na.rm = TRUE)
  }
  
  # Overall health index
  if ("health_index" %in% include_metrics) {
    # Calculate health index based on key metrics
    
    # Check for Lactobacillus dominance
    lacto_abundance <- sum(sample_data[which(taxa_info$Genus == "Lactobacillus")])
    lacto_dominance <- lacto_abundance > 0.5  # >50% Lactobacillus
    
    # Check for BV-associated bacteria
    bv_genera <- c("Gardnerella", "Atopobium", "Prevotella", "Sneathia", "Megasphaera")
    bv_abundance <- sum(sample_data[unlist(lapply(bv_genera, function(g) which(taxa_info$Genus == g)))])
    bv_presence <- bv_abundance > 0.2  # >20% BV-associated bacteria
    
    # Calculate health score (0-100)
    # Higher Lactobacillus, lower BV bacteria, and moderate diversity are preferred
    report$metrics$health_score <- 
      min(100, max(0,
                  50 * (lacto_abundance) +  # Up to 50 points for Lactobacillus
                  20 * (1 - bv_abundance) +  # Up to 20 points for low BV bacteria
                  15 * min(1, report$metrics$shannon / 2) +  # Up to 15 points for moderate diversity
                  15 * min(1, report$metrics$richness / 20)  # Up to 15 points for moderate richness
      ))
    
    # Classification based on score
    report$metrics$health_category <- cut(
      report$metrics$health_score,
      breaks = c(0, 25, 50, 75, 100),
      labels = c("Poor", "Fair", "Good", "Excellent"),
      include.lowest = TRUE
    )
    
    # Percentile compared to reference
    reference_scores <- sapply(1:ncol(reference_data), function(i) {
      sample_vec <- reference_data[, i]
      
      lacto_abund <- sum(sample_vec[which(taxa_info$Genus == "Lactobacillus")])
      bv_abund <- sum(sample_vec[unlist(lapply(bv_genera, function(g) which(taxa_info$Genus == g)))])
      shannon <- calculate_alpha_diversity(sample_vec, "shannon")
      richness <- calculate_alpha_diversity(sample_vec, "richness")
      
      min(100, max(0,
                  50 * (lacto_abund) +
                  20 * (1 - bv_abund) +
                  15 * min(1, shannon / 2) +
                  15 * min(1, richness / 20)
      ))
    })
    
    report$metrics$health_score_percentile <- mean(report$metrics$health_score >= reference_scores) * 100
    report$metrics$health_score_ref_mean <- mean(reference_scores)
    report$metrics$health_score_ref_sd <- sd(reference_scores)
  }
  
  # Generate plots
  report$plots <- list()
  
  # Composition plot (top 10 taxa)
  taxa_abundance <- sort(sample_data, decreasing = TRUE)
  top_taxa <- head(taxa_abundance, 10)
  other_abundance <- sum(taxa_abundance) - sum(top_taxa)
  
  # Create data for plotting
  plot_data <- data.frame(
    taxon = c(names(top_taxa), "Other"),
    abundance = c(top_taxa, other_abundance),
    stringsAsFactors = FALSE
  )
  
  # Add taxonomy information
  plot_data$genus <- taxa_info$Genus[match(plot_data$taxon, taxa_info$taxon_id)]
  plot_data$genus[plot_data$taxon == "Other"] <- "Other"
  
  # Replace taxon IDs with species names where available
  plot_data$label <- taxa_info$Species[match(plot_data$taxon, taxa_info$taxon_id)]
  plot_data$label[is.na(plot_data$label)] <- plot_data$genus[is.na(plot_data$label)]
  plot_data$label[plot_data$taxon == "Other"] <- "Other"
  
  # Sort by abundance
  plot_data <- plot_data[order(-plot_data$abundance), ]
  
  # Convert to factors for plotting
  plot_data$label <- factor(plot_data$label, levels = rev(plot_data$label))
  
  # Create composition plot
  composition_plot <- ggplot(plot_data, aes(x = "", y = abundance, fill = label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_viridis_d() +
    labs(
      title = "Microbiome Composition",
      subtitle = "Top 10 taxa",
      fill = "Taxon"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  report$plots$composition <- composition_plot
  
  # Diversity visualization
  if ("diversity" %in% include_metrics) {
    # Create data for density plot
    density_data <- data.frame(
      shannon = c(ref_shannon, report$metrics$shannon),
      sample = c(rep("Reference", length(ref_shannon)), "Current Sample"),
      stringsAsFactors = FALSE
    )
    
    # Create density plot
    diversity_plot <- ggplot(density_data, aes(x = shannon, fill = sample)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = report$metrics$shannon, linetype = "dashed", color = "red") +
      annotate("text", x = report$metrics$shannon, y = 0, 
               label = sprintf("%.2f (%d%%)", report$metrics$shannon, round(report$metrics$shannon_percentile)),
               hjust = -0.1, vjust = -0.5) +
      scale_fill_manual(values = c("Current Sample" = "red", "Reference" = "steelblue")) +
      labs(
        title = "Shannon Diversity",
        subtitle = "Compared to reference population",
        x = "Shannon Diversity Index",
        y = "Density",
        fill = ""
      ) +
      theme_minimal()
    
    report$plots$diversity <- diversity_plot
  }
  
  # Key taxa visualization
  if ("key_taxa" %in% include_metrics) {
    # Create data for bar plot
    key_taxa_data <- data.frame(
      genus = names(key_taxa_abundance)[!grepl("_percentile|_ref_", names(key_taxa_abundance))],
      abundance = unlist(key_taxa_abundance[!grepl("_percentile|_ref_", names(key_taxa_abundance))]),
      stringsAsFactors = FALSE
    )
    
    # Sort by abundance
    key_taxa_data <- key_taxa_data[order(-key_taxa_data$abundance), ]
    
    # Convert to factor for plotting
    key_taxa_data$genus <- factor(key_taxa_data$genus, levels = key_taxa_data$genus)
    
    # Create bar plot
    key_taxa_plot <- ggplot(key_taxa_data, aes(x = genus, y = abundance, fill = genus)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = sprintf("%.1f%%", abundance * 100)), 
                vjust = -0.5, size = 3) +
      scale_fill_viridis_d() +
      labs(
        title = "Key Taxa Abundance",
        x = "Genus",
        y = "Relative Abundance",
        fill = "Genus"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    report$plots$key_taxa <- key_taxa_plot
  }
  
  # Health score visualization
  if ("health_index" %in% include_metrics) {
    # Create gauge plot for health score
    gauge_data <- data.frame(
      value = report$metrics$health_score / 100,
      stringsAsFactors = FALSE
    )
    
    # Define breaks and colors
    breaks <- seq(0, 1, length.out = 101)
    colors <- colorRampPalette(c("red", "yellow", "green"))(100)
    
    # Create gauge plot
    health_gauge <- ggplot(gauge_data, aes(fill = value)) +
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "white") +
      annotate("text", x = 0.5, y = 0.65, label = "Microbiome Health Score",
               size = 5, fontface = "bold") +
      annotate("text", x = 0.5, y = 0.5, 
               label = paste0(round(report$metrics$health_score), "/100"),
               size = 8, fontface = "bold") +
      annotate("text", x = 0.5, y = 0.35, 
               label = paste0("Category: ", report$metrics$health_category),
               size = 5) +
      annotate("text", x = 0.5, y = 0.25, 
               label = paste0("Percentile: ", round(report$metrics$health_score_percentile), "%"),
               size = 4) +
      scale_fill_gradientn(colors = colors, limits = c(0, 1)) +
      coord_fixed() +
      theme_void() +
      theme(legend.position = "none") +
      
      # Add colored bar at the bottom to indicate score
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0.15, ymax = 0.16), fill = "gray80") +
      geom_rect(aes(xmin = 0, xmax = value, ymin = 0.15, ymax = 0.16), fill = "forestgreen") +
      
      # Add colored background based on category
      geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), 
                fill = case_when(
                  report$metrics$health_category == "Excellent" ~ "palegreen",
                  report$metrics$health_category == "Good" ~ "palegoldenrod",
                  report$metrics$health_category == "Fair" ~ "lightsalmon",
                  report$metrics$health_category == "Poor" ~ "lightpink",
                  TRUE ~ "white"
                ),
                alpha = 0.2)
    
    report$plots$health_gauge <- health_gauge
  }
  
  return(report)
}

# Function to create a combined report visual
create_report_visual <- function(report) {
  # Create a grid of plots
  title <- textGrob(paste("Microbiome Report for", report$sample_id),
                   gp = gpar(fontsize = 16, fontface = "bold"))
  
  subtitle <- textGrob(paste("Generated on", format(report$timestamp, "%B %d, %Y")),
                      gp = gpar(fontsize = 12))
  
  # Basic info
  info_text <- paste(
    "Sample ID:", report$sample_id, "\n",
    if (!is.null(report$metadata$patient_id)) paste("Patient ID:", report$metadata$patient_id, "\n") else "",
    if (!is.null(report$metadata$visit_type)) paste("Visit Type:", report$metadata$visit_type, "\n") else "",
    if (!is.null(report$metadata$visit_date)) paste("Visit Date:", format(report$metadata$visit_date, "%B %d, %Y"), "\n") else "",
    "\n",
    "MICROBIOME METRICS\n",
    "Shannon Diversity:", sprintf("%.2f", report$metrics$shannon), 
    sprintf("(Percentile: %d%%)", round(report$metrics$shannon_percentile)), "\n",
    "Species Richness:", sprintf("%d", report$metrics$richness),
    sprintf("(Percentile: %d%%)", round(report$metrics$richness_percentile)), "\n",
    "L/G Ratio:", ifelse(is.finite(report$metrics$lg_ratio), 
                        sprintf("%.2f", report$metrics$lg_ratio), "Infinity"),
    sprintf("(Percentile: %d%%)", round(report$metrics$lg_ratio_percentile)), "\n",
    "\n",
    "OVERALL HEALTH ASSESSMENT\n",
    "Health Score:", sprintf("%d/100", round(report$metrics$health_score)), "\n",
    "Category:", as.character(report$metrics$health_category), "\n",
    "Percentile:", sprintf("%d%%", round(report$metrics$health_score_percentile)), "\n"
  )
  
  info_grob <- textGrob(info_text, just = "left", x = 0.05,
                       gp = gpar(fontsize = 10))
  
  # Arrange plots
  plots <- list(
    composition = report$plots$composition,
    key_taxa = report$plots$key_taxa,
    diversity = report$plots$diversity,
    health = report$plots$health_gauge
  )
  
  # Remove any NULL plots
  plots <- Filter(Negate(is.null), plots)
  
  # Combine plots using patchwork
  if (length(plots) > 0) {
    combined_plots <- gridExtra::arrangeGrob(
      grobs = c(list(title, subtitle, info_grob), plots),
      layout_matrix = rbind(
        c(1, 1, 1, 1),
        c(2, 2, 2, 2),
        c(3, 3, 4, 4),
        c(3, 3, 4, 4),
        c(5, 5, 6, 6),
        c(5, 5, 6, 6)
      ),
      heights = c(0.5, 0.3, 1, 1, 1, 1)
    )
    
    return(combined_plots)
  } else {
    return(NULL)
  }
}

# Function to export the report
export_microbiome_report <- function(report, format = "png", output_dir = "inst/example_plots/reports/") {
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate visual report
  report_visual <- create_report_visual(report)
  
  if (is.null(report_visual)) {
    stop("Could not create report visual")
  }
  
  # Create filename
  filename <- file.path(output_dir, paste0("microbiome_report_", report$sample_id))
  
  # Export based on format
  if (format == "png") {
    png(paste0(filename, ".png"), width = 10, height = 12, units = "in", res = 300)
    grid.draw(report_visual)
    dev.off()
    return(paste0(filename, ".png"))
  } else if (format == "pdf") {
    pdf(paste0(filename, ".pdf"), width = 10, height = 12)
    grid.draw(report_visual)
    dev.off()
    return(paste0(filename, ".pdf"))
  } else if (format == "html") {
    # Create simple HTML report
    html_content <- paste0(
      "<!DOCTYPE html>\n",
      "<html>\n",
      "<head>\n",
      "  <title>Microbiome Report for ", report$sample_id, "</title>\n",
      "  <style>\n",
      "    body { font-family: Arial, sans-serif; margin: 20px; }\n",
      "    h1, h2 { color: #2c3e50; }\n",
      "    .container { max-width: 1000px; margin: 0 auto; }\n",
      "    .report-header { text-align: center; margin-bottom: 30px; }\n",
      "    .metrics { display: flex; flex-wrap: wrap; justify-content: space-between; margin-bottom: 30px; }\n",
      "    .metric-box { background-color: #f8f9fa; border-radius: 8px; padding: 15px; width: 45%; margin-bottom: 20px; }\n",
      "    .plots { display: flex; flex-wrap: wrap; justify-content: space-between; }\n",
      "    .plot-container { width: 45%; margin-bottom: 30px; }\n",
      "    .health-score { font-size: 24px; font-weight: bold; text-align: center; margin: 20px 0; }\n",
      "    .health-excellent { color: #27ae60; }\n",
      "    .health-good { color: #2ecc71; }\n",
      "    .health-fair { color: #f39c12; }\n",
      "    .health-poor { color: #e74c3c; }\n",
      "  </style>\n",
      "</head>\n",
      "<body>\n",
      "  <div class='container'>\n",
      "    <div class='report-header'>\n",
      "      <h1>Microbiome Report</h1>\n",
      "      <h2>", report$sample_id, "</h2>\n",
      "      <p>Generated on ", format(report$timestamp, "%B %d, %Y"), "</p>\n",
      "    </div>\n",
      "    \n",
      "    <div class='metrics'>\n",
      "      <div class='metric-box'>\n",
      "        <h3>Sample Information</h3>\n",
      "        <p><strong>Sample ID:</strong> ", report$sample_id, "</p>\n",
      if (!is.null(report$metadata$patient_id)) paste0("        <p><strong>Patient ID:</strong> ", report$metadata$patient_id, "</p>\n") else "",
      if (!is.null(report$metadata$visit_type)) paste0("        <p><strong>Visit Type:</strong> ", report$metadata$visit_type, "</p>\n") else "",
      if (!is.null(report$metadata$visit_date)) paste0("        <p><strong>Visit Date:</strong> ", format(report$metadata$visit_date, "%B %d, %Y"), "</p>\n") else "",
      "      </div>\n",
      "      \n",
      "      <div class='metric-box'>\n",
      "        <h3>Diversity Metrics</h3>\n",
      "        <p><strong>Shannon Diversity:</strong> ", sprintf("%.2f", report$metrics$shannon), 
      " (Percentile: ", round(report$metrics$shannon_percentile), "%)</p>\n",
      "        <p><strong>Species Richness:</strong> ", sprintf("%d", report$metrics$richness),
      " (Percentile: ", round(report$metrics$richness_percentile), "%)</p>\n",
      "        <p><strong>L/G Ratio:</strong> ", ifelse(is.finite(report$metrics$lg_ratio), 
                                                      sprintf("%.2f", report$metrics$lg_ratio), "Infinity"),
      " (Percentile: ", round(report$metrics$lg_ratio_percentile), "%)</p>\n",
      "      </div>\n",
      "    </div>\n",
      "    \n",
      "    <div class='health-score health-", tolower(as.character(report$metrics$health_category)), "'>\n",
      "      Health Score: ", round(report$metrics$health_score), "/100 (", as.character(report$metrics$health_category), ")\n",
      "    </div>\n",
      "    \n",
      "    <div class='plots'>\n",
      "      <div class='plot-container'>\n",
      "        <h3>Microbiome Composition</h3>\n",
      "        <img src='../images/composition_", report$sample_id, ".png' width='100%'>\n",
      "      </div>\n",
      "      \n",
      "      <div class='plot-container'>\n",
      "        <h3>Key Taxa Abundance</h3>\n",
      "        <img src='../images/key_taxa_", report$sample_id, ".png' width='100%'>\n",
      "      </div>\n",
      "      \n",
      "      <div class='plot-container'>\n",
      "        <h3>Diversity Distribution</h3>\n",
      "        <img src='../images/diversity_", report$sample_id, ".png' width='100%'>\n",
      "      </div>\n",
      "      \n",
      "      <div class='plot-container'>\n",
      "        <h3>Health Assessment</h3>\n",
      "        <img src='../images/health_gauge_", report$sample_id, ".png' width='100%'>\n",
      "      </div>\n",
      "    </div>\n",
      "  </div>\n",
      "</body>\n",
      "</html>"
    )
    
    # Save individual plots for HTML report
    ggsave(paste0(output_dir, "../images/composition_", report$sample_id, ".png"), 
           report$plots$composition, width = 6, height = 5, dpi = 150)
    
    ggsave(paste0(output_dir, "../images/key_taxa_", report$sample_id, ".png"), 
           report$plots$key_taxa, width = 6, height = 5, dpi = 150)
    
    ggsave(paste0(output_dir, "../images/diversity_", report$sample_id, ".png"), 
           report$plots$diversity, width = 6, height = 5, dpi = 150)
    
    ggsave(paste0(output_dir, "../images/health_gauge_", report$sample_id, ".png"), 
           report$plots$health_gauge, width = 6, height = 5, dpi = 150)
    
    # Write HTML file
    writeLines(html_content, paste0(filename, ".html"))
    return(paste0(filename, ".html"))
  } else {
    stop("Unsupported format. Use 'png', 'pdf', or 'html'.")
  }
}

# Generate report for a sample
# Get metadata to find a sample ID
metadata <- as.data.frame(colData(experimentData(fgt_exp)))
sample_id <- metadata$sample_id[1]  # Use first sample for example

# Generate report with all metrics
report <- generate_microbiome_report(
  fgt_exp, 
  sample_id = sample_id,
  include_metrics = c("diversity", "key_taxa", "health_index")
)

# Export as PNG
png_file <- export_microbiome_report(report, format = "png")

# Export as HTML
html_file <- export_microbiome_report(report, format = "html")

# Generate reports for multiple samples to show different profiles
if (ncol(fgt_exp) >= 3) {
  # Generate reports for additional samples
  for (i in 2:min(3, ncol(fgt_exp))) {
    sample_id <- metadata$sample_id[i]
    
    report <- generate_microbiome_report(
      fgt_exp, 
      sample_id = sample_id,
      include_metrics = c("diversity", "key_taxa", "health_index")
    )
    
    # Export
    export_microbiome_report(report, format = "png")
    export_microbiome_report(report, format = "html")
  }
}

# Print success message
cat("Personalized microbiome reports created successfully!\n")
cat("Reports saved to inst/example_plots/reports/\n")
cat("Additional plot images saved to inst/example_plots/images/\n")