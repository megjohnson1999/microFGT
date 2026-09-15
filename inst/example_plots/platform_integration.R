#' Platform Integration Visualization
#' 
#' This script demonstrates how microFGT integrates data from multiple platforms
#' (SpeciateIT, VALENCIA, and VIRGO) for comprehensive analysis.

# Create output directory if it doesn't exist
dir.create("inst/example_plots/images/integration", recursive = TRUE, showWarnings = FALSE)

# Load required libraries
library(microFGT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

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

# Generate mock data for platform integration visualization
generate_platform_data <- function() {
  set.seed(42)
  
  # Common sample set for integration
  n_samples <- 30
  sample_ids <- paste0("S", 1:n_samples)
  
  # Create common metadata structure
  metadata <- data.frame(
    sample_id = sample_ids,
    subject_id = rep(paste0("Subject_", 1:10), each = 3),
    timepoint = rep(1:3, 10),
    age = sample(21:45, n_samples, replace = TRUE),
    pH = round(runif(n_samples, 3.8, 5.5), 1),
    stringsAsFactors = FALSE
  )
  
  # Add Nugent score (clinical value)
  metadata$nugent_score <- round(runif(n_samples, 0, 10))
  metadata$bv_status <- ifelse(metadata$nugent_score < 4, "Negative", 
                               ifelse(metadata$nugent_score > 6, "Positive", "Intermediate"))
  
  # ---------------------------------------------------------------
  # 1. Generate SpeciateIT data (focused on taxonomic composition)
  # ---------------------------------------------------------------
  speciateit_taxa <- c(
    "Lactobacillus_crispatus", "Lactobacillus_iners", "Lactobacillus_gasseri", 
    "Lactobacillus_jensenii", "Gardnerella_vaginalis", "Atopobium_vaginae", 
    "Prevotella_bivia", "Sneathia_amnii", "Megasphaera_type_1", 
    "Mobiluncus_curtisii", "Streptococcus_agalactiae", "Escherichia_coli"
  )
  n_speciateit_taxa <- length(speciateit_taxa)
  
  # Generate counts based on BV status
  speciateit_counts <- matrix(0, nrow = n_speciateit_taxa, ncol = n_samples)
  rownames(speciateit_counts) <- speciateit_taxa
  colnames(speciateit_counts) <- sample_ids
  
  for (i in 1:n_samples) {
    bv_status <- metadata$bv_status[i]
    
    # Base low counts for all taxa
    speciateit_counts[, i] <- round(runif(n_speciateit_taxa, 0, 10))
    
    if (bv_status == "Negative") {
      # Healthy - dominated by Lactobacillus species
      lacto_taxa <- grep("Lactobacillus", speciateit_taxa)
      main_lacto <- sample(lacto_taxa, 1)
      speciateit_counts[main_lacto, i] <- round(runif(1, 8000, 15000))
      
      # Some low levels of other lactobacilli
      other_lacto <- setdiff(lacto_taxa, main_lacto)
      speciateit_counts[other_lacto, i] <- round(runif(length(other_lacto), 100, 1000))
      
    } else if (bv_status == "Positive") {
      # BV - diverse anaerobes
      bv_taxa <- c("Gardnerella_vaginalis", "Atopobium_vaginae", 
                    "Prevotella_bivia", "Sneathia_amnii", "Megasphaera_type_1", 
                    "Mobiluncus_curtisii")
      bv_idx <- match(bv_taxa, speciateit_taxa)
      bv_idx <- bv_idx[!is.na(bv_idx)]
      
      for (idx in bv_idx) {
        speciateit_counts[idx, i] <- round(runif(1, 1000, 8000))
      }
      
      # Some Lactobacillus iners may be present in BV
      l_iners_idx <- which(speciateit_taxa == "Lactobacillus_iners")
      speciateit_counts[l_iners_idx, i] <- round(runif(1, 500, 2000))
      
    } else {
      # Intermediate - mixed profile
      # Some Lactobacillus
      l_iners_idx <- which(speciateit_taxa == "Lactobacillus_iners")
      speciateit_counts[l_iners_idx, i] <- round(runif(1, 2000, 5000))
      
      # Some BV-associated taxa
      gv_idx <- which(speciateit_taxa == "Gardnerella_vaginalis")
      av_idx <- which(speciateit_taxa == "Atopobium_vaginae")
      speciateit_counts[gv_idx, i] <- round(runif(1, 1000, 4000))
      speciateit_counts[av_idx, i] <- round(runif(1, 500, 2000))
    }
  }
  
  # Create taxonomy information for SpeciateIT
  speciateit_taxonomy <- data.frame(
    taxon_id = speciateit_taxa,
    Kingdom = "Bacteria",
    Phylum = c(rep("Firmicutes", 4), "Actinobacteria", "Actinobacteria", 
               "Bacteroidetes", "Fusobacteria", "Firmicutes", 
               "Actinobacteria", "Firmicutes", "Proteobacteria"),
    Class = c(rep("Bacilli", 4), "Actinobacteria", "Coriobacteriia", 
              "Bacteroidia", "Fusobacteriia", "Negativicutes", 
              "Actinobacteria", "Bacilli", "Gammaproteobacteria"),
    Order = c(rep("Lactobacillales", 4), "Bifidobacteriales", "Coriobacteriales", 
              "Bacteroidales", "Fusobacteriales", "Veillonellales", 
              "Actinomycetales", "Lactobacillales", "Enterobacterales"),
    Family = c(rep("Lactobacillaceae", 4), "Bifidobacteriaceae", "Coriobacteriaceae", 
               "Prevotellaceae", "Leptotrichiaceae", "Veillonellaceae", 
               "Actinomycetaceae", "Streptococcaceae", "Enterobacteriaceae"),
    Genus = c(rep("Lactobacillus", 4), "Gardnerella", "Atopobium", 
              "Prevotella", "Sneathia", "Megasphaera", 
              "Mobiluncus", "Streptococcus", "Escherichia"),
    Species = speciateit_taxa,
    stringsAsFactors = FALSE
  )
  
  # Create SpeciateIT FGTExperiment
  speciateit_exp <- create_mock_fgt_experiment(
    counts = speciateit_counts,
    sample_metadata = metadata,
    taxonomy = speciateit_taxonomy,
    experiment_name = "Mock SpeciateIT Dataset"
  )
  
  # ---------------------------------------------------------------
  # 2. Generate VALENCIA data (focused on CST classification)
  # ---------------------------------------------------------------
  # VALENCIA assigns Community State Types
  cst_types <- c("CST I", "CST II", "CST III", "CST IV-A", "CST IV-B", "CST V")
  
  # Assign CSTs based on the SpeciateIT dominant taxa
  valencia_metadata <- metadata
  valencia_metadata$cst <- NA
  
  for (i in 1:n_samples) {
    sample_counts <- speciateit_counts[, i]
    dominant_taxon <- speciateit_taxa[which.max(sample_counts)]
    
    if (dominant_taxon == "Lactobacillus_crispatus") {
      valencia_metadata$cst[i] <- "CST I"
    } else if (dominant_taxon == "Lactobacillus_gasseri") {
      valencia_metadata$cst[i] <- "CST II"
    } else if (dominant_taxon == "Lactobacillus_iners") {
      valencia_metadata$cst[i] <- "CST III"
    } else if (dominant_taxon == "Lactobacillus_jensenii") {
      valencia_metadata$cst[i] <- "CST V"
    } else {
      # For diverse anaerobic communities
      if (runif(1) < 0.7) {
        valencia_metadata$cst[i] <- "CST IV-A"  # Dominated by G. vaginalis
      } else {
        valencia_metadata$cst[i] <- "CST IV-B"  # Other anaerobes
      }
    }
  }
  
  # Generate auxiliary data that might come from VALENCIA
  valencia_metadata$shannon_diversity <- runif(n_samples, 0, 3)
  valencia_metadata$simpson_diversity <- runif(n_samples, 0, 1)
  valencia_metadata$lactobacillus_rel_abundance <- numeric(n_samples)
  
  # Calculate Lactobacillus abundance
  for (i in 1:n_samples) {
    lacto_idx <- grep("Lactobacillus", rownames(speciateit_counts))
    total_abundance <- sum(speciateit_counts[, i])
    lacto_abundance <- sum(speciateit_counts[lacto_idx, i])
    valencia_metadata$lactobacillus_rel_abundance[i] <- lacto_abundance / total_abundance
  }
  
  # ---------------------------------------------------------------
  # 3. Generate VIRGO data (focused on functional profiles)
  # ---------------------------------------------------------------
  # VIRGO would have genes, pathways, etc.
  n_virgo_functions <- 15
  virgo_functions <- c(
    "Lactic_acid_production", "Hydrogen_peroxide_production", 
    "Glycogen_metabolism", "Biofilm_formation", "Bacteriocin_production",
    "Adhesion_to_mucus", "D-lactic_acid_production", "L-lactic_acid_production",
    "Sialidase_activity", "Cholesterol_metabolism", "Antimicrobial_resistance",
    "Epithelial_cell_cytotoxicity", "Immunomodulation", "Vitamin_production",
    "Short_chain_fatty_acid_production"
  )
  
  # Generate functional abundance based on taxa
  virgo_counts <- matrix(0, nrow = n_virgo_functions, ncol = n_samples)
  rownames(virgo_counts) <- virgo_functions
  colnames(virgo_counts) <- sample_ids
  
  for (i in 1:n_samples) {
    bv_status <- metadata$bv_status[i]
    
    # Base counts for all functions
    virgo_counts[, i] <- runif(n_virgo_functions, 0, 0.1)
    
    if (bv_status == "Negative") {
      # Lactobacillus dominated - high lactic acid, H2O2, etc.
      virgo_counts["Lactic_acid_production", i] <- runif(1, 0.7, 1.0)
      virgo_counts["Hydrogen_peroxide_production", i] <- runif(1, 0.5, 0.9)
      virgo_counts["Bacteriocin_production", i] <- runif(1, 0.4, 0.8)
      virgo_counts["Adhesion_to_mucus", i] <- runif(1, 0.5, 0.9)
      
      # L vs D lactic acid depends on the Lactobacillus species
      cst <- valencia_metadata$cst[i]
      if (cst == "CST I" || cst == "CST V") {
        virgo_counts["D-lactic_acid_production", i] <- runif(1, 0.7, 0.9)
        virgo_counts["L-lactic_acid_production", i] <- runif(1, 0.1, 0.3)
      } else {
        virgo_counts["D-lactic_acid_production", i] <- runif(1, 0.1, 0.3)
        virgo_counts["L-lactic_acid_production", i] <- runif(1, 0.7, 0.9)
      }
      
    } else if (bv_status == "Positive") {
      # BV - high sialidase, cytotoxicity, etc.
      virgo_counts["Sialidase_activity", i] <- runif(1, 0.6, 0.9)
      virgo_counts["Biofilm_formation", i] <- runif(1, 0.5, 0.8)
      virgo_counts["Epithelial_cell_cytotoxicity", i] <- runif(1, 0.4, 0.7)
      virgo_counts["Short_chain_fatty_acid_production", i] <- runif(1, 0.5, 0.8)
      virgo_counts["Antimicrobial_resistance", i] <- runif(1, 0.3, 0.6)
      
    } else {
      # Intermediate - mixed profile
      virgo_counts["Lactic_acid_production", i] <- runif(1, 0.3, 0.6)
      virgo_counts["Sialidase_activity", i] <- runif(1, 0.2, 0.5)
      virgo_counts["Biofilm_formation", i] <- runif(1, 0.3, 0.6)
      virgo_counts["Adhesion_to_mucus", i] <- runif(1, 0.2, 0.5)
    }
  }
  
  # Create VIRGO-specific metadata
  virgo_metadata <- metadata
  
  # Add metabolomic correlations (mock data)
  virgo_metadata$lactic_acid_mmol <- 10 * virgo_counts["Lactic_acid_production",] + rnorm(n_samples, 0, 0.5)
  virgo_metadata$acetate_mmol <- 8 * virgo_counts["Short_chain_fatty_acid_production",] + rnorm(n_samples, 0, 0.3)
  virgo_metadata$propionate_mmol <- 5 * virgo_counts["Short_chain_fatty_acid_production",] + rnorm(n_samples, 0, 0.2)
  virgo_metadata$butyrate_mmol <- 3 * virgo_counts["Short_chain_fatty_acid_production",] + rnorm(n_samples, 0, 0.1)
  
  # Function annotation
  virgo_taxonomy <- data.frame(
    function_id = virgo_functions,
    Category = c(rep("Metabolic", 3), "Community", "Defense", 
                "Adhesion", rep("Metabolic", 2), "Enzymatic", 
                "Metabolic", "Defense", "Virulence", "Host_interaction", 
                "Metabolic", "Metabolic"),
    Pathway = c("Fermentation", "Defense", "Carbohydrate_metabolism", 
                "Community_structure", "Defense_mechanism", "Adhesion_factor", 
                "Fermentation", "Fermentation", "Glycoside_hydrolase", 
                "Lipid_metabolism", "Drug_resistance", "Virulence_factor", 
                "Immune_modulation", "Micronutrients", "Fermentation"),
    stringsAsFactors = FALSE
  )
  
  # Create VIRGO experiment object
  virgo_exp <- create_mock_fgt_experiment(
    counts = virgo_counts,
    sample_metadata = virgo_metadata,
    taxonomy = virgo_taxonomy,
    experiment_name = "Mock VIRGO Dataset"
  )
  
  # Return all platform data
  return(list(
    speciateit = speciateit_exp,
    valencia = valencia_metadata,
    virgo = virgo_exp,
    shared_metadata = metadata
  ))
}

# Generate the integrated platform data
platform_data <- generate_platform_data()

# -----------------------------------------------------------------
# Visualization 1: Multi-platform Composition Plot
# -----------------------------------------------------------------
# Prepare relative abundance data
rel_abundances <- transformAbundance(platform_data$speciateit, "relative")$counts

# Focus on top taxa
top_taxa <- order(rowMeans(rel_abundances), decreasing = TRUE)[1:8]
top_taxa_names <- rownames(rel_abundances)[top_taxa]
other_abundance <- colSums(rel_abundances[-top_taxa,, drop = FALSE])

# Create plotting data frame
plot_data <- data.frame(
  sample_id = rep(colnames(rel_abundances), each = length(top_taxa) + 1),
  taxon = c(rep(top_taxa_names, ncol(rel_abundances)), rep("Other", ncol(rel_abundances))),
  abundance = c(as.vector(rel_abundances[top_taxa,]), other_abundance),
  stringsAsFactors = FALSE
)

# Add metadata
plot_data <- merge(plot_data, platform_data$valencia, by = "sample_id")

# Order samples by CST
plot_data$sample_id <- factor(plot_data$sample_id, 
                             levels = unique(platform_data$valencia$sample_id[order(platform_data$valencia$cst)]))

# Create the integrated stacked bar plot
multi_platform_plot <- ggplot(plot_data, aes(x = sample_id, y = abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ cst, scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "lightblue", color = "gray"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Taxonomic Composition by Community State Type",
    subtitle = "Integration of SpeciateIT composition data with VALENCIA CST classification",
    y = "Relative Abundance",
    fill = "Taxon"
  )

# Save the plot
ggsave("inst/example_plots/images/integration/multi_platform_integration.png", 
       multi_platform_plot, width = 12, height = 6, dpi = 300)

# -----------------------------------------------------------------
# Visualization 2: Correlation Heatmap between Function and Taxonomy
# -----------------------------------------------------------------
# Prepare data
taxa_abundance <- rel_abundances
function_abundance <- platform_data$virgo$counts

# Calculate correlation between top taxa and all functions
top_taxa_abundance <- taxa_abundance[top_taxa,]
correlation_matrix <- matrix(0, nrow = length(top_taxa), ncol = nrow(function_abundance))
rownames(correlation_matrix) <- rownames(taxa_abundance)[top_taxa]
colnames(correlation_matrix) <- rownames(function_abundance)

# Calculate correlations
for (i in 1:length(top_taxa)) {
  for (j in 1:nrow(function_abundance)) {
    correlation_matrix[i, j] <- cor(top_taxa_abundance[i,], function_abundance[j,], method = "spearman")
  }
}

# Plot correlation heatmap
correlation_heatmap <- pheatmap(
  correlation_matrix,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  breaks = seq(-1, 1, length.out = 101),
  border_color = "white",
  main = "Taxa-Function Correlation Heatmap\nIntegration of SpeciateIT with VIRGO data",
  angle_col = 45,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 14,
  cellheight = 14
)

# Save as PNG
png("inst/example_plots/images/integration/correlation_heatmap.png", 
    width = 10, height = 8, units = "in", res = 300)
correlation_heatmap
dev.off()

# -----------------------------------------------------------------
# Visualization 3: CST-Metabolites Relationship
# -----------------------------------------------------------------
# Prepare data for boxplots
cst_metabolite_data <- platform_data$valencia
cst_metabolite_data$lactic_acid_mmol <- platform_data$virgo$metadata$lactic_acid_mmol
cst_metabolite_data$acetate_mmol <- platform_data$virgo$metadata$acetate_mmol
cst_metabolite_data$propionate_mmol <- platform_data$virgo$metadata$propionate_mmol
cst_metabolite_data$butyrate_mmol <- platform_data$virgo$metadata$butyrate_mmol

# Reshape for plotting
cst_metabolite_long <- tidyr::pivot_longer(
  cst_metabolite_data,
  cols = c("lactic_acid_mmol", "acetate_mmol", "propionate_mmol", "butyrate_mmol"),
  names_to = "metabolite",
  values_to = "concentration"
)

# Format metabolite names
cst_metabolite_long$metabolite <- factor(cst_metabolite_long$metabolite,
                                      levels = c("lactic_acid_mmol", "acetate_mmol", "propionate_mmol", "butyrate_mmol"),
                                      labels = c("Lactic Acid", "Acetate", "Propionate", "Butyrate"))

# Create boxplot
cst_metabolite_plot <- ggplot(cst_metabolite_long, aes(x = cst, y = concentration, fill = cst)) +
  geom_boxplot() +
  facet_wrap(~ metabolite, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "lightblue", color = "gray"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Metabolite Concentrations by Community State Type",
    subtitle = "Integration of VALENCIA CST classification with VIRGO functional data",
    x = "Community State Type",
    y = "Concentration (mmol)"
  )

# Save the plot
ggsave("inst/example_plots/images/integration/valencia_cst_metabolites.png", 
       cst_metabolite_plot, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------
# Visualization 4: Multi-platform Profile Heatmap
# -----------------------------------------------------------------
# Prepare data for the integrated heatmap
# 1. Top taxa from SpeciateIT
top_taxa_data <- rel_abundances[top_taxa,]

# 2. Scale and prepare the VIRGO functional data
selected_functions <- c("Lactic_acid_production", "Hydrogen_peroxide_production", 
                        "Sialidase_activity", "Biofilm_formation", 
                        "Epithelial_cell_cytotoxicity", "Short_chain_fatty_acid_production")
function_data <- platform_data$virgo$counts[selected_functions,]

# 3. Prepare clinical data
clinical_data <- data.frame(
  pH = scale(platform_data$shared_metadata$pH),
  Nugent = scale(platform_data$shared_metadata$nugent_score)
)
rownames(clinical_data) <- platform_data$shared_metadata$sample_id

# Combine all data for the heatmap
combined_data <- rbind(
  top_taxa_data,
  function_data,
  t(clinical_data)
)

# Add some metadata for annotation
annotation_col <- data.frame(
  CST = platform_data$valencia$cst,
  BV_Status = platform_data$shared_metadata$bv_status
)
rownames(annotation_col) <- colnames(combined_data)

# Prepare row annotation to distinguish data types
annotation_row <- data.frame(
  DataType = c(
    rep("Taxonomy", nrow(top_taxa_data)),
    rep("Function", nrow(function_data)),
    rep("Clinical", nrow(t(clinical_data)))
  )
)
rownames(annotation_row) <- rownames(combined_data)

# Define annotation colors
ann_colors <- list(
  CST = c("CST I" = "#E41A1C", "CST II" = "#377EB8", "CST III" = "#4DAF4A", 
          "CST IV-A" = "#984EA3", "CST IV-B" = "#FF7F00", "CST V" = "#FFFF33"),
  BV_Status = c("Negative" = "#91CF60", "Intermediate" = "#FFFFBF", "Positive" = "#FC8D59"),
  DataType = c("Taxonomy" = "#8DD3C7", "Function" = "#BEBADA", "Clinical" = "#FDB462")
)

# Generate the heatmap
integration_heatmap <- pheatmap(
  combined_data,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  scale = "row",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cutree_cols = 6,
  border_color = NA,
  main = "Multi-Platform Integration Heatmap",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45
)

# Save as PNG
png("inst/example_plots/images/integration/integration_heatmap.png", 
    width = 12, height = 8, units = "in", res = 300)
integration_heatmap
dev.off()

# -----------------------------------------------------------------
# Visualization 5: SpeciateIT Composition by Platform
# -----------------------------------------------------------------
# Prepare data - same sample but across different platforms
# For the plot, we'll use the CST info from VALENCIA, and 
# function info from VIRGO to augment the SpeciateIT data

# Create a scatterplot with pH vs. lactic acid, colored by dominant taxa
scatter_data <- data.frame(
  sample_id = platform_data$shared_metadata$sample_id,
  pH = platform_data$shared_metadata$pH,
  lactic_acid = platform_data$virgo$metadata$lactic_acid_mmol,
  stringsAsFactors = FALSE
)

# Add dominant taxon information from SpeciateIT
dominant_taxa <- apply(rel_abundances, 2, function(x) rownames(rel_abundances)[which.max(x)])
scatter_data$dominant_taxon <- dominant_taxa

# Simplify the taxon names for better visualization
scatter_data$dominant_group <- "Other"
scatter_data$dominant_group[grep("Lactobacillus_crispatus", scatter_data$dominant_taxon)] <- "L. crispatus"
scatter_data$dominant_group[grep("Lactobacillus_iners", scatter_data$dominant_taxon)] <- "L. iners"
scatter_data$dominant_group[grep("Lactobacillus_", scatter_data$dominant_taxon)] <- "Other Lactobacillus"
scatter_data$dominant_group[grep("Gardnerella", scatter_data$dominant_taxon)] <- "Gardnerella"

# Add CST information from VALENCIA
scatter_data$cst <- platform_data$valencia$cst

# Create the scatter plot
ph_lactobacillus_plot <- ggplot(scatter_data, aes(x = pH, y = lactic_acid, color = dominant_group, shape = cst)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = 1:6) +
  theme_minimal() +
  labs(
    title = "pH vs. Lactic Acid Production",
    subtitle = "Integration of SpeciateIT taxonomy, VIRGO metabolic data, and VALENCIA CST classification",
    x = "Vaginal pH",
    y = "Lactic Acid (mmol)",
    color = "Dominant Taxa",
    shape = "CST"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Save the plot
ggsave("inst/example_plots/images/integration/ph_lactobacillus_correlation.png", 
       ph_lactobacillus_plot, width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------
# Visualization 6: Platform Comparison - Individual Taxa
# -----------------------------------------------------------------

# SpeciateIT - Focus on key taxa across multiple samples
speciateit_key_taxa <- c("Lactobacillus_crispatus", "Lactobacillus_iners", 
                         "Gardnerella_vaginalis", "Atopobium_vaginae")
speciateit_key_data <- rel_abundances[match(speciateit_key_taxa, rownames(rel_abundances)),]

# VIRGO - Focus on key functions
virgo_key_functions <- c("Lactic_acid_production", "Sialidase_activity")
virgo_key_data <- platform_data$virgo$counts[match(virgo_key_functions, rownames(platform_data$virgo$counts)),]

# Lactobacillus crispatus plot
# ---------------------------
l_crispatus_idx <- which(rownames(rel_abundances) == "Lactobacillus_crispatus")
l_crispatus_data <- data.frame(
  sample_id = colnames(rel_abundances),
  abundance = rel_abundances[l_crispatus_idx,],
  lactic_acid = platform_data$virgo$counts["Lactic_acid_production",],
  ph = platform_data$shared_metadata$pH,
  cst = platform_data$valencia$cst,
  stringsAsFactors = FALSE
)

# Sort by abundance
l_crispatus_data <- l_crispatus_data[order(l_crispatus_data$abundance, decreasing = TRUE),]
l_crispatus_data$sample_id <- factor(l_crispatus_data$sample_id, levels = l_crispatus_data$sample_id)

# Create plot
l_crispatus_plot <- ggplot(l_crispatus_data) +
  geom_col(aes(x = sample_id, y = abundance, fill = cst), width = 0.7) +
  geom_line(aes(x = as.integer(sample_id), y = lactic_acid), 
            color = "darkred", size = 1.5, group = 1) +
  geom_point(aes(x = as.integer(sample_id), y = lactic_acid), 
             color = "darkred", size = 3) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    name = "L. crispatus Abundance",
    sec.axis = sec_axis(~ ., name = "Lactic Acid Production")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = "Lactobacillus crispatus Abundance vs. Lactic Acid Production",
    subtitle = "Integration of SpeciateIT taxonomy with VIRGO functional data",
    fill = "VALENCIA CST"
  )

# Gardnerella vaginalis plot
# --------------------------
g_vaginalis_idx <- which(rownames(rel_abundances) == "Gardnerella_vaginalis")
g_vaginalis_data <- data.frame(
  sample_id = colnames(rel_abundances),
  abundance = rel_abundances[g_vaginalis_idx,],
  sialidase = platform_data$virgo$counts["Sialidase_activity",],
  biofilm = platform_data$virgo$counts["Biofilm_formation",],
  cst = platform_data$valencia$cst,
  stringsAsFactors = FALSE
)

# Sort by abundance
g_vaginalis_data <- g_vaginalis_data[order(g_vaginalis_data$abundance, decreasing = TRUE),]
g_vaginalis_data$sample_id <- factor(g_vaginalis_data$sample_id, levels = g_vaginalis_data$sample_id)

# Create plot
g_vaginalis_plot <- ggplot(g_vaginalis_data) +
  geom_col(aes(x = sample_id, y = abundance, fill = cst), width = 0.7) +
  geom_line(aes(x = as.integer(sample_id), y = sialidase), 
            color = "darkblue", size = 1.5, group = 1) +
  geom_point(aes(x = as.integer(sample_id), y = sialidase), 
             color = "darkblue", size = 3) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    name = "G. vaginalis Abundance",
    sec.axis = sec_axis(~ ., name = "Sialidase Activity")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = "Gardnerella vaginalis Abundance vs. Sialidase Activity",
    subtitle = "Integration of SpeciateIT taxonomy with VIRGO functional data",
    fill = "VALENCIA CST"
  )

# Combine the two plots
combined_taxa_plots <- l_crispatus_plot / g_vaginalis_plot + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Save the combined plot
ggsave("inst/example_plots/images/integration/virgo_gardnerella_lactobacillus.png", 
       combined_taxa_plots, width = 12, height = 10, dpi = 300)

# -----------------------------------------------------------------
# Visualization 7: CST Transition Profiles
# -----------------------------------------------------------------
# Create mock longitudinal data to show CST transitions over time
# Get unique subject IDs
subjects <- unique(platform_data$shared_metadata$subject_id)
n_subjects <- length(subjects)

# For each subject, get the timepoints
timepoint_data <- list()
for (subj in subjects) {
  # Get this subject's samples
  subject_samples <- platform_data$shared_metadata$sample_id[platform_data$shared_metadata$subject_id == subj]
  subject_timepoints <- platform_data$shared_metadata$timepoint[platform_data$shared_metadata$subject_id == subj]
  
  # Sort by timepoint
  ord <- order(subject_timepoints)
  subject_samples <- subject_samples[ord]
  subject_timepoints <- subject_timepoints[ord]
  
  # Get CST for each timepoint
  subject_cst <- platform_data$valencia$cst[match(subject_samples, platform_data$valencia$sample_id)]
  
  # Create data frame for this subject
  timepoint_data[[subj]] <- data.frame(
    subject_id = subj,
    timepoint = subject_timepoints,
    cst = subject_cst,
    sample_id = subject_samples,
    stringsAsFactors = FALSE
  )
}

# Combine all subjects
transition_data <- do.call(rbind, timepoint_data)

# Create temporal CST profile plot
cst_transition_plot <- ggplot(transition_data, aes(x = timepoint, y = subject_id, fill = cst)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(breaks = unique(transition_data$timepoint)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(hjust = 1)
  ) +
  labs(
    title = "CST Transitions Over Time",
    subtitle = "VALENCIA CST classification across longitudinal samples",
    x = "Timepoint",
    y = "Subject ID",
    fill = "Community State Type (CST)"
  )

# Save the plot
ggsave("inst/example_plots/images/integration/cst_transitions_over_time.png", 
       cst_transition_plot, width = 10, height = 8, dpi = 300)

# Print success message
cat("Platform integration plots created successfully!\n")
cat("Plots saved to inst/example_plots/images/integration/\n")