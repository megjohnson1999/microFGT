#' @title Example Data Generation and Loading for microFGT
#' @description
#' This file contains functions for generating, loading, and exporting example datasets
#' for the microFGT package. These functions are used for testing, demonstrations, and
#' tutorials.

#' @importFrom utils read.csv write.csv write.table
#' @importFrom methods is new existsClass
#' @importFrom stats runif rpois
#' @importFrom SummarizedExperiment SummarizedExperiment assays assayNames rowData colData
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @importFrom S4Vectors SimpleList

#' Generate example FGT microbiome data
#'
#' Creates synthetic microbiome data that simulates realistic female genital tract
#' (FGT) community compositions.
#'
#' @param n_samples Number of samples to generate
#' @param n_features Number of features (taxa) to generate
#' @param sample_groups Character vector of group names (e.g., c("Healthy", "BV"))
#' @param group_proportions Numeric vector of proportions for each group
#' @param community_types Character vector of CSTs to include
#' @param sequencing_depth Range of sequencing depths (c(min, max))
#' @param include_tree Whether to generate a phylogenetic tree
#' @param random_seed Optional seed for reproducibility
#' @param output_dir Directory to save files (if NULL, no files are saved)
#' @param format Output format ("FGTExperiment", "list")
#'
#' @return An FGTExperiment object or list of data components
#' @export
#'
#' @examples
#' # Generate a small example dataset
#' example_data <- generate_fgt_example_data(
#'   n_samples = 10,
#'   n_features = 50,
#'   community_types = c("CST-I", "CST-III", "CST-IV")
#' )
#'
#' # Generate data with specific parameters
#' custom_data <- generate_fgt_example_data(
#'   n_samples = 20, 
#'   n_features = 100,
#'   sample_groups = c("Healthy", "BV", "Intermediate"),
#'   group_proportions = c(0.5, 0.3, 0.2)
#' )
generate_fgt_example_data <- function(
  n_samples = 20,
  n_features = 100,
  sample_groups = c("Healthy", "BV"),
  group_proportions = c(0.7, 0.3),
  community_types = c("CST-I", "CST-III", "CST-IV"),
  sequencing_depth = c(5000, 20000),
  include_tree = FALSE,
  random_seed = NULL,
  output_dir = NULL,
  format = "FGTExperiment"
) {
  # Input validation
  if (length(sample_groups) != length(group_proportions)) {
    stop("sample_groups and group_proportions must have the same length")
  }
  
  if (abs(sum(group_proportions) - 1) > 0.001) {
    stop("group_proportions must sum to 1")
  }
  
  valid_csts <- c("CST-I", "CST-II", "CST-III", "CST-IV", "CST-V")
  invalid_csts <- setdiff(community_types, valid_csts)
  if (length(invalid_csts) > 0) {
    stop("Invalid community types: ", paste(invalid_csts, collapse = ", "), 
         ". Valid types are: ", paste(valid_csts, collapse = ", "))
  }
  
  # Set random seed if provided
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  # Use pre-defined CST constants
  # Constants are defined in R/constants/cst_constants.R
  cst_taxa <- list(
    "CST-I" = list(
      dominant = "Lactobacillus crispatus",
      subdominant = c("Lactobacillus jensenii", "Lactobacillus gasseri"),
      rare = c("Gardnerella vaginalis", "Atopobium vaginae")
    ),
    "CST-II" = list(
      dominant = "Lactobacillus gasseri",
      subdominant = c("Lactobacillus crispatus", "Lactobacillus jensenii"),
      rare = c("Gardnerella vaginalis", "Atopobium vaginae")
    ),
    "CST-III" = list(
      dominant = "Lactobacillus iners",
      subdominant = c("Gardnerella vaginalis"),
      rare = c("Prevotella bivia", "Atopobium vaginae", "Megasphaera")
    ),
    "CST-IV" = list(
      dominant = c("Gardnerella vaginalis", "Prevotella bivia", "Atopobium vaginae"),
      subdominant = c("Sneathia", "Megasphaera", "Mobiluncus"),
      rare = c("Lactobacillus iners", "Dialister", "Mycoplasma hominis")
    ),
    "CST-V" = list(
      dominant = "Lactobacillus jensenii",
      subdominant = c("Lactobacillus crispatus", "Lactobacillus gasseri"),
      rare = c("Gardnerella vaginalis", "Atopobium vaginae")
    )
  )

  # Note: In a full implementation, we would use the more detailed constants
  # from the constants file, e.g.:
  # cst_taxa <- list(
  #   "CST-I" = CST_I_TAXA,
  #   "CST-III" = CST_III_TAXA,
  #   "CST-IV" = CST_IV_TAXA,
  #   "CST-V" = CST_V_TAXA
  # )
  
  # Create count matrix
  counts <- matrix(0, nrow = n_features, ncol = n_samples)
  rownames(counts) <- paste0("Feature", 1:n_features)
  colnames(counts) <- paste0("Sample", 1:n_samples)
  
  # Create taxonomy data
  taxonomy <- data.frame(
    Kingdom = rep("Bacteria", n_features),
    Phylum = character(n_features),
    Class = character(n_features),
    Order = character(n_features),
    Family = character(n_features),
    Genus = character(n_features),
    Species = character(n_features),
    stringsAsFactors = FALSE
  )
  rownames(taxonomy) <- rownames(counts)
  
  # Define typical FGT taxa classifications
  common_taxa <- list(
    "Lactobacillus crispatus" = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus crispatus"),
    "Lactobacillus iners" = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus iners"),
    "Lactobacillus jensenii" = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus jensenii"),
    "Lactobacillus gasseri" = c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus gasseri"),
    "Gardnerella vaginalis" = c("Bacteria", "Actinobacteria", "Actinobacteria", "Bifidobacteriales", "Bifidobacteriaceae", "Gardnerella", "Gardnerella vaginalis"),
    "Atopobium vaginae" = c("Bacteria", "Actinobacteria", "Coriobacteriia", "Coriobacteriales", "Coriobacteriaceae", "Atopobium", "Atopobium vaginae"),
    "Prevotella bivia" = c("Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella bivia"),
    "Sneathia" = c("Bacteria", "Fusobacteria", "Fusobacteriia", "Fusobacteriales", "Leptotrichiaceae", "Sneathia", "Sneathia sp."),
    "Megasphaera" = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Megasphaera", "Megasphaera sp."),
    "Mobiluncus" = c("Bacteria", "Actinobacteria", "Actinobacteria", "Actinomycetales", "Actinomycetaceae", "Mobiluncus", "Mobiluncus sp."),
    "Dialister" = c("Bacteria", "Firmicutes", "Negativicutes", "Veillonellales", "Veillonellaceae", "Dialister", "Dialister sp."),
    "Mycoplasma hominis" = c("Bacteria", "Tenericutes", "Mollicutes", "Mycoplasmatales", "Mycoplasmataceae", "Mycoplasma", "Mycoplasma hominis")
  )
  
  # Assign taxonomies to features
  feature_assignments <- sample(
    names(common_taxa),
    n_features,
    replace = TRUE
  )
  
  for (i in 1:n_features) {
    taxon <- feature_assignments[i]
    if (taxon %in% names(common_taxa)) {
      taxonomy[i, 2:7] <- common_taxa[[taxon]][2:7]
    } else {
      # Fallback for unknown taxa
      taxonomy[i, "Phylum"] <- sample(c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria"), 1)
      taxonomy[i, "Genus"] <- sample(c("Lactobacillus", "Gardnerella", "Prevotella", "Atopobium"), 1)
      taxonomy[i, "Species"] <- paste0(taxonomy[i, "Genus"], " sp.")
    }
  }
  
  # Assign samples to groups based on proportions
  sample_groups_assignment <- character(n_samples)
  group_counts <- round(n_samples * group_proportions)
  
  # Adjust last group count to ensure sum is n_samples
  group_counts[length(group_counts)] <- n_samples - sum(group_counts[-length(group_counts)])
  
  idx <- 1
  for (g in 1:length(sample_groups)) {
    sample_groups_assignment[idx:(idx + group_counts[g] - 1)] <- sample_groups[g]
    idx <- idx + group_counts[g]
  }
  
  # Shuffle sample group assignments
  sample_groups_assignment <- sample(sample_groups_assignment)
  
  # Assign CSTs - each group tends to have specific CSTs
  cst_assignments <- character(n_samples)
  
  for (i in 1:n_samples) {
    if (sample_groups_assignment[i] == "Healthy") {
      cst_assignments[i] <- sample(c("CST-I", "CST-II", "CST-III", "CST-V"), 1, 
                                 prob = c(0.4, 0.1, 0.4, 0.1))
    } else if (sample_groups_assignment[i] == "BV") {
      cst_assignments[i] <- sample(c("CST-III", "CST-IV"), 1, prob = c(0.3, 0.7))
    } else {
      cst_assignments[i] <- sample(community_types, 1)
    }
    
    # Ensure selected CST is in the allowed list
    if (!cst_assignments[i] %in% community_types) {
      cst_assignments[i] <- sample(community_types, 1)
    }
  }
  
  # Generate abundance data based on CSTs
  for (s in 1:n_samples) {
    # Get assigned CST
    cst <- cst_assignments[s]
    
    # Generate depths
    depth <- round(runif(1, sequencing_depth[1], sequencing_depth[2]))
    
    # Generate abundances based on CST
    abundances <- numeric(n_features)
    
    # Find features matching dominant taxa for this CST
    dominant_taxa <- cst_taxa[[cst]]$dominant
    dominant_features <- which(feature_assignments %in% dominant_taxa)
    
    # If no exact matches, use the genus level
    if (length(dominant_features) == 0) {
      dominant_genera <- sapply(dominant_taxa, function(x) {
        strsplit(x, " ")[[1]][1]
      })
      dominant_features <- which(taxonomy$Genus %in% dominant_genera)
    }
    
    # Still no matches, use random features
    if (length(dominant_features) == 0) {
      dominant_features <- sample(1:n_features, min(3, n_features))
    }
    
    # Assign high abundance to dominant features
    abundances[dominant_features] <- runif(length(dominant_features), 0.6, 0.95)
    
    # Find features matching subdominant taxa
    subdominant_taxa <- cst_taxa[[cst]]$subdominant
    subdominant_features <- which(feature_assignments %in% subdominant_taxa)
    
    # If no exact matches, use the genus level
    if (length(subdominant_features) == 0) {
      subdominant_genera <- sapply(subdominant_taxa, function(x) {
        strsplit(x, " ")[[1]][1]
      })
      subdominant_features <- which(taxonomy$Genus %in% subdominant_genera)
    }
    
    # Get a few subdominant features
    subdominant_features <- setdiff(subdominant_features, dominant_features)
    if (length(subdominant_features) > 0) {
      subdominant_features <- sample(subdominant_features, 
                                    min(length(subdominant_features), 5))
      abundances[subdominant_features] <- runif(length(subdominant_features), 0.05, 0.3)
    }
    
    # Assign small abundance to remaining features
    remaining_features <- setdiff(1:n_features, c(dominant_features, subdominant_features))
    
    if (length(remaining_features) > 0) {
      # Only make some features present
      present_features <- sample(remaining_features, 
                              min(length(remaining_features), rpois(1, 10)))
      abundances[present_features] <- runif(length(present_features), 0.001, 0.05)
    }
    
    # Normalize to sum to 1
    abundances <- abundances / sum(abundances)
    
    # Convert to counts
    counts[, s] <- round(abundances * depth)
  }
  
  # Create metadata
  metadata <- data.frame(
    condition = sample_groups_assignment,
    community_state_type = cst_assignments,
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- colnames(counts)
  
  # Add clinical parameters based on group
  metadata$pH <- numeric(n_samples)
  metadata$Nugent_Score <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    if (metadata$condition[i] == "Healthy") {
      metadata$pH[i] <- runif(1, 3.8, 4.5)
      metadata$Nugent_Score[i] <- sample(0:3, 1)
    } else if (metadata$condition[i] == "BV") {
      metadata$pH[i] <- runif(1, 4.5, 6.0)
      metadata$Nugent_Score[i] <- sample(7:10, 1)
    } else {
      metadata$pH[i] <- runif(1, 4.3, 5.2)
      metadata$Nugent_Score[i] <- sample(4:6, 1)
    }
  }
  
  # Add categorical pH and BV diagnosis
  metadata$pH_category <- cut(metadata$pH, 
                            breaks = c(0, 4.5, 5.5, 10), 
                            labels = c("Normal", "Intermediate", "High"))
  
  metadata$BV_Nugent <- cut(metadata$Nugent_Score, 
                          breaks = c(-1, 3, 6, 10), 
                          labels = c("Negative", "Intermediate", "Positive"))
  
  # Create phylogenetic tree if requested
  tree <- NULL
  if (include_tree && requireNamespace("ape", quietly = TRUE)) {
    # Create a tree based on taxonomy
    dist_mat <- matrix(0, nrow = n_features, ncol = n_features)
    rownames(dist_mat) <- colnames(dist_mat) <- rownames(counts)
    
    for (i in 1:(n_features-1)) {
      for (j in (i+1):n_features) {
        # Calculate distance based on taxonomy
        same_levels <- 0
        for (level in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
          if (!is.na(taxonomy[i, level]) && !is.na(taxonomy[j, level]) && 
              taxonomy[i, level] == taxonomy[j, level]) {
            same_levels <- same_levels + 1
          } else {
            break
          }
        }
        
        # Convert to distance
        dist_mat[i, j] <- dist_mat[j, i] <- 1 - (same_levels / 6)
      }
    }
    
    # Create tree with neighbor joining algorithm
    tree <- ape::nj(dist_mat)
    tree <- ape::root(tree, 1)
  }
  
  # Save to files if output_dir is provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    saveRDS(counts, file.path(output_dir, "counts.rds"))
    saveRDS(taxonomy, file.path(output_dir, "taxonomy.rds"))
    saveRDS(metadata, file.path(output_dir, "metadata.rds"))
    
    if (!is.null(tree)) {
      saveRDS(tree, file.path(output_dir, "tree.rds"))
    }
  }
  
  # Create FGTExperiment or return as list
  if (format == "FGTExperiment") {
    fgt_exp <- FGTExperiment(
      assays = list(counts = counts),
      rowData = taxonomy,
      colData = metadata,
      rowTree = tree,
      experimentType = "amplicon"
    )
    
    return(fgt_exp)
  } else {
    # Return as list
    return(list(
      counts = counts,
      taxonomy = taxonomy,
      metadata = metadata,
      tree = tree
    ))
  }
}

#' Load example data
#'
#' Loads pre-built example data for testing and demonstration.
#'
#' @param size Size of dataset ("small", "medium", "large")
#' @param type Type of data ("amplicon", "metagenomic")
#' @param as_fgt_experiment Whether to return as FGTExperiment or list
#'
#' @return An FGTExperiment object or list of data components
#' @export
#'
#' @examples
#' # Load small amplicon dataset
#' example_data <- load_example_data("small", "amplicon")
#'
#' # Examine the data
#' example_data
load_example_data <- function(size = "small", type = "amplicon", as_fgt_experiment = TRUE) {
  # Validate inputs
  size <- match.arg(size, c("small", "medium", "large"))
  type <- match.arg(type, c("amplicon", "metagenomic"))
  
  # Construct file paths
  base_name <- paste0("microFGT_example_", size, "_", type)
  data_dir <- system.file("extdata", package = "microFGT")
  
  if (data_dir == "") {
    stop("Package 'microFGT' not properly installed or example data not found.")
  }
  
  # Define file paths
  count_file <- file.path(data_dir, paste0(base_name, "_counts.rds"))
  tax_file <- file.path(data_dir, paste0(base_name, "_taxonomy.rds"))
  meta_file <- file.path(data_dir, paste0(base_name, "_metadata.rds"))
  tree_file <- file.path(data_dir, paste0(base_name, "_tree.rds"))
  
  # Check if files exist
  if (!file.exists(count_file) || !file.exists(tax_file) || 
      !file.exists(meta_file)) {
    # Check if we have any example data
    available_files <- list.files(data_dir, pattern = "^microFGT_example_.+\\.rds$")
    
    if (length(available_files) == 0) {
      # No pre-built examples found, generate one on the fly
      message("No example data found. Generating a small example dataset.")
      
      # Set parameters based on size
      n_samples <- switch(size,
                        "small" = 10,
                        "medium" = 30,
                        "large" = 100)
      
      n_features <- switch(size,
                          "small" = 50,
                          "medium" = 150,
                          "large" = 500)
      
      # Generate example data
      result <- generate_fgt_example_data(
        n_samples = n_samples,
        n_features = n_features,
        format = if(as_fgt_experiment) "FGTExperiment" else "list"
      )
      
      return(result)
    } else {
      # Find what is available
      available_datasets <- unique(gsub("_(counts|taxonomy|metadata|tree)\\.rds$", "", available_files))
      stop("Requested example data (", base_name, ") not found. Available datasets: ", 
           paste(available_datasets, collapse = ", "))
    }
  }
  
  # Load data files
  counts <- readRDS(count_file)
  taxonomy <- readRDS(tax_file)
  metadata <- readRDS(meta_file)
  
  # Load tree if it exists
  tree <- NULL
  if (file.exists(tree_file)) {
    tree <- readRDS(tree_file)
  }
  
  # Return as list if not converting to FGTExperiment
  if (!as_fgt_experiment) {
    return(list(
      counts = counts,
      taxonomy = taxonomy,
      metadata = metadata,
      tree = tree
    ))
  }
  
  # Create FGTExperiment object
  fgt_exp <- FGTExperiment(
    assays = list(counts = counts),
    rowData = taxonomy,
    colData = metadata,
    rowTree = tree,
    experimentType = type
  )
  
  return(fgt_exp)
}

#' Export Example Data to Common Bioinformatics Formats
#'
#' This function exports microFGT example data to formats compatible with
#' common bioinformatics tools like phyloseq, DADA2, or QIIME2.
#'
#' @param data Example data (list or FGTExperiment object)
#' @param output_dir Directory to save the exported files
#' @param format Format to export to ("phyloseq", "dada2", "qiime2")
#' @param ... Additional parameters passed to format-specific export functions
#'
#' @return Invisibly returns a list of the exported file paths
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' example_data <- generate_fgt_example_data()
#'
#' # Export to DADA2 format
#' export_example_data(example_data, tempdir(), format = "dada2")
#'
#' # Export to QIIME2 format
#' export_example_data(example_data, tempdir(), format = "qiime2")
#' }
#' @export
export_example_data <- function(data, output_dir, format = c("phyloseq", "dada2", "qiime2"), ...) {
  format <- match.arg(format)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract components from FGTExperiment if needed
  if (methods::is(data, "FGTExperiment")) {
    counts <- SummarizedExperiment::assays(data)$counts
    taxonomy <- SummarizedExperiment::rowData(data)
    metadata <- SummarizedExperiment::colData(data)
    tree <- NULL
    if (methods::is(data, "TreeSummarizedExperiment")) {
      tree <- TreeSummarizedExperiment::rowTree(data)
    }
  } else if (is.list(data)) {
    counts <- data$counts
    taxonomy <- data$taxonomy
    metadata <- data$metadata
    tree <- data$tree
  } else {
    stop("data must be an FGTExperiment object or a list with counts, taxonomy, and metadata components")
  }
  
  # Export to the requested format
  if (format == "phyloseq") {
    return(export_to_phyloseq(counts, taxonomy, metadata, tree, output_dir, ...))
  } else if (format == "dada2") {
    return(export_to_dada2(counts, taxonomy, metadata, output_dir, ...))
  } else if (format == "qiime2") {
    return(export_to_qiime2(counts, taxonomy, metadata, tree, output_dir, ...))
  }
}

#' Export Example Data to phyloseq Format
#'
#' @param counts Count matrix
#' @param taxonomy Taxonomy table
#' @param metadata Sample metadata
#' @param tree Phylogenetic tree
#' @param output_dir Directory to save the exported files
#' @param base_name Base name for the output files
#'
#' @return Invisibly returns a list of the exported file paths
#'
#' @keywords internal
export_to_phyloseq <- function(counts, taxonomy, metadata, tree = NULL, 
                             output_dir, base_name = "phyloseq") {
  # Check if phyloseq package is available
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required for this export format")
  }
  
  # Create phyloseq object
  otu_table <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  tax_table <- phyloseq::tax_table(as.matrix(taxonomy))
  sample_data <- phyloseq::sample_data(metadata)
  
  # Create complete phyloseq object
  if (!is.null(tree) && requireNamespace("ape", quietly = TRUE)) {
    phy_tree <- phyloseq::phy_tree(tree)
    ps <- phyloseq::phyloseq(otu_table, tax_table, sample_data, phy_tree)
  } else {
    ps <- phyloseq::phyloseq(otu_table, tax_table, sample_data)
  }
  
  # Save to RDS file
  file_path <- file.path(output_dir, paste0(base_name, ".rds"))
  saveRDS(ps, file_path)
  
  # Return file path
  invisible(list(phyloseq = file_path))
}

#' Export Example Data to DADA2 Format
#'
#' @param counts Count matrix
#' @param taxonomy Taxonomy table
#' @param metadata Sample metadata
#' @param output_dir Directory to save the exported files
#' @param base_name Base name for the output files
#'
#' @return Invisibly returns a list of the exported file paths
#'
#' @keywords internal
export_to_dada2 <- function(counts, taxonomy, metadata, output_dir, 
                          base_name = "dada2") {
  # Convert taxonomy to DADA2 format if needed
  tax_cols <- colnames(taxonomy)
  
  # Check if taxonomy has the expected columns
  expected_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (!all(expected_cols %in% tax_cols)) {
    # Create a matrix with the expected format
    tax_matrix <- matrix(NA, nrow = nrow(taxonomy), ncol = length(expected_cols))
    colnames(tax_matrix) <- expected_cols
    rownames(tax_matrix) <- rownames(taxonomy)
    
    # Fill in the available levels
    for (col in intersect(tax_cols, expected_cols)) {
      tax_matrix[, col] <- taxonomy[, col]
    }
    
    # Convert to DADA2 format (semicolon-separated)
    tax_strings <- apply(tax_matrix, 1, function(row) {
      paste(sapply(seq_along(row), function(i) {
        level <- expected_cols[i]
        value <- row[i]
        if (is.na(value)) {
          paste0(level, ":NA")
        } else {
          paste0(level, ":", value)
        }
      }), collapse = ";")
    })
    
    tax_dada2 <- matrix(tax_strings, ncol = 1)
    colnames(tax_dada2) <- "Taxonomy"
    rownames(tax_dada2) <- rownames(taxonomy)
  } else {
    # Just format as semicolon-separated string
    tax_dada2 <- apply(taxonomy, 1, function(row) {
      paste(sapply(seq_along(row), function(i) {
        paste0(colnames(taxonomy)[i], ":", row[i])
      }), collapse = ";")
    })
    tax_dada2 <- matrix(tax_dada2, ncol = 1)
    colnames(tax_dada2) <- "Taxonomy"
    rownames(tax_dada2) <- rownames(taxonomy)
  }
  
  # Save files
  seqtab_file <- file.path(output_dir, paste0(base_name, "_seqtab.rds"))
  taxa_file <- file.path(output_dir, paste0(base_name, "_taxa.rds"))
  metadata_file <- file.path(output_dir, paste0(base_name, "_metadata.csv"))
  
  saveRDS(counts, seqtab_file)
  saveRDS(tax_dada2, taxa_file)
  utils::write.csv(metadata, metadata_file, row.names = TRUE)
  
  # Return file paths
  invisible(list(
    seqtab = seqtab_file,
    taxa = taxa_file,
    metadata = metadata_file
  ))
}

#' Export Example Data to QIIME2 Format
#'
#' @param counts Count matrix
#' @param taxonomy Taxonomy table
#' @param metadata Sample metadata
#' @param tree Phylogenetic tree
#' @param output_dir Directory to save the exported files
#' @param base_name Base name for the output files
#'
#' @return Invisibly returns a list of the exported file paths
#'
#' @keywords internal
export_to_qiime2 <- function(counts, taxonomy, metadata, tree = NULL, 
                           output_dir, base_name = "qiime2") {
  # Prepare feature table (biom format)
  feature_table_file <- file.path(output_dir, paste0(base_name, "_feature-table.tsv"))
  
  # Write feature table
  feature_table <- data.frame(
    "#OTU ID" = rownames(counts),
    counts,
    check.names = FALSE
  )
  utils::write.table(feature_table, file = feature_table_file, 
                   sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Prepare taxonomy file
  taxonomy_file <- file.path(output_dir, paste0(base_name, "_taxonomy.tsv"))
  
  # Format taxonomy for QIIME2
  tax_strings <- apply(taxonomy, 1, function(row) {
    paste(sapply(seq_along(row), function(i) {
      paste0(colnames(taxonomy)[i], "__", row[i])
    }), collapse = "; ")
  })
  
  tax_qiime <- data.frame(
    "Feature ID" = rownames(taxonomy),
    "Taxon" = tax_strings,
    "Confidence" = 1.0,
    check.names = FALSE
  )
  
  utils::write.table(tax_qiime, file = taxonomy_file, 
                   sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Prepare metadata file
  metadata_file <- file.path(output_dir, paste0(base_name, "_metadata.tsv"))
  
  # Format metadata for QIIME2
  metadata_qiime <- data.frame(
    "#SampleID" = rownames(metadata),
    metadata,
    check.names = FALSE
  )
  
  utils::write.table(metadata_qiime, file = metadata_file, 
                   sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Prepare tree file if provided
  tree_file <- NULL
  if (!is.null(tree) && requireNamespace("ape", quietly = TRUE)) {
    tree_file <- file.path(output_dir, paste0(base_name, "_tree.nwk"))
    ape::write.tree(tree, file = tree_file)
  }
  
  # Return file paths
  result <- list(
    feature_table = feature_table_file,
    taxonomy = taxonomy_file,
    metadata = metadata_file
  )
  
  if (!is.null(tree_file)) {
    result$tree <- tree_file
  }
  
  invisible(result)
}