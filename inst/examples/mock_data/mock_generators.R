#' Mock Data Generators for microFGT
#'
#' A collection of functions to generate realistic mock data for SpeciateIT, VIRGO, and VALENCIA
#' tools that can be used for testing and development purposes.
#'
#' @section Key Functions:
#' - `generate_mock_speciateit()`: Generate mock SpeciateIT taxonomic classifications
#' - `generate_mock_virgo()`: Generate mock VIRGO gene abundance data
#' - `generate_mock_valencia()`: Generate mock VALENCIA community state type classifications
#' - `generate_mock_fgt_dataset()`: Generate a complete set of mock data for all three tools
#'
#' @docType package
#' @name mock_data_generators
NULL

#' Generate mock SpeciateIT results
#'
#' Creates simulated speciateIT classification results for testing purposes.
#' The function generates data that mimics the format of the "MC_order7_results.txt"
#' file produced by speciateIT.
#'
#' @param n_sequences Number of sequences to generate
#' @param seq_id_prefix Prefix for sequence IDs
#' @param species_distribution Named vector of probabilities for different species. 
#'   Names should be species names, values should be probabilities that sum to 1.
#' @param confidence_range Range of confidence values (posterior probabilities) to generate
#' @param decisions_range Range of decision counts to generate
#' @param include_unclassified Proportion of sequences that should be unclassified (0-1)
#' @param seed Random seed for reproducibility
#'
#' @return A data frame with the columns: "Sequence ID", "Classification", 
#'   "posterior probability", and "number of Decisions"
#' @export
#'
#' @examples
#' \donttest{
#' # Generate default mock data
#' mock_data <- generate_mock_speciateit(n_sequences = 100)
#'
#' # Generate data with custom species distribution
#' species_dist <- c(
#'   "Lactobacillus crispatus" = 0.5,
#'   "Lactobacillus iners" = 0.3, 
#'   "Gardnerella vaginalis" = 0.2
#' )
#' mock_data <- generate_mock_speciateit(
#'   n_sequences = 100, 
#'   species_distribution = species_dist
#' )
#' }
generate_mock_speciateit <- function(
  n_sequences = 100,
  seq_id_prefix = "Seq",
  species_distribution = NULL,
  confidence_range = c(0.7, 1.0),
  decisions_range = c(5, 20),
  include_unclassified = 0.05,
  seed = NULL
) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Default species distribution for vaginal microbiome
  if (is.null(species_distribution)) {
    species_distribution <- c(
      "Lactobacillus crispatus" = 0.35,
      "Lactobacillus iners" = 0.25,
      "Lactobacillus gasseri" = 0.10,
      "Lactobacillus jensenii" = 0.05,
      "Gardnerella vaginalis" = 0.10,
      "Gardnerella leopoldii" = 0.02,
      "Gardnerella swidsinskii" = 0.02,
      "Gardnerella piotii" = 0.01,
      "Prevotella bivia" = 0.03,
      "Atopobium vaginae" = 0.03,
      "Sneathia amnii" = 0.02,
      "Megasphaera genomosp type 1" = 0.01,
      "Mobiluncus curtisii" = 0.01
    )
  }
  
  # Normalize species distribution
  species_distribution <- species_distribution / sum(species_distribution)
  species_names <- names(species_distribution)
  species_probs <- as.numeric(species_distribution)
  
  # Generate sequence IDs
  sequence_ids <- paste0(seq_id_prefix, "_", sprintf("%06d", 1:n_sequences))
  
  # Generate classifications
  if (include_unclassified > 0 && include_unclassified < 1) {
    # Determine how many sequences will be unclassified
    n_unclassified <- floor(n_sequences * include_unclassified)
    n_classified <- n_sequences - n_unclassified
    
    # For classified sequences
    classified_indices <- 1:n_classified
    classified_taxa <- sample(species_names, n_classified, replace = TRUE, prob = species_probs)
    
    # For unclassified sequences
    if (n_unclassified > 0) {  # Only do this if there are unclassified sequences
      unclassified_indices <- (n_classified + 1):n_sequences
      unclassified <- rep("Unclassified", n_unclassified)
      
      # Combine
      all_classifications <- character(n_sequences)
      all_classifications[classified_indices] <- classified_taxa
      all_classifications[unclassified_indices] <- unclassified
      
      # Shuffle to randomize positions
      rand_order <- sample(n_sequences)
      classifications <- all_classifications[rand_order]
      sequence_ids <- sequence_ids[rand_order]
    } else {
      # No unclassified sequences needed
      classifications <- classified_taxa
    }
  } else {
    # All sequences classified
    classifications <- sample(species_names, n_sequences, replace = TRUE, prob = species_probs)
  }
  
  # Generate confidence scores (posterior probabilities)
  confidence_scores <- runif(
    n_sequences, 
    min = confidence_range[1], 
    max = confidence_range[2]
  )
  
  # Set unclassified confidence to lower values
  if (include_unclassified > 0) {
    unclassified_indices <- which(classifications == "Unclassified")
    if (length(unclassified_indices) > 0) {
      confidence_scores[unclassified_indices] <- runif(
        length(unclassified_indices),
        min = 0.3,
        max = 0.6
      )
    }
  }
  
  # Generate decision counts
  decision_counts <- sample(
    decisions_range[1]:decisions_range[2],
    n_sequences,
    replace = TRUE
  )
  
  # Create data frame
  results <- data.frame(
    "Sequence ID" = sequence_ids,
    "Classification" = classifications,
    "posterior probability" = confidence_scores,
    "number of Decisions" = decision_counts,
    stringsAsFactors = FALSE,
    check.names = FALSE  # Preserve exact column names
  )
  
  return(results)
}

#' Write mock SpeciateIT results to file
#'
#' Generates mock SpeciateIT data and writes it to a file in the format
#' expected by the microFGT package.
#'
#' @param file_path Path to the output file
#' @param ... Additional parameters passed to `generate_mock_speciateit()`
#'
#' @return Path to the created file (invisibly)
#' @export
#'
#' @examples
#' \donttest{
#' # Generate and save default mock data
#' output_file <- write_mock_speciateit(tempfile(), n_sequences = 10)
#'
#' # Generate with custom parameters
#' output_file <- write_mock_speciateit(
#'   tempfile(),
#'   n_sequences = 20,
#'   include_unclassified = 0.1
#' )
#' }
write_mock_speciateit <- function(file_path, ...) {
  results <- generate_mock_speciateit(...)
  write.table(
    results, 
    file = file_path, 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
  )
  return(invisible(file_path))
}

#' Convert mock SpeciateIT results to FGTMicrobiome format
#'
#' Converts the output from `generate_mock_speciateit()` to a format compatible
#' with the `add_speciateit()` function.
#'
#' @param speciateit_results Data frame from `generate_mock_speciateit()`
#'
#' @return A list with an `assignments` component containing Feature and Species columns
#' @export
#'
#' @examples
#' \donttest{
#' # Generate mock data
#' mock_data <- generate_mock_speciateit(n_sequences = 20)
#'
#' # Convert to format for add_speciateit()
#' fgt_compatible <- convert_speciateit_to_fgt(mock_data)
#' }
convert_speciateit_to_fgt <- function(speciateit_results) {
  # Create assignments data frame
  assignments <- data.frame(
    Feature = speciateit_results[["Sequence ID"]],
    Species = speciateit_results[["Classification"]],
    Confidence = speciateit_results[["posterior probability"]],
    stringsAsFactors = FALSE
  )
  
  # Create result list
  result <- list(
    assignments = assignments,
    version = "mock"
  )
  
  return(result)
}

#' Generate mock VIRGO results
#'
#' Creates simulated VIRGO gene abundance results for testing purposes.
#' The function generates data that mimics the format of the VIRGO output files.
#'
#' @param n_genes Number of genes to include (VIRGO contains ~118K genes)
#' @param n_samples Number of samples to generate
#' @param sample_id_prefix Prefix for sample IDs
#' @param include_zeros Proportion of zero counts to include (sparsity)
#' @param abundance_distribution Distribution to use for non-zero counts ("lognormal" or "zipf")
#' @param abundance_params Parameters for the abundance distribution
#'   - For "lognormal": c(meanlog, sdlog)
#'   - For "zipf": c(alpha)
#' @param simulation_type Type of simulation to perform:
#'   - "realistic": Realistic vaginal microbiome with CST structure
#'   - "random": Random gene abundances without structure
#' @param cst_distribution Named vector of CST probabilities (only for realistic simulation)
#' @param seed Random seed for reproducibility
#'
#' @return A list with components:
#'   - `counts`: Matrix of gene counts (genes × samples)
#'   - `genes`: Data frame with gene information
#'   - `metadata`: Data frame with sample metadata
#' @export
#'
#' @examples
#' \donttest{
#' # Generate default mock VIRGO data
#' virgo_data <- generate_mock_virgo(n_genes = 1000, n_samples = 10)
#'
#' # Generate data with specific CST distribution
#' cst_dist <- c("CST I" = 0.5, "CST III" = 0.3, "CST IV" = 0.2)
#' virgo_data <- generate_mock_virgo(
#'   n_genes = 1000, 
#'   n_samples = 20,
#'   simulation_type = "realistic",
#'   cst_distribution = cst_dist
#' )
#' }
generate_mock_virgo <- function(
  n_genes = 1000,
  n_samples = 10,
  sample_id_prefix = "Sample",
  include_zeros = 0.7,
  abundance_distribution = "lognormal",
  abundance_params = c(1, 2),
  simulation_type = "realistic",
  cst_distribution = NULL,
  seed = NULL
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Default CST distribution
  if (is.null(cst_distribution)) {
    cst_distribution <- c(
      "CST I" = 0.35,    # L. crispatus dominated
      "CST II" = 0.10,   # L. gasseri dominated
      "CST III" = 0.25,  # L. iners dominated
      "CST V" = 0.10,    # L. jensenii dominated
      "CST IV-A" = 0.10, # G. vaginalis dominated
      "CST IV-B" = 0.05, # Diverse anaerobes
      "CST IV-C" = 0.05  # Mixed Lactobacillus and anaerobes
    )
  }
  
  # Normalize CST distribution
  cst_distribution <- cst_distribution / sum(cst_distribution)
  
  # Generate sample IDs
  sample_ids <- paste0(sample_id_prefix, "_", sprintf("%03d", 1:n_samples))
  
  # Generate gene IDs (VIRGO uses "V" followed by 7 digits)
  gene_ids <- paste0("V", sprintf("%07d", 1:n_genes))
  
  # Generate gene lengths (typically 100-5000bp)
  gene_lengths <- round(runif(n_genes, min = 100, max = 5000))
  
  # Create gene info data frame
  gene_info <- data.frame(
    gene_id = gene_ids,
    length = gene_lengths,
    stringsAsFactors = FALSE
  )
  
  # Assign CSTs to samples for realistic simulation
  if (simulation_type == "realistic") {
    # Assign CSTs based on distribution
    sample_csts <- sample(
      names(cst_distribution),
      n_samples,
      replace = TRUE,
      prob = as.numeric(cst_distribution)
    )
    
    # Define gene function groups (simplified)
    function_groups <- c(
      "L_crispatus_specific",
      "L_iners_specific",
      "L_gasseri_specific",
      "L_jensenii_specific",
      "G_vaginalis_specific",
      "Anaerobe_specific",
      "Housekeeping",
      "Mobile_elements",
      "Unknown",
      "CRISPR"
    )
    
    # Assign genes to function groups
    n_housekeeping <- floor(n_genes * 0.1)
    n_mobile <- floor(n_genes * 0.05)
    n_unknown <- floor(n_genes * 0.15)
    n_crispr <- floor(n_genes * 0.02)
    n_specific <- n_genes - n_housekeeping - n_mobile - n_unknown - n_crispr
    
    # Distribute specific genes among taxa
    n_l_crispatus <- floor(n_specific * 0.25)
    n_l_iners <- floor(n_specific * 0.2)
    n_l_gasseri <- floor(n_specific * 0.1)
    n_l_jensenii <- floor(n_specific * 0.1)
    n_g_vaginalis <- floor(n_specific * 0.15)
    n_anaerobe <- n_specific - n_l_crispatus - n_l_iners - n_l_gasseri - n_l_jensenii - n_g_vaginalis
    
    # Assign function groups
    gene_functions <- c(
      rep("L_crispatus_specific", n_l_crispatus),
      rep("L_iners_specific", n_l_iners),
      rep("L_gasseri_specific", n_l_gasseri),
      rep("L_jensenii_specific", n_l_jensenii),
      rep("G_vaginalis_specific", n_g_vaginalis),
      rep("Anaerobe_specific", n_anaerobe),
      rep("Housekeeping", n_housekeeping),
      rep("Mobile_elements", n_mobile),
      rep("Unknown", n_unknown),
      rep("CRISPR", n_crispr)
    )
    gene_info$function_group <- sample(gene_functions)
    
    # Create count matrix based on CSTs
    counts_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
    
    # Function to get non-zero counts
    get_nonzero_counts <- function(n, dist = abundance_distribution, params = abundance_params) {
      if (dist == "lognormal") {
        return(round(exp(rnorm(n, mean = params[1], sd = params[2]))))
      } else if (dist == "zipf") {
        # Zipf distribution
        ranks <- 1:n
        return(round(n * ranks^(-params[1])))
      } else {
        stop("Unsupported distribution type")
      }
    }
    
    # Generate counts for each sample based on its CST
    for (i in 1:n_samples) {
      cst <- sample_csts[i]
      
      # Set abundance based on CST
      if (cst == "CST I") {
        # L. crispatus dominated
        p_nonzero <- list(
          L_crispatus_specific = 0.9,
          L_iners_specific = 0.05,
          L_gasseri_specific = 0.05,
          L_jensenii_specific = 0.05,
          G_vaginalis_specific = 0.01,
          Anaerobe_specific = 0.01,
          Housekeeping = 0.8,
          Mobile_elements = 0.3,
          Unknown = 0.1,
          CRISPR = 0.2
        )
        
        abundance_factor <- list(
          L_crispatus_specific = 10,
          L_iners_specific = 0.1,
          L_gasseri_specific = 0.1,
          L_jensenii_specific = 0.1,
          G_vaginalis_specific = 0.01,
          Anaerobe_specific = 0.01,
          Housekeeping = 1,
          Mobile_elements = 1,
          Unknown = 0.5,
          CRISPR = 0.5
        )
      } else if (cst == "CST II") {
        # L. gasseri dominated
        p_nonzero <- list(
          L_crispatus_specific = 0.05,
          L_iners_specific = 0.05,
          L_gasseri_specific = 0.9,
          L_jensenii_specific = 0.05,
          G_vaginalis_specific = 0.01,
          Anaerobe_specific = 0.01,
          Housekeeping = 0.8,
          Mobile_elements = 0.3,
          Unknown = 0.1,
          CRISPR = 0.2
        )
        
        abundance_factor <- list(
          L_crispatus_specific = 0.1,
          L_iners_specific = 0.1,
          L_gasseri_specific = 10,
          L_jensenii_specific = 0.1,
          G_vaginalis_specific = 0.01,
          Anaerobe_specific = 0.01,
          Housekeeping = 1,
          Mobile_elements = 1,
          Unknown = 0.5,
          CRISPR = 0.5
        )
      } else if (cst == "CST III") {
        # L. iners dominated
        p_nonzero <- list(
          L_crispatus_specific = 0.05,
          L_iners_specific = 0.9,
          L_gasseri_specific = 0.05,
          L_jensenii_specific = 0.05,
          G_vaginalis_specific = 0.05,
          Anaerobe_specific = 0.05,
          Housekeeping = 0.8,
          Mobile_elements = 0.3,
          Unknown = 0.1,
          CRISPR = 0.2
        )
        
        abundance_factor <- list(
          L_crispatus_specific = 0.1,
          L_iners_specific = 10,
          L_gasseri_specific = 0.1,
          L_jensenii_specific = 0.1,
          G_vaginalis_specific = 0.1,
          Anaerobe_specific = 0.1,
          Housekeeping = 1,
          Mobile_elements = 1,
          Unknown = 0.5,
          CRISPR = 0.5
        )
      } else if (cst == "CST V") {
        # L. jensenii dominated
        p_nonzero <- list(
          L_crispatus_specific = 0.05,
          L_iners_specific = 0.05,
          L_gasseri_specific = 0.05,
          L_jensenii_specific = 0.9,
          G_vaginalis_specific = 0.01,
          Anaerobe_specific = 0.01,
          Housekeeping = 0.8,
          Mobile_elements = 0.3,
          Unknown = 0.1,
          CRISPR = 0.2
        )
        
        abundance_factor <- list(
          L_crispatus_specific = 0.1,
          L_iners_specific = 0.1,
          L_gasseri_specific = 0.1,
          L_jensenii_specific = 10,
          G_vaginalis_specific = 0.01,
          Anaerobe_specific = 0.01,
          Housekeeping = 1,
          Mobile_elements = 1,
          Unknown = 0.5,
          CRISPR = 0.5
        )
      } else if (cst == "CST IV-A") {
        # G. vaginalis dominated
        p_nonzero <- list(
          L_crispatus_specific = 0.01,
          L_iners_specific = 0.1,
          L_gasseri_specific = 0.01,
          L_jensenii_specific = 0.01,
          G_vaginalis_specific = 0.9,
          Anaerobe_specific = 0.3,
          Housekeeping = 0.8,
          Mobile_elements = 0.4,
          Unknown = 0.2,
          CRISPR = 0.3
        )
        
        abundance_factor <- list(
          L_crispatus_specific = 0.01,
          L_iners_specific = 0.5,
          L_gasseri_specific = 0.01,
          L_jensenii_specific = 0.01,
          G_vaginalis_specific = 10,
          Anaerobe_specific = 1,
          Housekeeping = 1,
          Mobile_elements = 1,
          Unknown = 0.5,
          CRISPR = 0.5
        )
      } else if (cst == "CST IV-B") {
        # Diverse anaerobes
        p_nonzero <- list(
          L_crispatus_specific = 0.01,
          L_iners_specific = 0.05,
          L_gasseri_specific = 0.01,
          L_jensenii_specific = 0.01,
          G_vaginalis_specific = 0.5,
          Anaerobe_specific = 0.8,
          Housekeeping = 0.8,
          Mobile_elements = 0.5,
          Unknown = 0.3,
          CRISPR = 0.4
        )
        
        abundance_factor <- list(
          L_crispatus_specific = 0.01,
          L_iners_specific = 0.1,
          L_gasseri_specific = 0.01,
          L_jensenii_specific = 0.01,
          G_vaginalis_specific = 2,
          Anaerobe_specific = 5,
          Housekeeping = 1,
          Mobile_elements = 1,
          Unknown = 0.5,
          CRISPR = 0.5
        )
      } else if (cst == "CST IV-C") {
        # Mixed Lactobacillus and anaerobes
        p_nonzero <- list(
          L_crispatus_specific = 0.2,
          L_iners_specific = 0.3,
          L_gasseri_specific = 0.1,
          L_jensenii_specific = 0.1,
          G_vaginalis_specific = 0.3,
          Anaerobe_specific = 0.4,
          Housekeeping = 0.8,
          Mobile_elements = 0.4,
          Unknown = 0.2,
          CRISPR = 0.3
        )
        
        abundance_factor <- list(
          L_crispatus_specific = 1,
          L_iners_specific = 2,
          L_gasseri_specific = 0.5,
          L_jensenii_specific = 0.5,
          G_vaginalis_specific = 2,
          Anaerobe_specific = 3,
          Housekeeping = 1,
          Mobile_elements = 1,
          Unknown = 0.5,
          CRISPR = 0.5
        )
      }
      
      # Generate counts for each function group
      for (group in names(p_nonzero)) {
        # Get genes in this function group
        genes_in_group <- which(gene_info$function_group == group)
        
        if (length(genes_in_group) > 0) {
          # Determine which genes have non-zero counts
          n_nonzero <- rbinom(1, length(genes_in_group), p_nonzero[[group]])
          if (n_nonzero > 0) {
            nonzero_indices <- sample(genes_in_group, n_nonzero)
            
            # Generate counts for non-zero genes
            raw_counts <- get_nonzero_counts(n_nonzero, abundance_distribution, abundance_params)
            
            # Apply abundance factor for this function group
            counts_matrix[nonzero_indices, i] <- round(raw_counts * abundance_factor[[group]])
          }
        }
      }
    }
  } else {
    # Random simulation
    # Create count matrix
    counts_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
    
    # For each gene-sample combination, decide if it's zero or non-zero
    for (i in 1:n_genes) {
      for (j in 1:n_samples) {
        if (runif(1) > include_zeros) {
          # Non-zero count
          if (abundance_distribution == "lognormal") {
            counts_matrix[i, j] <- round(exp(rnorm(1, mean = abundance_params[1], sd = abundance_params[2])))
          } else if (abundance_distribution == "zipf") {
            # Simple approximation for Zipf
            counts_matrix[i, j] <- round(10 * runif(1)^(-abundance_params[1]))
          }
        }
      }
    }
    
    # Assign random CSTs
    sample_csts <- sample(
      names(cst_distribution),
      n_samples,
      replace = TRUE,
      prob = as.numeric(cst_distribution)
    )
  }
  
  # Create metadata
  metadata <- data.frame(
    Sample = sample_ids,
    CST = sample_csts,
    Total_reads = colSums(counts_matrix),
    stringsAsFactors = FALSE
  )
  
  # Set row and column names on counts matrix
  rownames(counts_matrix) <- gene_ids
  colnames(counts_matrix) <- sample_ids
  
  # For VIRGO test.out format, create long format data frame
  test_out <- data.frame()
  for (j in 1:n_samples) {
    nonzero_indices <- which(counts_matrix[, j] > 0)
    if (length(nonzero_indices) > 0) {
      sample_data <- data.frame(
        gene_id = gene_ids[nonzero_indices],
        read_count = counts_matrix[nonzero_indices, j],
        gene_length = gene_lengths[nonzero_indices],
        Sample = sample_ids[j],
        stringsAsFactors = FALSE
      )
      test_out <- rbind(test_out, sample_data)
    } else {
      # Add at least one row for this sample with zero count
      sample_data <- data.frame(
        gene_id = gene_ids[1],
        read_count = 0,
        gene_length = gene_lengths[1],
        Sample = sample_ids[j],
        stringsAsFactors = FALSE
      )
      test_out <- rbind(test_out, sample_data)
    }
  }
  
  # Create results list
  results <- list(
    counts = counts_matrix,
    genes = gene_info,
    metadata = metadata,
    test_out = test_out
  )
  
  return(results)
}

#' Write mock VIRGO results to files
#'
#' Generates mock VIRGO data and writes it to files in the format
#' expected by the microFGT package.
#'
#' @param dir_path Directory path to write output files
#' @param base_name Base name for output files
#' @param ... Additional parameters passed to `generate_mock_virgo()`
#'
#' @return Path to the created directory (invisibly)
#' @export
#'
#' @examples
#' \donttest{
#' # Generate and save default mock VIRGO data
#' temp_dir <- tempdir()
#' output_dir <- write_mock_virgo(temp_dir, "mock_virgo", n_genes = 500, n_samples = 5)
#' }
write_mock_virgo <- function(dir_path, base_name, ...) {
  # Make directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Generate mock data
  virgo_data <- generate_mock_virgo(...)
  
  # Write test.out file (primary result)
  test_out_path <- file.path(dir_path, paste0(base_name, "_test.out"))
  write.table(
    virgo_data$test_out,
    file = test_out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Write genes metadata
  genes_path <- file.path(dir_path, paste0(base_name, "_genes.txt"))
  write.table(
    virgo_data$genes,
    file = genes_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Write sample metadata
  metadata_path <- file.path(dir_path, paste0(base_name, "_metadata.txt"))
  write.table(
    virgo_data$metadata,
    file = metadata_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Write count matrix as TSV
  counts_path <- file.path(dir_path, paste0(base_name, "_counts.tsv"))
  write.table(
    virgo_data$counts,
    file = counts_path,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
  )
  
  return(invisible(dir_path))
}

#' Convert mock VIRGO results to FGTMicrobiome format
#'
#' Converts the output from `generate_mock_virgo()` to a format compatible
#' with the `add_virgo()` function.
#'
#' @param virgo_results List from `generate_mock_virgo()`
#'
#' @return A list with components suitable for `add_virgo()`
#' @export
#'
#' @examples
#' \donttest{
#' # Generate mock data
#' virgo_data <- generate_mock_virgo(n_genes = 500, n_samples = 5)
#'
#' # Convert to format for add_virgo()
#' fgt_compatible <- convert_virgo_to_fgt(virgo_data)
#' }
convert_virgo_to_fgt <- function(virgo_results) {
  # Create a basic list structure as expected by add_virgo()
  result <- list(
    genes = virgo_results$counts,  # Gene abundance matrix
    version = "mock"
  )
  
  return(result)
}

#' Generate mock VALENCIA results
#'
#' Creates simulated VALENCIA community state type (CST) classifications
#' for testing purposes.
#'
#' @param n_samples Number of samples to generate
#' @param sample_id_prefix Prefix for sample IDs
#' @param cst_distribution Named vector of CST probabilities
#' @param include_scores Logical; whether to include CST probability scores
#' @param add_taxonomic_data Logical; whether to add taxonomic abundance data
#' @param seed Random seed for reproducibility
#'
#' @return A list with components:
#'   - `cst`: Data frame with CST assignments
#'   - `scores`: Matrix of CST probability scores (if include_scores=TRUE)
#'   - `abundance`: Data frame with taxonomic abundances (if add_taxonomic_data=TRUE)
#' @export
#'
#' @examples
#' \donttest{
#' # Generate default mock VALENCIA data
#' valencia_data <- generate_mock_valencia(n_samples = 20)
#'
#' # Generate data with specific CST distribution
#' cst_dist <- c("I" = 0.6, "III" = 0.3, "IV" = 0.1)
#' valencia_data <- generate_mock_valencia(
#'   n_samples = 10,
#'   cst_distribution = cst_dist,
#'   include_scores = TRUE,
#'   add_taxonomic_data = TRUE
#' )
#' }
generate_mock_valencia <- function(
  n_samples = 20,
  sample_id_prefix = "Sample",
  cst_distribution = NULL,
  include_scores = TRUE,
  add_taxonomic_data = TRUE,
  seed = NULL
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Default CST distribution
  if (is.null(cst_distribution)) {
    cst_distribution <- c(
      "I" = 0.35,    # L. crispatus dominated
      "II" = 0.10,   # L. gasseri dominated
      "III" = 0.25,  # L. iners dominated
      "V" = 0.10,    # L. jensenii dominated
      "IV-A" = 0.10, # G. vaginalis dominated
      "IV-B" = 0.05, # Diverse anaerobes
      "IV-C" = 0.05  # Mixed Lactobacillus and anaerobes
    )
  }
  
  # Normalize CST distribution
  cst_distribution <- cst_distribution / sum(cst_distribution)
  
  # Generate sample IDs
  sample_ids <- paste0(sample_id_prefix, "_", sprintf("%03d", 1:n_samples))
  
  # Generate CST assignments
  sample_csts <- sample(
    names(cst_distribution),
    n_samples,
    replace = TRUE,
    prob = as.numeric(cst_distribution)
  )
  
  # Create CST data frame
  cst_df <- data.frame(
    Sample = sample_ids,
    CST = sample_csts,
    stringsAsFactors = FALSE,
    row.names = sample_ids
  )
  
  # Create CST scores if requested
  if (include_scores) {
    # VALENCIA generates scores for 7 CSTs: I, II, III, IV-A, IV-B, IV-C, V
    all_csts <- c("I", "II", "III", "IV-A", "IV-B", "IV-C", "V")
    score_matrix <- matrix(0, nrow = n_samples, ncol = length(all_csts))
    colnames(score_matrix) <- paste0("CST-", all_csts)
    rownames(score_matrix) <- sample_ids
    
    # Generate scores for each sample
    for (i in 1:n_samples) {
      cst <- sample_csts[i]
      cst_col <- paste0("CST-", cst)
      
      # Generate main score
      main_score <- runif(1, min = 0.7, max = 0.99)
      
      # Distribute remaining probability
      remaining <- 1 - main_score
      other_cols <- setdiff(colnames(score_matrix), cst_col)
      other_scores <- tryCatch({
        rdirichlet(1, rep(1, length(other_cols)))
      }, error = function(e) {
        # Manual generation if rdirichlet not available
        raw_scores <- runif(length(other_cols))
        raw_scores / sum(raw_scores)
      })
      other_scores <- other_scores * remaining
      
      # Assign scores
      score_matrix[i, cst_col] <- main_score
      score_matrix[i, other_cols] <- other_scores
    }
  } else {
    score_matrix <- NULL
  }
  
  # Create taxonomic abundance data if requested
  if (add_taxonomic_data) {
    # Define common vaginal species
    common_taxa <- c(
      "Lactobacillus_crispatus",
      "Lactobacillus_iners",
      "Lactobacillus_gasseri",
      "Lactobacillus_jensenii",
      "Gardnerella_vaginalis",
      "Atopobium_vaginae",
      "Prevotella_bivia",
      "Sneathia_amnii",
      "Megasphaera_genomosp_type_1",
      "Mobiluncus_curtisii",
      "Other"
    )
    
    # Create abundance matrix
    abundance_matrix <- matrix(0, nrow = n_samples, ncol = length(common_taxa))
    colnames(abundance_matrix) <- common_taxa
    rownames(abundance_matrix) <- sample_ids
    
    # Generate abundances based on CST
    for (i in 1:n_samples) {
      cst <- sample_csts[i]
      
      # Create abundance profile based on CST
      if (cst == "I") {
        # L. crispatus dominated
        abundances <- c(
          runif(1, 0.8, 0.98),  # L. crispatus
          runif(1, 0, 0.05),    # L. iners
          runif(1, 0, 0.01),    # L. gasseri
          runif(1, 0, 0.01),    # L. jensenii
          runif(1, 0, 0.01),    # G. vaginalis
          runif(1, 0, 0.01),    # A. vaginae
          runif(1, 0, 0.01),    # P. bivia
          runif(1, 0, 0.01),    # S. amnii
          runif(1, 0, 0.01),    # Megasphaera
          runif(1, 0, 0.01),    # Mobiluncus
          0                     # Other (to be filled)
        )
      } else if (cst == "II") {
        # L. gasseri dominated
        abundances <- c(
          runif(1, 0, 0.05),    # L. crispatus
          runif(1, 0, 0.05),    # L. iners
          runif(1, 0.8, 0.98),  # L. gasseri
          runif(1, 0, 0.01),    # L. jensenii
          runif(1, 0, 0.01),    # G. vaginalis
          runif(1, 0, 0.01),    # A. vaginae
          runif(1, 0, 0.01),    # P. bivia
          runif(1, 0, 0.01),    # S. amnii
          runif(1, 0, 0.01),    # Megasphaera
          runif(1, 0, 0.01),    # Mobiluncus
          0                     # Other (to be filled)
        )
      } else if (cst == "III") {
        # L. iners dominated
        abundances <- c(
          runif(1, 0, 0.05),    # L. crispatus
          runif(1, 0.8, 0.98),  # L. iners
          runif(1, 0, 0.01),    # L. gasseri
          runif(1, 0, 0.01),    # L. jensenii
          runif(1, 0, 0.05),    # G. vaginalis
          runif(1, 0, 0.02),    # A. vaginae
          runif(1, 0, 0.02),    # P. bivia
          runif(1, 0, 0.01),    # S. amnii
          runif(1, 0, 0.01),    # Megasphaera
          runif(1, 0, 0.01),    # Mobiluncus
          0                     # Other (to be filled)
        )
      } else if (cst == "V") {
        # L. jensenii dominated
        abundances <- c(
          runif(1, 0, 0.05),    # L. crispatus
          runif(1, 0, 0.05),    # L. iners
          runif(1, 0, 0.01),    # L. gasseri
          runif(1, 0.8, 0.98),  # L. jensenii
          runif(1, 0, 0.01),    # G. vaginalis
          runif(1, 0, 0.01),    # A. vaginae
          runif(1, 0, 0.01),    # P. bivia
          runif(1, 0, 0.01),    # S. amnii
          runif(1, 0, 0.01),    # Megasphaera
          runif(1, 0, 0.01),    # Mobiluncus
          0                     # Other (to be filled)
        )
      } else if (cst == "IV-A") {
        # G. vaginalis dominated
        abundances <- c(
          runif(1, 0, 0.01),    # L. crispatus
          runif(1, 0, 0.1),     # L. iners
          runif(1, 0, 0.01),    # L. gasseri
          runif(1, 0, 0.01),    # L. jensenii
          runif(1, 0.5, 0.8),   # G. vaginalis
          runif(1, 0.05, 0.2),  # A. vaginae
          runif(1, 0.01, 0.1),  # P. bivia
          runif(1, 0, 0.05),    # S. amnii
          runif(1, 0, 0.05),    # Megasphaera
          runif(1, 0, 0.05),    # Mobiluncus
          0                     # Other (to be filled)
        )
      } else if (cst == "IV-B") {
        # Diverse anaerobes
        abundances <- c(
          runif(1, 0, 0.01),    # L. crispatus
          runif(1, 0, 0.05),    # L. iners
          runif(1, 0, 0.01),    # L. gasseri
          runif(1, 0, 0.01),    # L. jensenii
          runif(1, 0.1, 0.3),   # G. vaginalis
          runif(1, 0.1, 0.3),   # A. vaginae
          runif(1, 0.05, 0.2),  # P. bivia
          runif(1, 0.05, 0.2),  # S. amnii
          runif(1, 0.05, 0.2),  # Megasphaera
          runif(1, 0.01, 0.1),  # Mobiluncus
          0                     # Other (to be filled)
        )
      } else if (cst == "IV-C") {
        # Mixed Lactobacillus and anaerobes
        abundances <- c(
          runif(1, 0.05, 0.2),  # L. crispatus
          runif(1, 0.1, 0.3),   # L. iners
          runif(1, 0, 0.1),     # L. gasseri
          runif(1, 0, 0.1),     # L. jensenii
          runif(1, 0.1, 0.3),   # G. vaginalis
          runif(1, 0.05, 0.15), # A. vaginae
          runif(1, 0.01, 0.1),  # P. bivia
          runif(1, 0.01, 0.1),  # S. amnii
          runif(1, 0.01, 0.1),  # Megasphaera
          runif(1, 0, 0.05),    # Mobiluncus
          0                     # Other (to be filled)
        )
      }
      
      # Normalize to make sure they sum to 1, accounting for "Other"
      abundances[length(abundances)] <- max(0, 1 - sum(abundances[1:(length(abundances)-1)]))
      if (sum(abundances) < 1) {
        # In case of rounding errors
        abundances <- abundances / sum(abundances)
      }
      
      # Assign to matrix
      abundance_matrix[i, ] <- abundances
    }
    
    # Convert to data frame for output
    abundance_df <- as.data.frame(abundance_matrix)
    abundance_df$Sample <- sample_ids
    abundance_df$read_count <- sample(10000:1000000, n_samples, replace = TRUE)
  } else {
    abundance_df <- NULL
  }
  
  # Create results list
  results <- list(
    cst = cst_df
  )
  
  if (include_scores) {
    results$scores <- score_matrix
  }
  
  if (add_taxonomic_data) {
    results$abundance <- abundance_df
  }
  
  return(results)
}

#' Write mock VALENCIA results to files
#'
#' Generates mock VALENCIA data and writes it to files in the format
#' expected by the microFGT package.
#'
#' @param dir_path Directory path to write output files
#' @param base_name Base name for output files
#' @param ... Additional parameters passed to `generate_mock_valencia()`
#'
#' @return Path to the created directory (invisibly)
#' @export
#'
#' @examples
#' \donttest{
#' # Generate and save default mock VALENCIA data
#' temp_dir <- tempdir()
#' output_dir <- write_mock_valencia(temp_dir, "mock_valencia", n_samples = 15)
#' }
write_mock_valencia <- function(dir_path, base_name, ...) {
  # Make directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Generate mock data
  valencia_data <- generate_mock_valencia(...)
  
  # Write CST assignments
  cst_path <- file.path(dir_path, paste0(base_name, "_cst.csv"))
  write.csv(
    valencia_data$cst,
    file = cst_path,
    quote = FALSE,
    row.names = FALSE
  )
  
  # Write scores if included
  if (!is.null(valencia_data$scores)) {
    scores_path <- file.path(dir_path, paste0(base_name, "_scores.csv"))
    write.csv(
      valencia_data$scores,
      file = scores_path,
      quote = FALSE
    )
  }
  
  # Write taxonomic data if included
  if (!is.null(valencia_data$abundance)) {
    abundance_path <- file.path(dir_path, paste0(base_name, "_abundance.csv"))
    write.csv(
      valencia_data$abundance,
      file = abundance_path,
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  return(invisible(dir_path))
}

#' Convert mock VALENCIA results to FGTMicrobiome format
#'
#' Converts the output from `generate_mock_valencia()` to a format compatible
#' with the `add_valencia()` function.
#'
#' @param valencia_results List from `generate_mock_valencia()`
#'
#' @return A list with components suitable for `add_valencia()`
#' @export
#'
#' @examples
#' \donttest{
#' # Generate mock data
#' valencia_data <- generate_mock_valencia(n_samples = 20, include_scores = TRUE)
#'
#' # Convert to format for add_valencia()
#' fgt_compatible <- convert_valencia_to_fgt(valencia_data)
#' }
convert_valencia_to_fgt <- function(valencia_results) {
  # Create a basic list structure as expected by add_valencia()
  result <- list(
    cst = valencia_results$cst,
    version = "mock"
  )
  
  # Add scores if available
  if (!is.null(valencia_results$scores)) {
    result$scores <- valencia_results$scores
  }
  
  return(result)
}

#' Generate a complete set of mock data for all tools
#'
#' Creates a coordinated set of mock data for SpeciateIT, VIRGO, and VALENCIA
#' with consistent samples and community state types.
#'
#' @param n_samples Number of samples to generate
#' @param n_sequences Number of sequences for SpeciateIT data
#' @param n_genes Number of genes for VIRGO data
#' @param sample_id_prefix Prefix for sample IDs
#' @param cst_distribution Named vector of CST probabilities
#' @param seed Random seed for reproducibility
#' @param create_files Logical; whether to write files to disk
#' @param output_dir Directory to write files to (if create_files=TRUE)
#'
#' @return A list with components for each tool's data
#' @export
#'
#' @examples
#' \donttest{
#' # Generate a complete set of mock data
#' mock_data <- generate_mock_fgt_dataset(
#'   n_samples = 10,
#'   n_sequences = 100,
#'   n_genes = 500
#' )
#'
#' # Generate files on disk
#' temp_dir <- tempdir()
#' mock_data_with_files <- generate_mock_fgt_dataset(
#'   n_samples = 5,
#'   create_files = TRUE,
#'   output_dir = temp_dir
#' )
#' }
generate_mock_fgt_dataset <- function(
  n_samples = 10,
  n_sequences = 1000,
  n_genes = 1000,
  sample_id_prefix = "Sample",
  cst_distribution = NULL,
  seed = NULL,
  create_files = FALSE,
  output_dir = NULL
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Default CST distribution
  if (is.null(cst_distribution)) {
    cst_distribution <- c(
      "I" = 0.35,    # L. crispatus dominated
      "II" = 0.10,   # L. gasseri dominated
      "III" = 0.25,  # L. iners dominated
      "V" = 0.10,    # L. jensenii dominated
      "IV-A" = 0.10, # G. vaginalis dominated
      "IV-B" = 0.05, # Diverse anaerobes
      "IV-C" = 0.05  # Mixed Lactobacillus and anaerobes
    )
  }
  
  # Generate sample IDs
  sample_ids <- paste0(sample_id_prefix, "_", sprintf("%03d", 1:n_samples))
  
  # Generate CST assignments (to ensure consistency across tools)
  sample_csts <- sample(
    names(cst_distribution),
    n_samples,
    replace = TRUE,
    prob = as.numeric(cst_distribution)
  )
  
  # Generate VALENCIA data first (as it's the CST classifier)
  message("Generating mock VALENCIA data...")
  valencia_data <- generate_mock_valencia(
    n_samples = n_samples,
    sample_id_prefix = sample_id_prefix,
    cst_distribution = NULL,  # Use default distribution initially
    include_scores = TRUE,
    add_taxonomic_data = TRUE,
    seed = if (is.null(seed)) NULL else seed + 1
  )
  
  # Force the CST assignments to match our pre-generated ones
  valencia_data$cst$CST <- sample_csts
  
  # Generate VIRGO data with matching CSTs
  message("Generating mock VIRGO data...")
  virgo_cst_dist <- setNames(
    as.numeric(table(sample_csts)) / n_samples,
    unique(sample_csts)
  )

  virgo_data <- generate_mock_virgo(
    n_genes = n_genes,
    n_samples = n_samples,
    sample_id_prefix = sample_id_prefix,
    simulation_type = "random",  # Use random simulation for simplicity in tests
    seed = if (is.null(seed)) NULL else seed + 2
  )
  
  # Generate speciateIT data with taxonomic distribution matching CSTs
  message("Generating mock SpeciateIT data...")
  
  # Function to create species distribution based on CST
  create_species_dist <- function(cst) {
    if (cst == "I") {
      # L. crispatus dominated
      return(c(
        "Lactobacillus crispatus" = 0.8,
        "Lactobacillus iners" = 0.05,
        "Lactobacillus gasseri" = 0.02,
        "Lactobacillus jensenii" = 0.02,
        "Gardnerella vaginalis" = 0.02,
        "Atopobium vaginae" = 0.02,
        "Prevotella bivia" = 0.02,
        "Other" = 0.05
      ))
    } else if (cst == "II") {
      # L. gasseri dominated
      return(c(
        "Lactobacillus crispatus" = 0.05,
        "Lactobacillus iners" = 0.05,
        "Lactobacillus gasseri" = 0.8,
        "Lactobacillus jensenii" = 0.02,
        "Gardnerella vaginalis" = 0.02,
        "Atopobium vaginae" = 0.01,
        "Prevotella bivia" = 0.01,
        "Other" = 0.04
      ))
    } else if (cst == "III") {
      # L. iners dominated
      return(c(
        "Lactobacillus crispatus" = 0.05,
        "Lactobacillus iners" = 0.8,
        "Lactobacillus gasseri" = 0.02,
        "Lactobacillus jensenii" = 0.02,
        "Gardnerella vaginalis" = 0.05,
        "Atopobium vaginae" = 0.01,
        "Prevotella bivia" = 0.01,
        "Other" = 0.04
      ))
    } else if (cst == "V") {
      # L. jensenii dominated
      return(c(
        "Lactobacillus crispatus" = 0.05,
        "Lactobacillus iners" = 0.05,
        "Lactobacillus gasseri" = 0.02,
        "Lactobacillus jensenii" = 0.8,
        "Gardnerella vaginalis" = 0.02,
        "Atopobium vaginae" = 0.01,
        "Prevotella bivia" = 0.01,
        "Other" = 0.04
      ))
    } else if (cst == "IV-A") {
      # G. vaginalis dominated
      return(c(
        "Lactobacillus crispatus" = 0.01,
        "Lactobacillus iners" = 0.1,
        "Lactobacillus gasseri" = 0.01,
        "Lactobacillus jensenii" = 0.01,
        "Gardnerella vaginalis" = 0.6,
        "Gardnerella leopoldii" = 0.05,
        "Gardnerella swidsinskii" = 0.05,
        "Gardnerella piotii" = 0.02,
        "Atopobium vaginae" = 0.05,
        "Prevotella bivia" = 0.05,
        "Other" = 0.05
      ))
    } else if (cst == "IV-B") {
      # Diverse anaerobes
      return(c(
        "Lactobacillus crispatus" = 0.01,
        "Lactobacillus iners" = 0.05,
        "Lactobacillus gasseri" = 0.01,
        "Lactobacillus jensenii" = 0.01,
        "Gardnerella vaginalis" = 0.2,
        "Atopobium vaginae" = 0.2,
        "Prevotella bivia" = 0.15,
        "Sneathia amnii" = 0.15,
        "Megasphaera genomosp type 1" = 0.1,
        "Mobiluncus curtisii" = 0.05,
        "Other" = 0.07
      ))
    } else if (cst == "IV-C") {
      # Mixed Lactobacillus and anaerobes
      return(c(
        "Lactobacillus crispatus" = 0.15,
        "Lactobacillus iners" = 0.2,
        "Lactobacillus gasseri" = 0.05,
        "Lactobacillus jensenii" = 0.05,
        "Gardnerella vaginalis" = 0.15,
        "Atopobium vaginae" = 0.1,
        "Prevotella bivia" = 0.1,
        "Sneathia amnii" = 0.05,
        "Megasphaera genomosp type 1" = 0.05,
        "Mobiluncus curtisii" = 0.05,
        "Other" = 0.05
      ))
    } else {
      # Default balanced distribution
      return(c(
        "Lactobacillus crispatus" = 0.25,
        "Lactobacillus iners" = 0.25,
        "Lactobacillus gasseri" = 0.1,
        "Lactobacillus jensenii" = 0.1,
        "Gardnerella vaginalis" = 0.1,
        "Atopobium vaginae" = 0.05,
        "Prevotella bivia" = 0.05,
        "Other" = 0.1
      ))
    }
  }
  
  # Generate mock data for each sample
  speciateit_results <- list()
  sequences_per_sample <- ceiling(n_sequences / n_samples)
  
  for (i in 1:n_samples) {
    # Get CST for this sample
    cst <- sample_csts[i]
    
    # Create species distribution
    species_dist <- create_species_dist(cst)
    
    # Generate mock data for this sample
    sample_data <- generate_mock_speciateit(
      n_sequences = sequences_per_sample,
      seq_id_prefix = paste0(sample_id_prefix, "_", sprintf("%03d", i), "_Seq"),
      species_distribution = species_dist,
      include_unclassified = 0.02,
      seed = if (is.null(seed)) NULL else seed + i + 100
    )
    
    # Add to results
    speciateit_results[[i]] <- sample_data
  }
  
  # Combine all samples
  combined_speciateit <- do.call(rbind, speciateit_results)
  
  # Generate FGT-compatible formats
  message("Creating FGT-compatible formats")

  # SpeciateIT format
  fgt_speciateit <- tryCatch({
    convert_speciateit_to_fgt(combined_speciateit)
  }, error = function(e) {
    message("Warning: Could not convert SpeciateIT data: ", e$message)
    list(
      assignments = data.frame(
        Feature = combined_speciateit$`Sequence ID`,
        Species = combined_speciateit$Classification,
        stringsAsFactors = FALSE
      ),
      version = "mock"
    )
  })

  # VIRGO format
  fgt_virgo <- tryCatch({
    convert_virgo_to_fgt(virgo_data)
  }, error = function(e) {
    message("Warning: Could not convert VIRGO data: ", e$message)
    list(
      genes = virgo_data$counts,
      version = "mock"
    )
  })

  # VALENCIA format
  fgt_valencia <- tryCatch({
    convert_valencia_to_fgt(valencia_data)
  }, error = function(e) {
    message("Warning: Could not convert VALENCIA data: ", e$message)
    list(
      cst = valencia_data$cst,
      version = "mock"
    )
  })
  
  # Write files if requested
  if (create_files) {
    if (is.null(output_dir)) {
      output_dir <- tempdir()
      message("No output directory specified, using temporary directory: ", output_dir)
    }
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    message("Writing output files to: ", output_dir)
    
    # SpeciateIT files
    speciateit_dir <- file.path(output_dir, "speciateit")
    if (!dir.exists(speciateit_dir)) {
      dir.create(speciateit_dir)
    }
    
    speciateit_out <- file.path(speciateit_dir, "MC_order7_results.txt")
    write.table(
      combined_speciateit,
      file = speciateit_out,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    # VIRGO files
    virgo_dir <- file.path(output_dir, "virgo")
    write_mock_virgo(virgo_dir, "mock_virgo", 
                    n_genes = n_genes, 
                    n_samples = n_samples,
                    sample_id_prefix = sample_id_prefix,
                    simulation_type = "realistic",
                    cst_distribution = virgo_cst_dist,
                    seed = if (is.null(seed)) NULL else seed + 2)
    
    # VALENCIA files
    valencia_dir <- file.path(output_dir, "valencia")
    write_mock_valencia(valencia_dir, "mock_valencia",
                        n_samples = n_samples,
                        sample_id_prefix = sample_id_prefix,
                        include_scores = TRUE,
                        add_taxonomic_data = TRUE,
                        seed = if (is.null(seed)) NULL else seed + 1)
    
    # Override with our forced CST assignments
    cst_path <- file.path(valencia_dir, "mock_valencia_cst.csv")
    write.csv(
      data.frame(
        Sample = sample_ids,
        CST = sample_csts
      ),
      file = cst_path,
      quote = FALSE,
      row.names = FALSE
    )
    
    message("Files written successfully.")
  }
  
  # Return a comprehensive results object
  results <- list(
    speciateit = list(
      raw = combined_speciateit,
      fgt_compatible = fgt_speciateit
    ),
    virgo = list(
      raw = virgo_data,
      fgt_compatible = fgt_virgo
    ),
    valencia = list(
      raw = valencia_data,
      fgt_compatible = fgt_valencia
    ),
    metadata = data.frame(
      Sample = sample_ids,
      CST = sample_csts
    ),
    parameters = list(
      n_samples = n_samples,
      n_sequences = n_sequences,
      n_genes = n_genes,
      sample_id_prefix = sample_id_prefix,
      seed = seed
    )
  )
  
  if (create_files) {
    results$file_paths <- list(
      base_dir = output_dir,
      speciateit = speciateit_dir,
      virgo = virgo_dir,
      valencia = valencia_dir
    )
  }
  
  return(results)
}