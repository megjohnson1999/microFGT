# Mock data generator functions for microFGT package

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