# Direct example data loading script
# This can be copy-pasted into an R console

# Define a standalone load_example_data function that doesn't depend on exports
load_example_data_direct <- function(
  size = "small",
  type = "amplicon",
  as_fgt_experiment = TRUE
) {
  # Validate inputs
  size <- match.arg(size, c("small", "medium", "large"))
  type <- match.arg(type, c("amplicon", "metagenomics"))
  
  # Library checks
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required. Please install it with BiocManager::install('SummarizedExperiment')")
  }
  
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    stop("Package 'TreeSummarizedExperiment' is required. Please install it with BiocManager::install('TreeSummarizedExperiment')")
  }
  
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required. Please install it with BiocManager::install('S4Vectors')")
  }
  
  # Construct the file path
  base_name <- paste0("microFGT_example_", size, "_", type)
  
  # Check if we have pre-built examples in the package
  data_dir <- system.file("extdata", package = "microFGT")
  
  if (file.exists(file.path(data_dir, paste0(base_name, "_counts.rds")))) {
    message("Loading example data from the package...")
    
    # Load pre-built example
    counts <- readRDS(file.path(data_dir, paste0(base_name, "_counts.rds")))
    taxonomy <- readRDS(file.path(data_dir, paste0(base_name, "_taxonomy.rds")))
    metadata <- readRDS(file.path(data_dir, paste0(base_name, "_metadata.rds")))
    
    # Check if tree exists
    tree <- NULL
    tree_file <- file.path(data_dir, paste0(base_name, "_tree.rds"))
    if (file.exists(tree_file)) {
      tree <- readRDS(tree_file)
    }
  } else {
    # Generate example data on the fly
    message("Pre-built example not found. Generating example data on the fly...")
    
    # Set parameters based on size
    n_samples <- switch(size,
                       "small" = 10,
                       "medium" = 30,
                       "large" = 100)
    
    n_features <- switch(size,
                        "small" = 50,
                        "medium" = 150,
                        "large" = 500)
    
    # Create a simple count matrix
    set.seed(42)  # For reproducibility
    counts <- matrix(rpois(n_samples * n_features, 20), 
                    nrow = n_features, ncol = n_samples)
    rownames(counts) <- paste0("Feature", 1:n_features)
    colnames(counts) <- paste0("Sample", 1:n_samples)
    
    # Create simple taxonomy
    taxonomy <- data.frame(
      Kingdom = rep("Bacteria", n_features),
      Phylum = sample(c("Firmicutes", "Bacteroidetes", "Proteobacteria"), 
                     n_features, replace = TRUE),
      Genus = sample(c("Lactobacillus", "Gardnerella", "Prevotella", "Sneathia"), 
                    n_features, replace = TRUE),
      row.names = rownames(counts)
    )
    
    # Create simple metadata
    metadata <- data.frame(
      condition = sample(c("Healthy", "BV"), n_samples, replace = TRUE),
      pH = runif(n_samples, 3.8, 5.5),
      Nugent_Score = sample(0:10, n_samples, replace = TRUE),
      row.names = colnames(counts)
    )
    
    # No tree generated on the fly
    tree <- NULL
  }
  
  # Return as list or FGTExperiment
  if (as_fgt_experiment) {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE) && 
        requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
      
      # Create using the two-step constructor approach
      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowData = taxonomy,
        colData = metadata
      )
      
      tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(se, rowTree = tree)
      
      if (methods::existsClass("FGTExperiment")) {
        # Check if we can access the FGTExperiment constructor
        if (exists("FGTExperiment", mode = "function", envir = as.environment("package:microFGT"))) {
          # Use the constructor
          result <- FGTExperiment(tse, experimentType = type)
          
          # Ensure assay names are set correctly
          SummarizedExperiment::assayNames(result)[1] <- "counts"
          return(result)
        } else {
          # Use two-step method
          message("Using two-step constructor method...")
          result <- methods::new("FGTExperiment", 
                              tse,
                              experimentType = type,
                              fgtMetadata = S4Vectors::SimpleList())
          
          # Ensure assay names are set correctly
          SummarizedExperiment::assayNames(result)[1] <- "counts"
          return(result)
        }
      } else {
        message("FGTExperiment class not available. Returning TreeSummarizedExperiment object.")
        return(tse)
      }
    } else {
      message("Required packages not available. Returning data as list.")
    }
  }
  
  # Return as list
  return(list(
    counts = counts,
    taxonomy = taxonomy,
    metadata = metadata,
    tree = tree
  ))
}

# Usage example:
# Remember to load the microFGT package first
# library(microFGT)
# 
# # Then use this function
# example_data <- load_example_data_direct(size = "small", type = "amplicon")
# 
# # Examine the data
# if (methods::is(example_data, "FGTExperiment")) {
#   print(example_data)
#   print(SummarizedExperiment::colData(example_data)[, c("condition", "pH", "Nugent_Score")])
# } else {
#   str(example_data)
# }