#' Import microbiome count data from various formats
#'
#' A flexible import function that creates a SummarizedExperiment-like object from
#' different input formats including count matrices, data frames, tables, and files.
#'
#' @param counts Count data in any of the following formats:
#'   - Matrix (features as rows, samples as columns)
#'   - Data frame (features as rows, samples as columns)
#'   - File path to CSV/TSV file
#'   - File path to RDS file containing count data
#' @param taxonomy Taxonomic information in any of the following formats:
#'   - Data frame (features as rows, taxa as columns)
#'   - Matrix (features as rows, taxa as columns)
#'   - File path to CSV/TSV file
#'   - File path to RDS file containing taxonomy data
#'   - Character vector of taxonomic strings (will be parsed)
#' @param sample_data Sample metadata in any of the following formats:
#'   - Data frame (samples as rows, metadata as columns)
#'   - File path to CSV/TSV file
#'   - File path to RDS file containing sample metadata
#' @param tree Phylogenetic tree (optional) in any of the following formats:
#'   - phylo object from ape package
#'   - File path to Newick (.nwk, .tre) or RDS file containing tree
#' @param tax_sep Separator for taxonomic strings if taxonomy is a character vector
#' @param sample_col_name Name of the column in sample_data that matches sample names
#' @param feature_col_name Name of the column in taxonomy that matches feature names
#' @param as_FGTExperiment Logical; if TRUE and FGTExperiment class exists, 
#'   returns an FGTExperiment object
#'
#' @return A SummarizedExperiment or FGTExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Import from matrices
#' counts <- matrix(1:60, nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' 
#' taxonomy <- data.frame(
#'   Kingdom = rep("Bacteria", 10),
#'   Phylum = sample(c("Firmicutes", "Bacteroidetes", "Proteobacteria"), 10, replace = TRUE),
#'   Genus = sample(c("Lactobacillus", "Bacteroides", "Escherichia"), 10, replace = TRUE)
#' )
#' rownames(taxonomy) <- rownames(counts)
#'
#' metadata <- data.frame(
#'   Group = sample(c("Control", "Treatment"), 6, replace = TRUE),
#'   Age = sample(20:50, 6)
#' )
#' rownames(metadata) <- colnames(counts)
#'
#' se <- import_microbiome(counts, taxonomy, metadata)
#' }
import_microbiome <- function(counts, taxonomy = NULL, sample_data = NULL, 
                             tree = NULL, tax_sep = ";", 
                             sample_col_name = NULL, feature_col_name = NULL,
                             as_FGTExperiment = TRUE) {
  
  # Process count data
  counts_matrix <- process_count_data(counts)
  
  # Process taxonomy if provided
  if (!is.null(taxonomy)) {
    taxonomy_df <- process_taxonomy_data(taxonomy, tax_sep, feature_col_name)
    
    # Match taxonomy rows to count matrix rows
    if (nrow(taxonomy_df) > 0) {
      # If taxonomy rownames and count matrix rownames don't match
      if (!identical(rownames(taxonomy_df), rownames(counts_matrix))) {
        common_rows <- intersect(rownames(taxonomy_df), rownames(counts_matrix))
        if (length(common_rows) == 0) {
          warning("No matching rownames between counts and taxonomy, using row indices")
          # Use as many rows as available for both
          n_rows <- min(nrow(counts_matrix), nrow(taxonomy_df))
          counts_matrix <- counts_matrix[1:n_rows, , drop = FALSE]
          taxonomy_df <- taxonomy_df[1:n_rows, , drop = FALSE]
          rownames(taxonomy_df) <- rownames(counts_matrix)
        } else {
          # Subset to matching rows
          counts_matrix <- counts_matrix[common_rows, , drop = FALSE]
          taxonomy_df <- taxonomy_df[common_rows, , drop = FALSE]
        }
      }
    } else {
      # Create empty taxonomy data frame
      taxonomy_df <- data.frame(row.names = rownames(counts_matrix))
    }
  } else {
    # Create empty taxonomy data frame
    taxonomy_df <- data.frame(row.names = rownames(counts_matrix))
  }
  
  # Process sample metadata if provided
  if (!is.null(sample_data)) {
    sample_df <- process_sample_data(sample_data, sample_col_name)
    
    # Match sample metadata rows to count matrix columns
    if (nrow(sample_df) > 0) {
      if (!identical(rownames(sample_df), colnames(counts_matrix))) {
        common_cols <- intersect(rownames(sample_df), colnames(counts_matrix))
        if (length(common_cols) == 0) {
          warning("No matching rownames between counts and sample metadata, using column indices")
          # Use as many columns as available for both
          n_cols <- min(ncol(counts_matrix), nrow(sample_df))
          counts_matrix <- counts_matrix[, 1:n_cols, drop = FALSE]
          sample_df <- sample_df[1:n_cols, , drop = FALSE]
          rownames(sample_df) <- colnames(counts_matrix)
        } else {
          # Subset to matching columns
          counts_matrix <- counts_matrix[, common_cols, drop = FALSE]
          sample_df <- sample_df[common_cols, , drop = FALSE]
        }
      }
    } else {
      # Create empty sample metadata data frame
      sample_df <- data.frame(row.names = colnames(counts_matrix))
    }
  } else {
    # Create empty sample metadata data frame
    sample_df <- data.frame(row.names = colnames(counts_matrix))
  }
  
  # Process phylogenetic tree if provided
  if (!is.null(tree)) {
    phy_tree <- process_tree_data(tree)
    
    # Check if tree tip labels match count matrix rows
    if (!is.null(phy_tree) && !all(rownames(counts_matrix) %in% phy_tree$tip.label)) {
      warning("Tree tip labels do not match all feature names in the count matrix")
    }
  } else {
    phy_tree <- NULL
  }
  
  # Ensure rowData is a DataFrame
  if (requireNamespace("S4Vectors", quietly = TRUE)) {
    row_data <- S4Vectors::DataFrame(taxonomy_df)
    col_data <- S4Vectors::DataFrame(sample_df)
  } else {
    row_data <- taxonomy_df
    col_data <- sample_df
  }
  
  # Create the appropriate object based on the available packages and user preference
  if (as_FGTExperiment && exists("FGTExperiment", mode = "function")) {
    # Use FGTExperiment if available and requested
    if (!is.null(phy_tree) && requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
      # With phylogenetic tree
      result <- FGTExperiment(
        assays = list(counts = counts_matrix),
        rowData = row_data,
        colData = col_data,
        rowTree = phy_tree,
        experimentType = "amplicon"
      )
    } else {
      # Without phylogenetic tree
      result <- FGTExperiment(
        assays = list(counts = counts_matrix),
        rowData = row_data,
        colData = col_data,
        experimentType = "amplicon"
      )
    }
  } else if (!is.null(phy_tree) && requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    # Use TreeSummarizedExperiment if tree is provided and package is available
    result <- TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = row_data,
      colData = col_data,
      rowTree = phy_tree
    )
  } else if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    # Use SummarizedExperiment as fallback
    result <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_matrix),
      rowData = row_data,
      colData = col_data
    )
  } else {
    # If none of the above are available, return a list
    result <- list(
      counts = counts_matrix,
      taxonomy = taxonomy_df,
      sample_data = sample_df,
      tree = phy_tree
    )
    warning("SummarizedExperiment package not available; returning a list instead")
  }
  
  # Add import information to metadata
  if (is(result, "SummarizedExperiment")) {
    SummarizedExperiment::metadata(result)$import_info <- list(
      date = Sys.time(),
      n_features = nrow(counts_matrix),
      n_samples = ncol(counts_matrix)
    )
  }
  
  return(result)
}

#' Process count data input
#'
#' Internal function to handle various count data input formats
#'
#' @param counts Input count data
#' @return A matrix with features as rows and samples as columns
#' @keywords internal
process_count_data <- function(counts) {
  if (is.null(counts)) {
    stop("Count data is required")
  }
  
  # If input is a file path
  if (is.character(counts) && length(counts) == 1) {
    if (file.exists(counts)) {
      if (grepl("\\.rds$", counts, ignore.case = TRUE)) {
        # Load RDS file
        counts_data <- readRDS(counts)
      } else if (grepl("\\.(csv|tsv|txt)$", counts, ignore.case = TRUE)) {
        # Determine separator
        sep <- if (grepl("\\.tsv$", counts, ignore.case = TRUE)) "\t" else ","
        # Read CSV/TSV file
        counts_data <- read.table(counts, header = TRUE, sep = sep, 
                                 row.names = 1, check.names = FALSE)
      } else {
        stop("Unsupported file format for counts. Use CSV, TSV, or RDS.")
      }
    } else {
      stop("File not found: ", counts)
    }
  } else {
    counts_data <- counts
  }
  
  # Convert to matrix if data frame
  if (is.data.frame(counts_data)) {
    counts_matrix <- as.matrix(counts_data)
  } else if (is.matrix(counts_data)) {
    counts_matrix <- counts_data
  } else {
    stop("Counts must be a matrix, data frame, or file path")
  }
  
  # Check if features are in rows and samples in columns
  # This is a simple heuristic - in typical microbiome data, there are more features than samples
  if (ncol(counts_matrix) > nrow(counts_matrix)) {
    # More columns than rows suggests samples might be in columns, which is what we want
    # But let's check if we need to transpose (features in columns, samples in rows)
    
    # If all column names look like feature/OTU identifiers (e.g., "ASV1", "OTU123")
    if (!is.null(colnames(counts_matrix)) && 
        all(grepl("^(Feature|OTU|ASV|Seq|[A-Z0-9]{8})", colnames(counts_matrix)))) {
      message("Count matrix appears to have features as columns. Transposing...")
      counts_matrix <- t(counts_matrix)
    }
  }
  
  # Ensure we have row and column names
  if (is.null(rownames(counts_matrix))) {
    rownames(counts_matrix) <- paste0("Feature", 1:nrow(counts_matrix))
  }
  
  if (is.null(colnames(counts_matrix))) {
    colnames(counts_matrix) <- paste0("Sample", 1:ncol(counts_matrix))
  }
  
  return(counts_matrix)
}

#' Process taxonomy data input
#'
#' Internal function to handle various taxonomy data input formats
#'
#' @param taxonomy Input taxonomy data
#' @param tax_sep Separator for taxonomic strings if taxonomy is a character vector
#' @param feature_col_name Name of the column that contains feature IDs
#' @return A data frame with taxonomic information
#' @keywords internal
process_taxonomy_data <- function(taxonomy, tax_sep = ";", feature_col_name = NULL) {
  if (is.null(taxonomy)) {
    return(data.frame())
  }
  
  # If input is a file path
  if (is.character(taxonomy) && length(taxonomy) == 1) {
    if (file.exists(taxonomy)) {
      if (grepl("\\.rds$", taxonomy, ignore.case = TRUE)) {
        # Load RDS file
        taxonomy_data <- readRDS(taxonomy)
      } else if (grepl("\\.(csv|tsv|txt)$", taxonomy, ignore.case = TRUE)) {
        # Determine separator
        sep <- if (grepl("\\.tsv$", taxonomy, ignore.case = TRUE)) "\t" else ","
        # Read CSV/TSV file
        taxonomy_data <- read.table(taxonomy, header = TRUE, sep = sep, 
                                   row.names = 1, check.names = FALSE)
      } else {
        stop("Unsupported file format for taxonomy. Use CSV, TSV, or RDS.")
      }
    } else {
      stop("File not found: ", taxonomy)
    }
  } else {
    taxonomy_data <- taxonomy
  }
  
  # Handle character vector of taxonomic strings
  if (is.character(taxonomy_data) && length(taxonomy_data) > 0 && 
      !is.data.frame(taxonomy_data) && !is.matrix(taxonomy_data)) {
    # Parse taxonomic strings
    taxonomy_df <- parse_taxonomy_strings(taxonomy_data, tax_sep)
    return(taxonomy_df)
  }
  
  # Convert to data frame if matrix
  if (is.matrix(taxonomy_data)) {
    taxonomy_df <- as.data.frame(taxonomy_data, stringsAsFactors = FALSE)
  } else if (is.data.frame(taxonomy_data)) {
    taxonomy_df <- taxonomy_data
  } else {
    stop("Taxonomy must be a data frame, matrix, character vector, or file path")
  }
  
  # If a feature column name is provided, set it as rownames
  if (!is.null(feature_col_name) && feature_col_name %in% colnames(taxonomy_df)) {
    rownames(taxonomy_df) <- taxonomy_df[[feature_col_name]]
    taxonomy_df <- taxonomy_df[, colnames(taxonomy_df) != feature_col_name, drop = FALSE]
  }
  
  return(taxonomy_df)
}

#' Parse taxonomic strings into a data frame
#'
#' Internal function to parse taxonomic strings like "k__Bacteria;p__Firmicutes;g__Lactobacillus"
#'
#' @param tax_strings Character vector of taxonomic strings
#' @param tax_sep Separator for taxonomic strings
#' @return A data frame with taxonomic information
#' @keywords internal
parse_taxonomy_strings <- function(tax_strings, tax_sep = ";") {
  # Remove any empty strings
  tax_strings <- tax_strings[nzchar(tax_strings)]
  
  if (length(tax_strings) == 0) {
    return(data.frame())
  }
  
  # Split the strings
  tax_parts <- strsplit(tax_strings, tax_sep)
  
  # Determine the maximum number of levels
  max_levels <- max(sapply(tax_parts, length))
  
  # Try to detect rank prefixes (k__, p__, etc.)
  first_parts <- sapply(tax_parts, function(x) if(length(x) > 0) x[1] else "")
  has_prefix <- any(grepl("^[kpcofgs]__", first_parts))
  
  # Determine column names
  if (has_prefix) {
    # Extract prefixes and use as column names
    col_prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    col_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Map prefixes to full rank names
    rank_map <- setNames(col_names, col_prefixes)
    
    # Initialize the data frame
    tax_df <- data.frame(matrix("", nrow = length(tax_strings), ncol = length(col_names)))
    colnames(tax_df) <- col_names
    
    # Fill in the data frame
    for (i in seq_along(tax_strings)) {
      parts <- tax_parts[[i]]
      for (part in parts) {
        # Extract the prefix (first 3 characters)
        prefix <- substr(part, 1, 3)
        # If it's a known prefix, add to the data frame
        if (prefix %in% names(rank_map)) {
          tax_df[i, rank_map[prefix]] <- sub("^[kpcofgs]__", "", part)
        }
      }
    }
  } else {
    # Use standard rank names
    col_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    if (max_levels > length(col_names)) {
      extra_cols <- paste0("Rank", (length(col_names)+1):max_levels)
      col_names <- c(col_names, extra_cols)
    } else if (max_levels < length(col_names)) {
      col_names <- col_names[1:max_levels]
    }
    
    # Initialize the data frame
    tax_df <- data.frame(matrix("", nrow = length(tax_strings), ncol = length(col_names)))
    colnames(tax_df) <- col_names
    
    # Fill in the data frame
    for (i in seq_along(tax_strings)) {
      parts <- tax_parts[[i]]
      for (j in seq_along(parts)) {
        if (j <= ncol(tax_df)) {
          tax_df[i, j] <- parts[j]
        }
      }
    }
  }
  
  # Convert empty strings to NA
  tax_df[tax_df == ""] <- NA
  
  return(tax_df)
}

#' Process sample metadata input
#'
#' Internal function to handle various sample metadata input formats
#'
#' @param sample_data Input sample metadata
#' @param sample_col_name Name of the column that contains sample IDs
#' @return A data frame with sample metadata
#' @keywords internal
process_sample_data <- function(sample_data, sample_col_name = NULL) {
  if (is.null(sample_data)) {
    return(data.frame())
  }
  
  # If input is a file path
  if (is.character(sample_data) && length(sample_data) == 1) {
    if (file.exists(sample_data)) {
      if (grepl("\\.rds$", sample_data, ignore.case = TRUE)) {
        # Load RDS file
        sample_df <- readRDS(sample_data)
      } else if (grepl("\\.(csv|tsv|txt)$", sample_data, ignore.case = TRUE)) {
        # Determine separator
        sep <- if (grepl("\\.tsv$", sample_data, ignore.case = TRUE)) "\t" else ","
        # Read CSV/TSV file
        sample_df <- read.table(sample_data, header = TRUE, sep = sep, 
                             row.names = 1, check.names = FALSE)
      } else {
        stop("Unsupported file format for sample_data. Use CSV, TSV, or RDS.")
      }
    } else {
      stop("File not found: ", sample_data)
    }
  } else if (is.data.frame(sample_data)) {
    sample_df <- sample_data
  } else if (is.matrix(sample_data)) {
    sample_df <- as.data.frame(sample_data, stringsAsFactors = FALSE)
  } else {
    stop("Sample metadata must be a data frame, matrix, or file path")
  }
  
  # If a sample column name is provided, set it as rownames
  if (!is.null(sample_col_name) && sample_col_name %in% colnames(sample_df)) {
    rownames(sample_df) <- sample_df[[sample_col_name]]
    sample_df <- sample_df[, colnames(sample_df) != sample_col_name, drop = FALSE]
  }
  
  return(sample_df)
}

#' Process phylogenetic tree input
#'
#' Internal function to handle various phylogenetic tree input formats
#'
#' @param tree Input phylogenetic tree
#' @return A phylo object or NULL
#' @keywords internal
process_tree_data <- function(tree) {
  if (is.null(tree)) {
    return(NULL)
  }
  
  # Check if ape package is available
  if (!requireNamespace("ape", quietly = TRUE)) {
    warning("ape package is required for processing phylogenetic trees")
    return(NULL)
  }
  
  # If input is a file path
  if (is.character(tree) && length(tree) == 1) {
    if (file.exists(tree)) {
      if (grepl("\\.rds$", tree, ignore.case = TRUE)) {
        # Load RDS file
        phy_tree <- readRDS(tree)
      } else if (grepl("\\.(nwk|tre|tree)$", tree, ignore.case = TRUE)) {
        # Read Newick file
        phy_tree <- ape::read.tree(tree)
      } else {
        stop("Unsupported file format for tree. Use Newick (.nwk, .tre) or RDS.")
      }
    } else {
      stop("File not found: ", tree)
    }
  } else {
    phy_tree <- tree
  }
  
  # Verify it's a phylo object
  if (!inherits(phy_tree, "phylo")) {
    stop("Tree must be a phylo object or file path to a valid tree file")
  }
  
  return(phy_tree)
}

#' Import data from DADA2 results
#'
#' Creates a SummarizedExperiment or FGTExperiment object from DADA2 output.
#'
#' @param seqtab Sequence table from DADA2 (ASVs x samples)
#' @param taxa Taxonomy table from DADA2
#' @param metadata Sample metadata (optional)
#' @param refseq Reference sequences, e.g., DNAStringSet (optional)
#' @param as_FGTExperiment Logical; if TRUE and FGTExperiment exists, return FGTExperiment
#'
#' @return A SummarizedExperiment or FGTExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Load data from dada2 outputs
#' seqtab <- readRDS("seqtab.rds")
#' taxa <- readRDS("taxa.rds")
#' metadata <- read.csv("metadata.csv", row.names = 1)
#' 
#' # Create SummarizedExperiment object
#' se <- import_from_dada2(seqtab, taxa, metadata)
#' }
import_from_dada2 <- function(seqtab, taxa = NULL, metadata = NULL, refseq = NULL,
                            as_FGTExperiment = TRUE) {
  # Handle file path inputs
  if (is.character(seqtab) && file.exists(seqtab)) {
    seqtab <- readRDS(seqtab)
  }
  
  if (!is.null(taxa) && is.character(taxa) && file.exists(taxa)) {
    taxa <- readRDS(taxa)
  }
  
  # Process input data
  counts_matrix <- process_dada2_seqtab(seqtab)
  taxonomy_df <- process_dada2_taxa(taxa, rownames(counts_matrix))
  
  # Import using the general import function
  result <- import_microbiome(
    counts = counts_matrix,
    taxonomy = taxonomy_df,
    sample_data = metadata,
    as_FGTExperiment = as_FGTExperiment
  )
  
  # Add reference sequences if available
  if (!is.null(refseq) && methods::is(result, "FGTExperiment")) {
    if (requireNamespace("Biostrings", quietly = TRUE)) {
      if (!methods::is(refseq, "DNAStringSet")) {
        refseq <- Biostrings::DNAStringSet(refseq)
        names(refseq) <- rownames(counts_matrix)
      }
      
      # Store reference sequences
      fgt_meta <- S4Vectors::SimpleList(referenceSeq = refseq)
      fgtMetadata(result) <- fgt_meta
    }
  }
  
  # Add import information
  if (methods::is(result, "SummarizedExperiment")) {
    SummarizedExperiment::metadata(result)$import_info <- list(
      source = "dada2",
      date = Sys.time(),
      n_features = nrow(counts_matrix),
      n_samples = ncol(counts_matrix)
    )
  }
  
  return(result)
}

#' Process DADA2 sequence table
#'
#' Internal function to process DADA2 sequence tables
#'
#' @param seqtab DADA2 sequence table
#' @return A count matrix with features as rows and samples as columns
#' @keywords internal
process_dada2_seqtab <- function(seqtab) {
  # Ensure seqtab is a matrix
  if (!is.matrix(seqtab)) {
    seqtab <- as.matrix(seqtab)
  }
  
  # Create sequence names if they don't exist
  if (is.null(rownames(seqtab))) {
    if (is.null(colnames(seqtab))) {
      stop("Sequence table must have row or column names")
    }
    seqtab <- t(seqtab)  # Transpose to make sequences as rows
  }
  
  # Check if sequences are in columns (samples in rows)
  if (ncol(seqtab) > nrow(seqtab) && all(grepl("^[ACGT]+$", colnames(seqtab)))) {
    seqtab <- t(seqtab)  # Transpose to make sequences as rows
  }
  
  # Extract ASV sequences and create IDs
  asv_seqs <- rownames(seqtab)
  if (all(nchar(asv_seqs) > 20)) {
    # These are likely DNA sequences, create shorter IDs
    asv_ids <- paste0("ASV", seq_along(asv_seqs))
    
    # Add sequences to rowData
    seq_info <- data.frame(sequence = asv_seqs, row.names = asv_ids)
    rownames(seqtab) <- asv_ids
  }
  
  return(seqtab)
}

#' Process DADA2 taxonomy table
#'
#' Internal function to process DADA2 taxonomy tables
#'
#' @param taxa DADA2 taxonomy table
#' @param feature_ids Feature IDs to match
#' @return A data frame with taxonomic information
#' @keywords internal
process_dada2_taxa <- function(taxa, feature_ids) {
  if (is.null(taxa)) {
    return(data.frame(row.names = feature_ids))
  }
  
  # Convert to data frame if matrix
  if (is.matrix(taxa)) {
    taxa_df <- as.data.frame(taxa, stringsAsFactors = FALSE)
  } else if (is.data.frame(taxa)) {
    taxa_df <- taxa
  } else {
    stop("Taxa must be a data frame or matrix")
  }
  
  # For DADA2 format with a single Taxonomy column, parse it
  if (ncol(taxa_df) == 1 && colnames(taxa_df)[1] %in% c("Taxonomy", "taxonomy")) {
    # Parse the taxonomy strings
    taxa_strings <- taxa_df[[1]]
    parsed_taxa <- parse_taxonomy_strings(taxa_strings, tax_sep = ";")
    
    # Set the rownames
    if (!is.null(rownames(taxa_df))) {
      rownames(parsed_taxa) <- rownames(taxa_df)
    }
    
    taxa_df <- parsed_taxa
  }
  
  # Match rownames with feature_ids
  if (!all(feature_ids %in% rownames(taxa_df))) {
    # If not all features have taxonomy, create a new data frame
    out_df <- data.frame(row.names = feature_ids)
    common_rows <- intersect(feature_ids, rownames(taxa_df))
    
    if (length(common_rows) > 0) {
      # Copy data for matching rows
      for (col in colnames(taxa_df)) {
        out_df[common_rows, col] <- taxa_df[common_rows, col]
      }
    }
    
    taxa_df <- out_df
  } else if (!all(rownames(taxa_df) %in% feature_ids)) {
    # If taxonomy has extra rows, subset to matching ones
    taxa_df <- taxa_df[feature_ids, , drop = FALSE]
  }
  
  return(taxa_df)
}

#' Import from Biom format files
#'
#' Creates a SummarizedExperiment or FGTExperiment object from a BIOM format file.
#'
#' @param biom_file Path to BIOM format file
#' @param metadata_file Optional path to metadata file
#' @param as_FGTExperiment Logical; if TRUE and FGTExperiment exists, return FGTExperiment
#'
#' @return A SummarizedExperiment or FGTExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Import from BIOM file
#' se <- import_from_biom("table.biom", "metadata.tsv")
#' }
import_from_biom <- function(biom_file, metadata_file = NULL, as_FGTExperiment = TRUE) {
  # Check if biomformat package is available
  if (!requireNamespace("biomformat", quietly = TRUE)) {
    stop("biomformat package is required for importing BIOM files")
  }
  
  # Read BIOM file
  biom <- biomformat::read_biom(biom_file)
  
  # Extract counts
  counts_matrix <- biomformat::as.matrix(biom)
  
  # Extract taxonomy if available
  taxonomy_df <- NULL
  if ("taxonomy" %in% biomformat::observation_metadata_names(biom)) {
    # Get taxonomy from observation metadata
    obs_meta <- biomformat::observation_metadata(biom)
    taxonomy_vec <- sapply(obs_meta, function(x) {
      if ("taxonomy" %in% names(x)) {
        if (is.character(x$taxonomy)) {
          return(x$taxonomy)
        } else if (is.list(x$taxonomy) || is.vector(x$taxonomy)) {
          return(paste(x$taxonomy, collapse = ";"))
        }
      }
      return(NA)
    })
    
    # Parse taxonomy strings
    taxonomy_df <- parse_taxonomy_strings(taxonomy_vec, tax_sep = ";")
    rownames(taxonomy_df) <- rownames(counts_matrix)
  }
  
  # Read metadata file if provided
  sample_df <- NULL
  if (!is.null(metadata_file) && file.exists(metadata_file)) {
    # Determine separator
    sep <- if (grepl("\\.tsv$", metadata_file, ignore.case = TRUE)) "\t" else ","
    # Read metadata file
    sample_df <- read.table(metadata_file, header = TRUE, sep = sep,
                         row.names = 1, check.names = FALSE)
  }
  
  # Import using the general import function
  result <- import_microbiome(
    counts = counts_matrix,
    taxonomy = taxonomy_df,
    sample_data = sample_df,
    as_FGTExperiment = as_FGTExperiment
  )
  
  # Add import information
  if (methods::is(result, "SummarizedExperiment")) {
    SummarizedExperiment::metadata(result)$import_info <- list(
      source = "biom",
      biom_file = biom_file,
      metadata_file = metadata_file,
      date = Sys.time()
    )
  }
  
  return(result)
}

#' Import from phyloseq object
#'
#' Creates a SummarizedExperiment or FGTExperiment object from a phyloseq object.
#'
#' @param physeq A phyloseq object
#' @param assay_name Name for the count assay
#' @param as_FGTExperiment Logical; if TRUE and FGTExperiment exists, return FGTExperiment
#'
#' @return A SummarizedExperiment or FGTExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Load phyloseq object
#' library(phyloseq)
#' data(GlobalPatterns)
#' 
#' # Convert to SummarizedExperiment
#' se <- import_from_phyloseq(GlobalPatterns)
#' }
import_from_phyloseq <- function(physeq, assay_name = "counts", as_FGTExperiment = TRUE) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required for this function")
  }
  
  if (!methods::is(physeq, "phyloseq")) {
    stop("Input must be a phyloseq object")
  }
  
  # Extract components from phyloseq
  otu_table <- phyloseq::otu_table(physeq)
  tax_table <- tryCatch(phyloseq::tax_table(physeq), error = function(e) NULL)
  sample_data <- tryCatch(phyloseq::sample_data(physeq), error = function(e) NULL)
  phy_tree <- tryCatch(phyloseq::phy_tree(physeq), error = function(e) NULL)
  
  # Ensure OTU table is taxa x samples
  if (phyloseq::taxa_are_rows(otu_table) == FALSE) {
    otu_table <- t(otu_table)
  }
  
  # Convert to matrices/data frames
  counts_matrix <- as.matrix(otu_table)
  
  if (!is.null(tax_table)) {
    taxonomy_df <- as.data.frame(tax_table)
  } else {
    taxonomy_df <- NULL
  }
  
  if (!is.null(sample_data)) {
    sample_df <- as.data.frame(sample_data)
  } else {
    sample_df <- NULL
  }
  
  # Import using the general import function
  result <- import_microbiome(
    counts = counts_matrix,
    taxonomy = taxonomy_df,
    sample_data = sample_df,
    tree = phy_tree,
    as_FGTExperiment = as_FGTExperiment
  )
  
  # Add import information
  if (methods::is(result, "SummarizedExperiment")) {
    SummarizedExperiment::metadata(result)$import_info <- list(
      source = "phyloseq",
      date = Sys.time()
    )
  }
  
  return(result)
}