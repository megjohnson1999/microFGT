#' Get taxonomic ranks from a SummarizedExperiment-like object
#'
#' Identifies columns in the rowData that are likely taxonomic ranks
#'
#' @param x A SummarizedExperiment-like object with taxonomic information
#' @param standard_only Logical; if TRUE, returns only standard taxonomic ranks
#'
#' @return Character vector of available taxonomic ranks
#' @export
get_taxonomic_ranks <- function(x, standard_only = TRUE) {
  # Input validation
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required")
  }
  
  if (!methods::is(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment-like object")
  }
  
  # Get rowData column names
  rowdata <- SummarizedExperiment::rowData(x)
  
  if (ncol(rowdata) == 0) {
    return(character(0))
  }
  
  # Standard taxonomic ranks to look for
  standard_ranks <- c(
    "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    "kingdom", "phylum", "class", "order", "family", "genus", "species"
  )
  
  # Find which columns match standard ranks
  rank_cols <- colnames(rowdata)[colnames(rowdata) %in% standard_ranks]
  
  # If requested only standard ranks but none found, return empty vector
  if (standard_only && length(rank_cols) == 0) {
    return(character(0))
  }
  
  # If no standard ranks found or all columns requested, return all column names
  if (!standard_only || length(rank_cols) == 0) {
    return(colnames(rowdata))
  }
  
  return(rank_cols)
}

#' Aggregate taxa at a specific taxonomic rank
#'
#' Aggregates feature abundances at a specified taxonomic rank (e.g., Genus, Family).
#' Works with any SummarizedExperiment-like object containing taxonomic information.
#'
#' @param x A SummarizedExperiment-like object
#' @param rank Taxonomic rank to aggregate at (e.g., "Genus", "Family")
#' @param assay_name Name of the assay to use for aggregation
#' @param empty_label Label to use for features without classification at the specified rank
#' @param add_metadata Logical; whether to add aggregation metadata
#'
#' @return A SummarizedExperiment object with aggregated features
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data with taxonomy
#' counts <- matrix(sample(1:100, 60), nrow = 10, ncol = 6)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Sample", 1:6)
#' 
#' taxonomy <- data.frame(
#'   Phylum = sample(c("Firmicutes", "Bacteroidetes"), 10, replace = TRUE),
#'   Genus = sample(c("Lactobacillus", "Gardnerella"), 10, replace = TRUE),
#'   row.names = rownames(counts)
#' )
#' 
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = counts),
#'   rowData = taxonomy
#' )
#' 
#' # Aggregate counts at genus level
#' genus_level <- aggregate_taxa(se, rank = "Genus")
#' }
aggregate_taxa <- function(x, rank, assay_name = "counts", empty_label = "Unclassified", 
                          add_metadata = TRUE) {
  # Input validation
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required")
  }
  
  if (!methods::is(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment-like object")
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(x)) {
    stop(paste0("Assay '", assay_name, "' not found in the object"))
  }
  
  # Check if rank column exists in rowData
  rowdata <- SummarizedExperiment::rowData(x)
  if (!rank %in% colnames(rowdata)) {
    stop(paste0("Taxonomic rank '", rank, "' not found in rowData"))
  }
  
  # Extract count matrix and taxonomic information
  counts <- SummarizedExperiment::assays(x)[[assay_name]]
  taxa_info <- rowdata[[rank]]
  
  # Replace NA or empty values with Unclassified
  taxa_info[is.na(taxa_info) | taxa_info == ""] <- empty_label
  
  # Get unique taxa at the specified rank
  unique_taxa <- unique(taxa_info)
  
  # Create new count matrix for aggregated taxa
  agg_counts <- matrix(0, nrow = length(unique_taxa), ncol = ncol(counts))
  rownames(agg_counts) <- unique_taxa
  colnames(agg_counts) <- colnames(counts)
  
  # Aggregate counts by taxonomic rank
  for (i in seq_along(unique_taxa)) {
    taxon <- unique_taxa[i]
    rows_to_sum <- which(taxa_info == taxon)
    agg_counts[i, ] <- colSums(counts[rows_to_sum, , drop = FALSE])
  }
  
  # Create new rowData for aggregated taxa
  agg_rowdata <- S4Vectors::DataFrame(row.names = unique_taxa)
  agg_rowdata[[rank]] <- unique_taxa
  
  # Preserve all other assays with aggregation
  all_assays <- list()
  for (name in SummarizedExperiment::assayNames(x)) {
    current_assay <- SummarizedExperiment::assays(x)[[name]]
    agg_assay <- matrix(0, nrow = length(unique_taxa), ncol = ncol(current_assay))
    rownames(agg_assay) <- unique_taxa
    colnames(agg_assay) <- colnames(current_assay)
    
    for (i in seq_along(unique_taxa)) {
      taxon <- unique_taxa[i]
      rows_to_sum <- which(taxa_info == taxon)
      agg_assay[i, ] <- colSums(current_assay[rows_to_sum, , drop = FALSE])
    }
    all_assays[[name]] <- agg_assay
  }
  
  # Create new SummarizedExperiment with aggregated data
  if (methods::is(x, "TreeSummarizedExperiment") && 
      requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    # Create a TreeSummarizedExperiment if the input was one
    agg_se <- TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = all_assays,
      rowData = agg_rowdata,
      colData = SummarizedExperiment::colData(x)
    )
  } else if (methods::is(x, "FGTExperiment") && 
             exists("FGTExperiment", mode = "function")) {
    # Create an FGTExperiment if input was one and function exists
    agg_se <- FGTExperiment(
      assays = all_assays,
      rowData = agg_rowdata,
      colData = SummarizedExperiment::colData(x),
      experimentType = experimentType(x)
    )
  } else {
    # Default to regular SummarizedExperiment
    agg_se <- SummarizedExperiment::SummarizedExperiment(
      assays = all_assays,
      rowData = agg_rowdata,
      colData = SummarizedExperiment::colData(x)
    )
  }
  
  # Add aggregation information to metadata
  if (add_metadata) {
    agg_metadata <- list(
      original_feature_count = nrow(x),
      aggregated_feature_count = length(unique_taxa),
      rank = rank,
      assay_used = assay_name,
      date = Sys.time()
    )
    
    # Get existing metadata if any
    existing_metadata <- SummarizedExperiment::metadata(x)
    if (length(existing_metadata) > 0) {
      existing_metadata$aggregation_info <- agg_metadata
      SummarizedExperiment::metadata(agg_se) <- existing_metadata
    } else {
      SummarizedExperiment::metadata(agg_se) <- list(aggregation_info = agg_metadata)
    }
  }
  
  return(agg_se)
}

#' Normalize taxonomic names in a SummarizedExperiment-like object
#'
#' Cleans and normalizes taxonomic names in the rowData to ensure consistency.
#'
#' @param x A SummarizedExperiment-like object
#' @param ranks Character vector of taxonomic ranks to normalize
#' @param remove_prefixes Logical; whether to remove prefixes like "k__", "p__"
#' @param remove_confidence Logical; whether to remove confidence values in brackets
#' @param make_unique Logical; whether to make taxa unique by adding numeric suffixes
#'
#' @return SummarizedExperiment with normalized taxonomic names
#' @export
normalize_taxonomy <- function(x, ranks = NULL, remove_prefixes = TRUE, 
                             remove_confidence = TRUE, make_unique = TRUE) {
  # Input validation
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required")
  }
  
  if (!methods::is(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment-like object")
  }
  
  # If ranks not provided, get them from the object
  if (is.null(ranks)) {
    ranks <- get_taxonomic_ranks(x)
    if (length(ranks) == 0) {
      warning("No taxonomic ranks found in the object")
      return(x)
    }
  }
  
  # Get rowData
  rowdata <- SummarizedExperiment::rowData(x)
  
  # Process each taxonomic rank
  for (rank in ranks) {
    if (rank %in% colnames(rowdata)) {
      taxa <- rowdata[[rank]]
      
      # Skip if column is not character
      if (!is.character(taxa)) {
        next
      }
      
      # Remove standard prefixes if requested
      if (remove_prefixes) {
        taxa <- gsub("^[kpcofgs]__", "", taxa)
      }
      
      # Remove confidence values in brackets if requested
      if (remove_confidence) {
        taxa <- gsub("\\s*\\([^)]+\\)\\s*$", "", taxa)
      }
      
      # Replace empty or NA values with "Unclassified"
      taxa[is.na(taxa) | taxa == ""] <- paste0("Unclassified_", rank)
      
      # Make taxa unique if requested
      if (make_unique) {
        # Find duplicates
        dups <- duplicated(taxa)
        if (any(dups)) {
          # Create a table of taxa
          taxa_table <- table(taxa)
          duplicated_taxa <- names(taxa_table[taxa_table > 1])
          
          # Add suffix to duplicates
          for (dup_taxon in duplicated_taxa) {
            dup_indices <- which(taxa == dup_taxon)
            for (i in seq_along(dup_indices)[-1]) {  # Skip the first occurrence
              taxa[dup_indices[i]] <- paste0(dup_taxon, "_", i)
            }
          }
        }
      }
      
      # Update rowData
      rowdata[[rank]] <- taxa
    }
  }
  
  # Update object with modified rowData
  SummarizedExperiment::rowData(x) <- rowdata
  
  return(x)
}

#' Create taxonomic strings from individual ranks
#'
#' Creates full taxonomic strings (e.g., "p__Firmicutes;g__Lactobacillus") from individual rank columns.
#'
#' @param x A SummarizedExperiment-like object
#' @param ranks Character vector of taxonomic ranks to include (in order)
#' @param format Format string: "lineage" (semi-colon separated), "text" (human readable)
#' @param add_prefixes Logical; whether to add standard rank prefixes (k__, p__, etc.)
#' @param new_column Name for the new column to add to rowData
#'
#' @return SummarizedExperiment with taxonomic strings added to rowData
#' @export
create_tax_strings <- function(x, ranks = NULL, format = c("lineage", "text"), 
                            add_prefixes = TRUE, new_column = "taxonomy") {
  # Input validation
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required")
  }
  
  if (!methods::is(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment-like object")
  }
  
  format <- match.arg(format)
  
  # If ranks not provided, get standard ranks from the object
  if (is.null(ranks)) {
    standard_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    rowdata <- SummarizedExperiment::rowData(x)
    available_ranks <- colnames(rowdata)[tolower(colnames(rowdata)) %in% standard_ranks]
    
    if (length(available_ranks) == 0) {
      stop("No standard taxonomic ranks found in rowData")
    }
    
    # Sort ranks in standard order
    rank_order <- match(tolower(available_ranks), standard_ranks)
    ranks <- available_ranks[order(rank_order)]
  }
  
  # Get rowData
  rowdata <- SummarizedExperiment::rowData(x)
  
  # Check that all ranks exist
  missing_ranks <- ranks[!ranks %in% colnames(rowdata)]
  if (length(missing_ranks) > 0) {
    stop("The following ranks are missing from rowData: ", 
         paste(missing_ranks, collapse = ", "))
  }
  
  # Standard rank prefixes
  rank_prefixes <- c(
    kingdom = "k__", phylum = "p__", class = "c__", order = "o__", 
    family = "f__", genus = "g__", species = "s__"
  )
  
  # Create taxonomic strings
  tax_strings <- character(nrow(rowdata))
  
  for (i in seq_len(nrow(rowdata))) {
    rank_values <- sapply(ranks, function(r) {
      val <- rowdata[i, r]
      if (is.na(val) || val == "") {
        return(NA_character_)
      }
      
      if (add_prefixes) {
        prefix <- rank_prefixes[tolower(r)]
        if (!is.na(prefix)) {
          return(paste0(prefix, val))
        }
      }
      
      return(as.character(val))
    })
    
    # Remove NAs
    rank_values <- rank_values[!is.na(rank_values)]
    
    if (length(rank_values) == 0) {
      tax_strings[i] <- "Unclassified"
    } else if (format == "lineage") {
      tax_strings[i] <- paste(rank_values, collapse = ";")
    } else {
      # Text format
      if (add_prefixes) {
        # Remove prefixes for text format
        rank_values <- gsub("^[kpcofgs]__", "", rank_values)
      }
      tax_strings[i] <- paste(rank_values, collapse = " ")
    }
  }
  
  # Add to rowData
  rowdata[[new_column]] <- tax_strings
  SummarizedExperiment::rowData(x) <- rowdata
  
  return(x)
}

#' Parse taxonomic strings into separate rank columns
#'
#' Splits taxonomic strings (e.g., "p__Firmicutes;g__Lactobacillus") into individual rank columns.
#'
#' @param x A SummarizedExperiment-like object
#' @param tax_column Name of the column containing taxonomic strings
#' @param split Character to split taxonomy strings on (default ";")
#' @param remove_prefixes Logical; whether to remove rank prefixes (k__, p__, etc.)
#'
#' @return SummarizedExperiment with individual taxonomic ranks added to rowData
#' @export
parse_tax_strings <- function(x, tax_column = "taxonomy", split = ";", remove_prefixes = TRUE) {
  # Input validation
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required")
  }
  
  if (!methods::is(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment-like object")
  }
  
  # Get rowData
  rowdata <- SummarizedExperiment::rowData(x)
  
  # Check if tax_column exists
  if (!tax_column %in% colnames(rowdata)) {
    stop(paste0("Column '", tax_column, "' not found in rowData"))
  }
  
  # Get taxonomy strings
  tax_strings <- rowdata[[tax_column]]
  
  # Define standard ranks and their prefixes
  rank_prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  rank_names <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # Split the taxonomy strings
  tax_parts <- strsplit(tax_strings, split)
  
  # Create a matrix to store the parsed taxonomy
  max_parts <- max(sapply(tax_parts, length))
  tax_matrix <- matrix(NA_character_, nrow = length(tax_strings), ncol = max_parts)
  
  # Fill the matrix
  for (i in seq_along(tax_parts)) {
    parts <- tax_parts[[i]]
    tax_matrix[i, seq_along(parts)] <- parts
  }
  
  # Detect rank names from prefixes
  rank_names_detected <- character(ncol(tax_matrix))
  for (j in seq_len(ncol(tax_matrix))) {
    # Look at the first 3 characters of non-NA values in the column
    prefixes <- substr(tax_matrix[!is.na(tax_matrix[, j]), j], 1, 3)
    if (length(prefixes) > 0) {
      # Find the most common prefix
      most_common <- names(sort(table(prefixes), decreasing = TRUE))[1]
      
      # Match with standard prefixes
      idx <- match(most_common, rank_prefixes)
      if (!is.na(idx)) {
        rank_names_detected[j] <- rank_names[idx]
      } else {
        rank_names_detected[j] <- paste0("rank", j)
      }
    } else {
      rank_names_detected[j] <- paste0("rank", j)
    }
  }
  
  # Remove prefixes if requested
  if (remove_prefixes) {
    for (j in seq_len(ncol(tax_matrix))) {
      # Remove prefixes from non-NA values
      idx <- !is.na(tax_matrix[, j])
      tax_matrix[idx, j] <- gsub("^[kpcofgs]__", "", tax_matrix[idx, j])
    }
  }
  
  # Add parsed taxonomy to rowData
  for (j in seq_len(ncol(tax_matrix))) {
    rowdata[[rank_names_detected[j]]] <- tax_matrix[, j]
  }
  
  # Update object with modified rowData
  SummarizedExperiment::rowData(x) <- rowdata
  
  return(x)
}