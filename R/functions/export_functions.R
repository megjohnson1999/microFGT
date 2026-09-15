#' Export microbiome data to various formats
#'
#' A flexible export function that saves data from a SummarizedExperiment-like
#' object to various file formats.
#'
#' @param x A SummarizedExperiment or FGTExperiment object
#' @param output_dir Directory where files will be saved
#' @param format Format to export to: "csv", "rds", "biom", "qiime2", "phyloseq", "dada2"
#' @param assay_name Name of the assay to export
#' @param prefix Prefix for output filenames
#' @param compress Logical; whether to compress the output files
#' @param replace Logical; whether to replace existing files
#'
#' @return Invisibly returns a list of the exported file paths
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' se <- import_microbiome(counts_matrix, taxonomy_df, metadata_df)
#' 
#' # Export to CSV files
#' export_microbiome(se, "output_dir", format = "csv")
#' 
#' # Export to BIOM format
#' export_microbiome(se, "output_dir", format = "biom")
#' }
export_microbiome <- function(x, output_dir, format = c("csv", "rds", "biom", "qiime2", "phyloseq", "dada2"),
                             assay_name = "counts", prefix = NULL, compress = TRUE, replace = FALSE) {
  
  # Validate input
  if (!inherits(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object")
  }
  
  format <- match.arg(format)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set default prefix if not provided
  if (is.null(prefix)) {
    prefix <- ifelse(
      inherits(x, "FGTExperiment"), 
      paste0("fgt_export_", format_time(Sys.time())),
      paste0("microbiome_export_", format_time(Sys.time()))
    )
  }
  
  # Extract data from SummarizedExperiment
  if (!(assay_name %in% SummarizedExperiment::assayNames(x))) {
    stop("Assay '", assay_name, "' not found. Available assays: ",
         paste(SummarizedExperiment::assayNames(x), collapse = ", "))
  }
  
  counts <- SummarizedExperiment::assays(x)[[assay_name]]
  rowdata <- SummarizedExperiment::rowData(x)
  coldata <- SummarizedExperiment::colData(x)
  
  # Extract tree if available (for TreeSummarizedExperiment)
  tree <- NULL
  if (inherits(x, "TreeSummarizedExperiment") && requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    tree <- TreeSummarizedExperiment::rowTree(x)
  }
  
  # Call format-specific export function
  result <- switch(
    format,
    "csv" = export_to_csv(counts, rowdata, coldata, tree, output_dir, prefix, replace),
    "rds" = export_to_rds(counts, rowdata, coldata, tree, output_dir, prefix, compress, replace),
    "biom" = export_to_biom(counts, rowdata, coldata, output_dir, prefix, replace),
    "qiime2" = export_to_qiime2(counts, rowdata, coldata, tree, output_dir, prefix, replace),
    "phyloseq" = export_to_phyloseq(counts, rowdata, coldata, tree, output_dir, prefix, compress, replace),
    "dada2" = export_to_dada2(counts, rowdata, coldata, output_dir, prefix, compress, replace)
  )
  
  # Add export information
  result$info <- list(
    format = format,
    date = Sys.time(),
    output_dir = normalizePath(output_dir),
    prefix = prefix
  )
  
  return(invisible(result))
}

#' Format time for filenames
#'
#' @param time A POSIXct time object
#' @return A string with the formatted time
#' @keywords internal
format_time <- function(time) {
  format(time, "%Y%m%d_%H%M%S")
}

#' Export to CSV files
#'
#' @param counts Count matrix
#' @param rowdata Row data (taxonomy)
#' @param coldata Column data (sample metadata)
#' @param tree Phylogenetic tree (ignored for CSV export)
#' @param output_dir Output directory
#' @param prefix Prefix for output filenames
#' @param replace Whether to replace existing files
#'
#' @return A list of exported file paths
#' @keywords internal
export_to_csv <- function(counts, rowdata, coldata, tree = NULL, 
                         output_dir, prefix, replace = FALSE) {
  # Define file paths
  counts_file <- file.path(output_dir, paste0(prefix, "_counts.csv"))
  rowdata_file <- file.path(output_dir, paste0(prefix, "_taxonomy.csv"))
  coldata_file <- file.path(output_dir, paste0(prefix, "_metadata.csv"))
  
  # Check if files exist and should not be replaced
  if (!replace) {
    if (file.exists(counts_file)) {
      stop("File already exists: ", counts_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(rowdata_file)) {
      stop("File already exists: ", rowdata_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(coldata_file)) {
      stop("File already exists: ", coldata_file, ". Use replace = TRUE to overwrite.")
    }
  }
  
  # Write counts to CSV
  write.csv(counts, file = counts_file, row.names = TRUE)
  
  # Convert rowdata to data.frame if it's a DataFrame
  if (inherits(rowdata, "DataFrame")) {
    rowdata_df <- as.data.frame(rowdata)
  } else {
    rowdata_df <- rowdata
  }
  
  # Write taxonomy to CSV
  write.csv(rowdata_df, file = rowdata_file, row.names = TRUE)
  
  # Convert coldata to data.frame if it's a DataFrame
  if (inherits(coldata, "DataFrame")) {
    coldata_df <- as.data.frame(coldata)
  } else {
    coldata_df <- coldata
  }
  
  # Write metadata to CSV
  write.csv(coldata_df, file = coldata_file, row.names = TRUE)
  
  # Return file paths
  return(list(
    counts = counts_file,
    taxonomy = rowdata_file,
    metadata = coldata_file
  ))
}

#' Export to RDS files
#'
#' @param counts Count matrix
#' @param rowdata Row data (taxonomy)
#' @param coldata Column data (sample metadata)
#' @param tree Phylogenetic tree
#' @param output_dir Output directory
#' @param prefix Prefix for output filenames
#' @param compress Whether to compress the RDS files
#' @param replace Whether to replace existing files
#'
#' @return A list of exported file paths
#' @keywords internal
export_to_rds <- function(counts, rowdata, coldata, tree = NULL, 
                        output_dir, prefix, compress = TRUE, replace = FALSE) {
  # Define file paths
  counts_file <- file.path(output_dir, paste0(prefix, "_counts.rds"))
  rowdata_file <- file.path(output_dir, paste0(prefix, "_taxonomy.rds"))
  coldata_file <- file.path(output_dir, paste0(prefix, "_metadata.rds"))
  tree_file <- NULL
  if (!is.null(tree)) {
    tree_file <- file.path(output_dir, paste0(prefix, "_tree.rds"))
  }
  
  # Check if files exist and should not be replaced
  if (!replace) {
    if (file.exists(counts_file)) {
      stop("File already exists: ", counts_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(rowdata_file)) {
      stop("File already exists: ", rowdata_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(coldata_file)) {
      stop("File already exists: ", coldata_file, ". Use replace = TRUE to overwrite.")
    }
    if (!is.null(tree_file) && file.exists(tree_file)) {
      stop("File already exists: ", tree_file, ". Use replace = TRUE to overwrite.")
    }
  }
  
  # Write to RDS files
  saveRDS(counts, file = counts_file, compress = compress)
  saveRDS(rowdata, file = rowdata_file, compress = compress)
  saveRDS(coldata, file = coldata_file, compress = compress)
  
  # Save tree if available
  if (!is.null(tree)) {
    saveRDS(tree, file = tree_file, compress = compress)
  }
  
  # Return file paths
  result <- list(
    counts = counts_file,
    taxonomy = rowdata_file,
    metadata = coldata_file
  )
  
  if (!is.null(tree_file)) {
    result$tree <- tree_file
  }
  
  return(result)
}

#' Export to BIOM format
#'
#' @param counts Count matrix
#' @param rowdata Row data (taxonomy)
#' @param coldata Column data (sample metadata)
#' @param output_dir Output directory
#' @param prefix Prefix for output filenames
#' @param replace Whether to replace existing files
#'
#' @return A list of exported file paths
#' @keywords internal
export_to_biom <- function(counts, rowdata, coldata, output_dir, prefix, replace = FALSE) {
  # Check if biomformat package is available
  if (!requireNamespace("biomformat", quietly = TRUE)) {
    stop("biomformat package is required for BIOM format export")
  }
  
  # Define file paths
  biom_file <- file.path(output_dir, paste0(prefix, ".biom"))
  metadata_file <- file.path(output_dir, paste0(prefix, "_metadata.tsv"))
  
  # Check if files exist and should not be replaced
  if (!replace) {
    if (file.exists(biom_file)) {
      stop("File already exists: ", biom_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(metadata_file)) {
      stop("File already exists: ", metadata_file, ". Use replace = TRUE to overwrite.")
    }
  }
  
  # Prepare observation metadata (taxonomy)
  obs_meta <- NULL
  if (ncol(rowdata) > 0) {
    # Convert rowdata to list for BIOM format
    taxa_df <- as.data.frame(rowdata)
    
    # Check if we have standard taxonomy columns
    tax_cols <- intersect(colnames(taxa_df), 
                         c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                           "kingdom", "phylum", "class", "order", "family", "genus", "species"))
    
    if (length(tax_cols) > 0) {
      # Create taxonomy lists
      obs_meta <- lapply(seq_len(nrow(taxa_df)), function(i) {
        tax_values <- as.character(taxa_df[i, tax_cols])
        names(tax_values) <- tax_cols
        list(taxonomy = tax_values)
      })
      names(obs_meta) <- rownames(counts)
    } else {
      # No standard taxonomy columns, use all columns
      obs_meta <- lapply(seq_len(nrow(taxa_df)), function(i) {
        as.list(taxa_df[i, ])
      })
      names(obs_meta) <- rownames(counts)
    }
  }
  
  # Create BIOM object
  biom_obj <- biomformat::make_biom(
    data = counts,
    observation_metadata = obs_meta,
    sample_metadata = NULL # Sample metadata saved separately
  )
  
  # Convert coldata to data.frame if it's a DataFrame
  if (inherits(coldata, "DataFrame")) {
    sample_meta <- as.data.frame(coldata)
  } else {
    sample_meta <- coldata
  }
  
  # Write BIOM file
  biomformat::write_biom(biom_obj, biom_file)
  
  # Write sample metadata to TSV
  write.table(
    cbind("#SampleID" = rownames(sample_meta), sample_meta),
    file = metadata_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # Return file paths
  return(list(
    biom = biom_file,
    metadata = metadata_file
  ))
}

#' Export to QIIME2 format
#'
#' @param counts Count matrix
#' @param rowdata Row data (taxonomy)
#' @param coldata Column data (sample metadata)
#' @param tree Phylogenetic tree
#' @param output_dir Output directory
#' @param prefix Prefix for output filenames
#' @param replace Whether to replace existing files
#'
#' @return A list of exported file paths
#' @keywords internal
export_to_qiime2 <- function(counts, rowdata, coldata, tree = NULL, 
                           output_dir, prefix, replace = FALSE) {
  # Define file paths
  feature_table_file <- file.path(output_dir, paste0(prefix, "_feature-table.tsv"))
  taxonomy_file <- file.path(output_dir, paste0(prefix, "_taxonomy.tsv"))
  metadata_file <- file.path(output_dir, paste0(prefix, "_metadata.tsv"))
  tree_file <- NULL
  if (!is.null(tree) && requireNamespace("ape", quietly = TRUE)) {
    tree_file <- file.path(output_dir, paste0(prefix, "_tree.nwk"))
  }
  
  # Check if files exist and should not be replaced
  if (!replace) {
    if (file.exists(feature_table_file)) {
      stop("File already exists: ", feature_table_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(taxonomy_file)) {
      stop("File already exists: ", taxonomy_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(metadata_file)) {
      stop("File already exists: ", metadata_file, ". Use replace = TRUE to overwrite.")
    }
    if (!is.null(tree_file) && file.exists(tree_file)) {
      stop("File already exists: ", tree_file, ". Use replace = TRUE to overwrite.")
    }
  }
  
  # Write feature table with OTU/ASV IDs as first column
  feature_table <- cbind(
    "#OTU ID" = rownames(counts),
    as.data.frame(counts)
  )
  
  write.table(
    feature_table,
    file = feature_table_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # Format taxonomy for QIIME2
  taxa_df <- as.data.frame(rowdata)
  
  # Check if we have standard taxonomy columns
  tax_cols <- intersect(colnames(taxa_df), 
                       c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                         "kingdom", "phylum", "class", "order", "family", "genus", "species"))
  
  if (length(tax_cols) > 0) {
    # Create taxonomy strings in QIIME2 format (k__Kingdom; p__Phylum; etc)
    tax_strings <- apply(taxa_df[, tax_cols, drop = FALSE], 1, function(row) {
      tax_values <- sapply(seq_along(row), function(i) {
        if (is.na(row[i]) || row[i] == "") {
          return("")
        }
        
        prefix <- substr(tax_cols[i], 1, 1)
        if (prefix %in% c("K", "k", "P", "p", "C", "c", "O", "o", "F", "f", "G", "g", "S", "s")) {
          return(paste0(tolower(prefix), "__", row[i]))
        } else {
          return(row[i])
        }
      })
      paste(tax_values[tax_values != ""], collapse = "; ")
    })
  } else {
    # No standard taxonomy columns, use a generic approach
    tax_strings <- rep("k__Bacteria", nrow(taxa_df))
  }
  
  # Write taxonomy file
  taxonomy_df <- data.frame(
    "Feature ID" = rownames(counts),
    "Taxon" = tax_strings,
    "Confidence" = 1.0,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  write.table(
    taxonomy_df,
    file = taxonomy_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # Write metadata file
  sample_meta <- as.data.frame(coldata)
  
  write.table(
    cbind("#SampleID" = rownames(sample_meta), sample_meta),
    file = metadata_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # Write tree file if available
  if (!is.null(tree_file)) {
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
  
  return(result)
}

#' Export to phyloseq format
#'
#' @param counts Count matrix
#' @param rowdata Row data (taxonomy)
#' @param coldata Column data (sample metadata)
#' @param tree Phylogenetic tree
#' @param output_dir Output directory
#' @param prefix Prefix for output filenames
#' @param compress Whether to compress the output files
#' @param replace Whether to replace existing files
#'
#' @return A list of exported file paths
#' @keywords internal
export_to_phyloseq <- function(counts, rowdata, coldata, tree = NULL, 
                             output_dir, prefix, compress = TRUE, replace = FALSE) {
  # Check if phyloseq package is available
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq package is required for phyloseq format export")
  }
  
  # Define file path
  phyloseq_file <- file.path(output_dir, paste0(prefix, "_phyloseq.rds"))
  
  # Check if file exists and should not be replaced
  if (!replace && file.exists(phyloseq_file)) {
    stop("File already exists: ", phyloseq_file, ". Use replace = TRUE to overwrite.")
  }
  
  # Create phyloseq components
  otu_table <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  
  if (ncol(rowdata) > 0) {
    tax_table <- phyloseq::tax_table(as.matrix(as.data.frame(rowdata)))
  } else {
    tax_table <- NULL
  }
  
  if (ncol(coldata) > 0) {
    sample_data <- phyloseq::sample_data(as.data.frame(coldata))
  } else {
    sample_data <- NULL
  }
  
  # Create phyloseq object
  if (!is.null(tree) && !is.null(tax_table) && !is.null(sample_data)) {
    phy_tree <- phyloseq::phy_tree(tree)
    ps <- phyloseq::phyloseq(otu_table, tax_table, sample_data, phy_tree)
  } else if (!is.null(tax_table) && !is.null(sample_data)) {
    ps <- phyloseq::phyloseq(otu_table, tax_table, sample_data)
  } else if (!is.null(tax_table)) {
    ps <- phyloseq::phyloseq(otu_table, tax_table)
  } else if (!is.null(sample_data)) {
    ps <- phyloseq::phyloseq(otu_table, sample_data)
  } else {
    ps <- phyloseq::phyloseq(otu_table)
  }
  
  # Save phyloseq object
  saveRDS(ps, file = phyloseq_file, compress = compress)
  
  # Return file path
  return(list(
    phyloseq = phyloseq_file
  ))
}

#' Export to DADA2 format
#'
#' @param counts Count matrix
#' @param rowdata Row data (taxonomy)
#' @param coldata Column data (sample metadata)
#' @param output_dir Output directory
#' @param prefix Prefix for output filenames
#' @param compress Whether to compress the output files
#' @param replace Whether to replace existing files
#'
#' @return A list of exported file paths
#' @keywords internal
export_to_dada2 <- function(counts, rowdata, coldata, output_dir, prefix, 
                          compress = TRUE, replace = FALSE) {
  # Define file paths
  seqtab_file <- file.path(output_dir, paste0(prefix, "_seqtab.rds"))
  taxa_file <- file.path(output_dir, paste0(prefix, "_taxa.rds"))
  metadata_file <- file.path(output_dir, paste0(prefix, "_metadata.csv"))
  
  # Check if files exist and should not be replaced
  if (!replace) {
    if (file.exists(seqtab_file)) {
      stop("File already exists: ", seqtab_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(taxa_file)) {
      stop("File already exists: ", taxa_file, ". Use replace = TRUE to overwrite.")
    }
    if (file.exists(metadata_file)) {
      stop("File already exists: ", metadata_file, ". Use replace = TRUE to overwrite.")
    }
  }
  
  # Convert rowdata to DADA2 format if needed
  taxa_df <- as.data.frame(rowdata)
  
  # Check if we have standard taxonomy columns
  tax_cols <- intersect(colnames(taxa_df), 
                       c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                         "kingdom", "phylum", "class", "order", "family", "genus", "species"))
  
  if (length(tax_cols) > 0) {
    # Create taxonomy matrix
    tax_matrix <- matrix(NA, nrow = nrow(taxa_df), ncol = 7)
    colnames(tax_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    rownames(tax_matrix) <- rownames(taxa_df)
    
    # Map columns to standard taxonomy
    for (col in tax_cols) {
      std_col <- gsub("^([kpcofgs]).*$", "\\U\\1", col, perl = TRUE)
      std_col <- sub("^[a-z]$", toupper(std_col), std_col)
      if (std_col %in% colnames(tax_matrix)) {
        tax_matrix[, std_col] <- taxa_df[, col]
      }
    }
  } else {
    # No standard taxonomy columns, create a minimal one
    tax_matrix <- matrix("", nrow = nrow(taxa_df), ncol = 1)
    colnames(tax_matrix) <- "Taxonomy"
    rownames(tax_matrix) <- rownames(taxa_df)
  }
  
  # Save files
  saveRDS(counts, file = seqtab_file, compress = compress)
  saveRDS(tax_matrix, file = taxa_file, compress = compress)
  write.csv(as.data.frame(coldata), file = metadata_file, row.names = TRUE)
  
  # Return file paths
  return(list(
    seqtab = seqtab_file,
    taxa = taxa_file,
    metadata = metadata_file
  ))
}

#' Export to phyloseq object
#'
#' Converts a SummarizedExperiment-like object to a phyloseq object
#'
#' @param x A SummarizedExperiment-like object
#' @param assay_name Name of the assay to use
#'
#' @return A phyloseq object
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert SummarizedExperiment to phyloseq
#' ps <- to_phyloseq(se)
#' }
to_phyloseq <- function(x, assay_name = "counts") {
  # Check if phyloseq package is available
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq package is required for this function")
  }
  
  # Validate input
  if (!inherits(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object")
  }
  
  if (!(assay_name %in% SummarizedExperiment::assayNames(x))) {
    stop("Assay '", assay_name, "' not found. Available assays: ",
         paste(SummarizedExperiment::assayNames(x), collapse = ", "))
  }
  
  # Extract data
  counts <- SummarizedExperiment::assays(x)[[assay_name]]
  rowdata <- SummarizedExperiment::rowData(x)
  coldata <- SummarizedExperiment::colData(x)
  
  # Create phyloseq components
  otu_table <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
  
  if (ncol(rowdata) > 0) {
    tax_table <- phyloseq::tax_table(as.matrix(as.data.frame(rowdata)))
  } else {
    tax_table <- NULL
  }
  
  if (ncol(coldata) > 0) {
    sample_data <- phyloseq::sample_data(as.data.frame(coldata))
  } else {
    sample_data <- NULL
  }
  
  # Extract tree if available (for TreeSummarizedExperiment)
  tree <- NULL
  if (inherits(x, "TreeSummarizedExperiment") && requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    tree <- TreeSummarizedExperiment::rowTree(x)
    if (!is.null(tree)) {
      phy_tree <- phyloseq::phy_tree(tree)
    } else {
      phy_tree <- NULL
    }
  } else {
    phy_tree <- NULL
  }
  
  # Create phyloseq object
  if (!is.null(tax_table) && !is.null(sample_data) && !is.null(phy_tree)) {
    ps <- phyloseq::phyloseq(otu_table, tax_table, sample_data, phy_tree)
  } else if (!is.null(tax_table) && !is.null(sample_data)) {
    ps <- phyloseq::phyloseq(otu_table, tax_table, sample_data)
  } else if (!is.null(tax_table)) {
    ps <- phyloseq::phyloseq(otu_table, tax_table)
  } else if (!is.null(sample_data)) {
    ps <- phyloseq::phyloseq(otu_table, sample_data)
  } else {
    ps <- phyloseq::phyloseq(otu_table)
  }
  
  return(ps)
}

#' Convert to Biom object
#'
#' Converts a SummarizedExperiment-like object to a Biom-class object
#'
#' @param x A SummarizedExperiment-like object
#' @param assay_name Name of the assay to use
#'
#' @return A Biom-class object
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert SummarizedExperiment to Biom
#' biom_obj <- to_biom(se)
#' }
to_biom <- function(x, assay_name = "counts") {
  # Check if biomformat package is available
  if (!requireNamespace("biomformat", quietly = TRUE)) {
    stop("biomformat package is required for this function")
  }
  
  # Validate input
  if (!inherits(x, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object")
  }
  
  if (!(assay_name %in% SummarizedExperiment::assayNames(x))) {
    stop("Assay '", assay_name, "' not found. Available assays: ",
         paste(SummarizedExperiment::assayNames(x), collapse = ", "))
  }
  
  # Extract data
  counts <- SummarizedExperiment::assays(x)[[assay_name]]
  rowdata <- SummarizedExperiment::rowData(x)
  coldata <- SummarizedExperiment::colData(x)
  
  # Prepare observation metadata (taxonomy)
  obs_meta <- NULL
  if (ncol(rowdata) > 0) {
    # Convert rowdata to list for BIOM format
    taxa_df <- as.data.frame(rowdata)
    
    # Check if we have standard taxonomy columns
    tax_cols <- intersect(colnames(taxa_df), 
                         c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                           "kingdom", "phylum", "class", "order", "family", "genus", "species"))
    
    if (length(tax_cols) > 0) {
      # Create taxonomy lists
      obs_meta <- lapply(seq_len(nrow(taxa_df)), function(i) {
        tax_values <- as.character(taxa_df[i, tax_cols])
        names(tax_values) <- tax_cols
        list(taxonomy = tax_values)
      })
      names(obs_meta) <- rownames(counts)
    } else {
      # No standard taxonomy columns, use all columns
      obs_meta <- lapply(seq_len(nrow(taxa_df)), function(i) {
        as.list(taxa_df[i, ])
      })
      names(obs_meta) <- rownames(counts)
    }
  }
  
  # Prepare sample metadata
  sample_meta <- NULL
  if (ncol(coldata) > 0) {
    # Convert coldata to list for BIOM format
    sample_meta <- lapply(seq_len(nrow(coldata)), function(i) {
      as.list(as.data.frame(coldata)[i, ])
    })
    names(sample_meta) <- rownames(coldata)
  }
  
  # Create BIOM object
  biom_obj <- biomformat::make_biom(
    data = counts,
    observation_metadata = obs_meta,
    sample_metadata = sample_meta
  )
  
  return(biom_obj)
}