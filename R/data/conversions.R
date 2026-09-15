#' Convert phyloseq to FGTExperiment
#'
#' Converts a phyloseq object to an FGTExperiment object.
#'
#' @param physeq A phyloseq object
#' @param assay_name Name for the count assay
#'
#' @return An FGTExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Load phyloseq object
#' data(GlobalPatterns, package = "phyloseq")
#' 
#' # Convert to FGTExperiment
#' fgt_exp <- phyloseq_to_fgt(GlobalPatterns)
#' }
phyloseq_to_fgt <- function(physeq, assay_name = "counts") {
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
  
  # Create assays list
  assays_list <- list()
  assays_list[[assay_name]] <- as.matrix(otu_table)
  
  # Convert taxonomy table to rowData
  if (!is.null(tax_table)) {
    rowData_df <- as.data.frame(tax_table)
  } else {
    rowData_df <- data.frame(row.names = rownames(otu_table))
  }
  
  # Convert sample_data to colData
  if (!is.null(sample_data)) {
    colData_df <- as.data.frame(sample_data)
  } else {
    colData_df <- data.frame(row.names = colnames(otu_table))
  }
  
  # Create FGTExperiment
  fgt_exp <- FGTExperiment(
    assays = assays_list,
    rowData = rowData_df,
    colData = colData_df,
    rowTree = phy_tree,
    experimentType = "amplicon"
  )
  
  return(fgt_exp)
}

#' Import dada2 results into FGTExperiment
#'
#' Creates an FGTExperiment object from dada2 sequence table and taxonomy results.
#'
#' @param seqtab Sequence table from dada2 (or path to RDS file containing it)
#' @param taxa Taxonomy table from dada2 (or path to RDS file containing it)
#' @param sample_data Sample metadata data frame (or path to CSV/RDS file)
#' @param refseq Optional DNAStringSet of reference sequences
#'
#' @return An FGTExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Using data from dada2 workflow
#' seqtab <- readRDS("seqtab.rds")
#' taxa <- readRDS("taxa.rds")
#' metadata <- read.csv("metadata.csv", row.names = 1)
#' 
#' fgt_exp <- import_dada2(seqtab, taxa, metadata)
#' }
import_dada2 <- function(seqtab, taxa = NULL, sample_data = NULL, refseq = NULL) {
  # Handle file path inputs for seqtab
  if (is.character(seqtab) && file.exists(seqtab)) {
    seqtab <- readRDS(seqtab)
  }
  
  # Handle file path inputs for taxa
  if (!is.null(taxa) && is.character(taxa) && file.exists(taxa)) {
    taxa <- readRDS(taxa)
  }
  
  # Handle file path inputs for sample_data
  if (!is.null(sample_data) && is.character(sample_data) && file.exists(sample_data)) {
    ext <- tools::file_ext(sample_data)
    if (ext == "rds") {
      sample_data <- readRDS(sample_data)
    } else if (ext %in% c("csv", "txt", "tsv")) {
      sep <- ifelse(ext == "tsv", "\t", ",")
      sample_data <- read.table(sample_data, header = TRUE, sep = sep, 
                                row.names = 1, stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file format for sample_data. Use CSV, TSV, or RDS.")
    }
  }
  
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
  
  # Process taxonomy table
  if (!is.null(taxa)) {
    # Convert to data frame if matrix
    if (is.matrix(taxa)) {
      taxa <- as.data.frame(taxa)
    }
    
    # Ensure sequence IDs match
    if (!all(rownames(seqtab) %in% rownames(taxa))) {
      warning("Not all sequences in seqtab have taxonomy assignments")
      # Subset to matching sequences
      shared_seqs <- intersect(rownames(seqtab), rownames(taxa))
      seqtab <- seqtab[shared_seqs, , drop = FALSE]
      taxa <- taxa[shared_seqs, , drop = FALSE]
    }
    
    rowData_df <- taxa
  } else {
    # Create empty rowData
    rowData_df <- data.frame(row.names = rownames(seqtab))
  }
  
  # Add sequences to rowData
  rowData_df$sequence <- rownames(seqtab)
  
  # Process sample metadata
  if (!is.null(sample_data)) {
    if (!all(colnames(seqtab) %in% rownames(sample_data))) {
      warning("Not all samples in seqtab have metadata")
      # Subset to matching samples
      shared_samples <- intersect(colnames(seqtab), rownames(sample_data))
      seqtab <- seqtab[, shared_samples, drop = FALSE]
    }
    colData_df <- sample_data[colnames(seqtab), , drop = FALSE]
  } else {
    # Create empty colData
    colData_df <- data.frame(row.names = colnames(seqtab))
  }
  
  # Create reference sequence object if provided
  if (!is.null(refseq)) {
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop("Package 'Biostrings' is required for handling reference sequences")
    }
    
    if (!methods::is(refseq, "DNAStringSet")) {
      refseq <- Biostrings::DNAStringSet(refseq)
      names(refseq) <- rownames(seqtab)
    }
    
    # Store reference sequences in fgtMetadata
    fgt_meta <- S4Vectors::SimpleList(referenceSeq = refseq)
  } else {
    fgt_meta <- S4Vectors::SimpleList()
  }
  
  # Create FGTExperiment
  fgt_exp <- FGTExperiment(
    assays = list(counts = seqtab),
    rowData = rowData_df,
    colData = colData_df,
    experimentType = "amplicon"
  )
  
  # Add reference sequences if available
  if (!is.null(refseq)) {
    fgtMetadata(fgt_exp) <- fgt_meta
  }
  
  return(fgt_exp)
}

#' Create a multi-omic FGT experiment
#'
#' Creates a MultiAssayExperiment object containing multiple FGTExperiment objects
#' for integrated analysis of different data types.
#'
#' @param ... Named FGTExperiment objects
#' @param sample_id_column Name of the column in colData containing sample identifiers
#' @param sample_relationships Optional data frame linking samples across experiments
#'
#' @return A MultiAssayExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create multi-omic experiment from amplicon and metagenomic data
#' multi_exp <- create_multiomic_fgt(
#'   amplicon = amplicon_fgt_exp,
#'   metagenomic = metagenomic_fgt_exp,
#'   sample_id_column = "subject_id"
#' )
#' }
create_multiomic_fgt <- function(..., sample_id_column = NULL, sample_relationships = NULL) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("Package 'MultiAssayExperiment' is required for this function")
  }
  
  # Collect the named experiments
  experiments <- list(...)
  
  # Check if all are FGTExperiment objects
  if (!all(vapply(experiments, function(x) is(x, "FGTExperiment"), logical(1)))) {
    stop("All inputs must be FGTExperiment objects")
  }
  
  # Extract common sample metadata if no relationships provided
  if (is.null(sample_relationships)) {
    if (is.null(sample_id_column)) {
      # Use row names as sample IDs
      sample_map <- lapply(experiments, function(x) {
        data.frame(
          primary = rownames(colData(x)),
          colname = rownames(colData(x)),
          stringsAsFactors = FALSE
        )
      })
      
      # Create sample map for MultiAssayExperiment
      sample_map <- do.call(rbind, sample_map)
      rownames(sample_map) <- NULL
    } else {
      # Use specified column as sample IDs
      sample_map <- lapply(names(experiments), function(name) {
        x <- experiments[[name]]
        if (!sample_id_column %in% colnames(colData(x))) {
          stop(paste("Column", sample_id_column, "not found in colData of", name))
        }
        
        data.frame(
          primary = colData(x)[[sample_id_column]],
          colname = rownames(colData(x)),
          stringsAsFactors = FALSE
        )
      })
      
      # Add assay names and combine
      for (i in seq_along(sample_map)) {
        sample_map[[i]]$assay <- names(experiments)[i]
      }
      
      sample_map <- do.call(rbind, sample_map)
      rownames(sample_map) <- NULL
    }
  } else {
    # Use provided sample relationships
    sample_map <- sample_relationships
  }
  
  # Create colData for MultiAssayExperiment
  # Combine metadata from all experiments
  all_sample_ids <- unique(sample_map$primary)
  
  # Create MultiAssayExperiment
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experiments,
    colData = data.frame(row.names = all_sample_ids),
    sampleMap = sample_map
  )
  
  return(mae)
}