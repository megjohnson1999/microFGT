#' Import data from dada2
#'
#' Creates an FGTExperiment object from dada2 output.
#'
#' @param seqtab Sequence table from dada2 (ASVs x samples)
#' @param taxa Taxonomy table from dada2
#' @param metadata Sample metadata (optional)
#' @param tree Phylogenetic tree (optional)
#'
#' @return An FGTExperiment object
#' @export
#'
#' @examples
#' \dontrun{
#' # Load data from dada2 outputs
#' seqtab <- readRDS("seqtab.rds")
#' taxa <- readRDS("taxa.rds")
#' metadata <- read.csv("metadata.csv", row.names = 1)
#' 
#' # Create FGTExperiment object
#' fgt_exp <- import_dada2(seqtab, taxa, metadata)
#' }
import_dada2 <- function(seqtab, taxa, metadata = NULL, tree = NULL) {
  # Validate seqtab
  if (!is.matrix(seqtab) && !is.data.frame(seqtab)) {
    stop("seqtab must be a matrix or data frame")
  }
  
  # Validate taxa
  if (!is.matrix(taxa) && !is.data.frame(taxa)) {
    stop("taxa must be a matrix or data frame")
  }
  
  # Convert to data.frame if needed
  if (is.matrix(seqtab)) {
    seqtab <- as.data.frame(seqtab)
  }
  
  if (is.matrix(taxa)) {
    taxa <- as.data.frame(taxa)
  }
  
  # Extract ASV sequences and create IDs
  asv_seqs <- colnames(seqtab)
  if (all(nchar(asv_seqs) > 20)) {
    # These are likely DNA sequences, create shorter IDs
    asv_ids <- paste0("ASV", seq_along(asv_seqs))
  } else {
    # Use existing IDs
    asv_ids <- asv_seqs
  }
  
  # Transpose sequence table to have ASVs as rows
  count_matrix <- t(seqtab)
  rownames(count_matrix) <- asv_ids
  
  # Prepare taxonomy data
  taxa_df <- as.data.frame(taxa)
  
  # Set standard taxonomic rank names if not present
  if (is.null(colnames(taxa_df))) {
    colnames(taxa_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(taxa_df)]
  }
  
  # Ensure rownames match ASV IDs
  if (nrow(taxa_df) != length(asv_ids)) {
    warning("Number of taxa does not match number of ASVs. Taxa information may be incomplete.")
    # Subset to matching ASVs
    taxa_df <- taxa_df[1:min(nrow(taxa_df), length(asv_ids)), , drop = FALSE]
  }
  rownames(taxa_df) <- asv_ids[1:nrow(taxa_df)]
  
  # Add ASV sequences if they're actual sequences
  if (all(nchar(asv_seqs) > 20)) {
    taxa_df$ASVSeq <- asv_seqs
  }
  
  # Create FGTExperiment
  fgt_exp <- FGTExperiment(
    assays = list(counts = count_matrix),
    rowData = taxa_df,
    colData = metadata,
    rowTree = tree,
    experimentType = "amplicon"
  )
  
  # Add import information to metadata
  metadata(fgt_exp)$import_info <- list(
    source = "dada2",
    date = Sys.time(),
    n_features = nrow(count_matrix),
    n_samples = ncol(count_matrix)
  )
  
  return(fgt_exp)
}