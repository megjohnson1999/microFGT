#' Check if speciateIT is available
#'
#' @return Logical indicating if speciateIT is installed
#' @export
is_speciateit_available <- function() {
  cmd_output <- tryCatch(
    system("speciateit --version", intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character(0),
    warning = function(w) character(0)
  )
  
  return(length(cmd_output) > 0)
}

#' Get default speciateIT database path
#'
#' @return Path to the default speciateIT database
#' @export
default_speciateit_db <- function() {
  # Check if custom path is set in options
  custom_path <- getOption("microFGT.speciateit_db")
  if (!is.null(custom_path) && file.exists(custom_path)) {
    return(custom_path)
  }
  
  # Check standard locations
  std_paths <- c(
    file.path(system.file(package = "microFGT"), "extdata", "reference_dbs", "speciateit"),
    "~/speciateit_db",
    "/usr/local/share/speciateit_db",
    "/usr/share/speciateit_db"
  )
  
  for (path in std_paths) {
    if (file.exists(path)) {
      return(path)
    }
  }
  
  return(NULL)
}

#' Run speciateIT classification on sequence data
#'
#' @param fasta_file Path to FASTA file with sequences
#' @param output_file Path to output file (will be created)
#' @param reference_db Path to reference database
#' @param threads Number of threads to use
#' @param other_params Additional parameters for speciateIT
#'
#' @return Path to output file
#' @export
run_speciateit <- function(fasta_file, output_file = NULL, 
                          reference_db = default_speciateit_db(),
                          threads = 1, other_params = "") {
  
  # Check if speciateIT is installed
  if (!is_speciateit_available()) {
    stop("speciateIT is not installed or not in PATH")
  }
  
  # Check if reference database exists
  if (is.null(reference_db) || !file.exists(reference_db)) {
    stop("speciateIT reference database not found")
  }
  
  # Create output file name if not provided
  if (is.null(output_file)) {
    output_file <- paste0(tools::file_path_sans_ext(fasta_file), "_speciateit.tsv")
  }
  
  # Build command
  cmd <- paste0(
    "speciateit",
    " --input ", shQuote(fasta_file),
    " --output ", shQuote(output_file),
    " --refdb ", shQuote(reference_db),
    " --threads ", threads,
    " ", other_params
  )
  
  # Run command
  message("Running speciateIT: ", cmd)
  system_result <- system(cmd, intern = TRUE)
  
  # Check if output file was created
  if (!file.exists(output_file)) {
    stop("speciateIT failed to create output file")
  }
  
  return(output_file)
}

#' Write sequences to a FASTA file
#'
#' @param sequences Character vector of sequences or DNAStringSet
#' @param file Path to output FASTA file
#' @param names Names for the sequences (if NULL, uses names from sequences or generates them)
#'
#' @return Path to the created FASTA file
#' @keywords internal
write_fasta <- function(sequences, file, names = NULL) {
  if (requireNamespace("Biostrings", quietly = TRUE) && methods::is(sequences, "DNAStringSet")) {
    # Use Biostrings to write DNAStringSet
    Biostrings::writeXStringSet(sequences, file)
  } else {
    # Convert to character vector if needed
    sequences <- as.character(sequences)
    
    # Generate names if not provided
    if (is.null(names)) {
      if (!is.null(names(sequences))) {
        names <- names(sequences)
      } else {
        names <- paste0("seq", seq_along(sequences))
      }
    }
    
    # Create FASTA content
    fasta_lines <- character(length(sequences) * 2)
    line_idx <- 1
    
    for (i in seq_along(sequences)) {
      fasta_lines[line_idx] <- paste0(">", names[i])
      fasta_lines[line_idx + 1] <- sequences[i]
      line_idx <- line_idx + 2
    }
    
    # Write to file
    writeLines(fasta_lines, file)
  }
  
  return(file)
}

#' Apply speciateIT classification to an FGTExperiment
#'
#' Runs speciateIT taxonomic classifier on sequences in an FGTExperiment object
#' and updates the taxonomy information.
#'
#' @param fgt_exp An FGTExperiment object
#' @param sequence_col Column name in rowData containing sequences
#' @param reference_db Path to reference database
#' @param threads Number of threads to use
#'
#' @return FGTExperiment with updated taxonomy
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming fgt_exp has sequences in rowData$sequence
#' classified_exp <- classify_with_speciateit(fgt_exp)
#' }
classify_with_speciateit <- function(fgt_exp, sequence_col = "sequence",
                                    reference_db = default_speciateit_db(),
                                    threads = 1) {
  # Input validation
  if (!methods::is(fgt_exp, "FGTExperiment")) {
    stop("Input must be an FGTExperiment object")
  }
  
  # Check if speciateIT is available
  if (!is_speciateit_available()) {
    stop("speciateIT is not installed or not in PATH. Install it before using this function.")
  }
  
  # Extract sequences
  row_data <- SummarizedExperiment::rowData(fgt_exp)
  sequences <- row_data[[sequence_col]]
  
  if (is.null(sequences)) {
    stop("Sequence column not found in rowData. Specify the correct column name with sequence_col parameter.")
  }
  
  # Create temporary FASTA file
  temp_fasta <- tempfile(fileext = ".fasta")
  
  # Write sequences to FASTA file
  write_fasta(sequences, temp_fasta, names = rownames(row_data))
  
  # Run speciateIT
  temp_out <- tempfile(fileext = ".tsv")
  speciateit_out <- run_speciateit(
    fasta_file = temp_fasta,
    output_file = temp_out,
    reference_db = reference_db,
    threads = threads
  )
  
  # Read results
  tax_results <- read.delim(speciateit_out, stringsAsFactors = FALSE)
  
  # Match results with sequences in FGTExperiment
  # Assuming the first column in the speciateIT output is the sequence ID
  id_col <- colnames(tax_results)[1]
  tax_results <- tax_results[match(rownames(row_data), tax_results[[id_col]]), ]
  
  # Update rowData with taxonomic information
  # Assuming speciateIT output has standard taxonomy columns
  tax_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  for (col in tax_cols) {
    if (col %in% colnames(tax_results)) {
      row_data[[col]] <- tax_results[[col]]
    }
  }
  
  # Add confidence scores if available
  conf_cols <- grep("confidence$|score$", colnames(tax_results), value = TRUE)
  for (col in conf_cols) {
    row_data[[col]] <- tax_results[[col]]
  }
  
  # Update the FGTExperiment
  SummarizedExperiment::rowData(fgt_exp) <- row_data
  
  # Add classification info to metadata
  metadata(fgt_exp)$classification_info <- list(
    tool = "speciateIT",
    date = Sys.time(),
    database = reference_db,
    parameters = list(threads = threads)
  )
  
  # Clean up temporary files
  unlink(c(temp_fasta, temp_out))
  
  return(fgt_exp)
}