# Demonstration of using mock data generators for microFGT
# This script shows how to use the mock data generators for SpeciateIT, VIRGO, and VALENCIA

# Load the mock data generators
# When running from the package directory during development
if (file.exists("R/mock_data_generators.R")) {
  source("R/mock_data_generators.R")
} else {
  # In installed package, use the loadNamespace approach
  library(microFGT)
  tryCatch({
    # Attempt to access the functions directly
    if (!exists("generate_mock_speciateit")) {
      stop("Mock data generators not available in the loaded package")
    }
  }, error = function(e) {
    stop("Could not load mock data generators. Make sure you have the latest version of microFGT.")
  })
}

# Part 1: Generate individual mock datasets
# =========================================

# Example 1: Generate SpeciateIT mock data
# ----------------------------------------
cat("Example 1: Generating mock SpeciateIT data\n")
cat("------------------------------------------\n")

speciateit_data <- generate_mock_speciateit(
  n_sequences = 100,
  include_unclassified = 0.05,
  seed = 123
)

cat("Generated SpeciateIT data with", nrow(speciateit_data), "sequences\n")
cat("Sample of data:\n")
print(head(speciateit_data, 3))

# Species distribution
species_counts <- table(speciateit_data$Classification)
species_freq <- sort(prop.table(species_counts), decreasing = TRUE)
cat("\nSpecies distribution:\n")
print(head(species_freq, 5))

# Example 2: Generate VIRGO mock data
# -----------------------------------
cat("\n\nExample 2: Generating mock VIRGO data\n")
cat("-------------------------------------\n")

virgo_data <- generate_mock_virgo(
  n_genes = 500,
  n_samples = 5,
  simulation_type = "realistic",
  seed = 456
)

cat("Generated VIRGO data with", nrow(virgo_data$genes), "genes and", 
    ncol(virgo_data$counts), "samples\n")

# Display sample metadata
cat("\nSample metadata:\n")
print(virgo_data$metadata)

# Display counts distribution summary
counts_summary <- data.frame(
  Sample = colnames(virgo_data$counts),
  Total_Reads = colSums(virgo_data$counts),
  Non_Zero_Genes = colSums(virgo_data$counts > 0),
  Mean_Count = colMeans(virgo_data$counts),
  Sparsity = colSums(virgo_data$counts == 0) / nrow(virgo_data$counts)
)
cat("\nCounts summary:\n")
print(counts_summary)

# Example 3: Generate VALENCIA mock data
# --------------------------------------
cat("\n\nExample 3: Generating mock VALENCIA data\n")
cat("-----------------------------------------\n")

valencia_data <- generate_mock_valencia(
  n_samples = 20,
  include_scores = TRUE,
  add_taxonomic_data = TRUE,
  seed = 789
)

cat("Generated VALENCIA data with", nrow(valencia_data$cst), "samples\n")

# Display CST distribution
cst_dist <- table(valencia_data$cst$CST)
cat("\nCST distribution:\n")
print(cst_dist)

# Display score summary for a sample
sample_id <- valencia_data$cst$Sample[1]
cat("\nScore breakdown for sample", sample_id, ":\n")
scores <- valencia_data$scores[1, ]
scores_df <- data.frame(
  CST = sub("CST-", "", names(scores)),
  Score = as.numeric(scores)
)
scores_df <- scores_df[order(scores_df$Score, decreasing = TRUE), ]
print(scores_df)

# Display taxonomic data for the same sample
if (!is.null(valencia_data$abundance)) {
  cat("\nTaxonomic profile for sample", sample_id, ":\n")
  taxa_row <- valencia_data$abundance[valencia_data$abundance$Sample == sample_id, ]
  taxa_row <- taxa_row[, !colnames(taxa_row) %in% c("Sample", "read_count")]
  taxa_values <- sort(as.numeric(taxa_row), decreasing = TRUE)
  taxa_names <- names(sort(as.numeric(taxa_row), decreasing = TRUE))
  taxa_df <- data.frame(
    Taxon = taxa_names,
    Abundance = taxa_values
  )
  print(head(taxa_df, 5))
}

# Part 2: Convert to FGTMicrobiome compatible formats
# ==================================================
cat("\n\nPart 2: Converting to FGTMicrobiome compatible formats\n")
cat("======================================================\n")

# Convert SpeciateIT data
speciateit_fgt <- convert_speciateit_to_fgt(speciateit_data)
cat("Converted SpeciateIT data to FGTMicrobiome format\n")
cat("Structure of FGT-compatible SpeciateIT data:\n")
str(speciateit_fgt, max.level = 1)

# Convert VIRGO data
virgo_fgt <- convert_virgo_to_fgt(virgo_data)
cat("\nConverted VIRGO data to FGTMicrobiome format\n")
cat("Structure of FGT-compatible VIRGO data:\n")
str(virgo_fgt, max.level = 1)

# Convert VALENCIA data
valencia_fgt <- convert_valencia_to_fgt(valencia_data)
cat("\nConverted VALENCIA data to FGTMicrobiome format\n")
cat("Structure of FGT-compatible VALENCIA data:\n")
str(valencia_fgt, max.level = 1)

# Part 3: Generate a coordinated dataset for all tools
# ===================================================
cat("\n\nPart 3: Generating a coordinated dataset for all tools\n")
cat("====================================================\n")

# Create a temporary directory for output files
temp_dir <- tempfile()
dir.create(temp_dir)
cat("Created temporary directory for output files:", temp_dir, "\n")

# Generate the coordinated dataset
coordinated_data <- generate_mock_fgt_dataset(
  n_samples = 10,
  n_sequences = 200,
  n_genes = 300,
  create_files = TRUE,
  output_dir = temp_dir,
  seed = 101
)

cat("Generated coordinated dataset with", nrow(coordinated_data$metadata), "samples\n")

# Display metadata with CST assignments
cat("\nSample metadata with CST assignments:\n")
print(coordinated_data$metadata)

# Show output file structure
cat("\nOutput files created in:", temp_dir, "\n")
cat("Directory structure:\n")
dir_structure <- list.dirs(temp_dir, recursive = TRUE, full.names = FALSE)
for (dir in dir_structure) {
  if (dir == "") {
    cat("- [root directory]\n")
    files <- list.files(temp_dir, recursive = FALSE)
    if (length(files) > 0) {
      for (file in files) {
        cat("  - ", file, "\n", sep="")
      }
    }
  } else {
    cat("- ", dir, "\n", sep="")
    files <- list.files(file.path(temp_dir, dir), recursive = FALSE)
    if (length(files) > 0) {
      for (file in files) {
        cat("  - ", file, "\n", sep="")
      }
    }
  }
}

# Part 4: Working with the mock data in FGTMicrobiome
# ==================================================
cat("\n\nPart 4: Working with the mock data in FGTMicrobiome\n")
cat("==================================================\n")

cat("Example of how to use the mock data with FGTMicrobiome:\n\n")
cat("```r\n")
cat("library(microFGT)\n\n")
cat("# Load real data or create mock data\n")
cat("mock_data <- generate_mock_fgt_dataset(n_samples = 10, n_sequences = 200)\n\n")
cat("# Create a basic FGTMicrobiome object\n")
cat("# (In a real scenario, you would create this from your actual data)\n")
cat("counts <- matrix(sample(0:100, 200*10, replace=TRUE), nrow = 200, ncol = 10)\n")
cat("rownames(counts) <- mock_data$speciateit$raw$`Sequence ID`[1:200]\n")
cat("colnames(counts) <- mock_data$metadata$Sample\n\n")
cat("taxonomy <- data.frame(\n")
cat("  Feature = rownames(counts),\n")
cat("  Kingdom = rep(\"Bacteria\", 200),\n")
cat("  Phylum = rep(\"Firmicutes\", 200),\n")
cat("  stringsAsFactors = FALSE\n")
cat(")\n\n")
cat("metadata <- data.frame(\n")
cat("  Sample = mock_data$metadata$Sample,\n")
cat("  Subject = paste0(\"Subject_\", rep(1:5, each=2)),\n")
cat("  row.names = mock_data$metadata$Sample,\n")
cat("  stringsAsFactors = FALSE\n")
cat(")\n\n")
cat("# Create FGTMicrobiome object\n")
cat("fgt <- FGTMicrobiome(counts, taxonomy, metadata)\n\n")
cat("# Add mock tool results\n")
cat("fgt_with_speciate <- add_speciateit(fgt, mock_data$speciateit$fgt_compatible)\n")
cat("fgt_with_virgo <- add_virgo(fgt_with_speciate, mock_data$virgo$fgt_compatible)\n")
cat("fgt_with_valencia <- add_valencia(fgt_with_virgo, mock_data$valencia$fgt_compatible)\n\n")
cat("# The object now contains all three types of data\n")
cat("fgt_complete <- fgt_with_valencia\n")
cat("```\n")

# Clean up
cat("\nCleaning up temporary files...\n")
unlink(temp_dir, recursive = TRUE)
cat("Done.\n")