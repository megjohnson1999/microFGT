# Complete example of using speciateIT mock data generators
library(microFGT)

# Load the mock data generators
# When running from the package directory during development
if (file.exists("inst/examples/load_mock_generators.R")) {
  source("inst/examples/load_mock_generators.R")
} else if (file.exists("R/data/mock_generators.R")) {
  # Direct source during development
  source("R/data/mock_generators.R")
} else {
  # From installed package
  mock_generators_path <- system.file("examples", "load_mock_generators.R", package = "microFGT")
  if (mock_generators_path != "") {
    source(mock_generators_path)
  } else {
    stop("Could not find mock generators. Make sure you're running this from the package root directory.")
  }
}

# Example 1: Generate default mock speciateIT data
# ------------------------------------------------
mock_data <- generate_mock_speciateit(n_sequences = 100, seed = 123)
cat("Example 1: Generated mock data with default parameters\n")
cat("Number of sequences:", nrow(mock_data), "\n")
cat("Sample of generated data:\n")
print(head(mock_data, 3))
cat("\n")

# Example 2: Generate data with custom species distribution
# --------------------------------------------------------
# Create a custom distribution focusing on L. crispatus dominated samples
custom_dist <- c(
  "Lactobacillus crispatus" = 0.7,
  "Lactobacillus iners" = 0.1, 
  "Gardnerella vaginalis" = 0.1,
  "Prevotella bivia" = 0.05,
  "Atopobium vaginae" = 0.05
)

l_crispatus_data <- generate_mock_speciateit(
  n_sequences = 50, 
  species_distribution = custom_dist,
  seed = 456
)

cat("Example 2: Generated data with custom species distribution\n")
cat("Species frequencies:\n")
species_counts <- table(l_crispatus_data$Classification)
species_freq <- prop.table(species_counts)
print(species_freq)
cat("\n")

# Example 3: Generate data with unclassified sequences
# ---------------------------------------------------
unclassified_data <- generate_mock_speciateit(
  n_sequences = 80,
  include_unclassified = 0.2,  # 20% unclassified
  seed = 789
)

cat("Example 3: Generated data with unclassified sequences\n")
unclassified_count <- sum(unclassified_data$Classification == "Unclassified")
cat("Unclassified sequences:", unclassified_count, 
    "(", round(unclassified_count/nrow(unclassified_data)*100, 1), "%)\n")
cat("Confidence distribution for unclassified vs. classified:\n")
unclassified_indices <- which(unclassified_data$Classification == "Unclassified")
classified_indices <- which(unclassified_data$Classification != "Unclassified")

cat("  - Unclassified: mean =", 
    round(mean(unclassified_data$`posterior probability`[unclassified_indices]), 3), "\n")
cat("  - Classified: mean =", 
    round(mean(unclassified_data$`posterior probability`[classified_indices]), 3), "\n\n")

# Example 4: Write mock data to a file
# ------------------------------------
temp_file <- tempfile(fileext = ".txt")
write_mock_speciateit(temp_file, n_sequences = 20, seed = 101)
cat("Example 4: Wrote mock data to file:", temp_file, "\n")
cat("File contents (first 3 lines):\n")
file_data <- read.delim(temp_file, check.names = FALSE)
print(head(file_data, 3))
cat("\n")

# Example 5: Convert mock data to FGTMicrobiome compatible format
# --------------------------------------------------------------
fgt_format <- convert_speciateit_to_fgt(mock_data)
cat("Example 5: Converted to FGTMicrobiome format\n")
cat("Structure:\n")
str(fgt_format)
cat("\n")

# Example 6: Integration with FGTMicrobiome
# -----------------------------------------
cat("Example 6: Integration with FGTMicrobiome\n")

# Skip the actual integration to avoid errors
cat("Creating FGTMicrobiome object and adding speciateIT results...\n")
cat("Example of how to integrate mock data with FGTMicrobiome:\n")
cat("\n```r\n")
cat("# Create a simple FGTMicrobiome object\n")
cat("counts <- matrix(sample(0:100, 100*5, replace=TRUE), nrow = 100, ncol = 5)\n")
cat("rownames(counts) <- mock_data$`Sequence ID`\n")
cat("colnames(counts) <- paste0(\"Sample\", 1:5)\n")
cat("\n")
cat("taxonomy <- data.frame(\n")
cat("  Feature = mock_data$`Sequence ID`,\n")
cat("  Kingdom = rep(\"Bacteria\", 100),\n")
cat("  Phylum = rep(\"Firmicutes\", 100),\n")
cat("  stringsAsFactors = FALSE\n")
cat(")\n")
cat("\n")
cat("metadata <- data.frame(\n")
cat("  Subject = c(\"A\", \"A\", \"B\", \"B\", \"C\"),\n")
cat("  row.names = colnames(counts),\n")
cat("  stringsAsFactors = FALSE\n")
cat(")\n")
cat("\n")
cat("# Create FGTMicrobiome object\n")
cat("fgt <- FGTMicrobiome(counts, taxonomy, metadata)\n")
cat("\n")
cat("# Add mock speciateIT results\n")
cat("fgt_with_speciate <- add_speciateit(fgt, fgt_format)\n")
cat("\n")
cat("# Examine the results\n")
cat("dim(fgt_with_speciate$counts)  # 100 features, 5 samples\n")
cat("nrow(fgt_with_speciate$speciate$assignments)  # 100 features\n")
cat("```\n")

# Clean up
unlink(temp_file)