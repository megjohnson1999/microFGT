## Example of generating and using mock SpeciateIT data
library(microFGT)

# Generate mock SpeciateIT data with default settings
mock_data <- generate_mock_speciateit(n_sequences = 100)
head(mock_data)

# Create a custom species distribution
custom_distribution <- c(
  "Lactobacillus crispatus" = 0.6,
  "Lactobacillus iners" = 0.3,
  "Gardnerella vaginalis" = 0.1
)

# Generate data with custom distribution and other parameters
custom_mock <- generate_mock_speciateit(
  n_sequences = 50,
  species_distribution = custom_distribution,
  include_unclassified = 0.1,
  confidence_range = c(0.8, 1.0),
  seed = 123
)
head(custom_mock)

# Convert to FGTMicrobiome compatible format
fgt_compatible <- convert_speciateit_to_fgt(custom_mock)
str(fgt_compatible)

# Create a simple FGTMicrobiome object
counts <- matrix(sample(0:100, 50 * 5), nrow = 50, ncol = 5)
rownames(counts) <- custom_mock[["Sequence ID"]]
colnames(counts) <- paste0("Sample", 1:5)

taxonomy <- data.frame(
  Feature = rownames(counts),
  Kingdom = rep("Bacteria", 50),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  Subject = c("A", "A", "B", "B", "C"),
  row.names = colnames(counts),
  stringsAsFactors = FALSE
)

# Create FGTMicrobiome object
fgt <- FGTMicrobiome(counts, taxonomy, metadata)

# Add mock SpeciateIT results
fgt_with_speciate <- add_speciateit(fgt, fgt_compatible)

# View the results
str(fgt_with_speciate$speciate)

# Write mock data to a file
temp_file <- tempfile(fileext = ".txt")
write_mock_speciateit(temp_file, n_sequences = 20)
cat("Mock data written to:", temp_file, "\n")

# Read the file back
file_data <- read.delim(temp_file, check.names = FALSE)
head(file_data)

# Clean up
unlink(temp_file)