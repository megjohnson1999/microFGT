library(microFGT)

# ------------- IMPORTING DATA -------------

# Load example data that simulates realistic FGT microbiome profiles
fgt_exp <- load_example_data(size = "small", type = "amplicon")

# Save the example data to various formats
output_dir <- tempdir()
cat("Saving data to:", output_dir, "\n")

# Export to CSV format
csv_files <- export_microbiome(fgt_exp, output_dir, format = "csv", prefix = "microFGT_demo")
cat("CSV files exported to:", paste(unlist(csv_files[1:3]), collapse = ", "), "\n\n")

# Export to RDS format
rds_files <- export_microbiome(fgt_exp, output_dir, format = "rds", prefix = "microFGT_demo")
cat("RDS files exported to:", paste(unlist(rds_files[1:3]), collapse = ", "), "\n\n")

# ------------- REIMPORTING VARIOUS FORMATS -------------

# Import from CSV files
csv_se <- import_microbiome(
  counts = csv_files$counts,
  taxonomy = csv_files$taxonomy,
  sample_data = csv_files$metadata,
  as_FGTExperiment = FALSE
)
cat("Imported from CSV:", dim(csv_se)[1], "features,", dim(csv_se)[2], "samples\n")

# Import from RDS files
rds_se <- import_microbiome(
  counts = rds_files$counts,
  taxonomy = rds_files$taxonomy,
  sample_data = rds_files$metadata,
  as_FGTExperiment = FALSE
)
cat("Imported from RDS:", dim(rds_se)[1], "features,", dim(rds_se)[2], "samples\n\n")

# ------------- WORKING WITH TAXONOMIC DATA -------------

# Get available taxonomic ranks
ranks <- get_taxonomic_ranks(csv_se)
cat("Available taxonomic ranks:", paste(ranks, collapse = ", "), "\n")

# Normalize taxonomy
norm_se <- normalize_taxonomy(csv_se)

# Create taxonomic strings
tax_se <- create_tax_strings(norm_se, format = "lineage")
tax_strings <- rowData(tax_se)$taxonomy[1:3]
cat("Example taxonomic strings:\n", paste(tax_strings, collapse = "\n "), "\n\n")

# Export taxonomy strings to demonstrate reimporting
taxonomy_file <- file.path(output_dir, "tax_strings.txt")
writeLines(rowData(tax_se)$taxonomy, taxonomy_file)

# Reimport from taxonomy strings
tax_string_se <- import_microbiome(
  counts = csv_files$counts,
  taxonomy = readLines(taxonomy_file),
  sample_data = csv_files$metadata,
  as_FGTExperiment = FALSE
)
cat("Imported from taxonomy strings:", dim(tax_string_se)[1], "features,", dim(tax_string_se)[2], "samples\n\n")

# ------------- IMPORTING FROM DADA2 FORMAT -------------

# Export to DADA2 format
dada2_files <- export_microbiome(fgt_exp, output_dir, format = "dada2", prefix = "microFGT_demo_dada2")

# Import from DADA2 format
dada2_se <- import_from_dada2(
  seqtab = dada2_files$seqtab,
  taxa = dada2_files$taxa,
  metadata = dada2_files$metadata
)
cat("Imported from DADA2 format:", dim(dada2_se)[1], "features,", dim(dada2_se)[2], "samples\n\n")

# ------------- DEMONSTRATES CONVERSION FUNCTIONS -------------

# Only run if phyloseq is available
if (requireNamespace("phyloseq", quietly = TRUE)) {
  # Convert to phyloseq
  ps <- to_phyloseq(norm_se)
  cat("Converted to phyloseq object with", phyloseq::ntaxa(ps), "taxa and", 
      phyloseq::nsamples(ps), "samples\n")
  
  # Export to phyloseq format
  phyloseq_file <- export_microbiome(norm_se, output_dir, format = "phyloseq", 
                                    prefix = "microFGT_demo_phyloseq")
  cat("Exported to phyloseq format:", phyloseq_file$phyloseq, "\n")
  
  # Reimport from phyloseq
  ps_se <- import_from_phyloseq(ps)
  cat("Reimported from phyloseq:", dim(ps_se)[1], "features,", dim(ps_se)[2], "samples\n\n")
} else {
  cat("Phyloseq package not available, skipping phyloseq conversion examples\n\n")
}

# ------------- DEMONSTRATES FILE FORMAT DETECTION -------------

# Create a count matrix with some special properties
counts <- matrix(rpois(60, lambda = 20), nrow = 10, ncol = 6)
rownames(counts) <- paste0("ASV", 1:10)
colnames(counts) <- paste0("Sample", 1:6)

# Save to CSV with non-standard name
unusual_file <- file.path(output_dir, "unusual_format.data")
write.csv(counts, unusual_file)

# Import should still work due to format detection
unusual_se <- import_microbiome(unusual_file)
cat("Imported from unusual file format:", dim(unusual_se)[1], "features,", dim(unusual_se)[2], "samples\n\n")

# ------------- DEMONSTRATES MATRIX TRANSPOSITION -------------

# Create a transposed count matrix (samples in rows)
transposed_counts <- t(counts)
transposed_file <- file.path(output_dir, "transposed_counts.csv")
write.csv(transposed_counts, transposed_file)

# Import should detect and correct the orientation
transposed_se <- import_microbiome(transposed_file)
cat("Imported from transposed matrix:", dim(transposed_se)[1], "features,", dim(transposed_se)[2], "samples\n\n")

cat("Import/Export demo completed successfully!\n")