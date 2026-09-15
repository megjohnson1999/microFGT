test_that("import_microbiome correctly handles count matrices", {
  # Create test data
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- c("OTU1", "OTU2", "OTU3")
  colnames(counts) <- c("Sample1", "Sample2", "Sample3", "Sample4")
  
  # Import with just counts
  se <- import_microbiome(counts, as_FGTExperiment = FALSE)
  
  # Check dimensions
  expect_equal(dim(se), c(3, 4))
  
  # Check that count matrix is preserved
  expect_equal(SummarizedExperiment::assays(se)$counts, counts)
  
  # Check rownames and colnames
  expect_equal(rownames(se), rownames(counts))
  expect_equal(colnames(se), colnames(counts))
})

test_that("import_microbiome correctly handles taxonomy and metadata", {
  # Skip if S4Vectors is not available
  skip_if_not_installed("S4Vectors")
  
  # Create test data
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- c("OTU1", "OTU2", "OTU3")
  colnames(counts) <- c("Sample1", "Sample2", "Sample3", "Sample4")
  
  taxonomy <- data.frame(
    Kingdom = c("Bacteria", "Bacteria", "Bacteria"),
    Phylum = c("Firmicutes", "Bacteroidetes", "Proteobacteria"),
    Genus = c("Lactobacillus", "Bacteroides", "Escherichia"),
    row.names = rownames(counts)
  )
  
  metadata <- data.frame(
    Group = c("Control", "Control", "Treatment", "Treatment"),
    row.names = colnames(counts)
  )
  
  # Import with counts, taxonomy, and metadata
  se <- import_microbiome(counts, taxonomy, metadata, as_FGTExperiment = FALSE)
  
  # Check dimensions
  expect_equal(dim(se), c(3, 4))
  
  # Check that taxonomy is correctly stored
  row_data <- SummarizedExperiment::rowData(se)
  expect_equal(row_data$Kingdom, taxonomy$Kingdom)
  expect_equal(row_data$Phylum, taxonomy$Phylum)
  expect_equal(row_data$Genus, taxonomy$Genus)
  
  # Check that metadata is correctly stored
  col_data <- SummarizedExperiment::colData(se)
  expect_equal(col_data$Group, metadata$Group)
})

test_that("import_microbiome correctly handles taxonomic strings", {
  # Create test data
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- c("OTU1", "OTU2", "OTU3")
  colnames(counts) <- c("Sample1", "Sample2", "Sample3", "Sample4")
  
  # Create taxonomic strings
  tax_strings <- c(
    "k__Bacteria;p__Firmicutes;g__Lactobacillus",
    "k__Bacteria;p__Bacteroidetes;g__Bacteroides",
    "k__Bacteria;p__Proteobacteria;g__Escherichia"
  )
  names(tax_strings) <- rownames(counts)
  
  # Import with taxonomic strings
  se <- import_microbiome(counts, tax_strings, as_FGTExperiment = FALSE)
  
  # Check dimensions
  expect_equal(dim(se), c(3, 4))
  
  # Check that taxonomy is correctly parsed
  row_data <- SummarizedExperiment::rowData(se)
  expect_equal(as.character(row_data$Kingdom), rep("Bacteria", 3))
  expect_equal(as.character(row_data$Phylum), c("Firmicutes", "Bacteroidetes", "Proteobacteria"))
  expect_equal(as.character(row_data$Genus), c("Lactobacillus", "Bacteroides", "Escherichia"))
})

test_that("parse_taxonomy_strings correctly parses taxonomic strings", {
  # Test with standard prefixes
  tax_strings <- c(
    "k__Bacteria;p__Firmicutes;g__Lactobacillus;s__L. crispatus",
    "k__Bacteria;p__Bacteroidetes;g__Bacteroides"
  )
  
  parsed <- parse_taxonomy_strings(tax_strings)
  
  expect_equal(nrow(parsed), 2)
  expect_equal(ncol(parsed), 7)  # Standard taxonomy ranks
  expect_equal(as.character(parsed$Kingdom), c("Bacteria", "Bacteria"))
  expect_equal(as.character(parsed$Phylum), c("Firmicutes", "Bacteroidetes"))
  expect_equal(as.character(parsed$Genus), c("Lactobacillus", "Bacteroides"))
  expect_equal(as.character(parsed$Species), c("L. crispatus", NA))
  
  # Test without prefixes
  tax_strings2 <- c(
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;L. crispatus",
    "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides"
  )
  
  parsed2 <- parse_taxonomy_strings(tax_strings2)
  
  expect_equal(nrow(parsed2), 2)
  expect_equal(ncol(parsed2), 7)
  expect_equal(as.character(parsed2$Kingdom), c("Bacteria", "Bacteria"))
  expect_equal(as.character(parsed2$Phylum), c("Firmicutes", "Bacteroidetes"))
})

test_that("export_microbiome correctly exports to CSV", {
  # Skip the test if the required packages are not available
  skip_if_not_installed("SummarizedExperiment")
  
  # Create test data
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- c("OTU1", "OTU2", "OTU3")
  colnames(counts) <- c("Sample1", "Sample2", "Sample3", "Sample4")
  
  taxonomy <- data.frame(
    Kingdom = c("Bacteria", "Bacteria", "Bacteria"),
    Phylum = c("Firmicutes", "Bacteroidetes", "Proteobacteria"),
    Genus = c("Lactobacillus", "Bacteroides", "Escherichia"),
    row.names = rownames(counts)
  )
  
  metadata <- data.frame(
    Group = c("Control", "Control", "Treatment", "Treatment"),
    row.names = colnames(counts)
  )
  
  # Create a SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = taxonomy,
    colData = metadata
  )
  
  # Set up temporary directory
  temp_dir <- tempdir()
  
  # Export to CSV
  result <- export_microbiome(se, temp_dir, format = "csv", prefix = "test")
  
  # Check that files were created
  expect_true(file.exists(result$counts))
  expect_true(file.exists(result$taxonomy))
  expect_true(file.exists(result$metadata))
  
  # Read the exported files
  counts_csv <- read.csv(result$counts, row.names = 1)
  taxonomy_csv <- read.csv(result$taxonomy, row.names = 1)
  metadata_csv <- read.csv(result$metadata, row.names = 1)
  
  # Check that data was preserved
  expect_equal(as.matrix(counts_csv), counts)
  expect_equal(taxonomy_csv$Kingdom, taxonomy$Kingdom)
  expect_equal(taxonomy_csv$Phylum, taxonomy$Phylum)
  expect_equal(taxonomy_csv$Genus, taxonomy$Genus)
  expect_equal(metadata_csv$Group, metadata$Group)
  
  # Clean up
  unlink(result$counts)
  unlink(result$taxonomy)
  unlink(result$metadata)
})

test_that("export_microbiome correctly exports to RDS", {
  # Skip the test if the required packages are not available
  skip_if_not_installed("SummarizedExperiment")
  
  # Create test data
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts) <- c("OTU1", "OTU2", "OTU3")
  colnames(counts) <- c("Sample1", "Sample2", "Sample3", "Sample4")
  
  taxonomy <- data.frame(
    Kingdom = c("Bacteria", "Bacteria", "Bacteria"),
    Phylum = c("Firmicutes", "Bacteroidetes", "Proteobacteria"),
    Genus = c("Lactobacillus", "Bacteroides", "Escherichia"),
    row.names = rownames(counts)
  )
  
  metadata <- data.frame(
    Group = c("Control", "Control", "Treatment", "Treatment"),
    row.names = colnames(counts)
  )
  
  # Create a SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = taxonomy,
    colData = metadata
  )
  
  # Set up temporary directory
  temp_dir <- tempdir()
  
  # Export to RDS
  result <- export_microbiome(se, temp_dir, format = "rds", prefix = "test")
  
  # Check that files were created
  expect_true(file.exists(result$counts))
  expect_true(file.exists(result$taxonomy))
  expect_true(file.exists(result$metadata))
  
  # Read the exported files
  counts_rds <- readRDS(result$counts)
  taxonomy_rds <- readRDS(result$taxonomy)
  metadata_rds <- readRDS(result$metadata)
  
  # Check that data was preserved
  expect_equal(counts_rds, counts)
  expect_equal(as.data.frame(taxonomy_rds)$Kingdom, taxonomy$Kingdom)
  expect_equal(as.data.frame(taxonomy_rds)$Phylum, taxonomy$Phylum)
  expect_equal(as.data.frame(taxonomy_rds)$Genus, taxonomy$Genus)
  expect_equal(as.data.frame(metadata_rds)$Group, metadata$Group)
  
  # Clean up
  unlink(result$counts)
  unlink(result$taxonomy)
  unlink(result$metadata)
})

test_that("import_from_dada2 correctly handles DADA2 outputs", {
  # Skip the test if the required packages are not available
  skip_if_not_installed("SummarizedExperiment")
  
  # Create test data that mimics DADA2 outputs
  seqtab <- matrix(1:12, nrow = 4, ncol = 3)
  colnames(seqtab) <- c("Sample1", "Sample2", "Sample3")
  rownames(seqtab) <- c("ACGT", "ACGTA", "ACGTC", "ACGTG")
  
  taxa <- matrix(NA, nrow = 4, ncol = 7)
  colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(taxa) <- rownames(seqtab)
  
  taxa[, "Kingdom"] <- "Bacteria"
  taxa[, "Phylum"] <- c("Firmicutes", "Firmicutes", "Bacteroidetes", "Proteobacteria")
  taxa[, "Genus"] <- c("Lactobacillus", "Lactobacillus", "Bacteroides", "Escherichia")
  
  metadata <- data.frame(
    Group = c("Control", "Control", "Treatment"),
    row.names = colnames(seqtab)
  )
  
  # Import from DADA2
  se <- import_from_dada2(seqtab, taxa, metadata, as_FGTExperiment = FALSE)
  
  # Check dimensions
  expect_equal(dim(se), c(4, 3))
  
  # Check that data was preserved
  counts <- SummarizedExperiment::assays(se)$counts
  row_data <- SummarizedExperiment::rowData(se)
  col_data <- SummarizedExperiment::colData(se)
  
  expect_equal(dim(counts), c(4, 3))
  expect_equal(row_data$Kingdom, rep("Bacteria", 4))
  expect_equal(row_data$Phylum, c("Firmicutes", "Firmicutes", "Bacteroidetes", "Proteobacteria"))
  expect_equal(row_data$Genus, c("Lactobacillus", "Lactobacillus", "Bacteroides", "Escherichia"))
  expect_equal(col_data$Group, c("Control", "Control", "Treatment"))
})

test_that("import_microbiome correctly transpooses when needed", {
  # Create test data with features in columns (needs transpose)
  counts <- matrix(1:12, nrow = 4, ncol = 3)
  rownames(counts) <- c("Sample1", "Sample2", "Sample3", "Sample4")
  colnames(counts) <- c("Feature1", "Feature2", "Feature3")
  
  # Import with transpose needed
  se <- import_microbiome(counts, as_FGTExperiment = FALSE)
  
  # Check dimensions - should be transposed
  expect_equal(dim(se), c(3, 4))
  
  # Check that count matrix is transposed
  expect_equal(SummarizedExperiment::assays(se)$counts, t(counts))
  
  # Check rownames and colnames
  expect_equal(rownames(se), colnames(counts))
  expect_equal(colnames(se), rownames(counts))
})