# prototype/importers.R
# Importers rewritten against REAL tool output shapes (see real_fixtures/FORMATS.md).
# Each does only the FGT-specific reshape; the container is generic Bioconductor.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(TreeSummarizedExperiment)
})

# --- speciateIT --------------------------------------------------------------
# REAL shape: speciateIT classifies ASVs, not samples. `results_path`
# (MC_order7_results.txt) is one row per ASV: Sequence ID / Classification /
# posterior probability / number of Decisions. Sample identity is NOT here -- it
# lives in the ASV count table (the dada2 feature table: rows = sampleID,
# cols = ASVs). taxon x sample = join(ASV->Classification, ASV x sample counts).
# GLUE: recover the join speciateIT can't express, then aggregate ASV->taxon.
import_speciateit <- function(results_path, count_table_path,
                              unclassified_label = "Unclassified") {
  res <- read.delim(results_path, check.names = FALSE, stringsAsFactors = FALSE)
  asv2taxon <- setNames(res[["Classification"]], res[["Sequence ID"]])

  ct <- as.matrix(read.csv(count_table_path, row.names = 1, check.names = FALSE)) # samples x ASVs
  taxon <- asv2taxon[colnames(ct)]
  taxon[is.na(taxon)] <- unclassified_label

  counts <- rowsum(t(ct), group = taxon)   # ASVs x samples -> taxa x samples
  storage.mode(counts) <- "integer"
  taxonomy <- DataFrame(Classification = rownames(counts),
                        Genus = sub(" .*$", "", rownames(counts)),
                        row.names = rownames(counts))
  TreeSummarizedExperiment(assays = list(counts = counts), rowData = taxonomy)
}

# --- VIRGO -------------------------------------------------------------------
# REAL shape: one file PER SAMPLE, <sample>.out, NO header, 3 columns
# (geneID, read_count, gene_length). Sample = filename. Only nonzero genes appear.
# GLUE: stack per-sample files into a gene x sample matrix, zero-filling.
import_virgo <- function(dir, pattern = "\\.out$") {
  files   <- list.files(dir, pattern = pattern, full.names = TRUE)
  samples <- sub(pattern, "", basename(files))
  per <- lapply(files, function(f)
    read.delim(f, header = FALSE, stringsAsFactors = FALSE,
               col.names = c("gene_id", "read_count", "gene_length")))

  genes  <- sort(unique(unlist(lapply(per, `[[`, "gene_id"))))
  counts <- matrix(0L, length(genes), length(samples),
                   dimnames = list(genes, samples))
  glen   <- setNames(rep(NA_integer_, length(genes)), genes)
  for (i in seq_along(per)) {
    d <- per[[i]]
    counts[d$gene_id, i] <- as.integer(d$read_count)
    glen[d$gene_id]      <- as.integer(d$gene_length)
  }
  rd <- DataFrame(gene_id = genes, length = glen[genes], row.names = genes)
  SummarizedExperiment(assays = list(counts = counts), rowData = rd)
}

# --- VALENCIA ----------------------------------------------------------------
# REAL shape: ONE wide CSV = input (sampleID, read_count, taxa...) + appended
# columns: 13 <subCST>_sim, then subCST, score, CST (Valencia.py:125-135).
# GLUE: pull the per-sample labels out of the trailing columns, keyed by sampleID.
import_valencia <- function(csv_path, sample_col = "sampleID") {
  d <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  pick <- function(nm) if (nm %in% names(d)) d[[nm]] else NA
  DataFrame(CST = pick("CST"), subCST = pick("subCST"), score = pick("score"),
            row.names = d[[sample_col]])
}
