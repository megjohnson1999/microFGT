#!/usr/bin/env Rscript
# microFGT-owned glue: flatten a phyloseq .rds so the Python side can build an
# ASV-grain composition AnnData. We do NOT reimplement phyloseq -- we read the
# .rds and write its four slots as flat CSVs (the same orchestrate-an-R-tool
# stance as dada2_run.R).
#
# Usage:
#   Rscript phyloseq_export.R --rds <ps.rds> --outdir <dir>
#
# Writes into <dir> (all keyed so the Python side can align them):
#   counts.csv       rows = sample ids, cols = ASV ids (ASV1..ASVn)  -> counts
#   taxa_names.csv   rows = ASV ids, col 'sequence'                   -> the taxa_names
#                    (sequences ARE the taxa names; there is no refseq slot)
#   tax_table.csv    rows = ASV ids, cols = taxonomic ranks           -> tax_table
#   sample_data.csv  rows = sample ids, cols = all sample variables   -> sample_data

suppressPackageStartupMessages(library(phyloseq))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NA) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

rds <- getarg("--rds")
outdir <- getarg("--outdir")
if (is.na(rds) || is.na(outdir)) stop("both --rds and --outdir are required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ps <- readRDS(rds)
if (!methods::is(ps, "phyloseq")) {
  stop("Object in ", rds, " is not a phyloseq object (got ", class(ps)[1], ").")
}

# Sequences ARE the taxa names (no refseq slot). Assign stable ASV ids in taxa
# order so the huge sequences never have to be CSV column headers.
seqs <- taxa_names(ps)
asv_ids <- paste0("ASV", seq_along(seqs))

# counts: orient samples x ASVs regardless of how the otu_table is stored.
m <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) m <- t(m)          # -> samples x taxa
colnames(m) <- asv_ids
write.csv(as.data.frame(m, check.names = FALSE),
          file.path(outdir, "counts.csv"))

# taxa_names: ASV id -> sequence
taxa_df <- data.frame(sequence = seqs, row.names = asv_ids)
write.csv(taxa_df, file.path(outdir, "taxa_names.csv"))

# tax_table: ASV id -> ranks (Domain..Genus_Species). May be absent on a bare object.
tt_slot <- access(ps, "tax_table", errorIfNULL = FALSE)
if (!is.null(tt_slot)) {
  tt <- as.data.frame(as(tt_slot, "matrix"), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(tt) <- asv_ids
} else {
  tt <- data.frame(row.names = asv_ids)
}
write.csv(tt, file.path(outdir, "tax_table.csv"))

# sample_data: sample id -> all sample variables (kept whole, incl. CST/subCST/score).
sd_slot <- access(ps, "sam_data", errorIfNULL = FALSE)
if (!is.null(sd_slot)) {
  sd <- as(sd_slot, "data.frame")
} else {
  sd <- data.frame(row.names = sample_names(ps))
}
write.csv(sd, file.path(outdir, "sample_data.csv"))
