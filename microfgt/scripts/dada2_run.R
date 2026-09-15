#!/usr/bin/env Rscript
# microFGT-owned glue: orchestrate DADA2 (we do NOT reimplement it).
# Denoise primer-trimmed paired FASTQs into an ASV count table + representative
# sequences, and emit a per-position quality profile so the user can set truncation
# knowledgeably. Truncation/trimming are passed in as first-class args (region-aware
# defaults are chosen in Python and handed here); every value is overridable.
#
# Usage:
#   Rscript dada2_run.R --input <trimmed_dir> --asv-table <out.csv> \
#       --asv-seqs <out.fasta> --quality-profile <out.tsv> \
#       [--trunc-len F,R] [--trim-left F,R]

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))  # qa() lives here; dada2 imports but does not attach it

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NA) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
as_pair <- function(s, default) {
  if (is.na(s)) return(default)
  as.integer(strsplit(s, ",")[[1]])
}

input_dir <- getarg("--input")
asv_table_out <- getarg("--asv-table")
asv_seqs_out <- getarg("--asv-seqs")
qprofile_out <- getarg("--quality-profile")
truncLen <- as_pair(getarg("--trunc-len"), c(0, 0))     # 0 = no truncation (default)
trimLeft <- as_pair(getarg("--trim-left"), c(0, 0))

# Accept Illumina _R1/_R2; fall back to ENA/SRA _1/_2 (marker just before the extension),
# mirroring the Python discover_pairs so both front-ends pair reads identically.
fnFs <- sort(list.files(input_dir, pattern = "_R1.*fastq", full.names = TRUE))
fnRs <- sort(list.files(input_dir, pattern = "_R2.*fastq", full.names = TRUE))
if (length(fnFs) > 0) {
  sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
} else {
  fnFs <- sort(list.files(input_dir, pattern = "_1\\.(fastq|fq)(\\.gz)?$", full.names = TRUE))
  fnRs <- sort(list.files(input_dir, pattern = "_2\\.(fastq|fq)(\\.gz)?$", full.names = TRUE))
  sample.names <- sub("_1\\.(fastq|fq)(\\.gz)?$", "", basename(fnFs))
}

# Per-position quality profile (so truncation can be set from data, not guessed).
qp <- do.call(rbind, lapply(fnFs, function(f) {
  qa <- qa(f)[["perCycle"]]$quality
  data.frame(file = basename(f), cycle = qa$Cycle, mean_quality = qa$Score, count = qa$Count)
}))
write.table(qp, qprofile_out, sep = "\t", row.names = FALSE, quote = FALSE)

filtFs <- file.path(input_dir, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(input_dir, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))
filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = truncLen, trimLeft = trimLeft,
              maxEE = c(2, 2), rm.phix = TRUE, multithread = TRUE)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
mergers <- mergePairs(dada(filtFs, err = errF, multithread = TRUE), filtFs,
                      dada(filtRs, err = errR, multithread = TRUE), filtRs)

# Key the merged results by CLEAN sample id BEFORE building the sequence table, so the table's
# row names are sample ids (e.g. 'M8092') — not filt-file basenames — and the ASV table joins
# to other assays by SAMPLE, not filename (needed for co-assayed 16S+shotgun runs). Doing it
# here (rather than re-keying the table afterwards) also fixes the single-sample case: with one
# sample, mergePairs returns a bare, name-less data.frame, which otherwise makes a row-name-less
# 1-row table and crashes the CSV write below.
if (is.data.frame(mergers)) {
  mergers <- setNames(list(mergers), sample.names)   # n = 1: wrap as a one-element named list
} else {
  names(mergers) <- sample.names                     # n > 1: re-key the list to sample ids
}
seqtab <- removeBimeraDenovo(makeSequenceTable(mergers),
                             method = "consensus", multithread = TRUE)

# Emit ASV ids (ASV1..ASVn) + a count table (samples x ASVs) + a rep-seq FASTA,
# matching the shapes speciateIT + import_speciateit expect.
asv.ids <- paste0("ASV", seq_len(ncol(seqtab)))
asv.seqs <- colnames(seqtab)
counts <- seqtab
colnames(counts) <- asv.ids
write.csv(data.frame(sampleID = rownames(counts), counts, check.names = FALSE),
          asv_table_out, row.names = FALSE)
writeLines(as.vector(rbind(paste0(">", asv.ids), asv.seqs)), asv_seqs_out)
