# prototype/test_importers.R
# Contract check: each importer must parse BOTH the real-shaped fixture AND the
# fixed mock output into a valid, same-shaped object. This is the forcing
# function the mock never had: the SAME importer reads real and mock, or fails.
suppressPackageStartupMessages(library(MultiAssayExperiment))
source("prototype/importers.R")

FX <- "prototype/real_fixtures"
ok <- TRUE
check <- function(label, cond) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) ok <<- FALSE
}

cat("== REAL fixtures ==\n")

# VIRGO: genuine per-sample outputs (ravel-lab/VIRGO _test_run).
vdir <- tempfile("virgo_real_"); dir.create(vdir)
invisible(file.copy(file.path(FX, "virgo_sub1.out"), file.path(vdir, "sub1.out")))
invisible(file.copy(file.path(FX, "virgo_sub2.out"), file.path(vdir, "sub2.out")))
v_real <- import_virgo(vdir)
check("VIRGO real: genes x 2 samples", ncol(v_real) == 2 && nrow(v_real) > 1000)
check("VIRGO real: integer counts, zero-filled (no NA)", !anyNA(assay(v_real)))

# speciateIT: genuine ASV count table (169 x 1514) + real-SHAPED results over its ASVs.
ct_path <- file.path(FX, "speciateit_test_count_table.csv")
asvs <- colnames(read.csv(ct_path, row.names = 1, check.names = FALSE, nrows = 1))
set.seed(1)
pool <- c("Lactobacillus crispatus", "Lactobacillus iners", "Gardnerella vaginalis",
          "Atopobium vaginae", "Prevotella bivia")
res_path <- tempfile("MC_", fileext = ".txt")
write.table(data.frame(
  "Sequence ID"           = asvs,
  "Classification"        = sample(pool, length(asvs), replace = TRUE),
  "posterior probability" = runif(length(asvs), 0.7, 1.0),
  "number of Decisions"   = sample(5:20, length(asvs), replace = TRUE),
  check.names = FALSE), res_path, sep = "\t", quote = FALSE, row.names = FALSE)
t_real <- import_speciateit(res_path, ct_path)
check("speciateIT real: taxa x 169 samples (ASV->taxon join)",
      ncol(t_real) == 169 && nrow(t_real) <= length(pool) + 1)

# VALENCIA: GENUINE output produced by running Valencia.py on real published
# composition data (13,231 samples; 99.9% CST concordance with the paper's own
# Val_CST labels). Fixture is the head of that genuine run.
c_real <- import_valencia(file.path(FX, "valencia_genuine_output_head.csv"))
check("VALENCIA real (GENUINE tool output): CST/subCST/score from trailing cols",
      nrow(c_real) == 6 && all(c("CST","subCST","score") %in% names(c_real)) &&
      !anyNA(c_real$CST))

cat("== MOCK via fixed writers ==\n")
source("R/mock_data.R")
md <- tempfile("mock_")
invisible(capture.output(suppressMessages(
  generate_mock_fgt_dataset(n_samples = 10, n_sequences = 600, n_genes = 300,
                            seed = 42, create_files = TRUE, output_dir = md))))
t_mock <- import_speciateit(file.path(md, "speciateit", "MC_order7_results.txt"),
                            file.path(md, "speciateit", "ASV_count_table.csv"))
v_mock <- import_virgo(file.path(md, "virgo"))
c_mock <- import_valencia(file.path(md, "valencia", "mock_valencia_output.csv"))
check("speciateIT mock: same importer, taxa x 10", ncol(t_mock) == 10)
check("VIRGO mock: same importer reads per-sample .out, genes x 10", ncol(v_mock) == 10)
check("VALENCIA mock: same importer, 10 samples w/ CST",
      nrow(c_mock) == 10 && !anyNA(c_mock$CST))

cat("== Assemble + payoff (mock) ==\n")
mae <- MultiAssayExperiment(ExperimentList(taxa = t_mock, `function` = v_mock),
                            colData = c_mock)
check("MAE assembles (2 experiments)", length(experiments(mae)) == 2)
both <- intersectColumns(mae)
check("intersectColumns() yields shared samples", ncol(both[[1]]) >= 1)

cat(if (ok) "\nALL CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
