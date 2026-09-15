# prototype/demo.R
# End-to-end on REAL-shaped mock output: import -> assemble MAE -> payoff.
suppressPackageStartupMessages(library(MultiAssayExperiment))
source("prototype/importers.R")
source("R/mock_data.R")

# Generate coordinated, REAL-shaped mock files for all three tools.
md <- tempfile("mock_demo_")
invisible(capture.output(suppressMessages(
  generate_mock_fgt_dataset(n_samples = 12, n_sequences = 1200, n_genes = 400,
                            seed = 42, create_files = TRUE, output_dir = md))))

# Make metagenomics a SUBSET (amplicon broad, metagenomics subset): keep 8 of 12.
vfiles <- list.files(file.path(md, "virgo"), pattern = "\\.out$", full.names = TRUE)
invisible(file.remove(tail(sort(vfiles), 4)))

# --- step 3: import each tool's real-shaped output -------------------------
taxa <- import_speciateit(file.path(md, "speciateit", "MC_order7_results.txt"),
                          file.path(md, "speciateit", "ASV_count_table.csv"))
genes <- import_virgo(file.path(md, "virgo"))
cst  <- import_valencia(file.path(md, "valencia", "mock_valencia_output.csv"))
cat(sprintf("imported: taxa %dx%d | genes %dx%d | CST %d samples\n\n",
            nrow(taxa), ncol(taxa), nrow(genes), ncol(genes), nrow(cst)))

# --- step 4: assemble MultiAssayExperiment ---------------------------------
mae <- MultiAssayExperiment(ExperimentList(taxa = taxa, `function` = genes),
                            colData = cst)
print(mae)

# --- step 5: payoff --------------------------------------------------------
both <- intersectColumns(mae)
cat(sprintf("\nintersectColumns(): %d samples have BOTH modalities:\n  %s\n",
            ncol(both[[1]]), paste(colnames(both)[[1]], collapse = ", ")))

cst_both <- colData(both)$CST
genes_both <- assay(both[["function"]])
mean_detected <- function(cols) if (length(cols)) mean(colSums(genes_both[, cols, drop = FALSE] > 0)) else NA
iv <- colnames(both)[[2]][grepl("^IV", cst_both)]
i  <- colnames(both)[[2]][cst_both == "I"]
cat(sprintf("\nCST-IV* (n=%d) vs CST-I (n=%d): mean VIRGO genes detected/sample = %.0f vs %.0f\n",
            length(iv), length(i), mean_detected(iv), mean_detected(i)))
