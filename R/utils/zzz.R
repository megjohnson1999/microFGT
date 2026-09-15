#' @import methods
#' @importFrom SummarizedExperiment SummarizedExperiment assay assays assayNames
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
NULL

.onLoad <- function(libname, pkgname) {
  # Package initialization code
}

.onAttach <- function(libname, pkgname) {
  msg <- sprintf(
    "microFGT %s: Integrated Framework for Female Genital Tract Microbiome Analysis",
    utils::packageVersion("microFGT")
  )
  packageStartupMessage(msg)
}