#' microFGT: Integrated Analysis of Female Genital Tract Microbiome Data
#'
#' @description
#' An integrated framework for analyzing microbiome data from the female
#' genital tract (FGT). This package provides a streamlined data structure for FGT
#' microbiome analysis with specialized support for integrating outputs from speciateIT,
#' VALENCIA, and VIRGO. Focused on simplicity and usability, it offers a clean API
#' for transformation, analysis, and visualization of microbiome data.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{FGTMicrobiome}}: Create a new FGTMicrobiome object
#'   \item \code{\link{counts}}: Get or set count data
#'   \item \code{\link{taxonomy}}: Get or set taxonomic classifications
#'   \item \code{\link{metadata}}: Get or set sample metadata
#'   \item \code{\link{transform_abundance}}: Apply transformations to count data
#'   \item \code{\link{classify_with_speciateit}}: Classify with speciateIT
#'   \item \code{\link{generate_mock_fgt_dataset}}: Generate mock data for testing
#' }
#'
#' @section Mock Data Generators:
#' The package includes functions to generate realistic test data:
#' \itemize{
#'   \item \code{\link{generate_mock_speciateit}}: Generate mock SpeciateIT classifications
#'   \item \code{\link{generate_mock_virgo}}: Generate mock VIRGO gene abundances
#'   \item \code{\link{generate_mock_valencia}}: Generate mock VALENCIA CST data
#'   \item \code{\link{generate_mock_fgt_dataset}}: Generate coordinated data for all tools
#' }
#'
#' @name microFGT-package
#' @aliases microFGT
#' @useDynLib microFGT, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods is
#' @importFrom stats setNames dist median p.adjust prcomp wilcox.test
#' @importFrom utils head read.table
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @importFrom dplyr filter group_by mutate select summarize
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal theme element_text
NULL