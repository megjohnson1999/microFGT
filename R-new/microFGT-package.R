#' microFGT: Comprehensive Analysis of Female Genital Tract Microbiome Data
#'
#' @description
#' An integrated framework for analyzing microbiome data from the female
#' genital tract (FGT). This package unifies specialized tools within a cohesive
#' framework that handles both amplicon and metagenomic sequencing data.
#' Built on TreeSummarizedExperiment with efficient data structures,
#' it provides specialized functionality for FGT microbiome research.
#'
#' @details
#' The package provides:
#' \itemize{
#'   \item \code{FGTExperiment}: A specialized S4 class for FGT microbiome data
#'   \item Data transformation functions
#'   \item Import/export utilities
#'   \item Helper functions for common operations
#' }
#'
#' @docType package
#' @name microFGT-package
#' @aliases microFGT
#' @importFrom SummarizedExperiment SummarizedExperiment assay assays assayNames
#' @importFrom S4Vectors SimpleList DataFrame
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @author Scott Handley
NULL