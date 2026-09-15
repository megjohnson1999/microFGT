\name{microFGT-package}
\alias{microFGT-package}
\alias{microFGT}
\docType{package}
\title{
  Comprehensive Analysis of Female Genital Tract Microbiome Data
}
\description{
  An integrated framework for analyzing microbiome data from the female
  genital tract (FGT). This package unifies specialized tools (dada2, phyloseq,
  speciateIT, VALENCIA, and VIRGO) within a cohesive framework that handles both
  amplicon and metagenomic sequencing data. Built on TreeSummarizedExperiment
  with tidyverse design principles, it provides specialized functionality for
  FGT microbiome research including community state typing and multi-omics integration.
}
\details{
  The package provides a set of functions for manipulating, analyzing, and visualizing
  female genital tract microbiome data. It is built on the TreeSummarizedExperiment 
  data structure and integrates with the wider Bioconductor ecosystem.
  
  Key functions:
  \itemize{
    \item \code{FGTExperiment()}: Create a new FGTExperiment object
    \item \code{transform_abundance()}: Transform abundance data (relative, log, etc.)
    \item \code{filter_taxa()}: Filter taxa based on prevalence and abundance
    \item \code{aggregate_taxa()}: Aggregate taxa at a specific taxonomic rank
    \item \code{plot_taxa_composition()}: Plot taxonomic composition bar charts
    \item \code{plot_alpha_diversity()}: Plot alpha diversity metrics
  }
  
  Example data utilities:
  \itemize{
    \item \code{simple_load_example()}: Load pre-built example data
    \item \code{generate_fgt_example_data()}: Generate synthetic FGT microbiome data
    \item \code{export_example_data()}: Export data to various formats
  }
}
\author{
  Scott Handley

  Maintainer: Scott Handley <scott.handley@example.com>
}
\references{
  Development based on current literature in FGT microbiome research.
}
\seealso{
  \code{\link[SummarizedExperiment]{SummarizedExperiment}},
  \code{\link[TreeSummarizedExperiment]{TreeSummarizedExperiment}}
}
\examples{
  # Load package
  library(microFGT)
  
  # Load example data
  \dontrun{
  example_data <- simple_load_example()
  
  # View basic object structure
  example_data
  
  # Transform to relative abundance
  rel_data <- transform_abundance(example_data, type = "relative")
  }
}
\keyword{ package }