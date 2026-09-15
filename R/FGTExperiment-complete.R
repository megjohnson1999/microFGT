#' FGTExperiment Class - Compositional Implementation
#'
#' The FGTExperiment class provides a container for female genital tract microbiome data.
#' This implementation uses a compositional design where the class contains a 
#' TreeSummarizedExperiment rather than extending it.
#'
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment assay assays assayNames assayNames<- rowData colData
#' @importFrom methods is new setClass setMethod setGeneric setValidity show
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @import methods
#' 
#' @name FGTExperiment-class
#' @rdname FGTExperiment-class
#' @exportClass FGTExperiment
#' @export FGTExperiment
setClass("FGTExperiment",
         slots = list(
           experimentData = "TreeSummarizedExperiment",
           experimentType = "character",
           fgtMetadata = "SimpleList"
         ),
         prototype = list(
           experimentType = "amplicon",
           fgtMetadata = S4Vectors::SimpleList()
         ))

#' Validate FGTExperiment objects
#'
#' @param object FGTExperiment object to validate
#' @return TRUE if valid, otherwise error messages
setValidity("FGTExperiment", function(object) {
  msg <- NULL
  
  # Validate experimentType
  valid_types <- c("amplicon", "metagenomic", "integrated")
  if (length(object@experimentType) != 1 || 
      !object@experimentType %in% valid_types) {
    msg <- c(msg, paste("experimentType must be one of:",
                        paste(valid_types, collapse = ", ")))
  }
  
  # Validate fgtMetadata
  if (!is(object@fgtMetadata, "SimpleList")) {
    msg <- c(msg, "fgtMetadata must be a SimpleList object")
  }
  
  # Return TRUE if no errors, otherwise return error messages
  if (length(msg) == 0) TRUE else msg
})

#' Get the experiment type
#'
#' @param object An FGTExperiment object
#'
#' @return Character string indicating the experiment type
#' @export
setGeneric("experimentType", function(object) standardGeneric("experimentType"))

#' @rdname experimentType
#' @export
setMethod("experimentType", "FGTExperiment", function(object) {
  object@experimentType
})

#' Set the experiment type
#'
#' @param object An FGTExperiment object
#' @param value Character string with the new experiment type
#'
#' @return Updated FGTExperiment object
#' @export
setGeneric("experimentType<-", function(object, value) standardGeneric("experimentType<-"))

#' @rdname experimentType-set
#' @export
setMethod("experimentType<-", "FGTExperiment", function(object, value) {
  valid_types <- c("amplicon", "metagenomic", "integrated")
  if (!value %in% valid_types) {
    stop("experimentType must be one of: ", paste(valid_types, collapse = ", "))
  }
  object@experimentType <- value
  validObject(object)
  return(object)
})

#' Get the FGT metadata
#'
#' @param object An FGTExperiment object
#'
#' @return SimpleList with FGT metadata
#' @export
setGeneric("fgtMetadata", function(object) standardGeneric("fgtMetadata"))

#' @rdname fgtMetadata
#' @export
setMethod("fgtMetadata", "FGTExperiment", function(object) {
  object@fgtMetadata
})

#' Set the FGT metadata
#'
#' @param object An FGTExperiment object
#' @param value SimpleList with the new metadata
#'
#' @return Updated FGTExperiment object
#' @export
setGeneric("fgtMetadata<-", function(object, value) standardGeneric("fgtMetadata<-"))

#' @rdname fgtMetadata-set
#' @export
setMethod("fgtMetadata<-", "FGTExperiment", function(object, value) {
  if (!is(value, "SimpleList")) {
    value <- S4Vectors::SimpleList(value)
  }
  object@fgtMetadata <- value
  validObject(object)
  return(object)
})

#' Transform abundance data
#'
#' @param object An FGTExperiment object
#' @param type Type of transformation ("relative", "log", "clr", "presence")
#' @param assay_name Name of the assay to transform
#' @param pseudocount Pseudocount to add for log and clr transformations
#'
#' @return Transformed FGTExperiment object
#' @export
setGeneric("transformAbundance", function(object, type = "relative", assay_name = "counts", pseudocount = 1) 
  standardGeneric("transformAbundance"))

#' @rdname transformAbundance
#' @export
setMethod("transformAbundance", "FGTExperiment", function(object, type = "relative", assay_name = "counts", pseudocount = 1) {
  # Validate inputs
  valid_types <- c("relative", "log", "clr", "presence")
  if (!type %in% valid_types) {
    stop("type must be one of: ", paste(valid_types, collapse = ", "))
  }
  
  if (!assay_name %in% assayNames(object@experimentData)) {
    stop("assay_name '", assay_name, "' not found in assays")
  }
  
  # Get the count matrix
  counts <- assay(object@experimentData, assay_name)
  
  # Apply transformation
  if (type == "relative") {
    # Relative abundance (proportions)
    result <- t(t(counts) / colSums(counts))
    new_name <- paste0("relative_", assay_name)
  } else if (type == "log") {
    # Log transformation
    result <- log(counts + pseudocount)
    new_name <- paste0("log_", assay_name)
  } else if (type == "clr") {
    # Centered log-ratio transformation
    counts_adj <- counts + pseudocount
    # Calculate geometric mean for each sample
    geo_means <- apply(counts_adj, 2, function(x) exp(mean(log(x))))
    # Apply CLR transformation
    result <- log(t(t(counts_adj) / geo_means))
    new_name <- paste0("clr_", assay_name)
  } else if (type == "presence") {
    # Presence/absence (binary)
    result <- counts > 0
    mode(result) <- "numeric"  # Convert logical to numeric
    new_name <- paste0("presence_", assay_name)
  }
  
  # Add the transformed matrix as a new assay
  assays <- assays(object@experimentData)
  assays[[new_name]] <- result
  
  # Update the TreeSummarizedExperiment
  expData <- object@experimentData
  assays(expData) <- assays
  
  # Update the FGTExperiment
  object@experimentData <- expData
  
  return(object)
})

# Core accessor methods that proxy to the experimentData TreeSummarizedExperiment

#' @export
setMethod("dim", "FGTExperiment", function(x) {
  dim(x@experimentData)
})

#' @export
setMethod("dimnames", "FGTExperiment", function(x) {
  dimnames(x@experimentData)
})

#' @export
setMethod("rownames", "FGTExperiment", function(x) {
  rownames(x@experimentData)
})

#' @export
setMethod("colnames", "FGTExperiment", function(x) {
  colnames(x@experimentData)
})

#' @export
setMethod("assay", "FGTExperiment", function(x, i, ...) {
  assay(x@experimentData, i, ...)
})

#' @export
setMethod("assays", "FGTExperiment", function(x, ...) {
  assays(x@experimentData, ...)
})

#' @export
setMethod("assayNames", "FGTExperiment", function(x, ...) {
  assayNames(x@experimentData, ...)
})

#' @export
setMethod("rowData", "FGTExperiment", function(x, ...) {
  rowData(x@experimentData, ...)
})

#' @export
setMethod("colData", "FGTExperiment", function(x, ...) {
  colData(x@experimentData, ...)
})

#' @export
#' @export
setMethod("rowTree", "FGTExperiment", function(x) {
  TreeSummarizedExperiment::rowTree(x@experimentData)
})
#' @export
setMethod("show", "FGTExperiment", function(object) {
  cat("FGTExperiment object with", nrow(object), "features and", ncol(object), "samples
")
  cat("experimentType:", object@experimentType, "
")
  
  fgt_meta <- object@fgtMetadata
  if (length(fgt_meta) > 0) {
    cat("fgtMetadata: SimpleList with", length(fgt_meta), "elements
")
    meta_names <- names(fgt_meta)
    for (i in seq_along(fgt_meta)) {
      meta_item <- fgt_meta[[i]]
      cat("  $", meta_names[i], ": ", class(meta_item)[1],
          if (length(meta_item) == 1) {
            paste0(" [", as.character(meta_item), "]")
          } else {
            paste0(" [length: ", length(meta_item), "]")
          },
          "
", sep="")
    }
  } else {
    cat("fgtMetadata: Empty SimpleList
")
  }
  
  # Show assay information
  assay_names <- assayNames(object)
  if (length(assay_names) > 0) {
    cat("assays(", length(assay_names), "): ", paste(assay_names, collapse=", "), "
", sep="")
  }
  
  # Show dimensions
  cat("rownames(", nrow(object), "): ", paste(head(rownames(object), 3), collapse=", "), 
      if (nrow(object) > 3) "..." else "", "
", sep="")
  cat("colnames(", ncol(object), "): ", paste(head(colnames(object), 3), collapse=", "), 
      if (ncol(object) > 3) "..." else "", "
", sep="")
})

#' Create a new FGTExperiment object
#'
#' @param assays List of matrices or similar objects containing count data,
#'        or a SummarizedExperiment or TreeSummarizedExperiment object
#' @param rowData DataFrame of feature metadata
#' @param colData DataFrame of sample metadata
#' @param rowTree Phylogenetic tree of type phylo
#' @param experimentType Type of experiment ("amplicon", "metagenomic", "integrated")
#' @param fgtMetadata Additional metadata as a SimpleList
#' @param ... Additional arguments passed to TreeSummarizedExperiment
#'
#' @return A new FGTExperiment object
#' @export
FGTExperiment <- function(assays = list(), rowData = NULL, colData = NULL, rowTree = NULL,
                         experimentType = "amplicon", fgtMetadata = S4Vectors::SimpleList(), ...) {
  
  # Create the TreeSummarizedExperiment
  if (methods::is(assays, "TreeSummarizedExperiment")) {
    tse <- assays
  } else if (methods::is(assays, "SummarizedExperiment")) {
    tse <- tryCatch({
      TreeSummarizedExperiment::TreeSummarizedExperiment(
        SummarizedExperiment::assays(assays),
        rowData = SummarizedExperiment::rowData(assays),
        colData = SummarizedExperiment::colData(assays),
        rowTree = rowTree
      )
    }, error = function(e) {
      warning("Error converting SummarizedExperiment to TreeSummarizedExperiment: ", conditionMessage(e))
      # Create an empty TSE with the same dimensions
      tse_empty <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(counts = matrix(0, nrow=nrow(assays), ncol=ncol(assays))),
        rowData = SummarizedExperiment::rowData(assays),
        colData = SummarizedExperiment::colData(assays)
      )
      # Copy the assays over
      SummarizedExperiment::assays(tse_empty) <- SummarizedExperiment::assays(assays)
      tse_empty
    })
  } else {
    # Input validation for creating from scratch
    if (is.list(assays) && length(assays) > 0) {
      # Verify the first assay 
      first_assay <- assays[[1]]
      if (!is.matrix(first_assay) && !methods::is(first_assay, "Matrix")) {
        stop("The first element of 'assays' must be a matrix or Matrix object")
      }
      
      # Ensure row and column names are present
      if (is.null(rownames(first_assay))) {
        rownames(first_assay) <- paste0("Feature", seq_len(nrow(first_assay)))
        assays[[1]] <- first_assay
      }
      
      if (is.null(colnames(first_assay))) {
        colnames(first_assay) <- paste0("Sample", seq_len(ncol(first_assay)))
        assays[[1]] <- first_assay
      }
      
      # Prepare rowData if needed
      if (is.null(rowData)) {
        rowData <- S4Vectors::DataFrame(row.names = rownames(first_assay))
      } else if (!is.data.frame(rowData) && !methods::is(rowData, "DataFrame")) {
        rowData <- S4Vectors::DataFrame(rowData, row.names = rownames(first_assay))
      }
      
      # Prepare colData if needed
      if (is.null(colData)) {
        colData <- S4Vectors::DataFrame(row.names = colnames(first_assay))
      } else if (!is.data.frame(colData) && !methods::is(colData, "DataFrame")) {
        colData <- S4Vectors::DataFrame(colData, row.names = colnames(first_assay))
      }
      
      # Create TreeSummarizedExperiment
      tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = assays,
        rowData = rowData,
        colData = colData,
        rowTree = rowTree,
        ...
      )
    } else {
      stop("'assays' must be a non-empty list, SummarizedExperiment, or TreeSummarizedExperiment")
    }
  }
  
  # Convert fgtMetadata to SimpleList if necessary
  if (!methods::is(fgtMetadata, "SimpleList")) {
    fgtMetadata <- S4Vectors::SimpleList(fgtMetadata)
  }
  
  # Validate experimentType
  validTypes <- c("amplicon", "metagenomic", "integrated")
  if (!experimentType %in% validTypes) {
    stop("experimentType must be one of: ", paste(validTypes, collapse=", "))
  }
  
  # Create FGTExperiment object using new()
  fgtObj <- methods::new("FGTExperiment",
                       experimentData = tse,
                       experimentType = experimentType,
                       fgtMetadata = fgtMetadata)
  
  # Validate and return the object
  methods::validObject(fgtObj)
  return(fgtObj)
}

#' Migrate legacy FGTExperiment objects to new structure
#'
#' This function converts FGTExperiment objects created with versions prior to 
#' the compositional redesign (where FGTExperiment extends TreeSummarizedExperiment)
#' to the new compositional structure (where FGTExperiment contains a TreeSummarizedExperiment).
#'
#' @param object Legacy FGTExperiment object or a compatible object
#'
#' @return A new FGTExperiment object with the updated structure
#' @export
migrateToNewFGTExperiment <- function(object) {
  if (is(object, "FGTExperiment")) {
    # Check if this is already a new-style FGTExperiment
    if (methods::existsMethod("experimentData", "FGTExperiment") && 
        methods::existsSlot(object, "experimentData")) {
      # This appears to be a new-style object
      message("Object already has the new FGTExperiment structure.")
      return(object)
    }
    
    # This is an old-style object where FGTExperiment extends TreeSummarizedExperiment
    # Extract the key parts from the old object
    exp_type <- if (methods::existsSlot(object, "experimentType")) {
      object@experimentType 
    } else {
      "amplicon"
    }
    
    fgt_meta <- if (methods::existsSlot(object, "fgtMetadata")) {
      object@fgtMetadata
    } else {
      S4Vectors::SimpleList()
    }
    
    # Create a new object by treating the old object as a TreeSummarizedExperiment
    new_obj <- FGTExperiment(
      assays = object, # The object itself is passed as the TSE
      experimentType = exp_type,
      fgtMetadata = fgt_meta
    )
    
    return(new_obj)
  } else if (is(object, "TreeSummarizedExperiment") || is(object, "SummarizedExperiment")) {
    # For TreeSummarizedExperiment or SummarizedExperiment, convert directly
    return(FGTExperiment(object))
  } else {
    stop("Cannot convert object of class '", class(object), "' to FGTExperiment")
  }
}
