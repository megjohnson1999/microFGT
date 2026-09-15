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
#'
#' @examples
#' \dontrun{
#' # For a legacy FGTExperiment object 'old_fgt':
#' new_fgt <- migrateToNewFGTExperiment(old_fgt)
#' }
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
    # and adding the extracted values
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