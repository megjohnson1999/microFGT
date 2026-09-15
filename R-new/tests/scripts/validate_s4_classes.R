#' Validate S4 Class Registration
#'
#' Tests S4 class registration, method resolution, and object creation
#'
#' @param verbose Logical indicating whether to print detailed information
#'
#' @return Logical indicating if all checks passed
validate_s4_classes <- function(verbose = TRUE) {
  # Load required packages
  if (!requireNamespace("methods", quietly = TRUE)) {
    stop("Package 'methods' is required but not installed.")
  }
  
  # Try to load the package
  tryCatch({
    library(microFGT)
    if (verbose) cat("✓ Package loaded successfully\n")
  }, error = function(e) {
    if (verbose) cat("✗ Failed to load microFGT package:", conditionMessage(e), "\n")
    return(FALSE)
  })
  
  result <- TRUE
  
  # List of expected S4 classes - update as needed
  expected_classes <- c("FGTExperiment")
  
  # Check class registration
  if (verbose) cat("\nChecking S4 class registration:\n")
  for (class_name in expected_classes) {
    if (!methods::isClass(class_name)) {
      if (verbose) cat("✗ Class", class_name, "is not properly registered\n")
      result <- FALSE
    } else if (!methods::is(methods::getClass(class_name), "classRepresentation")) {
      if (verbose) cat("✗ Class", class_name, "is not a proper S4 class\n")
      result <- FALSE
    } else {
      if (verbose) cat("✓ Class", class_name, "is properly registered\n")
    }
  }
  
  # List of expected slots for each class
  expected_slots <- list(
    FGTExperiment = c("experimentType", "fgtMetadata")
  )
  
  # Check slot existence
  if (verbose) cat("\nChecking class slots:\n")
  for (class_name in expected_classes) {
    if (methods::isClass(class_name)) {
      slots <- methods::slotNames(methods::getClass(class_name))
      required_slots <- expected_slots[[class_name]]
      
      missing_slots <- setdiff(required_slots, slots)
      if (length(missing_slots) > 0) {
        if (verbose) cat("✗ Class", class_name, "is missing slots:", 
                         paste(missing_slots, collapse = ", "), "\n")
        result <- FALSE
      } else {
        if (verbose) cat("✓ Class", class_name, "has all required slots\n")
      }
    }
  }
  
  # List of expected methods
  expected_methods <- list(
    list(generic = "experimentType", signature = "FGTExperiment"),
    list(generic = "experimentType<-", signature = "FGTExperiment"),
    list(generic = "fgtMetadata", signature = "FGTExperiment"),
    list(generic = "fgtMetadata<-", signature = "FGTExperiment"),
    list(generic = "transformAbundance", signature = "FGTExperiment"),
    list(generic = "transformAbundance", signature = "SummarizedExperiment")
  )
  
  # Check method resolution
  if (verbose) cat("\nChecking method resolution:\n")
  for (method in expected_methods) {
    generic <- method$generic
    signature <- method$signature
    
    # Check if the generic function exists
    if (!methods::existsFunction(generic)) {
      if (verbose) cat("✗ Generic function", generic, "does not exist\n")
      result <- FALSE
      next
    }
    
    # Check if method resolves
    resolved_method <- tryCatch({
      methods::selectMethod(generic, signature, optional = TRUE)
    }, error = function(e) {
      if (verbose) cat("✗ Error selecting method", generic, "for", signature, ":", 
                      conditionMessage(e), "\n")
      result <- FALSE
      return(NULL)
    })
    
    if (is.null(resolved_method)) {
      if (verbose) cat("✗ Method", generic, "for", signature, "not found\n")
      result <- FALSE
    } else {
      if (verbose) cat("✓ Method", generic, "for", signature, "resolves correctly\n")
    }
  }
  
  # Test object creation
  if (verbose) cat("\nTesting object creation:\n")
  tryCatch({
    # Create test data
    counts <- matrix(1:12, nrow = 3, ncol = 4)
    rownames(counts) <- paste0("Gene", 1:3)
    colnames(counts) <- paste0("Sample", 1:4)
    
    # Create FGTExperiment object
    fgt <- FGTExperiment(assays = list(counts = counts))
    
    if (verbose) cat("✓ Successfully created FGTExperiment object\n")
    
    # Test if we can access slots through accessors
    exp_type <- experimentType(fgt)
    if (verbose) cat("✓ experimentType accessor returned:", exp_type, "\n")
    
    meta <- fgtMetadata(fgt)
    if (verbose) cat("✓ fgtMetadata accessor returned object of class:", class(meta)[1], "\n")
    
  }, error = function(e) {
    if (verbose) cat("✗ Error during object creation:", conditionMessage(e), "\n")
    result <- FALSE
  })
  
  # Return overall result
  if (verbose) {
    if (result) {
      cat("\n✓ All S4 class validation checks passed\n")
    } else {
      cat("\n✗ Some S4 class validation checks failed\n")
    }
  }
  
  return(result)
}

# If script is executed directly, run the validation
if (identical(environment(), globalenv())) {
  result <- validate_s4_classes(verbose = TRUE)
  
  # Set exit code based on result
  if (!result) {
    quit(status = 1)
  }
}