#' Check S4 Classes Registration and Methods
#'
#' Tests if S4 classes are properly registered, methods resolve correctly,
#' and objects can be created and manipulated.
#'
#' @param verbose Logical indicating whether to print detailed information
#'
#' @return Logical indicating if all checks passed
#' @export
check_s4_classes <- function(verbose = TRUE) {
  # Load required packages
  if (!requireNamespace("methods", quietly = TRUE)) {
    stop("Package 'methods' is required but not installed.")
  }
  
  result <- TRUE
  
  # Try to load the package
  pkg_loaded <- FALSE
  tryCatch({
    library(microFGT)
    pkg_loaded <- TRUE
    if (verbose) cat("✓ Package loaded successfully\n")
  }, error = function(e) {
    if (verbose) cat("✗ Failed to load microFGT package:", conditionMessage(e), "\n")
    result <- FALSE
  })
  
  # Only continue if package loaded
  if (!pkg_loaded) {
    return(result)
  }
  
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
    list(generic = "fgtMetadata<-", signature = "FGTExperiment")
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
  
  # Test object creation and manipulation
  if (verbose) cat("\nTesting object creation and manipulation:\n")
  tryCatch({
    # Create test data
    counts <- matrix(1:12, nrow = 3, ncol = 4)
    rownames(counts) <- paste0("Gene", 1:3)
    colnames(counts) <- paste0("Sample", 1:4)
    
    # Create FGTExperiment object
    fgt <- FGTExperiment(assays = list(counts = counts))
    
    if (verbose) cat("✓ Successfully created FGTExperiment object\n")
    
    # Test accessor
    exp_type <- experimentType(fgt)
    if (verbose) cat("✓ Default experimentType:", exp_type, "\n")
    
    # Test setter
    experimentType(fgt) <- "metagenomic"
    new_type <- experimentType(fgt)
    
    if (new_type != "metagenomic") {
      if (verbose) cat("✗ experimentType<- did not work correctly\n")
      result <- FALSE
    } else {
      if (verbose) cat("✓ experimentType<- works correctly\n")
    }
    
    # Test metadata
    fgtMetadata(fgt) <- S4Vectors::SimpleList(test = "value")
    if (!"test" %in% names(fgtMetadata(fgt))) {
      if (verbose) cat("✗ fgtMetadata<- did not work correctly\n")
      result <- FALSE
    } else {
      if (verbose) cat("✓ fgtMetadata<- works correctly\n")
    }
    
  }, error = function(e) {
    if (verbose) cat("✗ Error during object creation/manipulation:", conditionMessage(e), "\n")
    result <- FALSE
  })
  
  # Return overall result
  if (verbose) {
    if (result) {
      cat("\n✓ All S4 class checks passed\n")
    } else {
      cat("\n✗ Some S4 class checks failed\n")
    }
  }
  
  return(result)
}