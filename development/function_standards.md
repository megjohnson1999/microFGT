# microFGT Function Design Standards

This document outlines the standards for function design in the microFGT package, including parameter naming, error handling, return value patterns, and input validation.

## Parameter Naming Conventions

### General Rules

- Use `snake_case` for all parameter names
- Order parameters from most to least important
- Required parameters first, optional parameters with defaults later
- Use consistent parameter names across similar functions
- Use descriptive names that indicate the parameter's purpose
- Avoid abbreviations unless they are standard in the field

### Standard Parameter Names

| Parameter         | Type             | Description                                    |
|-------------------|------------------|------------------------------------------------|
| `fgt_exp`         | FGTExperiment    | FGTExperiment object to process                |
| `assay_name`      | character        | Name of assay to use (e.g., "counts", "tpm")   |
| `taxa_level`      | character        | Taxonomic level (e.g., "Phylum", "Genus")      |
| `min_prevalence`  | numeric          | Minimum prevalence threshold                   |
| `min_abundance`   | numeric          | Minimum abundance threshold                    |
| `method`          | character        | Method to use for calculation                  |
| `n_cores`         | integer          | Number of cores for parallel processing        |
| `verbose`         | logical          | Whether to print messages during execution     |
| `plot`            | logical          | Whether to generate a plot                     |
| `threshold`       | numeric          | Threshold value for filtering or detection     |

### Examples

#### Good Parameter Naming

```r
# Good example - consistent naming style and order
filter_taxa(fgt_exp, min_prevalence = 0.1, min_abundance = 5, assay_name = "counts")

# Good example - clarity in parameter purpose
transform_abundance(fgt_exp, method = "clr", pseudocount = 1, assay_name = "counts", new_assay_name = "clr_counts")
```

#### Bad Parameter Naming

```r
# Bad example - inconsistent naming style and unclear names
filter_taxa(x, prev = 0.1, abund = 5, assay = "counts")

# Bad example - unclear parameter purpose
transform_abundance(object, type = "clr", pc = 1, input = "counts", output = "clr_counts")
```

## Error Handling Standards

### Error Message Format

- Use the `format_error()` utility function for consistent error formatting
- Include context (function name) in all error messages
- Be specific about what went wrong
- Suggest how to fix the problem when possible

### Validation Approach

- Use `validate_param()` for parameter validation
- Check required parameters are not NULL
- Validate parameter values against expected types and ranges
- Provide specific error messages for each validation failure

### Error Handling Example

```r
filter_taxa <- function(fgt_exp, min_prevalence = 0.1, min_abundance = 0, assay_name = "counts") {
  # Validate FGTExperiment object
  if (!is_fgt_experiment(fgt_exp)) {
    stop(format_error("Expected an FGTExperiment object", "filter_taxa"))
  }
  
  # Validate numeric parameters
  if (!is.numeric(min_prevalence) || min_prevalence < 0 || min_prevalence > 1) {
    stop(format_error("min_prevalence must be a number between 0 and 1", "filter_taxa"))
  }
  
  if (!is.numeric(min_abundance) || min_abundance < 0) {
    stop(format_error("min_abundance must be a non-negative number", "filter_taxa"))
  }
  
  # Validate assay_name exists
  if (!assay_name %in% assayNames(fgt_exp)) {
    stop(format_error(
      sprintf("assay_name '%s' not found. Available assays: %s", 
              assay_name, paste(assayNames(fgt_exp), collapse = ", ")),
      "filter_taxa"))
  }
  
  # Function implementation...
}
```

## Return Value Standards

### General Rules

- Functions that modify FGTExperiment objects should return FGTExperiment objects
- Functions that extract data should return consistent types
- Document the return type and structure clearly
- Add metadata about operations performed to the returned object when possible
- Be consistent with return values for similar functions

### Return Value Patterns

| Function Type     | Return Type       | Metadata                                     |
|-------------------|-------------------|----------------------------------------------|
| Filters           | FGTExperiment     | Filter criteria in object metadata           |
| Transformations   | FGTExperiment     | New assay(s) with transformation info        |
| Extractors        | data.frame/matrix | Row and column names preserved               |
| Calculations      | list/vector/matrix| Named elements for results                   |
| Visualizations    | ggplot            | Return plot object, not just display         |

### Return Value Example

```r
transform_abundance <- function(fgt_exp, method = "clr", pseudocount = 1, 
                               assay_name = "counts", new_assay_name = NULL) {
  # Validation code...
  
  # Perform transformation
  counts <- assay(fgt_exp, assay_name)
  
  # Apply transformation logic...
  
  # Set default new_assay_name if not provided
  if (is.null(new_assay_name)) {
    new_assay_name <- paste0(method, "_", assay_name)
  }
  
  # Add transformed data as new assay
  result <- fgt_exp
  assay(result, new_assay_name) <- transformed_counts
  
  # Track transformation in metadata
  transformation_info <- list(
    original_assay = assay_name,
    method = method,
    pseudocount = pseudocount,
    date = Sys.Date()
  )
  
  # Store transformation info in fgtMetadata
  transformations <- fgtMetadata(result)$transformations
  if (is.null(transformations)) {
    transformations <- list()
  }
  transformations[[new_assay_name]] <- transformation_info
  fgtMetadata(result)$transformations <- transformations
  
  return(result)
}
```

## Input Validation Standards

### What to Validate

- Class/type of each parameter (e.g., is FGTExperiment, is numeric)
- Ranges for numeric parameters
- Valid options for character parameters with limited choices
- Existence of required data (e.g., assay names)
- Compatibility between parameters

### Validation Tools

- `validate_param()` - Utility function for standard validation
- `is_fgt_experiment()` - Check for FGTExperiment objects
- `check_package()` - Verify required packages are available
- Standard R checks: `is.numeric()`, `is.character()`, etc.

### Validation Example

```r
calculate_diversity <- function(fgt_exp, method = "shannon", assay_name = "counts", 
                              taxa_level = NULL, pseudocount = 0) {
  # Validate FGTExperiment
  if (!is_fgt_experiment(fgt_exp)) {
    stop(format_error("Expected an FGTExperiment object", "calculate_diversity"))
  }
  
  # Validate method
  method <- validate_param(method, 
                           c("shannon", "simpson", "invsimpson", "richness"), 
                           "method")
  
  # Validate assay_name
  if (!assay_name %in% assayNames(fgt_exp)) {
    stop(format_error(
      sprintf("assay_name '%s' not found. Available assays: %s", 
              assay_name, paste(assayNames(fgt_exp), collapse = ", ")),
      "calculate_diversity"))
  }
  
  # Check if required package is available for calculation
  if (!check_package("vegan", quietly = TRUE)) {
    stop(format_error(
      "Package 'vegan' is required for diversity calculations but is not installed",
      "calculate_diversity"))
  }
  
  # Function implementation...
}
```

## Function Documentation Standards

All functions should include complete roxygen2 documentation with:

1. Title and description
2. Parameter descriptions
3. Return value specification
4. Examples
5. References (if applicable)
6. Edge cases and limitations
7. See also section for related functions

### Documentation Example

```r
#' Calculate alpha diversity measures
#'
#' Calculates various alpha diversity indices from abundance data in an FGTExperiment object.
#'
#' @param fgt_exp An FGTExperiment object
#' @param method Diversity index to calculate, one of: "shannon", "simpson", "invsimpson", "richness"
#' @param assay_name Name of the assay to use for calculations
#' @param taxa_level Optional taxonomic level to aggregate to before calculation
#' @param pseudocount Small value to add to counts to avoid zero issues (default: 0)
#'
#' @return A data.frame with samples as rows and the calculated diversity measure as column
#'
#' @details
#' This function calculates alpha diversity using the vegan package. The "shannon" method 
#' calculates Shannon's diversity index, "simpson" calculates Simpson's index, "invsimpson"
#' calculates the inverse Simpson index, and "richness" calculates species richness.
#'
#' @examples
#' # Load example data
#' data(example_fgt)
#'
#' # Calculate Shannon diversity
#' shannon_div <- calculate_diversity(example_fgt, method = "shannon")
#'
#' # Calculate richness at Genus level
#' genus_richness <- calculate_diversity(example_fgt, 
#'                                     method = "richness",
#'                                     taxa_level = "Genus")
#'
#' @seealso \code{\link{calculate_beta_diversity}} for between-sample diversity measures
#'
#' @export
calculate_diversity <- function(fgt_exp, method = "shannon", assay_name = "counts", 
                              taxa_level = NULL, pseudocount = 0) {
  # Function implementation...
}
```

## Implementation Approach

To standardize existing functions:

1. Identify groups of related functions
2. Create templates for each function type
3. Update functions to follow standards
4. Add comprehensive testing for edge cases
5. Update documentation to match implementation

New functions should be designed following these standards from the start.