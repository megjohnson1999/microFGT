#' Helper functions for internal package operations
#'
#' These functions are internal utilities used by the package and not exported.
#'
#' @keywords internal

#' Check if a required package is installed
#'
#' @param package_name Name of the package to check
#' @param min_version Minimum version required (optional)
#' @param quietly Suppress messages if TRUE
#'
#' @return Logical indicating if package is available
#' @keywords internal
check_package <- function(package_name, min_version = NULL, quietly = FALSE) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (!quietly) {
      message("Package '", package_name, "' is required but not installed.")
    }
    return(FALSE)
  }
  
  if (!is.null(min_version)) {
    pkg_version <- utils::packageVersion(package_name)
    if (pkg_version < min_version) {
      if (!quietly) {
        message("Package '", package_name, "' version ", pkg_version, 
                " is installed, but version >= ", min_version, " is required.")
      }
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Validate input parameters
#'
#' @param param Value to validate
#' @param valid_values Vector of valid values or function that returns TRUE/FALSE
#' @param param_name Name of parameter for error messages
#' @param allow_null Whether NULL is a valid value
#'
#' @return The validated parameter value
#' @keywords internal
validate_param <- function(param, valid_values, param_name, allow_null = FALSE) {
  if (is.null(param)) {
    if (allow_null) {
      return(NULL)
    } else {
      stop(param_name, " cannot be NULL", call. = FALSE)
    }
  }
  
  if (is.function(valid_values)) {
    is_valid <- valid_values(param)
    if (!is_valid) {
      stop(param_name, " has invalid value: ", param, call. = FALSE)
    }
  } else if (is.character(valid_values)) {
    if (!param %in% valid_values) {
      stop(param_name, " must be one of: ", paste(valid_values, collapse = ", "), 
           ", but got: ", param, call. = FALSE)
    }
  }
  
  return(param)
}

#' Format error messages consistently
#'
#' @param msg Error message text
#' @param context Context in which the error occurred
#'
#' @return Formatted error message
#' @keywords internal
format_error <- function(msg, context = NULL) {
  if (!is.null(context)) {
    return(paste0("[", context, "] ", msg))
  }
  return(msg)
}