#!/bin/bash
# Pre-commit hook for microFGT package

# Print colored output
green() { echo -e "\033[0;32m$1\033[0m"; }
yellow() { echo -e "\033[0;33m$1\033[0m"; }
red() { echo -e "\033[0;31m$1\033[0m"; }
blue() { echo -e "\033[0;34m$1\033[0m"; }

# Capture the start time
start_time=$(date +%s)

# Show header
blue "========================================================"
blue "              microFGT Package Pre-commit Hook           "
blue "========================================================"
yellow "Running pre-commit checks - this may take a minute..."

# Initialize error flag
errors_found=0

# Function to run a check and report status
run_check() {
  local name="$1"
  local command="$2"
  
  yellow "Running check: $name"
  if eval "$command"; then
    green "✓ $name check passed"
    return 0
  else
    red "✗ $name check failed"
    errors_found=1
    return 1
  fi
}

# Step 1: Run roxygen2 to update documentation
run_check "Roxygen documentation update" "Rscript -e 'roxygen2::roxygenize()'" || {
  yellow "Tip: Make sure roxygen2 is installed and all documentation is properly formatted"
}

# Step 2: Validate NAMESPACE
run_check "NAMESPACE validation" "Rscript tests/scripts/validate_namespace.R" || {
  yellow "Tip: Check that all exports are properly declared in roxygen comments"
}

# Step 3: Validate S4 classes
run_check "S4 class validation" "Rscript tests/scripts/validate_s4_classes.R" || {
  yellow "Tip: Make sure all S4 classes are properly defined and methods are registered"
}

# Step 4: Run unit tests
run_check "Run tests" "Rscript -e 'testthat::test_dir(\"tests/testthat\")'" || {
  yellow "Tip: Fix failing tests before committing"
}

# Step 5: Check package style (optional)
if command -v Rscript >/dev/null 2>&1 && Rscript -e "find.package('styler', quiet = TRUE)" >/dev/null 2>&1; then
  run_check "Code style check" "Rscript -e 'styler::style_pkg()'" || {
    yellow "Tip: Consider fixing style issues"
  }
fi

# Calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))

# Report result
blue "========================================================"
if [ $errors_found -eq 0 ]; then
  green "✓ All pre-commit checks passed! (completed in ${duration}s)"
  blue "========================================================"
  exit 0
else
  red "✗ Pre-commit checks failed (completed in ${duration}s)"
  yellow "Fix the issues above before committing your changes."
  blue "========================================================"
  exit 1
fi