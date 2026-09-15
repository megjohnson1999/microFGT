# microFGT Implementation: Next Steps

Based on our assessment of Phase 2 of the engineering plan, this document outlines the specific next steps for the microFGT package development, focusing on the highest priority items.

## Priority 1: S4 Class Design

### Task 1.1: Consolidate Class Definitions
- Merge duplicated class definitions from `R/core/FGTExperiment-class.R` and `R/core/data_structures.R`
- Keep only one comprehensive implementation
- Ensure all methods are properly documented in a single file

### Task 1.2: Implement Validation Methods
- Add `setValidity()` for `FGTExperiment` class
- Validate slot types and values:
  - `experimentType` must be one of "amplicon", "metagenomic", or "integrated"
  - `fgtMetadata` must be a SimpleList with valid contents
- Add helpful error messages for invalid class construction

### Task 1.3: Add Coercion Methods
- Implement `as()` and `coerce()` methods for:
  - `FGTExperiment` to/from `TreeSummarizedExperiment`
  - `FGTExperiment` to/from `SummarizedExperiment`
  - `FGTExperiment` to/from `phyloseq` (if phyloseq is available)
- Add thorough testing for coercion methods

### Task 1.4: Standardize Constructor
- Keep only the most comprehensive constructor implementation
- Enhance error handling and input validation
- Document constructor parameters extensively
- Add examples covering different construction scenarios

## Priority 2: Function Design Standardization

### Task 2.1: Establish Parameter Naming Guidelines
- Create a document defining standard parameter names:
  - Use `fgt_exp` consistently for FGTExperiment objects
  - Use snake_case for all parameter names
  - Define standard names for common parameters
  - Document parameter ordering (required then optional)

### Task 2.2: Implement Consistent Error Handling
- Create standardized error handling utilities in `R/utils/helpers.R`:
  - Enhance `format_error()` function
  - Add `validate_input()` function
  - Add `check_required_param()` function
- Update existing functions to use these utilities

### Task 2.3: Standardize Return Values
- Define return value patterns:
  - Functions modifying FGTExperiment should return FGTExperiment
  - Functions extracting data should return consistent types
  - Add operation metadata to track transformations
- Document expected return types clearly

### Task 2.4: Implement Input Validation
- Add comprehensive validation to:
  - All core data manipulation functions
  - Visualization functions
  - Analysis functions
- Use standardized validation utilities

## Priority 3: Dependency Management

### Task 3.1: Audit Dependencies
- Create a document listing each dependency with:
  - Purpose in the package
  - Required functions/methods
  - Whether it's essential or optional
  - Minimum version required

### Task 3.2: Optimize Import Structure
- Move non-essential packages from Imports to Suggests:
  - Review BiocParallel usage
  - Assess tidyverse components
  - Check visualization dependencies
- Update DESCRIPTION file accordingly

### Task 3.3: Enhance Dependency Handling
- Improve `check_package()` and `require_packages()` functions:
  - Add version checking
  - Provide more helpful installation instructions
  - Implement caching to avoid repeated checks
- Apply consistently throughout the package

### Task 3.4: Implement Fallbacks
- Add fallback functionality when optional packages are missing:
  - Provide basic alternatives for visualization
  - Add graceful degradation for advanced features
  - Clearly document limitations when packages are missing

## Implementation Schedule

### Week 1: S4 Class Design
- Days 1-2: Task 1.1 (Consolidate Class Definitions)
- Days 3-4: Task 1.2 (Implement Validation Methods)
- Day 5: Task 1.3 (Add Coercion Methods)

### Week 2: Function Design & Constructor
- Days 1-2: Task 1.4 (Standardize Constructor)
- Day 3: Task 2.1 (Establish Parameter Naming Guidelines)
- Days 4-5: Task 2.2 (Implement Consistent Error Handling)

### Week 3: Continue Function Design & Dependencies
- Days 1-2: Task 2.3 & 2.4 (Standardize Returns & Implement Input Validation)
- Day 3: Task 3.1 (Audit Dependencies)
- Days 4-5: Tasks 3.2 - 3.4 (Optimize Imports, Enhance Handling, Implement Fallbacks)

## Success Criteria

1. All duplicated class code is consolidated
2. All functions have consistent parameter names and validation
3. Dependencies are properly documented and managed
4. All code passes R CMD check without warnings

Upon completion of these tasks, we will have a solid foundation for moving to Phase 3 of the engineering plan, focusing on the comprehensive testing framework.