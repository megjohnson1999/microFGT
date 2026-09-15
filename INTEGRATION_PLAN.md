# R-rebuild Integration Plan

## Overview

This document outlines the plan for integrating the R-rebuild code with the main R codebase.

## Key Changes

1. **Class Structure Change**: The R-rebuild version uses a compositional approach where FGTExperiment *contains* a TreeSummarizedExperiment rather than extending it. This is a significant architectural change that will require careful migration.

2. **Improved Method Implementations**: The rebuild contains enhanced transformation methods and better accessor handling.

## Integration Steps

### 1. Class Definition

- [ ] Create a backup of current `R/core/FGTExperiment-class.R`
- [ ] Replace with R-rebuild's compositional design from `R-rebuild/R/AllClasses.R`
- [ ] Update NAMESPACE exports

### 2. Generic Definitions

- [ ] Integrate generics from `R-rebuild/R/AllGenerics.R` into `R/AllGenerics.R`
- [ ] Remove duplicates and ensure proper documentation

### 3. Constructor Function

- [ ] Replace existing constructor in `R/core/FGTExperiment-class.R` with the rebuild version from `R-rebuild/R/FGTExperiment-constructor.R`
- [ ] Update documentation to reflect the new structure

### 4. Methods Implementation

- [ ] Integrate methods from `R-rebuild/R/methods-FGTExperiment.R`
- [ ] Add missing coercion methods
- [ ] Ensure compatibility with existing code that may expect direct inheritance

### 5. Transform Functions

- [ ] Integrate transformation methods from `R-rebuild/R/transform.R`
- [ ] Ensure backward compatibility with existing transform functions

### 6. Package Initialization

- [ ] Update `R/utils/zzz.R` with the improvements from `R-rebuild/R/zzz.R`

### 7. Migration Support

- [ ] Create migration functions to convert between old and new class structures if needed
- [ ] Add deprecation warnings for functions that have changed significantly

## Testing Strategy

1. After each integration step, run tests to validate functionality
2. Create specific tests for the compositional class structure
3. Create tests that verify backward compatibility with code expecting inheritance
4. Test all transformation functions to ensure consistent behavior

## Documentation Updates

- [ ] Update class documentation to reflect the new structure
- [ ] Update vignettes to show the proper usage of the new class
- [ ] Add migration guide for users of the previous version

## Timeline

1. Class definition and generics: 1 day
2. Constructor and basic methods: 1 day
3. Transform functions: 1 day
4. Testing and debugging: 2 days
5. Documentation updates: 1 day

Total estimated time: 6 days