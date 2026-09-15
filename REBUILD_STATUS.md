# microFGT Package Rebuild Status

## Overview

This document tracks the status of the microFGT package rebuild. The package has undergone a complete restructuring to improve code organization, implement proper S4 class design, and enhance functionality.

## Current Status

- Complete rebuild with proper S4 class structure
- Implementation of mock data generators for SpeciateIT, VIRGO, and VALENCIA
- Updated CI/CD workflows and package validation
- Improved transformation function implementation
- **NEW**: Integrated R-rebuild code with compositional class design

## Main Components

1. **Core Infrastructure**
   - S4 Class Structure: FGTExperiment class
   - Constructor and accessor methods
   - Data structure validation
   - Improved composition-based architecture

2. **Mock Data Generators**
   - SpeciateIT mock data generator
   - VIRGO mock data generator
   - VALENCIA mock data generator
   - Parameterized generation capabilities

3. **Functional Components**
   - Diversity functions (alpha and beta diversity)
   - Taxonomic analysis functions
   - Import/Export functions
   - Transformation functions

## TODO

- [ ] Complete test coverage for all new functions
- [ ] Finalize documentation for all public functions
- [ ] Create comprehensive vignettes with usage examples
- [ ] Add additional validation checks for input data
- [ ] Implement any remaining functionality from the planning documents
- [ ] Run full test suite with new compositional class structure
- [ ] Update existing functions to work with new class design

## Recently Completed

- [x] Force push updated codebase to remote repository
- [x] Clean up repository files and update .gitignore
- [x] Create mock data generators for all three tools
- [x] Implement proper S4 class structure following best practices
- [x] Integrate R-rebuild code into main codebase
- [x] Create migration function for converting legacy FGTExperiment objects
- [x] Implement compositional class design (containing not extending TSE)
- [x] Add transformAbundance generic and methods

## Architectural Changes

### New Compositional Class Design

The FGTExperiment class has been redesigned to use composition rather than inheritance:

**Previous Design (Inheritance):**
- FGTExperiment extends TreeSummarizedExperiment
- Direct access to TSE slots and methods

**New Design (Composition):**
- FGTExperiment contains a TreeSummarizedExperiment
- Better encapsulation and separation of concerns
- More robust method dispatch
- Simplified validation

This architectural change provides better flexibility and maintainability while still preserving the core functionality and API.

## Directory Structure

The package follows standard R package conventions with the following organization:

- **R/**: R source code
  - **core/**: Core class definitions and infrastructure
  - **data/**: Data manipulation and import/export
  - **constants/**: Package constants
  - **utils/**: Utility functions
  - **visualization/**: Plotting and visualization

- **src/**: C++ source code for performance-critical functions

- **inst/examples/**: Example scripts

- **tests/testthat/**: Unit tests

- **development/**: Development documentation and planning files

## Documentation

Comprehensive documentation is available in:
- Function documentation (via roxygen)
- Vignettes (in the vignettes/ directory)
- Development plans (in the development/ directory)