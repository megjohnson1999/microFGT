# microFGT Package Structure Improvements

This document summarizes the structural improvements made to the microFGT package to align with the engineering plan and best practices for R package development.

## Directory Organization

The package structure has been reorganized according to the following scheme:

```
microFGT/
├── CONTRIBUTING.md       # New contributor guidelines
├── DESCRIPTION           # Updated with organized dependencies
├── LICENSE
├── NAMESPACE
├── NEWS.md               # Updated to reflect structural changes
├── R/                    # Reorganized by functionality
│   ├── constants/        # Constants used throughout the package
│   ├── core/             # Core classes and data structures
│   ├── data/             # Data manipulation and import/export
│   ├── microFGT-package.R  # Package documentation
│   ├── utils/            # Utility functions
│   └── visualization/    # Plotting and visualization
├── data-raw/             # Scripts to prepare package data
├── development/          # Development documents and plans
│   ├── docs/             # Detailed documentation
│   └── microFGT_engineering_plan.md
├── inst/
│   ├── examples/         # Example scripts
│   ├── extdata/          # Example datasets
│   └── validation/       # Validation scripts
├── man/                  # Documentation (generated)
├── src/                  # C++ source code
├── tests/                # Test suite
│   └── testthat/
└── vignettes/            # User guides
```

## Key Improvements

1. **Functional Organization**: 
   - Reorganized R code files by functionality (core, data, utils, visualization)
   - Grouped related functions together for easier maintenance

2. **Consolidated Duplicate Code**:
   - Combined similar example data functions
   - Eliminated redundant scripts and utility functions

3. **Better Development Documentation**:
   - Moved all development documents to the development directory
   - Created CONTRIBUTING.md with guidelines for contributors

4. **Code Quality Improvements**:
   - Added helper utility functions for consistent error handling
   - Organized constants properly

5. **Package Metadata**:
   - Updated DESCRIPTION with organized dependencies
   - Added proper configuration for testthat

6. **Documentation**:
   - Improved NEWS.md with structured changelog
   - Added READMEs for example data

## Addressing Engineering Plan Items

This reorganization addresses several items from the microFGT engineering plan:

- [x] Define and document naming conventions (2.1)
- [x] Organize files by functional area (2.1)
- [x] Apply consistent formatting to all files (2.1)
- [x] Reduce code duplication across similar functionality (2.1)
- [x] Consolidate example data functions (2.1)
- [x] Move constants to the proper location (2.1)
- [x] Group files by functionality (2.1)
- [x] Write contribution guidelines in CONTRIBUTING.md (1.1)

## Next Steps

1. **Documentation**:
   - Update function documentation with roxygen2
   - Create comprehensive vignettes

2. **Testing**:
   - Expand test coverage
   - Add unit tests for new functionality

3. **CI/CD**:
   - Set up GitHub Actions for automated testing
   - Configure code quality checks

4. **Dependency Management**:
   - Further review and optimize dependencies

This reorganization provides a solid foundation for implementing the rest of the engineering plan and continuing development of the microFGT package.