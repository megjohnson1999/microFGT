# microFGT Package Refactoring Summary

This document summarizes the refactoring work done to clean up the codebase and reduce duplication.

## Consolidation of Example Data Functions

The following duplicate files were consolidated into a single comprehensive file:

- `R/example_data.R`: Now contains all example data generation, loading, and exporting functionality
- `R/constants/cst_constants.R`: Contains constants used for generating realistic FGT microbiome data

### Files Removed (Duplicates)
- `R/generate_example_data.R`
- `R/load_example_data.R`
- `R/export_example_data.R`
- `R/example_data_export.R`
- `R/simple_load_example.R`
- `R/export_functions.R`
- `R/example_data/constants.R`

### Functions Consolidated

All example data functions are now in `R/example_data.R`:

1. **`generate_fgt_example_data()`**: Creates synthetic microbiome data with realistic FGT community state types
2. **`load_example_data()`**: Loads pre-built example datasets from the package
3. **`export_example_data()`**: Exports data to common bioinformatics formats (phyloseq, DADA2, QIIME2)
4. Helper functions for exporting to different formats:
   - `export_to_phyloseq()`
   - `export_to_dada2()`
   - `export_to_qiime2()`

### Constants Reorganization

Constants are now properly organized in `R/constants/cst_constants.R`, including:
- CST (Community State Type) taxonomic compositions
- Clinical parameter distributions
- Taxonomic ranks
- CST descriptions

## Example Data Generation Script

Added a script to generate example datasets for the package:
- `data-raw/prepare_example_data.R`

This script creates example datasets in `inst/extdata/` that can be loaded with `load_example_data()`.

## Benefits of Refactoring

1. **Reduced duplication**: Eliminated multiple files with similar or identical functions
2. **Better organization**: Grouped related functionality together
3. **Consistent API**: Now has a clear, consistent API for working with example data
4. **Maintainability**: Easier to maintain a single file than multiple scattered files
5. **Documentation**: All related functions are documented together

## Next Steps

1. **Update vignette**: Add documentation on how to use the example data functions
2. **Improve constants integration**: More fully integrate the constants into the example data generation
3. **Generate pre-built examples**: Run the `prepare_example_data.R` script to generate example datasets 
4. **Add example data tests**: Create tests to verify the example data functions work correctly