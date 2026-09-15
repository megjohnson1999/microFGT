# Taxonomic Functions in the Hybrid Approach

## Overview

The taxonomic functions in microFGT follow the hybrid approach by implementing standalone functions that work with any SummarizedExperiment-like object, while providing enhanced functionality for FGTExperiment objects when available.

## Key Features

1. **Universal Compatibility**: All functions work with standard SummarizedExperiment objects, TreeSummarizedExperiment objects, and FGTExperiment objects
2. **Preservation of Object Type**: Functions return the same type of object that was provided as input
3. **Minimal Dependencies**: Only requires SummarizedExperiment and S4Vectors
4. **Comprehensive Taxonomic Operations**: Supports aggregation, normalization, and formatting

## Core Functions

### `get_taxonomic_ranks()`

Identifies taxonomic rank columns in a SummarizedExperiment object's rowData.

```r
ranks <- get_taxonomic_ranks(se)
```

### `aggregate_taxa()`

Aggregates feature abundances at a specified taxonomic rank (e.g., genus, family).

```r
genus_level <- aggregate_taxa(se, rank = "genus")
```

### `normalize_taxonomy()`

Cleans and standardizes taxonomy names in rowData.

```r
clean_se <- normalize_taxonomy(se, remove_prefixes = TRUE)
```

### `create_tax_strings()`

Creates formatted taxonomic strings from individual rank columns.

```r
se_with_tax <- create_tax_strings(se, format = "lineage")
```

### `parse_tax_strings()`

Parses taxonomy strings into individual rank columns.

```r
parsed_se <- parse_tax_strings(se, tax_column = "taxonomy") 
```

## Integration with SummarizedExperiment

These functions integrate with SummarizedExperiment by:

1. **Operating on rowData**: Using standard rowData accessors
2. **Preserving Object Integrity**: Maintaining the original object type
3. **Working with Assays**: Supporting operations on any assay in the object

## Example Workflow

A typical taxonomic analysis workflow:

```r
# Start with a SummarizedExperiment or FGTExperiment object
se <- SummarizedExperiment(...)

# Clean up taxonomy
se <- normalize_taxonomy(se)

# Aggregate at genus level
genus_level <- aggregate_taxa(se, rank = "genus")

# Create formatted taxonomy strings for visualization
genus_level <- create_tax_strings(genus_level, format = "text")

# Plot composition at genus level
plot_composition(genus_level, tax_column = "taxonomy")
```

## Benefits of the Hybrid Approach

1. **Reusable Components**: Functions can be used with any SummarizedExperiment-like object
2. **Reduced Complexity**: No need to extend the S4 class system
3. **Easier Testing**: Functions can be tested in isolation
4. **Better Encapsulation**: Each function has a single, well-defined purpose
5. **Compatibility**: Works well with both Bioconductor classes and custom code