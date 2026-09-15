# Import/Export Functions in the Hybrid Approach

## Overview

The import/export functions in microFGT follow the hybrid approach design pattern, allowing seamless data interchange between different microbiome analysis tools while maintaining compatibility with both SummarizedExperiment and FGTExperiment objects.

## Design Goals

1. **Universal Compatibility**: Functions work with any standard microbiome format and SummarizedExperiment-like objects
2. **Format Flexibility**: Support for common formats used in microbiome research (BIOM, QIIME2, phyloseq, DADA2)
3. **Graceful Degradation**: Functions adapt to available packages and object types
4. **Type Preservation**: Original object types are maintained or enhanced when possible

## Core Import Functions

### `import_microbiome()`

The primary import function supporting multiple data formats:

```r
# Import from matrices or data frames
se <- import_microbiome(counts_matrix, taxonomy_df, metadata_df)

# Import from files
se <- import_microbiome("counts.csv", "taxonomy.csv", "metadata.csv")

# Import with taxonomic strings
se <- import_microbiome(counts_matrix, tax_strings)
```

### Format-Specific Import Functions

- `import_from_dada2()`: Import from DADA2 sequence table and taxonomy
- `import_from_biom()`: Import from BIOM format files
- `import_from_phyloseq()`: Import from phyloseq objects

## Core Export Functions

### `export_microbiome()`

The primary export function supporting multiple output formats:

```r
# Export to CSV files
export_microbiome(se, "output_dir", format = "csv")

# Export to RDS files
export_microbiome(se, "output_dir", format = "rds")

# Export to BIOM format
export_microbiome(se, "output_dir", format = "biom")

# Export to QIIME2 format
export_microbiome(se, "output_dir", format = "qiime2") 

# Export to phyloseq format
export_microbiome(se, "output_dir", format = "phyloseq")

# Export to DADA2 format
export_microbiome(se, "output_dir", format = "dada2")
```

### Conversion Functions

- `to_phyloseq()`: Convert to phyloseq object
- `to_biom()`: Convert to BIOM object

## Implementation Features

1. **Automatic Detection**:
   - Input format detection (matrix, data.frame, file path)
   - Matrix orientation detection (features in rows vs. columns)
   - Taxonomy format detection (data.frame vs. strings)

2. **Flexible Taxonomy Handling**:
   - Support for standard taxonomic ranks
   - Parsing of taxonomic strings with different formats
   - Handling of non-standard taxonomic classifications

3. **Multiple Output Formats**:
   - Standard formats: CSV, RDS, BIOM
   - Bioinformatics formats: QIIME2, phyloseq, DADA2
   - Customizable output naming and compression options

## Integration with SummarizedExperiment

These functions seamlessly integrate with SummarizedExperiment by:

1. **Using Standard Accessors**: assays(), rowData(), colData(), metadata()
2. **Preserving Metadata**: Maintaining information about data provenance
3. **Supporting Extensions**: Working with TreeSummarizedExperiment and FGTExperiment

## Example Workflow

A typical import/export workflow:

```r
# Import data from DADA2
se <- import_from_dada2("seqtab.rds", "taxa.rds", "metadata.csv")

# Perform analysis
se <- transform_abundance(se, method = "relative")
se <- filter_features(se, min_prevalence = 0.1)
genus_level <- aggregate_taxa(se, rank = "Genus")

# Export results to multiple formats
export_microbiome(genus_level, "results/qiime2", format = "qiime2")
export_microbiome(genus_level, "results/r", format = "rds")

# Convert to phyloseq for additional analysis
ps <- to_phyloseq(genus_level)
```

## Benefits of the Hybrid Approach

1. **Interoperability**: Seamless exchange with other tools and packages
2. **Flexibility**: Works with any SummarizedExperiment-like object
3. **Graceful Fallbacks**: Adapts to available packages
4. **Extended Functionality**: Provides specialized handling for FGTExperiment
5. **Format Consistency**: Ensures data integrity across conversions