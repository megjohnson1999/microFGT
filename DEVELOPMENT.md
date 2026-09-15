# microFGT Development Guide

## Project Status

microFGT is currently in **alpha stage development**. This means:

- Core data structures are implemented but may have bugs
- Basic functionality is coded but not fully tested
- APIs may change significantly before release
- Documentation is incomplete

## Known Issues and Workarounds

### FGTExperiment Constructor Issue

**Issue**: The `FGTExperiment()` constructor may fail with "invalid rownames length" error.

**Workaround**:
```r
# Create FGTExperiment objects using this two-step approach:
library(microFGT)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)

# Step 1: Create a SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = your_count_matrix),
  rowData = your_feature_metadata,  # Optional
  colData = your_sample_metadata    # Optional
)

# Step 2: Convert to TreeSummarizedExperiment
tse <- TreeSummarizedExperiment(se)

# Step 3: Create FGTExperiment directly
fgt <- methods::new("FGTExperiment", 
                  tse,
                  experimentType = "amplicon",
                  fgtMetadata = S4Vectors::SimpleList())

# Step 4: Ensure the assay is properly named
assayNames(fgt)[1] <- "counts"
```

### Function Parameter Mismatches

**Issue**: Function parameter names may not match documentation.

**Correct Parameter Names**:
- `transform_abundance()`: Use `type` instead of `method`
- `filter_taxa()`: Include `assay_name = "counts"`
- `aggregate_taxa()`: Include `assay_name = "counts"`
- `plot_alpha_diversity()`: Use `metrics` and `group_var` instead of `index` and `by`

## Development Roadmap

### Current Phase: Core Functionality

- Fix constructor issues
- Implement proper test suite
- Ensure basic data manipulation works
- Add error handling and validation

### Next Phase: Tool Integration

- Implement speciateIT integration
- Add VALENCIA community state typing
- Integrate VIRGO metagenomic pipeline

### Final Phase: Analysis Module Development

- Taxonomic analysis functions
- Community typing analysis
- Multi-omics integration

## Testing Approach

- Unit tests for individual functions
- Integration tests for workflows
- Manual validation scripts in `inst/validation/`

## Development Guidelines

1. **Error Handling**: All functions should validate inputs and provide informative error messages
2. **Documentation**: All functions need roxygen documentation with examples
3. **Testing**: New features should include tests
4. **Performance**: Consider using Rcpp for performance-critical functions
5. **Compatibility**: Maintain compatibility with Bioconductor ecosystem

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for your changes
4. Submit a pull request

## Resources

- [TreeSummarizedExperiment Documentation](https://bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html)
- [R Package Development Guide](https://r-pkgs.org/)
- [Bioconductor Guidelines](https://bioconductor.org/developers/package-guidelines/)