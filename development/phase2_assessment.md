# microFGT Engineering Plan: Phase 2 Assessment

This document provides an analysis of the current state of Phase 2 (Code Quality & Architecture) of the microFGT engineering plan, identifying what has been accomplished and what needs to be addressed next.

## 2.1 Code Style & Organization ✅

This section has been largely completed:

- ✅ Defined and documented naming conventions
- ✅ Organized files by functional area
- ✅ Applied consistent formatting to all files
- ✅ Reviewed and refactored internal helper functions
- ✅ Reduced code duplication across similar functionality
- ✅ Consolidated example data functions
- ✅ Moved constants to dedicated location
- ✅ Grouped files by functionality (core, data, utils, visualization)

## 2.2 S4 Class Design ⏳

This section requires focused attention:

### Current State
- `FGTExperiment` class extends `TreeSummarizedExperiment`
- Basic slots defined: `experimentType` and `fgtMetadata`
- Simple accessor/mutator methods implemented
- Custom show method created

### Issues Identified
- ⚠️ Duplicated class definitions across multiple files
- ⚠️ No formal validation methods implemented
- ⚠️ Missing coercion methods for compatibility
- ⚠️ Multiple constructor implementations with inconsistent behavior
- ⚠️ Limited documentation of class hierarchy

### Next Steps
1. **Consolidate class definitions** into a single file
2. **Implement proper validation** method with `setValidity()`
3. **Add coercion methods** for compatible classes
4. **Standardize constructors** - maintain one with comprehensive error handling
5. **Enhance documentation** with inheritance diagram and examples

## 2.3 Function Design ⏳

This section needs significant work:

### Parameter Naming Inconsistencies
- Parameter naming varies across similar functions:
  - `fgt_exp` vs. `x` vs. `physeq` for similar parameters
  - Inconsistent parameter ordering
  - Mixed use of camelCase and snake_case

### Return Value Issues
- Inconsistent return value documentation
- Varying return types for similar operations
- Inconsistent metadata handling for operation tracking

### Error Handling Gaps
- Multiple validation approaches used
- Error message formatting varies
- Utility function `format_error()` exists but is rarely used

### Input Validation Problems
- Different validation implementations across functions
- Inconsistent NULL value handling
- Multiple approaches to required package checking

### Next Steps
1. **Standardize parameter naming** across all functions
2. **Establish consistent return value patterns**
3. **Implement uniform error handling** using utility functions
4. **Enforce input validation** in all functions

## 2.4 Dependency Management ⏳

This section requires thorough review:

### Current State
- 13 packages in Imports
- 9 packages in Suggests
- Basic dependency checking through helper functions
- Some handling of optional packages

### Issues Identified
- ⚠️ No clear justification documented for each dependency
- ⚠️ Some imports may be better suited as suggestions
- ⚠️ Inconsistent fallback mechanisms for optional dependencies
- ⚠️ Minimal version requirements specified

### Next Steps
1. **Audit all dependencies** and document their necessity
2. **Minimize core dependencies** by moving non-essential packages to Suggests
3. **Implement better fallbacks** for optional functionality
4. **Specify minimum versions** for all dependencies
5. **Standardize dependency checking** throughout the codebase

## Implementation Priorities

Based on this assessment, we recommend the following implementation priorities:

### Immediate (1-2 weeks)
1. **S4 Class Design**
   - Consolidate class definitions
   - Implement validation methods
   - Standardize constructors

2. **Function Design**
   - Create parameter naming guidelines
   - Implement standard error handling utilities
   - Begin updating key functions

### Short-term (2-4 weeks)
1. **Dependency Management**
   - Audit and document dependencies
   - Move non-essential packages to Suggests
   - Enhance fallback mechanisms

2. **Complete Function Standardization**
   - Apply consistent patterns to all functions
   - Improve return value consistency
   - Enhance validation across all functions

## Measuring Success

Progress in Phase 2 can be measured by:

1. **S4 Classes**
   - Number of classes with proper validation
   - Percentage of methods with complete documentation
   - Reduction in duplicated code

2. **Function Design**
   - Percentage of functions following naming conventions
   - Error handling consistency score
   - Parameter validation coverage

3. **Dependency Management**
   - Core dependency count reduction
   - Documentation coverage for dependencies
   - Version specification completeness

## Conclusion

While significant progress has been made in code organization (2.1), there is still substantial work needed in S4 class design (2.2), function design (2.3), and dependency management (2.4). Focusing on these areas next will ensure the package has a solid foundation before moving on to Phase 3 (Testing Framework).