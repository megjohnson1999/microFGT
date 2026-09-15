# CI/CD Implementation for microFGT

This document outlines the continuous integration and continuous deployment (CI/CD) pipeline set up for the microFGT package.

## Overview

The CI/CD pipeline is implemented using GitHub Actions and consists of several workflows that handle different aspects of package testing, documentation, and validation.

## Workflows

### 1. R CMD Check

**File:** `.github/workflows/R-CMD-check.yaml`

This workflow runs the standard R package checks (`R CMD check`) across multiple platforms and R versions:

- Windows (R release)
- macOS (R release)
- Ubuntu (R release)
- Ubuntu (R devel)
- Ubuntu (R 4.0)

The workflow is triggered on:
- Push to main and develop branches
- Pull requests to main and develop branches
- Weekly schedule (Monday at 4:00 UTC)

### 2. Test Coverage

**File:** `.github/workflows/test-coverage.yaml`

This workflow calculates test coverage and reports it to [codecov.io](https://codecov.io/). It helps track how much of the codebase is covered by tests.

The workflow is triggered on:
- Push to main and develop branches
- Pull requests to main and develop branches
- Weekly schedule (Monday at 4:00 UTC)

### 3. Package Documentation

**File:** `.github/workflows/pkgdown.yaml`

This workflow builds and deploys the package documentation website using pkgdown. The site is deployed to GitHub Pages.

The workflow is triggered on:
- Push to main branch
- Release tags
- Manual trigger

### 4. Dependency Review

**File:** `.github/workflows/dependency-review.yaml`

This workflow reviews dependencies for security vulnerabilities and updates. It checks for high severity issues in dependencies when pull requests are opened.

## Configuration Files

### codecov.yml

Configures the code coverage reporting settings for codecov.io:
- Disables PR comments
- Sets thresholds for project and patch coverage
- Makes coverage checks informational

### _pkgdown.yml

Configures the package documentation website:
- Sets theme and appearance
- Defines navigation structure
- Organizes function reference into categories

## GitHub Templates

### Pull Request Template

**File:** `.github/pull_request_template.md`

Standardizes pull request submissions with:
- Description field
- Type of change checklist
- Testing information
- Implementation details

### Issue Templates

1. **Bug Report**: `.github/ISSUE_TEMPLATE/bug_report.yml`
   - Collects version information
   - Requests reproducible examples
   - Captures error messages

2. **Feature Request**: `.github/ISSUE_TEMPLATE/feature_request.yml`
   - Captures feature description
   - Documents use cases
   - Requests implementation suggestions

### CODEOWNERS

**File:** `.github/CODEOWNERS`

Defines who is responsible for reviewing changes to different parts of the codebase:
- Core implementation files
- Documentation
- Testing
- CI and build configuration

## Usage for Developers

### Local Testing

Before pushing changes, you can run the same checks locally:

```r
# Run R CMD check
devtools::check()

# Calculate test coverage
covr::package_coverage()

# Build documentation
pkgdown::build_site()
```

### Pull Request Process

1. Make changes on a feature branch
2. Run local checks
3. Open a pull request using the template
4. Review CI workflow results
5. Address any issues
6. Request review from appropriate CODEOWNERS

### Release Process

1. Ensure all tests pass on develop branch
2. Merge develop into main
3. Tag a new release
4. The pkgdown workflow will automatically deploy updated documentation

## Troubleshooting

### Common Issues

#### Multiple Workflow Instances Running

If multiple instances of the same workflow start running simultaneously after a push, there are several potential causes:

1. **Duplicate Workflow Files**:
   - GitHub Actions recognizes both `.yml` and `.yaml` extensions
   - Having both (e.g., `R-CMD-check.yml` and `R-CMD-check.yaml`) will cause duplicate runs
   - Solution: Remove duplicate workflow files

   ```bash
   git rm .github/workflows/duplicate-file.yml
   git commit -m "Remove duplicate workflow file"
   git push
   ```

2. **Workflow Concurrency Issues**:
   - GitHub Actions may process workflow files in parallel
   - Modifying workflow files can sometimes trigger the workflow itself
   - Solution: Add concurrency group to limit simultaneous runs

   ```yaml
   # Add this to your workflow
   concurrency:
     group: ${{ github.workflow }}-${{ github.ref }}
     cancel-in-progress: true
   ```

3. **Path Filtering**:
   - Use `paths-ignore` to skip workflow runs for irrelevant changes
   - Particularly useful for documentation and CI configuration changes

   ```yaml
   on:
     push:
       branches: [main, develop]
       paths-ignore:
         - '**.md'
         - '.github/workflows/**'
   ```

4. **Self-Referential Triggers**:
   - Changes to workflow files can trigger the workflows themselves
   - Create a separate workflow for CI/CD file validation
   - Use path filtering to prevent cycles

## Future Improvements

Planned enhancements to the CI/CD pipeline:

1. Add performance benchmarking tests
2. Implement dependency scanning and updates
3. Set up deployment to test CRAN-like repositories
4. Create automated issue reproduction tests
5. Add installation validation workflows