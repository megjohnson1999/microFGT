# Contributing to microFGT

Thank you for your interest in contributing to microFGT! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

Please be respectful and considerate of others when participating in this project. We expect all contributors to adhere to professional standards of communication and behavior.

## How to Contribute

### Reporting Bugs

If you find a bug, please create an issue on GitHub with the following information:

1. A clear, descriptive title
2. Steps to reproduce the issue
3. Expected and actual behavior
4. R version, microFGT version, and operating system
5. Any relevant error messages or screenshots

### Suggesting Features

We welcome feature suggestions! Please create an issue with:

1. A clear description of the feature
2. The rationale behind the feature request
3. Any relevant scientific references or use cases

### Pull Requests

We welcome code contributions through pull requests. To submit a PR:

1. Fork the repository
2. Create a new branch for your feature/fix
3. Make your changes, following the coding standards
4. Ensure tests pass and add new tests if needed
5. Submit a pull request with a clear description of the changes

## Development Workflow

We follow a branching model with:

- `main`: Stable release branch
- `develop`: Development branch for integrating features
- Feature branches: Named as `feature/short-description`
- Bugfix branches: Named as `fix/short-description`

## Coding Standards

Please follow these standards when contributing code:

1. Code Style:
   - Follow tidyverse style guide
   - Use meaningful variable and function names
   - Include roxygen2 documentation for all functions

2. Testing:
   - Add tests for new functions using testthat
   - Ensure all tests pass before submitting PRs

3. Commits:
   - Use informative commit messages
   - Reference issue numbers when appropriate

## Package Structure

The package is organized by functionality:

- `R/core/`: Core class definitions and constructors
- `R/data/`: Data manipulation and import/export
- `R/utils/`: Utility functions
- `R/visualization/`: Plotting and visualization
- `R/constants/`: Package constants 

## Documentation

Please document your code with:

1. Roxygen2 comments for function documentation
2. Examples for all exported functions
3. Clear in-line comments for complex code

## Getting Help

If you need assistance, you can:

- Open an issue with your question
- Contact the maintainers directly

Thank you for contributing to microFGT!