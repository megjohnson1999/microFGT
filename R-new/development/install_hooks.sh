#!/bin/bash
# Script to install git hooks for microFGT package

# Print colored output
green() { echo -e "\033[0;32m$1\033[0m"; }
yellow() { echo -e "\033[0;33m$1\033[0m"; }
red() { echo -e "\033[0;31m$1\033[0m"; }
blue() { echo -e "\033[0;34m$1\033[0m"; }

# Show header
blue "========================================================"
blue "       Installing Git Hooks for microFGT Package         "
blue "========================================================"

# Create hooks directory if it doesn't exist
if [ ! -d ".git/hooks" ]; then
  mkdir -p .git/hooks
  green "Created .git/hooks directory"
fi

# Copy the pre-commit hook to .git/hooks
cp development/pre-commit-hook.sh .git/hooks/pre-commit
if [ $? -ne 0 ]; then
  red "Failed to copy pre-commit hook"
  exit 1
fi

# Make the pre-commit hook executable
chmod +x .git/hooks/pre-commit
if [ $? -ne 0 ]; then
  red "Failed to make pre-commit hook executable"
  exit 1
fi

green "✓ Pre-commit hook installed successfully"
yellow "The hook will run these checks before each commit:"
echo "  - Roxygen documentation update"
echo "  - NAMESPACE validation"
echo "  - S4 class validation"
echo "  - Run unit tests"
echo "  - Code style check (if styler is installed)"

yellow "Note: To bypass the pre-commit hook temporarily, use:"
echo "  git commit --no-verify"

blue "========================================================"