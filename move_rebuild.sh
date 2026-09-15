#!/bin/bash

# This script moves the rebuilt package files to the main directory

# Define directories
REBUILD_DIR="../rebuild"
MAIN_DIR="."

# Create a backup of the current code
echo "Creating backup of current code..."
BACKUP_DIR="../microFGT_backup_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"
cp -r * "$BACKUP_DIR/"
echo "Backup created at $BACKUP_DIR"

# Copy new files from rebuild to current directory
echo "Copying files from rebuild directory..."

# First, clear out files we're going to replace, preserving specified ones
for file in R/* src/* man/* tests/* vignettes/*; do
  # Skip non-existent files (if glob doesn't match)
  [ -e "$file" ] || continue
  
  # Skip the files we want to keep
  if [[ "$file" == "microFGT_logo.png" || 
        "$file" == "microFGT_development_prompt.md" ||
        "$file" == "microFGT_architecture.md" || 
        "$file" == "microFGT_implementation_guide.md" ]]; then
    echo "Preserving $file"
    continue
  fi
  
  rm -f "$file"
done

# Copy new files
echo "Copying new files..."
cp -r "$REBUILD_DIR/R" "$MAIN_DIR/"
cp -r "$REBUILD_DIR/inst" "$MAIN_DIR/"
cp -r "$REBUILD_DIR/man" "$MAIN_DIR/"
cp -r "$REBUILD_DIR/tests" "$MAIN_DIR/"
cp -r "$REBUILD_DIR/vignettes" "$MAIN_DIR/"
cp -r "$REBUILD_DIR/src" "$MAIN_DIR/" 2>/dev/null || true  # May not exist

# Copy root files
cp "$REBUILD_DIR/DESCRIPTION" "$MAIN_DIR/"
cp "$REBUILD_DIR/NAMESPACE" "$MAIN_DIR/"
cp "$REBUILD_DIR/LICENSE" "$MAIN_DIR/"
cp "$REBUILD_DIR/README.md" "$MAIN_DIR/"

echo "Files copied successfully!"
echo ""
echo "Next steps:"
echo "1. Review the changes with 'git status'"
echo "2. Commit the changes to git"
echo "3. Push to GitHub"
echo "4. Install the updated package with 'devtools::install()'"

exit 0