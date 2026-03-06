#!/bin/bash 
set -euo pipefail

# SCRIPT NAME: mk-branch.sh
# PURPOSE: This script creates a new git branch, sets it to track the current branch, and pushes it to the remote repository.
#          It also checks for uncommited changes before proceeding 
# AUTHOR: Matthew Hirschberger
# DATE: 06 Mar 2026
# REVISION: 1.0
# ----------------------------------------------------------------------

# Check for uncommited changes.
# -n is the string test operator which checks if the output string is non-empty
if [[ -n $(git status --porcelain) ]]; then
    echo "Error: You have uncommited changes. Please commit or stash them before creating a new branch."
    exit 1
fi

# Get the current branch name
current_branch=$(git rev-parse --abbrev-ref HEAD)

# Prompt for new branch name and read input
read -p -r "Enter new branch name: " new_branch

# Create the new branch and set it to track the current branch
git checkout -b "$new_branch" "$current_branch"

# Push the new branch to the remote repository
git push -u origin "$new_branch"