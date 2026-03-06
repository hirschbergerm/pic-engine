#!/bin/bash
set -euo pipefail

# SCRIPT NAME: prune-deleted.sh
# PURPOSE: This script prunes deleted branches from the local git repository and then fetches updates
#          from the remote repository to ensure the local branch list is up to date.
# AUTHOR: Matthew Hirschberger
# DATE: 06 Mar 2026
# REVISION: 1.0
# ----------------------------------------------------------------------

# Fetch and prune the remote-tracking branches to remove any whose remote has been deleted
git fetch --prune

# get current branch
current_branch=$(git symbolic-ref --short HEAD)

# Make list of local branches whose upstream is gone
gone_branches=$(git branch --format '%(refname:short) %(upstream:track)' \
| awk '$2=="[gone]" {print $1}')

# Check if we're on any of these branches and throw a warning if so
if echo "$gone_branches" | grep -qx "$current_branch"; then
    echo "WARNING: You are currently on branch '$current_branch' which has been deleted from the remote."
    echo "Please switch to a different branch before running this script to avoid issues."
    exit 1
fi

# Delete all the other gone branches safely
if [ -n "gone_branches" ]; then # use string check operator -n to check if the list is empty
    echo "$gone_branches" | xargs -r git branch -d
else 
    echo "No local branches with deleted upstreams found. No branches pruned."
fi
