#!/bin/bash
set -euo pipefail

# SCRIPT NAME: prune-deleted.sh
# PURPOSE: This script prunes deleted branches from the local git repository and then fetches updates
#          from the remote repository to ensure the local branch list is up to date.
# AUTHOR: Matthew Hirschberger
# DATE: 06 Mar 2026
# REVISION: 1.0
# ----------------------------------------------------------------------

git branch -vv | grep ':gone]' | awk '{print $1}' | xargs git branch -d