#!/usr/bin/env bash
# Runs every verification script for the Symmetry paper.
set -e

PY=$(command -v python3 || command -v python) || {
  echo "Python was not found on PATH. Install Python 3.9 or later and try again." >&2
  exit 1
}

echo "============================================================"
echo " 1/2  verify_all.py"
echo "============================================================"
"$PY" verify_all.py

echo
echo "============================================================"
echo " 2/2  verify_dedup_independent.py"
echo "============================================================"
"$PY" verify_dedup_independent.py

echo
echo "All verification scripts completed successfully."
