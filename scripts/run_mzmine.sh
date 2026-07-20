#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

# Edit these paths as needed.
# Leave any variable empty ("") to use the default defined in mzmine_analysis.py.
# On Windows Git Bash, use forward slashes, e.g. E:/Nhi/Code/Python/MassCube.

PYTHON_CMD="${PYTHON_CMD:-python}"

PROJECT_DIR="sample_data/mzmine_sample"
INPUT_CSV=""
MZML_DIR="sample_data/mzML/100"
BLANK_MZML=""
OUTPUT_ROOT=""

cmd=("$PYTHON_CMD" "scripts/mzmine_analysis.py")

[[ -n "$PROJECT_DIR" ]] && cmd+=(--project-dir "$PROJECT_DIR")
[[ -n "$INPUT_CSV" ]] && cmd+=(--input-csv "$INPUT_CSV")
[[ -n "$MZML_DIR" ]] && cmd+=(--mzml-dir "$MZML_DIR")
[[ -n "$BLANK_MZML" ]] && cmd+=(--blank-mzml "$BLANK_MZML")
[[ -n "$OUTPUT_ROOT" ]] && cmd+=(--output-root "$OUTPUT_ROOT")

echo "Running from: $REPO_ROOT"
echo "Command:"
printf ' %q' "${cmd[@]}"
echo

"${cmd[@]}"