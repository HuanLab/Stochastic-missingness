#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

# Edit these paths as needed.
# Leave any variable empty ("") to use the default defined in masscube_analysis.py.
# On Windows Git Bash, use forward slashes, e.g. E:/Code/Python/MassCube.

PYTHON_CMD="${PYTHON_CMD:-python}"

PROJECT_DIR="sample_data/masscube_sample"
INPUT_CSV=""
SINGLE_TXT_DIR=""
SINGLE_CSV_DIR=""
MZML_ROOT="sample_data/mzML"
BLANK_MZML=""
OUTPUT_ROOT=""

SKIP_TXT_CONVERSION=false
OVERWRITE_CONVERTED_CSV=false

cmd=("$PYTHON_CMD" "scripts/masscube_analysis.py")

[[ -n "$PROJECT_DIR" ]] && cmd+=(--project-dir "$PROJECT_DIR")
[[ -n "$INPUT_CSV" ]] && cmd+=(--input-csv "$INPUT_CSV")
[[ -n "$SINGLE_TXT_DIR" ]] && cmd+=(--single-txt-dir "$SINGLE_TXT_DIR")
[[ -n "$SINGLE_CSV_DIR" ]] && cmd+=(--single-csv-dir "$SINGLE_CSV_DIR")
[[ -n "$MZML_ROOT" ]] && cmd+=(--mzml-root "$MZML_ROOT")
[[ -n "$BLANK_MZML" ]] && cmd+=(--blank-mzml "$BLANK_MZML")
[[ -n "$OUTPUT_ROOT" ]] && cmd+=(--output-root "$OUTPUT_ROOT")

[[ "$SKIP_TXT_CONVERSION" == true ]] && cmd+=(--skip-txt-conversion)
[[ "$OVERWRITE_CONVERTED_CSV" == true ]] && cmd+=(--overwrite-converted-csv)

echo "Running from: $REPO_ROOT"
echo "Command:"
printf ' %q' "${cmd[@]}"
echo

"${cmd[@]}"