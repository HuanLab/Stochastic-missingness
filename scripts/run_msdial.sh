#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

# Edit these paths as needed.
# Leave any variable empty ("") to use the default defined in msdial_analysis.py.
# On Windows Git Bash, use forward slashes, e.g. E:/Code/Python/MassCube.

PYTHON_CMD="${PYTHON_CMD:-python}"

PROJECT_DIR="sample_data/msdial_sample"
BLANK_INCLUDED_CSV=""
NOBLANK_CSV=""
INPUT_CSV=""
SINGLE_FILES_DIR=""
SINGLE_FILES_CSV_DIR=""
MZML_DIR="sample_data/mzML/100"
BLANK_MZML=""
OUTPUT_ROOT=""

SKIP_TXT_CONVERSION=false
OVERWRITE_CONVERTED_CSV=false

cmd=("$PYTHON_CMD" "scripts/msdial_analysis.py")

[[ -n "$PROJECT_DIR" ]] && cmd+=(--project-dir "$PROJECT_DIR")
[[ -n "$BLANK_INCLUDED_CSV" ]] && cmd+=(--blank-included-csv "$BLANK_INCLUDED_CSV")
[[ -n "$NOBLANK_CSV" ]] && cmd+=(--noblank-csv "$NOBLANK_CSV")
[[ -n "$INPUT_CSV" ]] && cmd+=(--input-csv "$INPUT_CSV")
[[ -n "$SINGLE_FILES_DIR" ]] && cmd+=(--single-files-dir "$SINGLE_FILES_DIR")
[[ -n "$SINGLE_FILES_CSV_DIR" ]] && cmd+=(--single-files-csv-dir "$SINGLE_FILES_CSV_DIR")
[[ -n "$MZML_DIR" ]] && cmd+=(--mzml-dir "$MZML_DIR")
[[ -n "$BLANK_MZML" ]] && cmd+=(--blank-mzml "$BLANK_MZML")
[[ -n "$OUTPUT_ROOT" ]] && cmd+=(--output-root "$OUTPUT_ROOT")

[[ "$SKIP_TXT_CONVERSION" == true ]] && cmd+=(--skip-txt-conversion)
[[ "$OVERWRITE_CONVERTED_CSV" == true ]] && cmd+=(--overwrite-converted-csv)

echo "Running from: $REPO_ROOT"
echo "Command:"
printf ' %q' "${cmd[@]}"
echo

"${cmd[@]}"