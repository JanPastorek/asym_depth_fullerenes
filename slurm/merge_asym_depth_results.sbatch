#!/bin/bash
#SBATCH --job-name=merge_asym_depth
#SBATCH --output=merge_asym_depth_%j.out
#SBATCH --error=merge_asym_depth_%j.err
#SBATCH --partition=CPU
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail

# Optional on TUKE PERUN: activate scratch if available.
if [ -f .activate_scratch ]; then
    source .activate_scratch
fi

INPUT_GLOB=${INPUT_GLOB:-results/asym_depth_parts/*.csv}
OUTPUT_CSV=${OUTPUT_CSV:-results/ipr_fullerenes_depths_to120.csv}
JULIA_BIN=${JULIA_BIN:-julia}
JULIA_MODULE=${JULIA_MODULE:-julia/1.11.2}

if ! command -v "$JULIA_BIN" >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load "$JULIA_MODULE"
    fi
fi

if ! command -v "$JULIA_BIN" >/dev/null 2>&1; then
    echo "ERROR: Julia executable '$JULIA_BIN' not found. Set JULIA_BIN or JULIA_MODULE."
    exit 127
fi

"$JULIA_BIN" --project=. scripts/merge_asym_depth_results.jl "$INPUT_GLOB" "$OUTPUT_CSV"

echo "Merged results saved to $OUTPUT_CSV"