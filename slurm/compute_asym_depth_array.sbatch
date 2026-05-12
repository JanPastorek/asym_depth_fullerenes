#!/bin/bash
#SBATCH --job-name=ipr_asym_depth
#SBATCH --output=asym_depth_%A_%a.out
#SBATCH --error=asym_depth_%A_%a.err
#SBATCH --partition=CPU
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-9

set -euo pipefail

# Optional on TUKE PERUN: activate scratch if available.
if [ -f .activate_scratch ]; then
    source .activate_scratch
fi

N_MIN=${N_MIN:-102}
N_MAX=${N_MAX:-120}
N_STEP=${N_STEP:-2}
MAX_DEPTH=${MAX_DEPTH:-4}
RESULTS_DIR=${RESULTS_DIR:-results/asym_depth_parts}
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

mkdir -p "$RESULTS_DIR"

idx=${SLURM_ARRAY_TASK_ID}
n_c=$((N_MIN + idx * N_STEP))

if [ "$n_c" -gt "$N_MAX" ]; then
    echo "Skipping task $idx because n_c=$n_c > N_MAX=$N_MAX"
    exit 0
fi

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

out_csv=$(printf "%s/ipr_depth_n%03d.csv" "$RESULTS_DIR" "$n_c")

echo "Task $idx processing n_c=$n_c with JULIA_NUM_THREADS=$JULIA_NUM_THREADS"

"$JULIA_BIN" --project=. scripts/calculate_fullerene_asym_depth_hpc.jl \
    "$n_c" "$n_c" "$out_csv" "$MAX_DEPTH"

echo "Done n_c=$n_c -> $out_csv"