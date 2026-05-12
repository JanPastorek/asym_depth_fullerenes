#!/bin/bash
set -euo pipefail

N_MIN=${N_MIN:-102}
N_MAX=${N_MAX:-120}
N_STEP=${N_STEP:-2}
PARTITION=${PARTITION:-CPU}

if [ "$N_MIN" -gt "$N_MAX" ]; then
    echo "N_MIN must be <= N_MAX"
    exit 1
fi

if [ "$N_STEP" -le 0 ]; then
    echo "N_STEP must be > 0"
    exit 1
fi

task_count=$(( (N_MAX - N_MIN) / N_STEP + 1 ))
array_spec="0-$((task_count - 1))"

echo "Submitting array with N_MIN=$N_MIN N_MAX=$N_MAX N_STEP=$N_STEP PARTITION=$PARTITION (array $array_spec)"

array_job_id=$(
    sbatch \
    --partition="$PARTITION" \
        --array="$array_spec" \
        --export=ALL,N_MIN="$N_MIN",N_MAX="$N_MAX",N_STEP="$N_STEP" \
        slurm/compute_asym_depth_array.sbatch | awk '{print $4}'
)

echo "Array job id: $array_job_id"

merge_job_id=$(
    sbatch \
    --partition="$PARTITION" \
        --dependency=afterok:"$array_job_id" \
        slurm/merge_asym_depth_results.sbatch | awk '{print $4}'
)

echo "Merge job id: $merge_job_id"