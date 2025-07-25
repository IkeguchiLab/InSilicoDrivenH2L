#!/bin/bash

total_jobs=1200
batch_size=100
threshold=300

current_job=0

submit_jobs() {
    for ((i=0; i<batch_size; i++)); do
        if [ $current_job -ge $total_jobs ]; then
            echo "All jobs have been submitted."
            exit 0
        fi
        echo "Submitting jobscript_transition_${current_job}"
        sbatch jobscript_transition_${current_job}
        if [ $? -ne 0 ]; then
            echo "Failed to submit jobscript_transition_${current_job}" >&2
            exit 1
        fi
        ((current_job++))
    done
}

while true; do
    running_jobs=$(squeue | wc -l)

    echo "Current running jobs: $running_jobs"

    if [ $running_jobs -lt $threshold ]; then
        echo "Running jobs below threshold. Submitting new jobs..."
        submit_jobs
    else
        echo "Running jobs above threshold. Waiting..."
    fi

    sleep 60
done
