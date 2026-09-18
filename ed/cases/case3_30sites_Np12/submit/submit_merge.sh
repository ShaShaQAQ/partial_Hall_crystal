#!/bin/bash
#SBATCH --job-name=ED_Np12_merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=amd_256
#SBATCH --time=04:00:00
#SBATCH --output=output/merge_%j.out
#SBATCH --error=output/merge_%j.err

export JULIA_NUM_THREADS=16
export OMP_NUM_THREADS=1

JULIA=~/sc71394/shajy/software/env_soft/julia_1.12.4/julia-1.12.4/bin/julia

echo "========================================"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURMD_NODENAME"
echo "Start:     $(date)"
$JULIA --version
echo "========================================"

cd $SLURM_SUBMIT_DIR

$JULIA --check-bounds=no --threads=16 merge.jl

EXIT_CODE=$?
echo "End: $(date)   Exit: $EXIT_CODE"
exit $EXIT_CODE
