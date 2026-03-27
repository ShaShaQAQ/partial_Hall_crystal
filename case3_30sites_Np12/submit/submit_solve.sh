#!/bin/bash
#SBATCH --job-name=ED_Np12_solve
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=amd_256
#SBATCH --time=240:00:00
#SBATCH --output=output/solve_%j_k${1}-${2}.out
#SBATCH --error=output/solve_%j_k${1}-${2}.err

export JULIA_NUM_THREADS=64
export OMP_NUM_THREADS=1
export SECTOR_START=$1
export SECTOR_END=$2

JULIA=~/sc71394/shajy/software/env_soft/julia_1.12.4/julia-1.12.4/bin/julia

echo "========================================"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURMD_NODENAME"
echo "Sectors:   $SECTOR_START – $SECTOR_END"
echo "Start:     $(date)"
$JULIA --version
echo "========================================"

mkdir -p output
cd $SLURM_SUBMIT_DIR

$JULIA --check-bounds=no --threads=64 main.jl

EXIT_CODE=$?
echo "End: $(date)   Exit: $EXIT_CODE"
exit $EXIT_CODE
