#!/bin/bash
#SBATCH --job-name=SF_Np12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=amd_256
#SBATCH --time=480:00:00
#SBATCH --output=output_sf/sf_%j_phi${1}-${2}.out
#SBATCH --error=output_sf/sf_%j_phi${1}-${2}.err

export JULIA_NUM_THREADS=64
export OMP_NUM_THREADS=1
export PHI_START=$1
export PHI_END=$2

JULIA=~/sc71394/shajy/software/env_soft/julia_1.12.4/julia-1.12.4/bin/julia

echo "========================================"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURMD_NODENAME"
echo "PHI:       $PHI_START – $PHI_END"
echo "Start:     $(date)"
$JULIA --version
echo "========================================"

mkdir -p output_sf
cd $SLURM_SUBMIT_DIR

$JULIA --check-bounds=no --threads=64 spectral_flow.jl

EXIT_CODE=$?
echo "End: $(date)   Exit: $EXIT_CODE"
exit $EXIT_CODE
