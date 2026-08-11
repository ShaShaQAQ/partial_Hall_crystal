#!/bin/bash
set -euo pipefail

cd /home/public/shajy/partial_Hall_crystal

CONFIG_DIR=dmrg/warm_configs
mkdir -p "$CONFIG_DIR" dmrg/logs

THREADS=${THREADS:-24}
CYCLES=${CYCLES:-3}
MAX_SWEEPS=${MAX_SWEEPS:-40}
RETRY_SWEEPS=${RETRY_SWEEPS:-40}
MAX_RETRIES=${MAX_RETRIES:-2}
MIN_SWEEPS=${MIN_SWEEPS:-4}
STABLE_SWEEPS=${STABLE_SWEEPS:-3}
ENERGY_TOL=${ENERGY_TOL:-1e-6}
DENSITY_TOL=${DENSITY_TOL:-2e-5}
TRUNCERR_TOL=${TRUNCERR_TOL:-5e-8}
REQUIRE_CONVERGED=${REQUIRE_CONVERGED:-true}
DRY_RUN=${DRY_RUN:-0}

write_config() {
  local name="$1"
  local outroot="$2"
  local direction="$3"
  local steps_per_2pi="$4"
  local seed="$5"
  mkdir -p "$outroot/logs"
  cat > "$CONFIG_DIR/$name.env" <<EOF
OUTROOT=$outroot
DIRECTION=$direction
STEPS_PER_2PI=$steps_per_2pi
CYCLES=$CYCLES
THREADS=$THREADS
SEED=$seed
MAX_SWEEPS=$MAX_SWEEPS
RETRY_SWEEPS=$RETRY_SWEEPS
MAX_RETRIES=$MAX_RETRIES
MIN_SWEEPS=$MIN_SWEEPS
STABLE_SWEEPS=$STABLE_SWEEPS
ENERGY_TOL=$ENERGY_TOL
DENSITY_TOL=$DENSITY_TOL
TRUNCERR_TOL=$TRUNCERR_TOL
REQUIRE_CONVERGED=$REQUIRE_CONVERGED
EOF
}

submit_job() {
  local name="$1"
  if [ "$DRY_RUN" = "1" ]; then
    echo "qsub -N $name dmrg/jobs/run_flux_lx15_warm.pbs"
    printf "DRYRUN_%s\n" "$name"
  else
    qsub -N "$name" dmrg/jobs/run_flux_lx15_warm.pbs
  fi
}

write_config warm_fwd8a dmrg/flux_Lx15_warm_fwd_pi8_seed1234 forward 16 1234
write_config warm_fwd8b dmrg/flux_Lx15_warm_fwd_pi8_seed4321 forward 16 4321
write_config warm_bwd8a dmrg/flux_Lx15_warm_bwd_pi8_seed2468 backward 16 2468
write_config warm_fwd12 dmrg/flux_Lx15_warm_fwd_pi12_seed1234 forward 24 1234

echo "Submitting 4 warm-start flux trajectories"
echo "THREADS=$THREADS CYCLES=$CYCLES MAX_SWEEPS=$MAX_SWEEPS RETRY_SWEEPS=$RETRY_SWEEPS MAX_RETRIES=$MAX_RETRIES"
jobs=()
for name in warm_fwd8a warm_fwd8b warm_bwd8a warm_fwd12; do
  jobid=$(submit_job "$name" | tail -n 1)
  echo "$name -> $jobid"
  jobs+=("$jobid")
done
echo "Submitted jobs: ${jobs[*]}"
