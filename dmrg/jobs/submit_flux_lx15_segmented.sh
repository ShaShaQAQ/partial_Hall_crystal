#!/bin/bash
set -euo pipefail

cd /home/public/shajy/partial_Hall_crystal

CONFIG_DIR=dmrg/segment_configs
mkdir -p "$CONFIG_DIR" dmrg/logs

RUN_TAG=${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}
JOB_TAG=${JOB_TAG:-$(date +%H%M%S)}
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-dmrg/flux_Lx15_pump_prod_${RUN_TAG}}
MANIFEST="$CAMPAIGN_ROOT/segment_jobs.tsv"

LX=${LX:-15}
LY=${LY:-6}
NP=${NP:-15}
THREADS=${THREADS:-24}
CYCLES=${CYCLES:-3}
MAIN_STEPS_PER_2PI=${MAIN_STEPS_PER_2PI:-24}
MAIN_SEGMENT_LEN=${MAIN_SEGMENT_LEN:-3}
BRANCH_SEGMENT_LEN=${BRANCH_SEGMENT_LEN:-3}

MAX_SWEEPS=${MAX_SWEEPS:-40}
RETRY_SWEEPS=${RETRY_SWEEPS:-40}
MAX_RETRIES=${MAX_RETRIES:-2}
MIN_SWEEPS=${MIN_SWEEPS:-4}
STABLE_SWEEPS=${STABLE_SWEEPS:-3}
ENERGY_TOL=${ENERGY_TOL:-1e-6}
DENSITY_TOL=${DENSITY_TOL:-2e-5}
TRUNCERR_TOL=${TRUNCERR_TOL:-5e-8}
REQUIRE_CONVERGED=${REQUIRE_CONVERGED:-true}
TARGET_MAXDIM=${TARGET_MAXDIM:-2000}
MAXDIM=${MAXDIM:-200,400,800,1200,1600,2000}
CUTOFF=${CUTOFF:-1e-9}

SEED_FORWARD=${SEED_FORWARD:-1234}
SEED_BACKWARD=${SEED_BACKWARD:-2468}
SEED_BRANCH48=${SEED_BRANCH48:-1357}
SEED_BRANCH72=${SEED_BRANCH72:-9753}

SUBMIT_FORWARD=${SUBMIT_FORWARD:-1}
SUBMIT_BACKWARD=${SUBMIT_BACKWARD:-1}
SUBMIT_BRANCH=${SUBMIT_BRANCH:-1}
DRY_RUN=${DRY_RUN:-0}

MAIN_OUT=${MAIN_OUT:-$CAMPAIGN_ROOT/forward_pi12}
BACKWARD_OUT=${BACKWARD_OUT:-$CAMPAIGN_ROOT/backward_pi12}
BRANCH48_OUT=${BRANCH48_OUT:-$CAMPAIGN_ROOT/branch_from_fwd17_pi24}
BRANCH72_OUT=${BRANCH72_OUT:-$CAMPAIGN_ROOT/branch_from_fwd17_pi36}

MAIN_TOTAL_STEPS=$((CYCLES * MAIN_STEPS_PER_2PI))
BRANCH_BASE_STEP=${BRANCH_BASE_STEP:-17}
BRANCH_START_STEP=${BRANCH_START_STEP:-18}
BRANCH48_END_STEP=${BRANCH48_END_STEP:-24}
BRANCH72_END_STEP=${BRANCH72_END_STEP:-26}

mkdir -p "$CAMPAIGN_ROOT"
printf "# label\tjobname\tjobid\toutroot\tdirection\tsteps_per_2pi\tstart_step\tend_step\tdependency\tcheckpoint_in\tphi_start\tdphi\n" > "$MANIFEST"

float_eval() {
  awk "BEGIN{printf \"%.17g\", $*}"
}

write_config() {
  local name="$1"
  local outroot="$2"
  local direction="$3"
  local steps_per_2pi="$4"
  local start_step="$5"
  local end_step="$6"
  local seed="$7"
  local checkpoint_in="$8"
  local phi_start="$9"
  local dphi="${10}"

  mkdir -p "$outroot/logs" "$outroot/checkpoints"
  cat > "$CONFIG_DIR/$name.env" <<EOF
LX=$LX
LY=$LY
NP=$NP
OUTROOT=$outroot
DIRECTION=$direction
STEPS_PER_2PI=$steps_per_2pi
CYCLES=$CYCLES
START_STEP=$start_step
END_STEP=$end_step
CHECKPOINT_IN=$checkpoint_in
PHI_START=$phi_start
DPHI=$dphi
THREADS=$THREADS
SEED=$seed
T1=1.0
T3=0.2
V1=1.0
V2=0.0
V3=0.0
MAX_SWEEPS=$MAX_SWEEPS
RETRY_SWEEPS=$RETRY_SWEEPS
MAX_RETRIES=$MAX_RETRIES
MIN_SWEEPS=$MIN_SWEEPS
STABLE_SWEEPS=$STABLE_SWEEPS
ENERGY_TOL=$ENERGY_TOL
DENSITY_TOL=$DENSITY_TOL
TRUNCERR_TOL=$TRUNCERR_TOL
REQUIRE_CONVERGED=$REQUIRE_CONVERGED
TARGET_MAXDIM=$TARGET_MAXDIM
MAXDIM=$MAXDIM
CUTOFF=$CUTOFF
EOF
}

submit_job() {
  local name="$1"
  local dependency="$2"
  if [ "$DRY_RUN" = "1" ]; then
    if [ -n "$dependency" ]; then
      echo "qsub -W depend=afterok:$dependency -N $name dmrg/jobs/run_flux_lx15_segment.pbs" >&2
    else
      echo "qsub -N $name dmrg/jobs/run_flux_lx15_segment.pbs" >&2
    fi
    printf "DRYRUN_%s\n" "$name"
    return
  fi
  if [ -n "$dependency" ]; then
    qsub -W depend=afterok:"$dependency" -N "$name" dmrg/jobs/run_flux_lx15_segment.pbs
  else
    qsub -N "$name" dmrg/jobs/run_flux_lx15_segment.pbs
  fi
}

LAST_CHAIN_JOB=""
MAIN_STEP17_JOB=""

submit_chain() {
  local label="$1"
  local prefix="$2"
  local outroot="$3"
  local direction="$4"
  local steps_per_2pi="$5"
  local chain_start="$6"
  local chain_end="$7"
  local segment_len="$8"
  local seed="$9"
  local first_dependency="${10}"
  local first_checkpoint="${11}"
  local chain_phi_start="${12}"
  local chain_dphi="${13}"
  local track_step="${14}"

  local dependency="$first_dependency"
  local checkpoint_in="$first_checkpoint"
  [ -n "$checkpoint_in" ] || checkpoint_in=auto

  local start="$chain_start"
  while [ "$start" -le "$chain_end" ]; do
    local end=$((start + segment_len - 1))
    if [ "$end" -gt "$chain_end" ]; then
      end="$chain_end"
    fi

    local name
    printf -v name "%s%s%03d" "$prefix" "$JOB_TAG" "$start"
    local phi_start=""
    local dphi=""
    if [ -n "$chain_phi_start" ] && [ -n "$chain_dphi" ]; then
      phi_start=$(float_eval "$chain_phi_start + $chain_dphi * ($start - $chain_start)")
      dphi="$chain_dphi"
    fi

    write_config "$name" "$outroot" "$direction" "$steps_per_2pi" "$start" "$end" "$seed" "$checkpoint_in" "$phi_start" "$dphi"
    local jobid
    jobid=$(submit_job "$name" "$dependency" | tail -n 1)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$label" "$name" "$jobid" "$outroot" "$direction" "$steps_per_2pi" "$start" "$end" "${dependency:-none}" "$checkpoint_in" "${phi_start:-auto}" "${dphi:-auto}" >> "$MANIFEST"
    echo "$label step $start-$end -> $jobid"

    if [ -n "$track_step" ] && [ "$start" -le "$track_step" ] && [ "$end" -ge "$track_step" ]; then
      MAIN_STEP17_JOB="$jobid"
    fi

    dependency="$jobid"
    checkpoint_in=auto
    start=$((end + 1))
  done
  LAST_CHAIN_JOB="$dependency"
}

echo "Submitting segmented Lx=15 flux pump campaign"
echo "CAMPAIGN_ROOT=$CAMPAIGN_ROOT"
echo "JOB_TAG=$JOB_TAG THREADS=$THREADS MAIN_STEPS_PER_2PI=$MAIN_STEPS_PER_2PI MAIN_SEGMENT_LEN=$MAIN_SEGMENT_LEN"
echo "Convergence: max_sweeps=$MAX_SWEEPS retry_sweeps=$RETRY_SWEEPS max_retries=$MAX_RETRIES energy_tol=$ENERGY_TOL density_tol=$DENSITY_TOL truncerr_tol=$TRUNCERR_TOL"

if [ "$SUBMIT_FORWARD" = "1" ]; then
  submit_chain "forward_main" "f" "$MAIN_OUT" "forward" "$MAIN_STEPS_PER_2PI" 0 "$MAIN_TOTAL_STEPS" "$MAIN_SEGMENT_LEN" "$SEED_FORWARD" "" "auto" "" "" "$BRANCH_BASE_STEP"
fi

if [ "$SUBMIT_BACKWARD" = "1" ]; then
  submit_chain "backward_check" "b" "$BACKWARD_OUT" "backward" "$MAIN_STEPS_PER_2PI" 0 "$MAIN_TOTAL_STEPS" "$MAIN_SEGMENT_LEN" "$SEED_BACKWARD" "" "auto" "" "" ""
fi

if [ "$SUBMIT_BRANCH" = "1" ]; then
  if [ -z "$MAIN_STEP17_JOB" ]; then
    echo "Cannot submit branch-window chains because the forward chain did not track step $BRANCH_BASE_STEP"
    exit 2
  fi
  base_phi=$(float_eval "2 * atan2(0,-1) * $BRANCH_BASE_STEP / $MAIN_STEPS_PER_2PI")
  dphi48=$(float_eval "2 * atan2(0,-1) / 48")
  phi_start48=$(float_eval "$base_phi + $dphi48")
  dphi72=$(float_eval "2 * atan2(0,-1) / 72")
  phi_start72=$(float_eval "$base_phi + $dphi72")
  branch_checkpoint="$MAIN_OUT/checkpoints/$(printf "state_%03d.jls" "$BRANCH_BASE_STEP")"

  submit_chain "branch_pi24" "p" "$BRANCH48_OUT" "forward" 48 "$BRANCH_START_STEP" "$BRANCH48_END_STEP" "$BRANCH_SEGMENT_LEN" "$SEED_BRANCH48" "$MAIN_STEP17_JOB" "$branch_checkpoint" "$phi_start48" "$dphi48" ""
  submit_chain "branch_pi36" "q" "$BRANCH72_OUT" "forward" 72 "$BRANCH_START_STEP" "$BRANCH72_END_STEP" "$BRANCH_SEGMENT_LEN" "$SEED_BRANCH72" "$MAIN_STEP17_JOB" "$branch_checkpoint" "$phi_start72" "$dphi72" ""
fi

echo "Manifest: $MANIFEST"
