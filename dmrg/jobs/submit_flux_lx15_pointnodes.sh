#!/bin/bash
set -euo pipefail

cd /home/public/shajy/partial_Hall_crystal

STEPS=${STEPS:-12}
LANES=${LANES:-4}
THREADS=${THREADS:-24}
OUTROOT=${OUTROOT:-dmrg/flux_Lx15_chi2000_pointnodes_steps${STEPS}}
DRY_RUN=${DRY_RUN:-0}

if [ "$LANES" -lt 1 ] || [ "$LANES" -gt 4 ]; then
  echo "LANES must be between 1 and 4, got $LANES"
  exit 2
fi

if [ "$THREADS" -gt 24 ]; then
  echo "THREADS=$THREADS exceeds node ppn=24"
  exit 2
fi

mkdir -p "$OUTROOT/logs" dmrg/logs
cat > "$OUTROOT/job_config.env" <<EOF
STEPS=$STEPS
THREADS=$THREADS
OUTROOT=$OUTROOT
EOF

submit_point() {
  local step="$1"
  local dep="${2:-}"
  local name
  name=$(printf "phc_phi_%03d" "$step")
  if [ "$DRY_RUN" = "1" ]; then
    if [ -n "$dep" ]; then
      echo "qsub -N $name -W depend=afterok:$dep dmrg/jobs/submit_flux_lx15_point.pbs"
    else
      echo "qsub -N $name dmrg/jobs/submit_flux_lx15_point.pbs"
    fi
    printf "DRYRUN_%03d\n" "$step"
    return 0
  fi
  if [ -n "$dep" ]; then
    qsub -N "$name" -W "depend=afterok:$dep" dmrg/jobs/submit_flux_lx15_point.pbs
  else
    qsub -N "$name" dmrg/jobs/submit_flux_lx15_point.pbs
  fi
}

submit_merge() {
  local deps="$1"
  if [ "$DRY_RUN" = "1" ]; then
    if [ -n "$deps" ]; then
      echo "qsub -N phc_phi_merge -W depend=afterok:$deps dmrg/jobs/merge_flux_lx15_pointnodes.pbs"
    else
      echo "qsub -N phc_phi_merge dmrg/jobs/merge_flux_lx15_pointnodes.pbs"
    fi
    printf "DRYRUN_MERGE\n"
    return 0
  fi
  if [ -n "$deps" ]; then
    qsub -N phc_phi_merge -W "depend=afterok:$deps" dmrg/jobs/merge_flux_lx15_pointnodes.pbs
  else
    qsub -N phc_phi_merge dmrg/jobs/merge_flux_lx15_pointnodes.pbs
  fi
}

echo "Submitting Lx=15 flux pump as $LANES dependency chains"
echo "OUTROOT=$OUTROOT STEPS=$STEPS THREADS=$THREADS"

all_jobs=()
last_jobs=()

for lane in $(seq 0 $((LANES - 1))); do
  dep=""
  step="$lane"
  while [ "$step" -le "$STEPS" ]; do
    point_dir=$(printf "%s/phi_%03d" "$OUTROOT" "$step")
    if [ -s "$point_dir/summary.dat" ] && [ -s "$point_dir/density.dat" ]; then
      echo "step $step already complete; skipping"
    else
      jobid=$(submit_point "$step" "$dep")
      jobid=$(echo "$jobid" | tail -n 1)
      echo "step $step -> $jobid"
      dep="$jobid"
      all_jobs+=("$jobid")
    fi
    step=$((step + LANES))
  done
  if [ -n "$dep" ]; then
    last_jobs+=("$dep")
  fi
done

if [ "${#last_jobs[@]}" -eq 0 ]; then
  echo "All flux points already complete; submitting merge without point dependencies"
  merge_job=$(submit_merge "")
else
  dep_list=$(IFS=:; echo "${last_jobs[*]}")
  merge_job=$(submit_merge "$dep_list")
fi
merge_job=$(echo "$merge_job" | tail -n 1)

echo "Submitted ${#all_jobs[@]} point jobs"
echo "Point jobs: ${all_jobs[*]:-none}"
echo "Merge job: $merge_job"
