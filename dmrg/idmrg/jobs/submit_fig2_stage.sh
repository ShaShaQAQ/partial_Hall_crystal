#!/usr/bin/env bash
set -euo pipefail

export PATH=/home/shajy/bin:/usr/local/bin:/usr/bin:/bin:/opt/bin:/opt/pbs/bin

results_base=/home/public/shajy/codex/results/fqahc-fig2
results_prefix=/home/public/shajy/codex/results/fqahc-fig2/
W003_REPO=${W003_REPO:-/home/public/shajy/codex/partial_Hall_crystal-idmrg-benchmark}
JULIA_BIN=${JULIA_BIN:-/home/public/shajy/codex/runtime/julia-1.12.5/bin/julia}
JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH:-/home/public/shajy/codex/depots/idmrg-julia-1.12.5}
FIG2_MANIFEST=${FIG2_MANIFEST:-$W003_REPO/dmrg/idmrg/benchmarks/fqahc_fig2.toml}
FIG2_STAGE=${FIG2_STAGE:-pilot}
FIG2_OUTPUT=${FIG2_OUTPUT:-$results_base/$FIG2_STAGE}
FIG2_DIMENSIONS=${FIG2_DIMENSIONS:-32,64,128}
FIG2_FLUX_UNITS=${FIG2_FLUX_UNITS:-0,0.5,1,1.5,2,2.5,3}
FIG2_THREADS=${FIG2_THREADS:-24}
FIG2_DEPENDENCY=${FIG2_DEPENDENCY:-}
DRY_RUN=${DRY_RUN:-0}

[[ "$FIG2_STAGE" =~ ^[A-Za-z0-9._-]+$ ]] || {
  echo "FIG2_STAGE contains an unsafe character" >&2
  exit 2
}
case "$FIG2_THREADS" in
  4|12|24) ;;
  *)
    echo "FIG2_THREADS must be one of 4, 12, or 24" >&2
    exit 2
    ;;
esac
[[ -z "$FIG2_DEPENDENCY" || "$FIG2_DEPENDENCY" =~ ^[0-9]+([.]w003)?$ ]] || {
  echo "FIG2_DEPENDENCY is not a PBS job ID" >&2
  exit 2
}

output=$(/usr/bin/realpath -m "$FIG2_OUTPUT")
case "$output" in
  "$results_prefix"*) ;;
  *)
    echo "FIG2_OUTPUT must be below $results_prefix" >&2
    exit 2
    ;;
esac
[[ "$output" != "$results_base" ]] || {
  echo "FIG2_OUTPUT must not be the shared results root" >&2
  exit 2
}

IFS=',' read -r -a dimensions <<< "$FIG2_DIMENSIONS"
((${#dimensions[@]} > 0)) || {
  echo "FIG2_DIMENSIONS must not be empty" >&2
  exit 2
}
max_dimension=0
for dimension in "${dimensions[@]}"; do
  [[ "$dimension" =~ ^[0-9]+$ ]] && ((dimension > 0)) || {
    echo "FIG2_DIMENSIONS must contain positive integers" >&2
    exit 2
  }
  ((dimension > max_dimension)) && max_dimension=$dimension
done

if [[ -z ${FIG2_WALLTIME:-} ]]; then
  if ((max_dimension <= 128)); then
    FIG2_WALLTIME=12:00:00
  elif ((max_dimension <= 256)); then
    FIG2_WALLTIME=36:00:00
  elif ((max_dimension <= 1000)); then
    FIG2_WALLTIME=72:00:00
  else
    FIG2_WALLTIME=120:00:00
  fi
fi
[[ "$FIG2_WALLTIME" =~ ^[0-9]{2,3}:[0-5][0-9]:[0-5][0-9]$ ]] || {
  echo "FIG2_WALLTIME must use HH:MM:SS" >&2
  exit 2
}

runner="$W003_REPO/dmrg/idmrg/jobs/run_fig2_stage.pbs"
test -f "$runner"
test -f "$FIG2_MANIFEST"

job_name="f2_${FIG2_STAGE}_d${max_dimension}"
job_name=${job_name//[^A-Za-z0-9_]/_}
job_name=${job_name:0:15}

if [[ "$DRY_RUN" == 1 ]]; then
  printf 'FIG2_STAGE=%s\n' "$FIG2_STAGE"
  printf 'FIG2_OUTPUT=%s\n' "$output"
  printf 'FIG2_DIMENSIONS=%s\n' "$FIG2_DIMENSIONS"
  printf 'FIG2_FLUX_UNITS=%s\n' "$FIG2_FLUX_UNITS"
  printf 'FIG2_THREADS=%s\n' "$FIG2_THREADS"
  printf 'FIG2_WALLTIME=%s\n' "$FIG2_WALLTIME"
  if [[ -n "$FIG2_DEPENDENCY" ]]; then
    printf 'qsub -N %s -l walltime=%s -W depend=afterok:%s -v FIG2_JOB_CONFIG=<generated-on-submit> %s\n' \
      "$job_name" "$FIG2_WALLTIME" "$FIG2_DEPENDENCY" "$runner"
  else
    printf 'qsub -N %s -l walltime=%s -v FIG2_JOB_CONFIG=<generated-on-submit> %s\n' \
      "$job_name" "$FIG2_WALLTIME" "$runner"
  fi
  exit 0
fi

config_directory="$results_base/job_configs"
mkdir -p "$config_directory"
umask 077
config=$(mktemp "$config_directory/${FIG2_STAGE}.XXXXXX.env")
{
  printf 'W003_REPO=%q\n' "$W003_REPO"
  printf 'JULIA_BIN=%q\n' "$JULIA_BIN"
  printf 'JULIA_DEPOT_PATH=%q\n' "$JULIA_DEPOT_PATH"
  printf 'FIG2_MANIFEST=%q\n' "$FIG2_MANIFEST"
  printf 'FIG2_STAGE=%q\n' "$FIG2_STAGE"
  printf 'FIG2_OUTPUT=%q\n' "$output"
  printf 'FIG2_DIMENSIONS=%q\n' "$FIG2_DIMENSIONS"
  printf 'FIG2_FLUX_UNITS=%q\n' "$FIG2_FLUX_UNITS"
  printf 'FIG2_THREADS=%q\n' "$FIG2_THREADS"
  printf 'FIG2_WALLTIME=%q\n' "$FIG2_WALLTIME"
} > "$config"
chmod 0444 "$config"

qsub_args=(
  -N "$job_name"
  -l "walltime=$FIG2_WALLTIME"
  -v "FIG2_JOB_CONFIG=$config"
)
if [[ -n "$FIG2_DEPENDENCY" ]]; then
  qsub_args+=(-W "depend=afterok:$FIG2_DEPENDENCY")
fi

job_id=$(qsub "${qsub_args[@]}" "$runner")
printf 'job_id=%s\n' "$job_id"
printf 'job_config=%s\n' "$config"
printf 'output=%s\n' "$output"
