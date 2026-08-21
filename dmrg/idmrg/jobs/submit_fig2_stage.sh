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
QSUB_BIN=${QSUB_BIN:-qsub}
DRY_RUN=${DRY_RUN:-0}
job_directory=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
job_contract="$job_directory/fig2_job_contract.sh"
test -f "$job_contract"
# shellcheck disable=SC1090
source "$job_contract"

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

max_dimension=$(fig2_max_dimension_csv "$FIG2_DIMENSIONS")

if [[ -z ${FIG2_WALLTIME:-} ]]; then
  FIG2_WALLTIME=$(fig2_max_walltime "$max_dimension")
fi
fig2_validate_walltime "$max_dimension" "$FIG2_WALLTIME"

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
  printf 'FIG2_JOB_MANIFEST=<copied-on-submit>\n'
  printf 'FIG2_JOB_LAUNCHER=<generated-on-submit>\n'
  if [[ -n "$FIG2_DEPENDENCY" ]]; then
    printf 'qsub -N %s -l walltime=%s -W depend=afterok:%s <generated-on-submit>\n' \
      "$job_name" "$FIG2_WALLTIME" "$FIG2_DEPENDENCY"
  else
    printf 'qsub -N %s -l walltime=%s <generated-on-submit>\n' \
      "$job_name" "$FIG2_WALLTIME"
  fi
  exit 0
fi

config_directory="$results_base/job_configs"
mkdir -p "$config_directory"
umask 077
config=
launcher=
job_manifest=
cleanup_unsubmitted_files() {
  for submission_file in "$launcher" "$config" "$job_manifest"; do
    case "$submission_file" in
      "$config_directory"/"$FIG2_STAGE".*.env|\
      "$config_directory"/"$FIG2_STAGE".*.pbs|\
      "$config_directory"/"$FIG2_STAGE".*.manifest.toml)
        rm -f -- "$submission_file"
        ;;
    esac
  done
}
trap cleanup_unsubmitted_files EXIT

job_manifest=$(mktemp "$config_directory/${FIG2_STAGE}.XXXXXX.manifest.toml")
cp -- "$FIG2_MANIFEST" "$job_manifest"
manifest_sha256=$(sha256sum "$job_manifest" | awk '{print $1}')

config=$(mktemp "$config_directory/${FIG2_STAGE}.XXXXXX.env")
{
  printf 'W003_REPO=%q\n' "$W003_REPO"
  printf 'JULIA_BIN=%q\n' "$JULIA_BIN"
  printf 'JULIA_DEPOT_PATH=%q\n' "$JULIA_DEPOT_PATH"
  printf 'FIG2_MANIFEST=%q\n' "$job_manifest"
  printf 'FIG2_MANIFEST_SHA256=%q\n' "$manifest_sha256"
  printf 'FIG2_STAGE=%q\n' "$FIG2_STAGE"
  printf 'FIG2_OUTPUT=%q\n' "$output"
  printf 'FIG2_DIMENSIONS=%q\n' "$FIG2_DIMENSIONS"
  printf 'FIG2_FLUX_UNITS=%q\n' "$FIG2_FLUX_UNITS"
  printf 'FIG2_THREADS=%q\n' "$FIG2_THREADS"
  printf 'FIG2_WALLTIME=%q\n' "$FIG2_WALLTIME"
} > "$config"

launcher=$(mktemp "$config_directory/${FIG2_STAGE}.XXXXXX.pbs")
sed -n '1,/^$/p' "$runner" > "$launcher"
{
  printf 'export FIG2_JOB_CONFIG=%q\n' "$config"
  printf 'export FIG2_JOB_LAUNCHER=%q\n' "$launcher"
  printf 'exec %q\n' "$runner"
} >> "$launcher"
chmod 0444 "$job_manifest" "$config" "$launcher"

qsub_args=(
  -N "$job_name"
  -l "walltime=$FIG2_WALLTIME"
)
if [[ -n "$FIG2_DEPENDENCY" ]]; then
  qsub_args+=(-W "depend=afterok:$FIG2_DEPENDENCY")
fi

if [[ "$QSUB_BIN" == */* ]]; then
  test -x "$QSUB_BIN"
else
  command -v "$QSUB_BIN" >/dev/null
fi
job_id=$("$QSUB_BIN" "${qsub_args[@]}" "$launcher")
trap - EXIT
printf 'job_id=%s\n' "$job_id"
printf 'job_config=%s\n' "$config"
printf 'job_launcher=%s\n' "$launcher"
printf 'job_manifest=%s\n' "$job_manifest"
printf 'output=%s\n' "$output"
