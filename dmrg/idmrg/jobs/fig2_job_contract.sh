#!/usr/bin/env bash

fig2_max_dimension_csv() {
  local csv=${1:-}
  [[ "$csv" =~ ^[0-9]+(,[0-9]+)*$ ]] || {
    echo "FIG2_DIMENSIONS must contain comma-separated positive integers" >&2
    return 2
  }

  local dimension
  local numeric
  local maximum=0
  local dimensions=()
  IFS=',' read -r -a dimensions <<< "$csv"
  for dimension in "${dimensions[@]}"; do
    numeric=$((10#$dimension))
    ((numeric > 0)) || {
      echo "FIG2_DIMENSIONS must contain comma-separated positive integers" >&2
      return 2
    }
    ((numeric > maximum)) && maximum=$numeric
  done
  printf '%d\n' "$maximum"
}

fig2_max_walltime() {
  local max_dimension=${1:-}
  [[ "$max_dimension" =~ ^[0-9]+$ ]] || {
    echo "Fig. 2 maximum dimension must be a positive integer" >&2
    return 2
  }
  max_dimension=$((10#$max_dimension))
  ((max_dimension > 0)) || {
    echo "Fig. 2 maximum dimension must be a positive integer" >&2
    return 2
  }

  if ((max_dimension <= 128)); then
    printf '12:00:00\n'
  elif ((max_dimension <= 256)); then
    printf '36:00:00\n'
  elif ((max_dimension <= 1000)); then
    printf '72:00:00\n'
  else
    printf '120:00:00\n'
  fi
}

fig2_walltime_seconds() {
  local walltime=${1:-}
  [[ "$walltime" =~ ^([0-9]{2,3}):([0-5][0-9]):([0-5][0-9])$ ]] || {
    echo "FIG2_WALLTIME must use HH:MM:SS" >&2
    return 2
  }
  local hours=${BASH_REMATCH[1]}
  local minutes=${BASH_REMATCH[2]}
  local seconds=${BASH_REMATCH[3]}
  printf '%d\n' "$((10#$hours * 3600 + 10#$minutes * 60 + 10#$seconds))"
}

fig2_validate_walltime() {
  local max_dimension=${1:-}
  local walltime=${2:-}
  local cap
  local requested_seconds
  local cap_seconds
  cap=$(fig2_max_walltime "$max_dimension") || return
  requested_seconds=$(fig2_walltime_seconds "$walltime") || return
  cap_seconds=$(fig2_walltime_seconds "$cap") || return
  ((requested_seconds <= cap_seconds)) || {
    echo "FIG2_WALLTIME=$walltime exceeds $cap cap for max dimension $max_dimension" >&2
    return 2
  }
}
