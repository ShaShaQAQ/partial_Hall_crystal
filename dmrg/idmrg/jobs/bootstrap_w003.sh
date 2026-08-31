#!/usr/bin/env bash
set -euo pipefail

runtime_root=/home/public/shajy/codex/runtime
julia_bin="$runtime_root/julia-1.12.5/bin/julia"
depot_root=/home/public/shajy/codex/depots/idmrg-julia-1.12.5
gitconfig_path="$depot_root/gitconfig"
julia_archive=julia-1.12.5-linux-x86_64.tar.gz
julia_url="https://mirrors.tuna.tsinghua.edu.cn/julia-releases/bin/linux/x64/1.12/$julia_archive"
julia_sha256=41b84d727e4e96fbf3ed9e92fa195d773d247b9097f73fad688f8b699758bae7

if [[ ${1:-} == --check ]]; then
  [[ -x "$julia_bin" ]]
  [[ -f "$gitconfig_path" ]]
  [[ "$($julia_bin --startup-file=no -e 'print(VERSION)')" == 1.12.5 ]]
  GIT_CONFIG_GLOBAL="$gitconfig_path" JULIA_PKG_USE_CLI_GIT=true \
    JULIA_DEPOT_PATH="$depot_root" "$julia_bin" --startup-file=no \
      --project=dmrg/idmrg -e '
      using Pkg
      Pkg.instantiate()
      using ITensors, ITensorMPS, ITensorInfiniteMPS
      import JLD2, MatrixAlgebraKit, MPSKit, TensorKit, TensorKitTensors
      VERSION == v"1.12.5" || error("wrong Julia version")
    '
  exit 0
fi

mkdir -p "$runtime_root" "$depot_root"
git config --file "$gitconfig_path" url.git@github.com:.insteadOf https://github.com/
if [[ ! -x "$julia_bin" ]]; then
  install_root="$runtime_root/julia-1.12.5"
  [[ ! -e "$install_root" ]] || {
    echo "refusing to replace incomplete runtime: $install_root" >&2
    exit 1
  }
  temporary_root=$(mktemp -d /home/public/shajy/codex/julia-bootstrap.XXXXXX)
  trap 'rm -rf "$temporary_root"' EXIT
  curl -fL "$julia_url" -o "$temporary_root/$julia_archive"
  printf '%s  %s\n' "$julia_sha256" "$temporary_root/$julia_archive" |
    sha256sum --check --status
  tar -xzf "$temporary_root/$julia_archive" -C "$temporary_root"
  mv "$temporary_root/julia-1.12.5" "$install_root"
fi

[[ "$($julia_bin --startup-file=no -e 'print(VERSION)')" == 1.12.5 ]]
GIT_CONFIG_GLOBAL="$gitconfig_path" JULIA_PKG_USE_CLI_GIT=true \
  JULIA_DEPOT_PATH="$depot_root" "$julia_bin" --startup-file=no \
    --project=dmrg/idmrg -e '
    using Pkg
    Pkg.instantiate()
    Pkg.precompile()
  '

exec "$0" --check
