# DMRG module

This directory contains ITensor/ITensorMPS DMRG workflows for the partial Hall crystal project.

## Main entry points

- `DMRGFullModel.jl`: shared cylinder lattice, Hamiltonian, DMRG, observables, and convergence logic.
- `run_dmrg.jl`: single flux or ground-state DMRG run.
- `run_flux_pump_warm.jl`: sequential warm-start flux insertion trajectory.
- `run_flux_pump_segment.jl`: checkpointed flux insertion segments for PBS jobs.
- `run_flux_point.jl` and `merge_flux_pump.jl`: independent flux-point workflow and post-merge utilities.
- `plot_charge_pump.py`: dependency-free SVG/CSV plotting for charge pump data.
- `analyze_ground_state.jl`: one-file finite-cylinder ground-state diagnostics
  from `summary.dat`, `density.dat`, `connected_density.dat`, `green.dat`,
  `convergence.dat`, and optionally a saved MPS checkpoint.

Example ground-state diagnostic pass:

```bash
julia --project=. dmrg/analyze_ground_state.jl \
  --input=dmrg/output_prl2014_tpl_Lx15_Np36_V10_2_2_unbiased_20260819/stage3_continue_chi3000 \
  --outdir=dmrg/output_prl2014_tpl_Lx15_Np36_V10_2_2_unbiased_20260819/stage3_continue_chi3000/diagnostics \
  --top_peaks=12
```

Add `--checkpoint=auto --entanglement_bonds=center` for a central-cut
entanglement entropy/spectrum check, or `--entanglement_bonds=all` when a full
profile is worth the extra checkpoint SVD cost.

## Cluster jobs

PBS submission templates are under `dmrg/jobs/`. The production Lx=15 charge pump workflow used:

- `jobs/submit_flux_lx15_segmented.sh`
- `jobs/run_flux_lx15_segment.pbs`

The scripts are written for the W003 PBS environment and assume 24 CPU cores per node.

## Data policy

Large or generated DMRG outputs are intentionally ignored by git, including checkpoints, logs, density profiles, cumulative charge files, pumping data, and plotted artifacts. Keep production results on the cluster or archive them separately.

## Infinite cylinder

The [`idmrg/`](idmrg/) directory is a separate-project infinite-cylinder
workflow. Its numerical algorithm is VUMPS on an infinite MPS, not traditional
growth iDMRG. It reuses the finite real-space Hamiltonian and transverse twist
convention but has its own pinned experimental ITensorInfiniteMPS environment,
commands, convergence criteria, observables, and HDF5 checkpoint contract.

See [`idmrg/README.md`](idmrg/README.md) for installation, exact CLI examples,
finite/infinite parameter mapping, raw flux-branch policy, smoke status, and
known backend limitations. Always invoke it with `--project=dmrg/idmrg`; do not
add its unregistered backend to the finite root project.

## Quick checks

From the repository root:

```bash
julia --project=. dmrg/test_convergence_logic.jl
julia --project=. dmrg/test_flux_merge.jl
julia --project=. dmrg/test_flux_segment_checkpoint.jl
julia --project=. dmrg/test_ground_state_diagnostics.jl
python3 -m py_compile dmrg/plot_charge_pump.py
```
