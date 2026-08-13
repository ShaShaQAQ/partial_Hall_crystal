# Infinite-cylinder VUMPS workflow

This directory contains the infinite-cylinder companion to the finite DMRG
workflow in [`dmrg/`](../). Despite the historical directory name `idmrg`, the
algorithm is **VUMPS on an infinite MPS**, not traditional growth-based iDMRG.
Programs and output metadata therefore use the name `VUMPS`.

The implementation uses the official but experimental and unregistered
[`ITensorInfiniteMPS.jl`](https://github.com/ITensor/ITensorInfiniteMPS.jl)
backend. It is isolated from the finite-DMRG environment and pinned in
`Manifest.toml` to:

- ITensorInfiniteMPS 0.2.2, commit
  `765f2777703bc1138b009adbed1b97bde1973402`;
- ITensors 0.9.30;
- ITensorMPS 0.3.45;
- KrylovKit 0.10.4;
- HDF5 0.17.3.

This stack is work in progress. Checkpoints and the adapter to the backend's
iteration interface should not be assumed compatible with another commit or
package resolution.

## Environment and tests

Run commands below from the repository root. Instantiate only the isolated
environment:

```bash
julia --project=dmrg/idmrg -e 'using Pkg; Pkg.instantiate()'
```

The fast test suite can be run directly or through `Pkg.test()`:

```bash
julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl
julia --project=dmrg/idmrg -e 'using Pkg; Pkg.test()'
```

Neither command changes the root finite-DMRG environment.

## Geometry, Hamiltonian, and filling

The MPS is infinite and untwisted along `x`. It repeats a finite reference cell
of `x_period` rings. Each ring has `Ly` physical sites and is periodic along
`y`. The site order is

```text
i(x,y) = x*Ly + y + 1,    0 <= x < x_period, 0 <= y < Ly.
```

The infinite and finite programs reuse the full two-band real-space
Hamiltonian: the hopping matrices determined by `t1` and `t3`, and density
interactions `V1`, `V2`, and `V3`. A hopping path that winds `w_y` times around
the periodic direction receives the same twist factor as finite DMRG,

```text
exp(i*w_y*phi_y).
```

Thus `phi_y` is the total transverse twist in radians, not `phi_y/(2*pi)`.
The Hamiltonian is `2*pi` periodic.

There are two physical sites in each microscopic two-band unit cell. Therefore

```text
sites_per_cell      = x_period*Ly
unit_cells_per_cell = x_period*Ly/2
Np_cell             = filling_num*unit_cells_per_cell/filling_den.
```

`Ly` must be positive and even, `x_period` must be positive, and `Np_cell` must
be an integer. Invalid cells are rejected before sites or the Hamiltonian are
built. At filling `1/3`:

| Geometry | Physical sites | Microscopic two-site cells | `Np_cell` |
|---|---:|---:|---:|
| `Ly=6, x_period=1` | 6 | 3 | 1 |
| `Ly=6, x_period=3` | 18 | 9 | 3 |
| `Ly=2, x_period=3` | 6 | 3 | 1 |

`Ly=2, x_period=1` is invalid because it would contain one third of a particle.
The `Ly=2` valid cell is useful only as a backend smoke geometry; transverse
wrapping is geometrically degenerate there.

The default product occupation pattern is deterministic. For example, the
default occupied sites are `[1]` for `Ly=6, x_period=1` and `[1,7,13]` for
`Ly=6, x_period=3`. A single point or the initial warm state of a scan can use
an explicit product state with `--occupied_sites=...`; the number of distinct
sites must equal `Np_cell`. A flux scan can add further cold product candidates
with `--cold_patterns`; semicolons separate candidates and commas separate
occupied sites within a candidate.

`--seed` (default `0`) controls deterministic `Random.Xoshiro` streams for
neutral-transfer and mixed-transfer eigensolver initial tensors. Scan streams
are derived separately for every point, candidate, and observable, so adding a
fidelity solve does not silently change a later transfer solve. The seed and
occupied pattern are stored in `summary.toml`. This makes runs reproducible at
the input and RNG-stream level; floating-point eigensolvers may still differ
across Julia/package versions, hardware, or thread configurations, all of
which must be treated as provenance.

## Single-point commands

Every argument uses `--key=value` syntax. This minimal backend check is the
smallest valid `1/3` cell:

```bash
julia --project=dmrg/idmrg dmrg/idmrg/bin/run_vumps.jl \
  --Ly=2 --x_period=3 --filling_num=1 --filling_den=3 \
  --t1=1.0 --t3=0.2 --V1=1.0 --V2=0.0 --V3=0.0 \
  --phi_y=0.0 --maxdim=1,4,8 --cutoff=1e-8 \
  --vumps_tol=1e-5 --maxiter=50 --transfer_tol=1e-8 \
  --output=dmrg/idmrg/output/minimal_Ly2_xp3
```

For production studies at `Ly=6`, compare both admissible translation cells.
The following matched commands are starting templates, not evidence that bond
dimension 512 is sufficient or that either state is the production ground
state:

```bash
julia --threads=24 --project=dmrg/idmrg dmrg/idmrg/bin/run_vumps.jl \
  --Ly=6 --x_period=1 --filling_num=1 --filling_den=3 \
  --t1=1.0 --t3=0.2 --V1=1.0 --V2=0.0 --V3=0.0 \
  --phi_y=0.0 --maxdim=1,4,16,64,256,512 --cutoff=1e-10 \
  --vumps_tol=1e-7 --energy_tol=1e-6 --mismatch_tol=1e-6 \
  --stability=3 --maxiter=100 --transfer_tol=1e-10 \
  --threads=24 --output=dmrg/idmrg/output/Ly6_xp1

julia --threads=24 --project=dmrg/idmrg dmrg/idmrg/bin/run_vumps.jl \
  --Ly=6 --x_period=3 --filling_num=1 --filling_den=3 \
  --t1=1.0 --t3=0.2 --V1=1.0 --V2=0.0 --V3=0.0 \
  --phi_y=0.0 --maxdim=1,4,16,64,256,512 --cutoff=1e-10 \
  --vumps_tol=1e-7 --energy_tol=1e-6 --mismatch_tol=1e-6 \
  --stability=3 --maxiter=100 --transfer_tol=1e-10 \
  --threads=24 --output=dmrg/idmrg/output/Ly6_xp3
```

The `x_period=1` ansatz enforces one-ring translation symmetry and can suppress
crystalline order. The `x_period=3` result has three physical transverse cuts
and can expose period-three density structure. Compare energies using the same
normalization, convergence tolerances, and adequately converged bond
dimensions; do not compare the raw per-cell energies of unequal cells.

Useful optional single-point keys are `--checkpoint`, `--occupied_sites`,
`--seed`, `--transfer_neigs`, `--imaginary_tol`, and
`--allow_nonconverged`. The parser recognizes `--load` so an attempted restart
gets a precise diagnostic, but both single-point and flux-scan optimization
reject it before configuring threads or invoking any state/backend callback;
see the reproduced pinned-backend limitation below. Checkpoints are supported
for audit and observable-only reload validation, not optimization restart. By
default the checkpoint is `<output>/state.h5`, `seed=0`, `transfer_neigs=4`,
`stability=2`, and both energy tolerances are `10*vumps_tol`. The default
behavior exits nonzero if VUMPS or a required observable is invalid. Setting
`--allow_nonconverged=true` permits a zero exit but does not change validity
fields or fabricate observables.

## Raw flux scan

`run_flux_scan.jl` uses `--phi_start`, `--phi_stop`, and `--phi_steps`; it does
not accept the single-point `--phi_y` or `--checkpoint` keys. This example
scans one flux quantum on the larger `Ly=6` cell, keeps two explicit cold
patterns in addition to the warm candidate, and selects the adiabatic branch
after the first point:

```bash
julia --threads=24 --project=dmrg/idmrg dmrg/idmrg/bin/run_flux_scan.jl \
  --Ly=6 --x_period=3 --filling_num=1 --filling_den=3 \
  --t1=1.0 --t3=0.2 --V1=1.0 --V2=0.0 --V3=0.0 \
  --phi_start=0.0 --phi_stop=6.283185307179586 --phi_steps=73 \
  --maxdim=1,4,16,64,256,512 --cutoff=1e-10 \
  --vumps_tol=1e-7 --energy_tol=1e-6 --mismatch_tol=1e-6 \
  --stability=3 --maxiter=100 --transfer_tol=1e-10 --threads=24 \
  --branch_mode=adiabatic --cold_patterns=1,7,13\;2,8,14 \
  --output=dmrg/idmrg/output/flux_Ly6_xp3
```

The backslash before the semicolon protects it from the shell. At every flux
point, `candidate_warm` and all requested `candidate_cold_NN` directories are
retained. `ground` selects the lowest finite energy per site among converged
candidates; `adiabatic` selects the highest valid mixed-transfer fidelity among
converged candidates. A failure of an ancillary observable such as the neutral
transfer spectrum remains visible in that candidate's validity fields but does
not change energy/fidelity eligibility. The first point is necessarily selected
by energy because there is no previous state.

The scan constructs one site-index set and reuses those exact indices for all
warm and cold candidates in memory. That identity is required by mixed
transfer matrices. Each candidate still receives its own checkpoint and text
files. `scan_summary.tsv` records the selected raw trajectory, while
`branch_events.tsv` records raw adjacent-point diagnostics and thresholds.
Every branch diagnostic uses `cut_x=1`, stored explicitly. The summary also
stores the selected energy difference from the lowest converged finite energy
at the same flux and both fidelity forms,

```text
fidelity_raw             = fidelity_cell
fidelity_density_x_raw   = -log(fidelity_cell)/x_period.
```

For an exactly zero valid fidelity, the density is retained as `Inf` and
`fidelity_density_x_divergent=true`; it is not confused with a failed solve.

`scan_metadata.toml` records the flux grid, geometry, model, optimization,
branch mode, fixed cut, thresholds, patterns, seed, and runtime provenance even
for a one-point scan with no adjacent event. `sector_gauge.tsv` compares every
pair of explicit cold candidates at the same flux and fixed cut under the best
integer raw-QN shift. The shift is evidence only and is never applied to either
candidate's stored sector weights or polarization.

## Finite/infinite mapping

The production finite campaign summarized in `dmrg/report/dmrg_summary.tex`
used the following parameters. Infinite entries describe the two `Ly=6`
comparison cells above.

| Quantity | Finite `Lx=15` DMRG | Infinite VUMPS |
|---|---|---|
| `Lx` | 15 | omitted; `x` is infinite |
| `Ly` | 6 | 6 |
| Size | `Ns=90` | `sites_per_cell=6` (`xp1`) or 18 (`xp3`) |
| Two-site cells | `Nuc=45` | 3 (`xp1`) or 9 (`xp3`) per MPS cell |
| Particles at `1/3` | total `Np=15` | `Np_cell=1` (`xp1`) or 3 (`xp3`) |
| Boundaries | open `x`, periodic/twisted `y` | infinite `x`, same periodic/twisted `y` |
| Model | `t1=1, t3=.2, V1=1, V2=V3=0` | identical full real-space model |
| Bond control | sweep `maxdim=200,400,800,1200,1600,2000` | staged subspace-expansion caps, e.g. `1,4,...,512` |
| Truncation | `cutoff=1e-9` | expansion `cutoff`; template `1e-10` |
| Convergence | `dE<=1e-6`, density change `<=2e-5`, truncation error `<=5e-8`, 3 stable sweeps | canonical residual, `dE`, and left/right energy mismatch; template `1e-7,1e-6,1e-6`, 3 stable iterations |
| Energy | finite total energy | four energy densities below |
| Pump observable | cumulative density transfer from an open edge | raw Schmidt-charge polarization at a fixed cut |
| State comparison | finite-state overlap | mixed-transfer fidelity per infinite MPS cell |

If `e_cell` is the Hamiltonian contribution anchored in one repeating MPS
cell, output stores all four normalizations:

```text
energy.per_cell      = e_cell
energy.per_x         = e_cell/x_period
energy.per_unit_cell = e_cell/unit_cells_per_cell
energy.per_site      = e_cell/sites_per_cell.
```

Use `per_x`, `per_unit_cell`, or `per_site` consistently when comparing
`x_period=1` and `x_period=3`. A finite total energy is not directly comparable
to any of these; use a separately established finite bulk energy estimator.

## Optimization and convergence records

At each entry in `--maxdim`, the runner first attempts QN-conserving subspace
expansion and then performs VUMPS iterations. A schedule entry is an expansion
cap, not a guarantee that every QN bond reaches that number. Expansion is
accepted only when at least one actual bond satisfies
`after_dimension[i] > before_dimension[i]`; an unchanged cell is an error.
QN constraints commonly make growth nonuniform. `expansion.tsv` records the
stage, requested cap, complete comma-separated before/after bond-dimension
vectors, actual progress flag, and elapsed time. `convergence.tsv` separately
records `maxlinkdim` after each VUMPS iteration.

`convergence.tsv` distinguishes three independent conditions:

- `precision_error=max(eps_left,eps_right)` is the VUMPS canonical residual and
  must be less than `vumps_tol`;
- `delta_energy` is the change in the average of left and right cell-energy
  estimates;
- `energy_mismatch=abs(energy_left-energy_right)` must be less than
  the corresponding mismatch threshold after normalization.

The energy columns remain raw totals for the complete MPS reference cell.
Convergence divides `delta_energy` and `energy_mismatch` by
`sites_per_cell` before comparing them with `energy_tol` and `mismatch_tol`.
The thresholds are therefore energy-per-site tolerances and remain comparable
between `x_period=1` and `x_period=3`; the output is not silently rescaled.

All three must hold for `stability` consecutive iterations at every stage.
The `eps_left`/`eps_right` values are not Krylov eigensolver residuals. Transfer
matrix eigenpair residuals appear separately in `transfer_spectrum.tsv` and
are checked against `transfer_tol`.

Before returning either a converged or nonconverged result, the implementation
runs the official `ITensorMPS.orthogonalize` canonicalization and verifies cell
size, exact site identities, link dimensions, QN flux, and center equations.
Only a validated canonical state can be written to HDF5.

The final Hamiltonian expectation first sums its complex per-cell estimates,
then requires the total imaginary component to be no larger than
`imaginary_tol`. A materially complex expectation is invalid; it is not made
real by discarding its imaginary part.

## Observables

A valid single-point output includes:

- site-resolved density for every `(x,y)` in the reference cell and the total
  density of every `x` ring;
- charge-resolved entanglement spectra at every physical cut after a complete
  `Ly` ring (one cut for `xp1`, three for `xp3`);
- singular values, probabilities, entanglement energies, raw Schmidt QNs,
  converted physical excess charges, sector weights, entropy, and raw
  Schmidt-charge polarization;
- the neutral transfer spectrum and its explicit Krylov residuals;
- normalized mixed-transfer fidelity per complete MPS cell during scans.

Centered QNs use `charge_scale=sites_per_cell/gcd(sites_per_cell,Np_cell)`.
For a consistently oriented left Schmidt index,
`physical_excess_charge=raw_charge/charge_scale`. The charge is an excess
relative to the cell background; a raw QN must not be labeled as a particle
number.

The correlation length uses the two leading eigenvalues found in the **neutral
QN channel** of a transfer matrix spanning one complete MPS cell:

```text
xi_cell = -1/log(abs(lambda2/lambda1))
xi_x    = x_period*xi_cell.
```

It is not a maximization over all charge sectors. The result is invalid unless
two eigenpairs converge with acceptable residuals. Mixed fidelity likewise
comes from a normalized dominant mixed-transfer eigenvalue and is reported per
MPS cell. For comparisons across different `x_period`, a useful derived
quantity is `-log(fidelity_cell)/x_period`; the current scan tables store the
raw per-cell fidelity and this derived per-ring density with separate validity
fields.

## Raw charge and branch policy

The finite `Lx=15` cumulative charge is raw data. In particular, the refined
finite trajectory jumps from raw `cumulative_mid=0.275054113` at
`Phi_y/(2*pi)=0.805556` to `-0.059271363` at `0.819444`. The infinite workflow
uses raw Schmidt-charge polarization and does not add `1/3` or perform modulo
unwrapping.

An integer raw-QN shift that best aligns two cold candidates is a comparison
diagnostic only. It is written as sector evidence and is never applied to raw
Schmidt QNs, sector weights, polarization, or pump data. A branch-switch
interpretation requires combined mixed-fidelity, Schmidt-sector,
charge-resolved entanglement-spectrum, and retained candidate-energy evidence;
no single jump flag relabels the data.

## Output and checkpoint contract

A single-point directory contains:

```text
summary.toml
convergence.tsv
expansion.tsv
density.tsv
ring_density.tsv
entanglement_spectrum.tsv
schmidt_sectors.tsv
transfer_spectrum.tsv
state.h5
```

`summary.toml` is the scalar contract. It records `algorithm="VUMPS"`, the
complete geometry and filling, `phi_y`, charge scale and configuration
signature, model and optimization settings, backend commit, convergence,
validity, energy normalizations, entropy/polarization summaries, and neutral
correlation length. Dependency provenance contains the Julia,
ITensorInfiniteMPS, ITensorMPS, ITensors, KrylovKit, and HDF5 versions plus the
pinned backend commit. The TSV files are the analysis contract and retain raw
numbers with explicit `valid`/`converged` fields. A nonconverged or invalid run
still writes its raw text and canonical checkpoint before the default nonzero
exit.

A flux-scan root contains `scan_metadata.toml`, `scan_summary.tsv`,
`branch_events.tsv`, `sector_gauge.tsv`, and one directory per retained
candidate at every flux. All tables preserve raw quantities and attach
independent validity fields; no derived sector alignment rewrites them.

HDF5 format `infinite_cylinder_vumps_v1` stores the canonical infinite state
and these compatibility attributes exactly:

```text
format, backend_commit,
Ly, x_period, filling_num, filling_den, phi_y,
sites_per_cell, particles_per_cell, charge_scale,
configuration_signature.
```

Observable-only load validates the pinned backend, derived arithmetic,
configuration signature, canonical tensors, and exact equality of `Ly`,
`x_period`, filling, and `phi_y`. Model couplings and optimization settings are
in `summary.toml`, not HDF5 compatibility metadata, so an audit must inspect
both files. Because `phi_y` is strict, a checkpoint at one flux cannot even be
validated as the state of a different flux. Flux-scan warm states instead stay
in memory and reuse the exact same site indices while the Hamiltonian twist
changes.

## Smoke status and known limitations

The opt-in smoke tests were run with Julia 1.12.5 and the pinned backend. The
ignored summaries currently report:

| Smoke | Status | Final energy per site | Neutral `xi_cell` | Charge scale |
|---|---|---:|---:|---:|
| `Ly=2, x_period=3` | converged and valid | -0.3953571006 | 0.1360013040 | 6 |
| `Ly=6, x_period=1` | converged and valid | -0.2539604335 | 0.2779948039 | 6 |

Both used `maxdim=1,4,8`, demonstrated actual QN bond growth, produced valid
charge-resolved and neutral-transfer data, and passed HDF5 reload comparisons.
They are data-contract tests, not production phase results. Run them locally
with:

```bash
IDMRG_RUN_SMOKE=true julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl
```

The 24-core W003 PBS wrapper runs the minimal smoke first and the `Ly=6` smoke
only after it succeeds:

```bash
qsub dmrg/idmrg/jobs/run_smoke.pbs
```

Known pinned-backend limitations and guards are:

- QN subspace expansion may grow only a subset of bonds. The runner checks
  actual per-bond progress and errors if no bond grows; a requested `maxdim`
  remains a cap rather than an attained dimension.
- Final canonicalization uses the official ITensorMPS routine and is validated
  before checkpointing.
- With the pinned backend, loading the fresh canonical `Ly=6` smoke checkpoint
  and requesting further QN subspace expansion reproduces an upstream failure:
  direct-sum link indices have opposite directions (`Out` versus `In`). Even a
  no-expansion VUMPS restart reaches final re-canonicalization and reproduces a
  related direction mismatch. `load_checkpoint` remains validated for
  observable/reload checks, but both CLIs reject `--load` optimization before
  any backend callback at this pinned commit. Start a fresh deterministic or
  explicit cold product candidate instead, and retain all candidates for
  comparison.
- VUMPS can converge to metastable branches, `x_period=1` can suppress broken
  translation symmetry, the reported correlation length covers only the
  neutral channel, and checkpoint representation compatibility is tied to the
  pinned commit.
