# Infinite-Cylinder VUMPS Design

**Date:** 2026-08-11

**Status:** Approved design

## Purpose

Add an infinite-cylinder tensor-network workflow to the DMRG branch without
replacing or destabilizing the validated finite-cylinder DMRG workflow. The new
workflow represents an infinite \(x\) direction with a finite repeating MPS unit
cell and keeps the existing periodic, flux-twisted \(y\) direction.

The numerical backend is VUMPS acting on an infinite MPS. The feature directory
is named `dmrg/idmrg` to match the project goal, but programs and output metadata
must identify the algorithm as VUMPS. This design does not claim to implement a
traditional growth-based iDMRG algorithm.

The first implementation must:

- reuse the full two-band real-space Hamiltonian used by finite DMRG;
- preserve the finite-cylinder raw-data policy;
- compute energy densities, a neutral-channel correlation length,
  charge-resolved entanglement spectra, and Schmidt-charge diagnostics;
- track flux branches using infinite-MPS fidelity and entanglement data;
- run a small-width backend smoke test and a `Ly=6` target smoke test.

## Existing Baseline

The finite implementation currently lives in `dmrg/DMRGFullModel.jl`. It
combines cylinder geometry, real-space Hamiltonian generation, finite MPO
construction, finite DMRG control, convergence observers, observables, and text
output utilities. Flux workflows use warm-started finite MPS checkpoints and
measure cumulative density transfer from an open edge.

At design time the root environment contains:

- Julia 1.12.5;
- ITensors 0.9.30;
- ITensorMPS 0.4.1;
- KrylovKit 0.10.3.

The three existing DMRG test programs pass with 19 total assertions. These
tests form a non-regression gate for the shared-model extraction.

The accepted finite `Lx=15` data policy is documented in
`dmrg/report/dmrg_summary.tex`: the raw charge jump between
`Phi/(2pi)=0.805556` and `0.819444` is retained without adding `1/3`. A sector
interpretation requires overlap, Schmidt-sector, or entanglement-spectrum
evidence.

## Backend Decision and Dependency Isolation

Standard ITensors and ITensorMPS do not expose iDMRG or VUMPS. The official
`ITensor/ITensorInfiniteMPS.jl` repository provides infinite MPS/MPO types,
VUMPS/TDVP, transfer matrices, QN unit cells, and subspace expansion. The
package is unregistered and explicitly marked as work in progress.

The audited repository state is commit
`765f2777703bc1138b009adbed1b97bde1973402`, package version 0.2.2. Its
compatibility bounds accept ITensors 0.9 but require ITensorMPS 0.3. An isolated
resolution on Julia 1.12.5 loaded successfully with ITensors 0.9.30 and
ITensorMPS 0.3.45.

The infinite workflow therefore has its own `dmrg/idmrg/Project.toml` and
`dmrg/idmrg/Manifest.toml`. The manifest pins the official repository commit
and the compatible ITensor versions. The root project and manifest do not gain
the infinite-MPS dependency. Checkpoints record the pinned backend commit
because their representation is not assumed to remain compatible with a future
backend revision.

All infinite commands use the isolated environment:

```bash
julia --project=dmrg/idmrg dmrg/idmrg/bin/run_vumps.jl ...
```

## Architecture and File Boundaries

The planned structure is:

```text
dmrg/
├── CylinderModelCore.jl
├── DMRGFullModel.jl
└── idmrg/
    ├── Project.toml
    ├── Manifest.toml
    ├── README.md
    ├── src/
    │   ├── InfiniteCylinderDMRG.jl
    │   ├── Configuration.jl
    │   ├── Hamiltonian.jl
    │   ├── VUMPSRunner.jl
    │   ├── Observables.jl
    │   ├── BranchTracking.jl
    │   └── Output.jl
    ├── bin/
    │   ├── run_vumps.jl
    │   └── run_flux_scan.jl
    ├── test/
    │   ├── runtests.jl
    │   ├── test_configuration.jl
    │   ├── test_hamiltonian.jl
    │   ├── test_observables.jl
    │   └── test_smoke.jl
    └── jobs/
        └── run_smoke.pbs
```

### Shared Model Core

`dmrg/CylinderModelCore.jl` owns only model data and transformations shared by
finite and infinite representations:

- `CylinderModelParams`;
- neighbor displacement tables;
- reciprocal-vector and two-band Fourier-transform utilities;
- the hopping matrices \(t_{\alpha\beta}(\Delta x,\Delta y)\);
- the periodic-\(y\) winding and twist-phase convention.

It does not depend on MPS, MPO, DMRG observers, VUMPS, checkpointing, or CLI
code. `DMRGFullModel.jl` imports and re-exports the existing public names so
finite callers retain their current API. Finite site-indexed term types,
finite-MPO construction, DMRG control, and finite observables remain in
`DMRGFullModel.jl`.

### Infinite Modules

`Configuration.jl` defines validated model, unit-cell, optimization, and flux
configuration types. It owns coordinate mapping and filling arithmetic.

`Hamiltonian.jl` converts shared real-space hopping and interaction data into
translation-normalized operator orbits and then into an
`InfiniteSum{MPO}`. It does not run VUMPS.

`VUMPSRunner.jl` is the only source file allowed to call experimental or
non-public ITensorInfiniteMPS optimization APIs. It owns product-state
construction, VUMPS iterations, subspace expansion, warm starts, and
convergence records. Pinning and isolating the backend makes this boundary
replaceable.

`Observables.jl` computes density, energy normalization, transfer spectra,
neutral-channel correlation length, entanglement spectra, and Schmidt charge
distributions.

`BranchTracking.jl` computes normalized mixed-transfer fidelity, compares
successive flux points, and labels possible branch events without modifying
the measured data.

`Output.jl` writes text diagnostics and HDF5 checkpoints with complete run
metadata and validity flags.

## Unit Cell and Filling

The MPS ordering remains the finite ordering: sites are ordered by \(x\), then
by \(y\). Within the reference MPS cell,

\[
  i(x,y)=x L_y+y+1,
  \qquad 0\le x<x_{\mathrm{period}},\quad 0\le y<L_y.
\]

For band filling \(\nu=p/q\),

\[
  N_{\mathrm{sites,cell}}=L_y x_{\mathrm{period}},
  \qquad
  N_{\mathrm{uc,cell}}=\frac{L_y x_{\mathrm{period}}}{2},
\]

\[
  N_{\mathrm{p,cell}}
  =\frac{p}{q}N_{\mathrm{uc,cell}}.
\]

`Ly` must be positive and even, `x_period` must be positive, and
`Np_cell` must be an integer. Invalid combinations fail before site-index or
Hamiltonian construction. The physical-site density is recorded separately as

\[
  \bar n=\frac{p}{2q}.
\]

For the production filling \(\nu=1/3\):

- `Ly=6, x_period=1` gives 6 sites, 3 two-band cells, and 1 particle;
- `Ly=6, x_period=3` gives 18 sites, 9 two-band cells, and 3 particles;
- `Ly=2, x_period=1` is invalid;
- `Ly=2, x_period=3` gives 6 sites, 3 two-band cells, and 1 particle.

The implementation supports configurable `x_period`. The `Ly=6` smoke test
uses `x_period=1`, while production studies must compare `x_period=1` and
`x_period=3`. Cell-resolved densities are never silently folded to a smaller
period.

QN-conserving infinite site indices use the centered-charge representation
created by `infsiteinds("Fermion", ...; initstate)`. A link QN measures charge
relative to the cell's mean background with an integer scale factor. Output
therefore stores both the raw QN and the converted physical excess charge. Raw
QN values are not labeled as particle numbers.

For \(N=N_{\mathrm{sites,cell}}\), \(N_p=N_{\mathrm{p,cell}}\), and
\(g=\gcd(N,N_p)\), the centered-charge scale is

\[
  s_Q=\frac{N}{g}.
\]

The empty and occupied site charges are respectively \(-N_p/g\) and
\(s_Q-N_p/g\). For the consistently oriented left Schmidt index, a raw sector
charge \(Q_{\mathrm{raw}}\) is converted to physical excess charge by

\[
  Q_{\mathrm{physical}}=\frac{Q_{\mathrm{raw}}}{s_Q}.
\]

The sign and normalization are covered by a product-state test at every
physical transverse cut.

## Infinite Hamiltonian Construction

The existing Fourier convention yields 11 nonzero displacement matrices for
the production hopping parameters. Their support is

\[
  \Delta x\in\{-2,-1,0,1,2\},
  \qquad
  \Delta y\in\{-1,0,1\}.
\]

The infinite builder represents a hopping path by source and target transverse
coordinates, its \(x\)-cell displacement, complex amplitude, and \(y\)-winding
number. The \(x\) direction is infinite and untwisted. The \(y\) direction is
periodic, and each path receives exactly the finite-model phase

\[
  \exp(i w_y \Phi_y).
\]

Each operator term is translated by whole MPS cells until its leftmost 1D
support site lies in the reference cell. Translation-equivalent terms are
deduplicated. Directed hopping terms retain both \(c_i^\dagger c_j\) and their
Hermitian conjugates. Density interactions are deduplicated as unordered
pairs. On-site hopping contributions are emitted as number operators.

The builder verifies Hermiticity after accumulation. It also exposes a
canonical term representation for tests and diagnostics before converting to
an ITensor `OpSum`. Fermionic Jordan-Wigner strings are generated by ITensor
from fermionic operators and are not implemented manually.

The twist is \(2\pi\)-periodic at the Hamiltonian level. `Phi_y` and
`Phi_y+2pi` must generate identical canonical terms within the numerical
tolerance inherited from the hopping transform.

## Finite/Infinite Hamiltonian Parity

The principal model regression test tiles infinite translation orbits over a
finite \(x\) window, removes terms crossing the open \(x\) boundary, and
compares the result with finite `build_model_terms`. It compares site pairs,
operator direction, interaction shell, complex amplitude, and winding phase.

Coverage includes:

- zero twist;
- a generic twist such as `0.37pi`;
- `2pi` twist;
- `x_period=1` and `x_period=3`;
- hopping and all enabled interaction shells.

The existing 19 finite-DMRG assertions must continue to pass after the shared
core extraction.

## Optimization Workflow

The initial infinite state is a QN-conserving product state with exactly
`Np_cell` occupied sites. The CLI accepts deterministic initial patterns and
seeds so different charge or crystal patterns can be tested reproducibly.

At a fixed Hamiltonian and target bond dimension, optimization alternates:

1. VUMPS iterations at the current bond dimension;
2. a measured subspace-expansion step;
3. VUMPS iterations at the expanded dimension.

The public backend VUMPS routine discards its precision residual on return.
Because the dependency commit is pinned, `VUMPSRunner.jl` wraps the backend's
iteration interface and records at every iteration:

- stage and iteration;
- link dimensions;
- left and right energy estimates;
- their mismatch;
- maximum left and right canonical precision errors;
- total VUMPS precision error;
- wall time.

The primary convergence condition is

\[
  \epsilon_{\mathrm{VUMPS}}
  =\max(\epsilon_L,\epsilon_R)<\texttt{vumps_tol}.
\]

The default convergence gate requires this condition for two consecutive
iterations, an energy-density change below `10*vumps_tol`, and a left/right
energy mismatch below `10*vumps_tol`. Each tolerance and the consecutive-count
requirement are recorded in output metadata. A
subspace-expansion stage records every link dimension before and after the
call. If the state is below the requested target dimension and no link grows,
the stage fails explicitly. This prevents the backend's known long-range QN
zero-expansion behavior from passing silently.

Nonconverged states may be checkpointed for diagnosis, but the default program
exits nonzero and marks all derived final results as nonconverged.

## Energy and Local Observables

Let the sum of the Hamiltonian contributions anchored in one repeating MPS
cell be \(e_{\mathrm{cell}}\). Output includes all of:

\[
  e_{\mathrm{cell}},\qquad
  e_x=\frac{e_{\mathrm{cell}}}{x_{\mathrm{period}}},
\]

\[
  e_{\mathrm{two-band\ uc}}
  =\frac{e_{\mathrm{cell}}}{N_{\mathrm{uc,cell}}},
  \qquad
  e_{\mathrm{site}}
  =\frac{e_{\mathrm{cell}}}{N_{\mathrm{sites,cell}}}.
\]

These values are not compared directly with a finite cylinder's total energy.
Finite/infinite comparisons use a finite-cylinder bulk slope or central local
energy after the implementation has established those estimators.

Density output retains every \((x,y)\) value in the MPS cell and the total
density of every \(x\) ring. This makes translation breaking within
`x_period=3` visible.

## Entanglement Cuts and Spectrum

Physical cuts are bonds perpendicular to \(x\), after each complete ring of
`Ly` sites. There is one such cut for `x_period=1` and three potentially
inequivalent cuts for `x_period=3`. Branch comparisons always use the same cut
label.

For a mixed-canonical center matrix \(C_b\), singular values define

\[
  p_\alpha=\frac{s_\alpha^2}{\sum_\beta s_\beta^2},
  \qquad
  \epsilon_\alpha=-\log p_\alpha,
\]

\[
  S_b=-\sum_\alpha p_\alpha\log p_\alpha.
\]

Each spectrum row contains cut, level, singular value, normalized probability,
entanglement energy, raw QN, and physical excess charge. A sector table stores

\[
  P_b(Q)=\sum_{\alpha\in Q}p_\alpha.
\]

The first implementation resolves U(1) charge but not transverse momentum.
No output is described as a `ky`-resolved entanglement spectrum.

## Transfer Spectrum and Correlation Length

The transfer matrix spans one complete repeating MPS cell. From its two
largest neutral-channel eigenvalues,

\[
  \xi_{\mathrm{cell}}
  =-\frac{1}{\log|\lambda_2/\lambda_1|},
  \qquad
  \xi_x=x_{\mathrm{period}}\xi_{\mathrm{cell}}.
\]

Output also records eigenvalue phases, eigensolver residuals, and the number of
converged eigenpairs. Because a QN eigensolve does not automatically search
every charge channel, the primary value is named
`correlation_length_neutral`. It is not presented as the maximum over all QN
sectors.

If fewer than two eigenpairs converge or the dominant eigenvalue is invalid,
the correlation length is marked invalid and the smoke test fails.

## Infinite-Cylinder Charge Pump

An infinite cylinder has no physical left or right edge. The finite cumulative
edge charge is therefore not copied into the infinite workflow. At a fixed
transverse cut, define the raw Schmidt polarization

\[
  Q_{\mathrm{Schmidt}}(\Phi_y)
  =\sum_\alpha p_\alpha(\Phi_y)q_\alpha^{\mathrm{physical}}.
\]

The raw pump diagnostic is

\[
  \Delta Q_{\mathrm{Schmidt}}(\Phi_y)
  =Q_{\mathrm{Schmidt}}(\Phi_y)-Q_{\mathrm{Schmidt}}(0).
\]

All flux points reuse the same site indices and cut convention. The workflow
stores raw QNs, converted charges, sector weights, and the QN scaling rule. It
does not add `1/3`, apply a hidden modulo operation, or unwrap a discontinuity.

Independent cold starts can use different QN gauge origins. A detected relative
charge shift is recorded as a separate diagnostic and is not applied to the
raw pump curve.

## Infinite-MPS Fidelity and Branch Tracking

The full overlap of distinct infinite MPS generally vanishes. Successive states
are compared using the dominant eigenvalue of their mixed transfer matrix. The
normalized fidelity per MPS cell is

\[
  f_{\mathrm{cell}}
  =\frac{|\lambda_{12}|}
  {\sqrt{|\lambda_{11}\lambda_{22}|}},
  \qquad 0\le f_{\mathrm{cell}}\le1.
\]

The reported fidelity density per \(x\) ring is

\[
  -\frac{\log f_{\mathrm{cell}}}{x_{\mathrm{period}}}.
\]

Mixed-transfer nonconvergence produces an invalid diagnostic, not a fabricated
number.

A flux point may retain a warm-start candidate, independent product-state
candidates, and refinement candidates loaded from an earlier checkpoint. Two
selection modes are explicit:

- `ground` selects the converged candidate with lowest energy density;
- `adiabatic` selects the converged candidate with greatest fidelity to the
  previously selected state while reporting its energy difference.

Every candidate is written to disk. Branch-event labels combine evidence from
mixed-transfer fidelity, raw Schmidt-polarization jumps, Schmidt-sector weight
changes, charge-resolved spectrum rearrangements, and candidate energy
crossings. Labels never modify observables.

For reproducible first-pass flags, every adjacent pair is written with its raw
metrics and these configurable defaults:

- fidelity flag when \(1-f_{\mathrm{cell}}>10^{-3}\);
- charge flag when the adjacent raw Schmidt-polarization change exceeds `0.1`;
- sector flag when the Bhattacharyya distance
  \(1-\sum_Q\sqrt{P_i(Q)P_{i+1}(Q)}\) exceeds `0.1`;
- spectrum flag when a charge-matched level among the lowest 20 levels changes
  by more than `0.5` in entanglement energy;
- energy-crossing flag when the energy ordering of two retained candidates
  reverses.

Cold-start sector tables are also compared under discrete raw-QN shifts to
report the shift that maximizes their Bhattacharyya overlap. That best shift is
stored as gauge/sector evidence but is never applied to the raw polarization.
Threshold values and the number of compared levels are recorded in scan
metadata.

## Output Contract

A single-point run writes:

```text
output/run_name/
├── summary.toml
├── convergence.tsv
├── density.tsv
├── entanglement_spectrum.tsv
├── schmidt_sectors.tsv
├── transfer_spectrum.tsv
└── state.h5
```

`summary.toml` records the full physical configuration, unit-cell arithmetic,
algorithm name, dependency versions and commit, optimization settings,
convergence state, energy normalizations, entropy, neutral correlation length,
and validity flags.

A flux scan writes:

```text
output/scan_name/
├── scan_summary.tsv
├── branch_events.tsv
├── phi_000/
│   ├── candidate_warm/
│   └── candidate_cold_01/
└── ...
```

Text tables are the analysis contract. HDF5 stores the infinite state and the
metadata required to validate a restart. A checkpoint loaded with incompatible
geometry, filling, site indices, or backend commit is rejected.

## Test and Smoke-Test Strategy

### Fast Tests

Fast tests cover configuration arithmetic, invalid filling combinations,
coordinate translation, canonical term deduplication, Hermiticity, twist
periodicity, finite/infinite term parity, and energy-normalization arithmetic.
They do not require VUMPS convergence.

### Minimal Backend Smoke

The first complete backend smoke uses:

```text
Ly=2
x_period=3
filling=1/3
t1=1.0
t3=0.2
V1=1.0
V2=0.0
V3=0.0
Phi_y=0
maxdim schedule = 1,4,8
cutoff = 1e-8
vumps_tol = 1e-5
max VUMPS iterations per stage = 50
transfer_tol = 1e-8
```

It must demonstrate an actual bond expansion, reach the configured smoke
residual, produce finite energy densities, converge at least two neutral
transfer eigenpairs, produce a finite positive correlation length, normalize
the entanglement probabilities to absolute error below `1e-10`, expose nonempty
Schmidt sectors, and reload an HDF5 checkpoint that reproduces scalar
observables to relative error below `1e-10`.

`Ly=2` is only a backend and data-contract smoke test. Periodic wrapping makes
it a narrow, geometrically degenerate cylinder, so its observables are not used
as evidence for the `Ly=6` phase.

### `Ly=6` Target Smoke

The target smoke uses:

```text
Ly=6
x_period=1
Np_cell=1
t1=1.0
t3=0.2
V1=1.0
V2=0.0
V3=0.0
Phi_y=0
maxdim schedule = 1,4,8
cutoff = 1e-8
vumps_tol = 1e-4
max VUMPS iterations per stage = 50
transfer_tol = 1e-8
```

It validates the real circumference and complete output path. Its small bond
dimension is not presented as a production ground state. The run must still
meet its declared residual and output energy density, neutral correlation
length, charge-resolved entanglement levels, and Schmidt-sector diagnostics. A
backend failure to expand or converge is a failed acceptance criterion, not a
successful partial result.

## Documentation Contract

`dmrg/idmrg/README.md` explains installation, pinned dependencies, local and
PBS commands, output files, and the finite/infinite mapping:

| Finite DMRG | Infinite VUMPS |
|---|---|
| finite `Lx` | no `Lx`; infinite translation |
| `Ly` | the same `Ly` |
| `Ns=Lx*Ly` | `sites_per_cell=x_period*Ly` |
| total `Np` | integer `Np_cell` |
| open \(x\) | infinite \(x\) |
| periodic/twisted \(y\) | identical twist convention |
| total energy | energy densities |
| cumulative edge charge | Schmidt-charge polarization |
| finite-state overlap | fidelity per infinite cell |
| sweeps and cutoff | VUMPS residual and expansion stages |

The README repeats that the finite `Lx=15` jump is raw data and describes how
the new infinite diagnostics can support or reject a branch-switch
interpretation without manually shifting the charge.

## Risks and Explicit Non-Goals

The implementation accepts these risks:

- ITensorInfiniteMPS is experimental, unregistered, and pinned to a commit;
- the VUMPS adapter uses a pinned internal iteration interface to retain
  residuals;
- ITensorMPS 0.3.45 is isolated from the finite root environment;
- QN-conserving long-range subspace expansion has a documented open issue;
- VUMPS can converge to metastable branches;
- `x_period=1` can suppress crystalline order;
- neutral correlation length is not a full sector maximization;
- finite edge charge and infinite Schmidt polarization are related but
  different observables;
- HDF5 checkpoints are tied to the pinned backend representation.

The first implementation does not include:

- a growth-based iDMRG algorithm;
- transverse-momentum-resolved entanglement spectra;
- automatic `1/3` or modulo unwrapping;
- a claim of production convergence at bond dimension 2000;
- GPU or distributed-memory support;
- automatic selection of the physically correct topological sector from a
  single candidate.

## Acceptance Criteria

The feature is complete only when:

1. the finite DMRG regression tests still pass;
2. finite/infinite canonical Hamiltonian terms agree at the tested twists;
3. the isolated backend environment instantiates reproducibly;
4. the minimal backend smoke passes every validity check;
5. the `Ly=6` target smoke converges to its declared tolerance and writes all
   required observables;
6. checkpoint reload reproduces the saved diagnostics;
7. README documents parameter mapping, dependency risks, output semantics, and
   raw branch-jump policy.
