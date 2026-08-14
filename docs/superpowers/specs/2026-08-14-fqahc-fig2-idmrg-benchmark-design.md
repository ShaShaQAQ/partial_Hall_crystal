# FQAHC Fig. 2 iDMRG Benchmark Design

**Date:** 2026-08-14

**Status:** Approved design

## Purpose

Reproduce the infinite-cylinder results in Figs. 2(d)-2(f) of Hongyu Lu,
Han-Qing Wu, Bin-Bin Chen, Wang Yao, and Zi Yang Meng, "Generic integer and
fractional quantum anomalous Hall crystals from interaction-driven band
folding," *Physical Review B* **113**, 245105 (2026), DOI
`10.1103/mh3v-fb9w`.

The primary benchmark is the fermionic fractional quantum anomalous Hall
crystal at

```text
Ny = 6 original two-orbital unit cells around the cylinder
Nx = 3 original two-orbital unit cells in the infinite-MPS repeat
nu = 7/9, nu_star = 1/3
t = 1, t_prime = 0.2
V1 = 10, V2 = V3 = 2
maximum published bond dimension D = 3000
flux interval phi_y = 0 ... 6*pi
```

The benchmark succeeds only when a converged calculation reproduces all of
the following:

1. one electron pumped after inserting `6*pi` flux, within absolute tolerance
   `0.05`;
2. a shift of the entanglement spectrum by one measured charge sector;
3. low-lying momentum-resolved entanglement counting `1, 1, 2, 3, 5`;
4. a pump curve consistent with a digitization of Fig. 2(d), within the
   digitization uncertainty;
5. conclusions stable under the final bond-dimension progression through
   `D = 1000, 2000, 3000` unless a documented earlier plateau makes the last
   stage numerically redundant.

The paper does not publish a raw numerical energy table. Energy density is a
convergence and branch diagnostic, not an independently quoted paper target.

## Source Evidence

The paper defines a two-site unit cell with primitive vectors

```text
a1 = (1, 0)
a2 = (0, sqrt(3))
```

and states that iDMRG uses up to `Ny = 6` unit cells along the periodic `a2`
direction, mainly `Nx = 3` for the iDMRG unit cell, charge conservation,
`D <= 3000`, and truncation errors of order `1e-6`. The calculation acts on the
unprojected two-band Hamiltonian.

For the target filling, the physical density is `nbar = nu/2 = 7/18`. The
paper-exact MPS repeat therefore contains

```text
2 * Nx * Ny = 36 physical fermion sites
Nx * Ny = 18 original two-orbital unit cells
(7/9) * 18 = 14 fermions
```

The Zotero record used to verify the paper and its supplement is item
`A7PSUI68`. The main PDF attachment is `9ZH9F49U`; the supplemental attachment
is `2ZGP5FE8`.

## Non-goals

- Do not use the existing finite `Lx=15`, `Np=36/45=4/5` topological-pinball
  jobs as a substitute for the `nu=7/9` Fig. 2 benchmark.
- Do not project the Hamiltonian into a Chern band or folded miniband.
- Do not modify the validated finite-DMRG behavior to accommodate the paper
  geometry.
- Do not call a tilted or helical cylinder a paper-exact result.
- Do not manually unwrap, offset, modulo-reduce, or relabel raw pump data.
- Do not require production computation on the user's Mac.

## Execution Boundary

The Mac is used only for reading, editing, Git commits, and Git transport. No
Julia unit test, smoke test, tensor contraction, or benchmark calculation is
run there after approval of this design.

All code tests and numerical calculations run on W003. The existing
`/home/public/shajy/partial_Hall_crystal` directory is a dirty production
workspace on `master` and must not be checked out, cleaned, or repurposed. A
separate W003 clone or worktree is used for this benchmark.

Large or tensor-intensive work is submitted through PBS to the `cmt` queue.
The observed compute-node contract is one node, 24 CPU cores, and about 93 GB
RAM. Existing user jobs must not be cancelled, renamed, or overwritten.

## Geometry Model

### Existing geometry

The existing infinite-cylinder implementation labels `Ly` physical triangular
lattice sites and groups two adjacent transverse sites into a two-orbital unit
cell. Its transverse repeat is the sheared vector `(1, sqrt(3))`. Thus the old
`Ly=6` contains only three two-orbital cells and does not represent the paper's
`Ny=6` cylinder.

### Paper geometry

Add a first-class `paper_straight` cylinder embedding whose periodic vector is
`Ny * (0, sqrt(3))`. Its site coordinate is

```text
(x_cell, y_cell, orbital),
x_cell in Z,
y_cell in 0:(Ny-1),
orbital in (A, B).
```

The physical position uses the paper basis, with the B orbital displaced by a
triangular-lattice nearest-neighbor vector. The reference MPS cell orders all
`2*Ny` sites of an `x_cell` ring contiguously and repeats `x_period=Nx` rings.

The existing sheared geometry remains available and is the default for old
CLI calls and old checkpoints. Paper benchmark commands must name
`paper_straight` and `Ny=6` explicitly. Output metadata stores both the
geometry name and the derived physical-site circumference.

### Hamiltonian construction

The paper geometry constructs the full real-space hopping orbits in the paper
unit-cell basis. The construction must reconstruct the repository's two-band
Bloch Hamiltonian at arbitrary momenta and must reproduce the hopping phases
in Eq. (1) and Fig. 1(a) of the paper. The `t_prime=0.2` hopping and all three
density-interaction shells are retained.

The transverse twist multiplies a hopping by

```text
exp(i * winding_y * phi_y)
```

where `winding_y` is computed around the straight paper circumference. The
geometry implementation, not an after-the-fact phase patch, owns winding
calculation.

Before VUMPS is allowed to run, automated tests verify:

- reconstruction of `H(k)` at deterministic generic momenta;
- exact geometry arithmetic for 36 sites, 18 unit cells, and 14 particles;
- Hermiticity of translation-normalized hopping terms;
- independently enumerated neighbor-shell and finite-window parity;
- straight-cylinder seam connectivity;
- equality of Hamiltonian terms at `phi_y` and `phi_y + 2*pi` to numerical
  tolerance;
- incompatibility between paper-straight and legacy checkpoint signatures.

## Configuration and Compatibility

The generic infinite-cylinder configuration gains an explicit geometry tag
and an unambiguous circumference in original two-orbital unit cells. Existing
`--Ly` calls remain valid for the legacy geometry. Paper benchmark entry points
use `--Ny=6` or a versioned benchmark TOML and reject simultaneous ambiguous
`Ly`/`Ny` input.

The HDF5 checkpoint format is versioned when geometry becomes part of state
compatibility. Loading an old checkpoint assigns the legacy geometry only
when all old-format invariants pass. A paper-straight checkpoint cannot be
loaded as a legacy state, or vice versa.

Every scalar and table output records:

- geometry and coordinate convention;
- `Ny`, `x_period`, sites, two-orbital cells, particles, and physical density;
- model couplings and flux;
- code commit, Julia version, package versions, and pinned backend commit;
- PBS job ID, node, requested threads, actual thread settings, and timing.

## Ground-state Candidate Strategy

The target state contains a period-three CDW background plus doped particles
in the interaction-generated topological miniband. A single evenly spaced
product state is not a sufficient initialization strategy.

At zero flux, generate a bounded deterministic candidate set containing:

- the three translations of the period-three CDW background;
- symmetry-related placements of the two particles above the `nu=2/3`
  background;
- the existing deterministic filling-compatible pattern as a control.

All patterns contain exactly 14 particles in 36 sites. Candidate generation is
unit tested; the benchmark manifest records every pattern. At `phi=0`, the
selected ground candidate is the converged state with the lowest finite energy
per physical site. Raw states and outputs for all candidates are retained.

At later flux values, retain a warm candidate and a bounded cold-candidate
set. The adiabatic branch is selected by valid normalized mixed-MPS fidelity
among converged candidates. Candidate eligibility is not invalidated by an
unrelated optional observable failure.

## Flux and Pump Data

The same site-index identities are reused across all flux values in a process.
The primary scan runs from `0` through `6*pi`, inclusive.

Each point stores:

- raw Schmidt polarization at the fixed physical transverse cut;
- raw Schmidt QNs, physical excess charges, and sector weights;
- charge-resolved entanglement levels;
- mixed-state fidelity to the prior selected state;
- all candidate energies and selection evidence;
- raw sector-alignment diagnostics without applying them to raw data.

The paper-comparison pump is derived only from the selected adiabatic state and
measured Schmidt-sector flow. A cumulative sector shift may be reported only
when an overlap-supported entanglement-spectrum mapping identifies it. The raw
polarization and raw sectors remain unchanged alongside the derived curve.

The sign convention is fixed once from the paper's flux orientation, cut
orientation, and hopping phase. A sign difference may be explained by an
explicit orientation map, but data are never silently multiplied by `-1`.

## Momentum-resolved Entanglement Spectrum

Charge QNs alone do not reproduce Fig. 2(f). Add a transverse-translation
diagnostic for the `paper_straight` cylinder.

For one original unit-cell translation around the cylinder:

1. construct the fermionic cyclic-translation operator, including fermion
   permutation signs;
2. build its action on the half-infinite Schmidt space through a
   symmetry-twisted transfer problem;
3. form the virtual transverse-translation operator in each charge sector;
4. validate unitarity, `Ty^Ny = I`, the transfer residual, and commutation with
   the Schmidt density matrix;
5. diagonalize only within charge and numerically degenerate Schmidt blocks;
6. quantize valid eigenphases to `ky = 2*pi*m/Ny` and attach a residual to every
   label.

If any validation fails, momentum labels are invalid and the code must not
manufacture a counting sequence. Tests compare the fermionic translation
operator against exact small-ring Fock-space translations and exercise
degenerate, nondegenerate, symmetry-broken, and invalid cases.

The accepted final spectrum must show the low-lying counting
`1, 1, 2, 3, 5` under a documented energy-window and momentum-sector rule that
is fixed before inspecting the production result.

## Checkpoint and Resume Contract

Every candidate at every flux point is atomically checkpointed as a canonical
infinite MPS. Text/TOML/TSV data are written before the driver advances to the
next point. A scan ledger records the last fully committed point and selected
candidate.

The pinned backend currently reproduces a QN link-direction failure when a
loaded state is used for further optimization. Production segmentation is
forbidden until a regression test demonstrates this sequence in two Julia
processes:

1. optimize and save a state;
2. terminate the first process;
3. load and validate the state in a new process;
4. run a no-expansion VUMPS iteration;
5. run a controlled expansion stage;
6. reproduce pre-restart energy, density, entanglement spectrum, and Schmidt
   sectors within stated tolerances;
7. advance to the next flux using the loaded site's exact index identities.

The restart failure is debugged at the link-direction source. A loader that
merely suppresses the error, discards QNs, or rewrites physical state data is
not acceptable.

Until this gate passes, only single-job pilot scans are allowed. Their
checkpoints remain valuable for audit and observable reload, but are not
advertised as optimization-resumable.

## W003 Runtime and Repository Isolation

Create a dedicated W003 clone/worktree from the synchronized `DMRG` branch.
Do not operate on the dirty production `master` tree. Output is written to a
separate benchmark-data root outside the Git worktree.

The locally validated runtime is Julia 1.12.5, whereas W003 currently defaults
to Julia 1.10.10. Install or select a user-owned Julia 1.12.5 runtime on W003
and invoke it by an explicit path recorded in PBS metadata. Use a dedicated
task depot variable rather than modifying a global system environment. The
pinned project must instantiate to the committed manifest and report the exact
dependency provenance before tests run.

## PBS Resource Policy

All Julia tests run on W003. Fast non-tensor tests use short PBS jobs with only
the cores they need. Tensor/VUMPS work uses the `cmt` queue and one complete
24-core node:

```text
#PBS -q cmt
#PBS -l nodes=1:ppn=24
JULIA_NUM_THREADS=24
OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
```

ITensor threaded block-sparse contractions are enabled. Nested BLAS/OMP
threading is disabled. Each stage records elapsed time, per-iteration timing,
maximum resident memory when available, and scheduler resource usage. A short
thread-utilization pilot checks that the 24-core configuration is actually
beneficial; allocation size is not treated as proof of parallel efficiency.

Initial walltime caps are:

| Stage | Bond-dimension cap | Flux coverage | Walltime cap |
|---|---:|---|---:|
| environment and fast tests | none | none | 2 h |
| Hamiltonian/geometry tensor smoke | small | selected points | 4 h |
| sparse-flux pilot | 64-128 | 7 points | 12 h |
| convergence pilot | 256 | sparse/dense decision | 36 h |
| intermediate scan | 512-1000 | full coarse scan | 72 h |
| production | 2000-3000 | full/refined scan | 120 h |

Walltimes are adjusted using measured pilot timings, not guesswork. Jobs use
unique names and output roots. After restart is validated, segments are linked
with PBS `afterok` dependencies and never write to the same candidate
directory concurrently.

## Staged Numerical Gates

### Gate 1: geometry

All paper-geometry Hamiltonian and configuration tests pass on W003. No VUMPS
result can waive this gate.

### Gate 2: state validity

At zero flux, a small-bond state must have:

- exact centered-QN charge conservation;
- valid canonical equations;
- finite consistent left/right energy density;
- recorded bond growth and expansion diagnostics;
- a recognizable period-three CDW density pattern among translated candidates.

Low bond dimension is not expected to reproduce the final topological data.

### Gate 3: restart

The cross-process optimization-restart contract passes before segmented
production is submitted.

### Gate 4: sparse flux

Run the inclusive grid

```text
phi_y/(2*pi) = 0, 0.5, 1, 1.5, 2, 2.5, 3.
```

This gate checks branch continuity, fidelity, sector flow, observable validity,
and runtime. It does not require quantization at low bond dimension.

### Gate 5: convergence and refinement

Increase bond dimensions through `256`, `512`, `1000`, `2000`, and `3000`.
Use a denser flux grid only after sparse-branch behavior is understood. At the
final stages, require stable endpoint pump, unchanged low-lying ES counting,
and stable energy density within predeclared tolerances.

## Acceptance Report

The benchmark root contains a machine-readable acceptance table and a human
readable note. It includes:

- the paper citation, figure, and digitized reference with uncertainty;
- exact geometry and filling arithmetic;
- every PBS job ID and dependency;
- code/runtime/dependency provenance;
- convergence versus bond dimension;
- raw and tracked pump tables;
- charge- and momentum-resolved entanglement spectra;
- CDW density diagnostics;
- checkpoint/restart validation;
- pass/fail status for every acceptance criterion.

An apparent match from a nonconverged state, an invalid momentum assignment,
or a manually shifted pump is reported as a failure.

## Debugging Policy

When a gate fails, investigate in this order:

1. paper-versus-code geometry and unit-cell arithmetic;
2. hopping phases, interaction shells, winding, and flux orientation;
3. QN/site-index and checkpoint compatibility;
4. product-state/CDW candidate initialization;
5. VUMPS residual, energy mismatch, and bond growth;
6. mixed-fidelity branch selection and Schmidt-sector flow;
7. transverse-translation and momentum-label residuals.

Each iteration states one root-cause hypothesis and changes one variable.
Production fixes require a failing regression test first and must pass the
full W003 test suite before a larger PBS stage is submitted. After three failed
fix attempts in different components, stop and review the architecture rather
than stacking a fourth patch.

## Approved Decisions

The user approved:

- the Fig. 2(d)-2(f) `nu=7/9`, `V1=10`, `V2=V3=2` target;
- pump, charge-sector, momentum-counting, curve, and bond-dimension acceptance
  criteria;
- paper-straight geometry isolated from finite DMRG;
- raw-data and adiabatic branch policy;
- W003-only testing and computation;
- small-scale-first progression;
- `cmt` single-node 24-core tensor jobs with staged walltimes;
- per-candidate wavefunction/checkpoint preservation and a restart gate.
