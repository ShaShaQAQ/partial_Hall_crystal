# MPSKit iDMRG Backend Design

**Date:** 2026-08-31

**Status:** Approved continuation design

## Purpose

Add a second, isolated infinite-MPS backend for the paper-straight Fig. 2
benchmark. The existing ITensor/ITensorInfiniteMPS implementation remains in
place as a diagnostic and regression reference, but it is not accepted as the
production solver for the 36-site paper unit cell.

The new backend converts the already-validated, backend-independent term lists
from `build_infinite_model_terms` into a fermionic MPSKit Hamiltonian, grows a
state with `IDMRG2`, refines it with VUMPS, and accepts convergence only from a
fresh MPSKit Galerkin residual. It must preserve the current checkpoint,
momentum-entanglement, raw-pump, provenance, and PBS audit requirements.

## Evidence for the Backend Change

The pinned ITensorInfiniteMPS backend is version 0.2.2 at commit
`765f2777703bc1138b009adbed1b97bde1973402`. Its upstream `main` branch is still
that commit and therefore offers no upgrade path.

The W003 paper-straight diagnostics showed that the 36-site `InfiniteSum{MPO}`
update is not a tolerance or sweep-origin problem:

- the retained parallel VUMPS update has a near-minus-one fixed-point mode and
  a residual two-cycle;
- cyclic sweep origins, reverse sweeps, parallel TDVP, and damped parallel
  updates did not produce a lower consistent residual;
- the final consistent residual after the `alpha=0.5` damped run was
  `1.158707978e-3`, compared with `1.144636e-3` for its source state, even though
  the reported energy decreased monotonically;
- upstream tests primarily cover one- through four-site cells, and the sole
  six-site FQHE test explicitly states that VUMPS cannot give the right result.

This is an architectural failure gate: energy monotonicity is not convergence,
and no further scalar damping, origin, or tolerance variants are eligible for
production evidence.

MPSKit commit `811ecf6c06c1f7c1bc656da61abcd679effcd428` provides maintained
`IDMRG2` and VUMPS implementations, recomputed environments, a VUMPS Galerkin
residual, U(1)-symmetric graded fermionic spaces, and long-range Jordan MPOs.
TensorKitTensors commit `3755705a1c44a3d5e32086e7d89b2c561b268cb1`
provides the spinless-fermion number and two-site hopping tensors used by the
adapter. Both commits are pinned because their relevant APIs are still moving.

The immediate MPSKit successor `4b579944066e9b4bf814205a99f9c24014f76885` changes
TimerOutputs support from 0.5 to 1 and rewrites timer plumbing, but does not
change the IDMRG2 or VUMPS numerical updates. ITensors 0.9.30 requires
TimerOutputs 0.5, so that successor cannot resolve in the legacy regression
environment. A source diff over IDMRG, VUMPS, environments, and gauge code
confirms that `811ecf6c06c1f7c1bc656da61abcd679effcd428` is the latest commit containing all
required numerical fixes while retaining TimerOutputs 0.5 compatibility.
BlockTensorKit 0.3.16 is pinned at
`230cec77c9c7527817d2216b9c6de87f6d8bdda8` because the W003 registry does not
yet contain that release.

## Scope and Isolation

The implementation is additive and lives under `dmrg/idmrg`:

```text
dmrg/idmrg/
├── src/
│   ├── MPSKitBackend.jl
│   ├── MPSKitHamiltonian.jl
│   ├── MPSKitSolver.jl
│   ├── MPSKitCheckpoint.jl
│   ├── MPSKitObservables.jl
│   └── MPSKitFig2Adapter.jl
├── test/
│   ├── test_mpskit_hamiltonian.jl
│   ├── test_mpskit_solver.jl
│   ├── test_mpskit_checkpoint.jl
│   ├── test_mpskit_observables.jl
│   └── test_mpskit_fig2_adapter.jl
└── jobs/
    └── run_mpskit_gate.pbs
```

The existing files `Hamiltonian.jl`, `Configuration.jl`, `PaperGeometry.jl`,
and `InitialCandidates.jl` remain the source of model truth. The MPSKit adapter
must not reimplement geometry, hopping transforms, interaction shells, flux
winding, or candidate generation. The finite-DMRG files and user-owned
untracked files remain untouched.

Qualified imports are mandatory:

```julia
import MPSKit
import TensorKit
import TensorKitTensors
import JLD2
```

This avoids collisions with the existing `InfiniteMPS`, MPO, and tensor names.

## Dependency Contract

The isolated `dmrg/idmrg/Project.toml` and `Manifest.toml` retain the pinned
ITensor stack for regression coverage and add exact source pins for MPSKit and
TensorKitTensors. `TensorKit` and `JLD2` are direct dependencies because the
adapter imports them directly. Dependency resolution and all Julia execution
occur on W003, never on the Mac.

Every result records:

- repository commit and dirty-state audit;
- Julia version and thread count;
- MPSKit, TensorKit, TensorKitTensors, MatrixAlgebraKit, JLD2, and legacy
  ITensor backend versions/commits;
- benchmark manifest bytes and SHA-256;
- PBS job ID, resource request, wall time, CPU time, and maximum RSS;
- solver stages, requested/achieved bond dimensions, energy, and residuals.

## Uniform Centered Charge Representation

The paper cell has 36 spinless-fermion sites and 14 particles. A site-specific
application of `MPSKit.add_physical_charge` would encode the background charge
on 14 chosen sites. Although physically gauge-equivalent, that representation
would make the local symmetric spaces site-dependent and obscure the
six-fold transverse translation needed by the momentum-resolved entanglement
spectrum. It is therefore not used.

Instead, the adapter extends the centered-charge convention already used by
the ITensor backend. For `N=36`, `Np=14`, and `g=gcd(N,Np)=2`, the raw charge
scale is `sQ=N/g=18`. Every physical site uses the same graded symmetric space

In TensorKit notation this is a `Vect` over the product sector
`ProductSector{Tuple{FermionParity,U1Irrep}}`, with sector multiplicities
`(0,-7)=>1` and `(1,11)=>1`; the implementation keeps the names qualified
through a small constructor helper.

where an empty site has raw charge `-7` and an occupied site has raw charge
`+11`. A 36-site cell with 14 occupied sites has total raw charge zero:

```text
14 * 11 + 22 * (-7) = 0.
```

The physical excess charge reported on a Schmidt bond is
`Q_physical = Q_raw / 18`, with a single documented bond orientation. All 36
physical spaces are identical, so transverse translation is represented in
one common symmetric Hilbert space.

The adapter symmetrizes the TensorKitTensors dense graded fermion operators
onto this custom centered space. Tests must prove that this preserves the dense
matrix elements of `n`, `c†_i c_j`, and their Hermitian combinations.

## Hamiltonian Adapter

`build_infinite_model_terms(config, params)` returns directed hopping terms and
undirected density interactions. The MPSKit adapter consumes those terms
without changing their amplitudes.

For an ordered pair `i < j`, let `a` be the coefficient of
`c†_i c_j`; Hermiticity supplies `conj(a)` for `c†_j c_i`. TensorKitTensors
defines `f_min_f_plus = -adjoint(f_plus_f_min)`, so the pair is represented as

```julia
a * f_plus_f_min - conj(a) * f_min_f_plus
```

on sites `(i,j)`. On-site hopping is `a*n_i`, and an interaction is
`V*(n_i ⊗ n_j)`. Each Hermitian directed pair is emitted once. MPSKit's
`InfiniteMPOHamiltonian` inserts identity/braiding tensors across nonadjacent
sites and across the repeating-cell boundary.

The adapter exposes both finite and infinite constructors. The finite
constructor is a test oracle: for two through six sites, converting the MPO to
a dense matrix must agree element by element with an independently assembled
Jordan-Wigner Fock-space Hamiltonian. Required cases include nonadjacent
hopping, complex hopping, boundary-crossing terms, density interactions, and
both centered and uncentered charge spaces. The test also compares the
Hamiltonian at `phi` and `phi+2pi` and checks the complete 36-site term count and
Hermiticity without densifying that cell.

## Initial States and Filling

An initial candidate is a 36-entry occupation pattern with exactly 14 occupied
sites. Its one-dimensional virtual spaces carry the cumulative centered raw
charge. At site `i`, the unique product tensor maps the preceding cumulative
sector and the site's physical raw charge to the next cumulative sector. The
final cumulative sector must return to zero.

Tests cover:

- exact 14-particle arithmetic and zero total raw charge;
- deterministic reconstruction of every configured Fig. 2 product candidate;
- rejection of wrong length, wrong particle count, or a nonclosing charge path;
- identical physical spaces at all sites;
- a unit-cell transverse translation that preserves the centered charge
  representation.

## Solver Contract

The solver has two distinct stages:

1. `MPSKit.IDMRG2` grows/truncates the bond dimension with an explicit
   MatrixAlgebraKit truncation strategy. Its returned error is recorded as an
   IDMRG sweep diagnostic only; it is not the production convergence gate.
2. `MPSKit.VUMPS` refines the state. Its returned `calc_galerkin` norm is the
   primary residual. Environments are recomputed for the final state and the
   Galerkin residual is evaluated again before acceptance.

A state is converged only when all of the following hold:

- the independently recomputed final Galerkin residual is finite and no larger
  than the manifest tolerance;
- the requested bond-dimension stage was reached, subject only to explicitly
  recorded symmetry-block rank limits;
- the energy density is finite and its imaginary part is below tolerance;
- a restart-and-refine process reproduces the energy and residual gates;
- required entanglement and transfer observables are valid.

IDMRG2's change-in-center-tensor norm, VUMPS log text, energy decrease, or wall
time can never substitute for the Galerkin gate.

The first numerical decision gate is the full 36-site paper model at `D=32`.
It must start from MPSKit-native candidates and reach a materially lower,
nonoscillatory Galerkin residual than the failed ITensor source residual. No
`D>=64` production sequence starts until this gate is audited.

## Checkpoint and Restart

MPSKit checkpoints use a new versioned payload while retaining the workflow
filename `state.h5` for compatibility with the existing directory contract.
JLD2 may write the payload regardless of filename extension. The checkpoint
contains the MPS tensors, centered-charge metadata, configuration signature,
solver/backend provenance, last completed stage, convergence history SHA-256,
and an atomic-completion marker.

Checkpoint support is accepted only after a two-process W003 test:

1. process A saves a nontrivial symmetric MPS and exits;
2. process B loads it in a fresh Julia process, validates metadata and tensor
   spaces, recomputes energy/residual, performs at least one refinement step,
   and saves a successor checkpoint;
3. a negative test rejects a legacy-ITensor checkpoint, wrong configuration,
   truncated file, or backend-pin mismatch with a specific error.

If JLD2 cannot round-trip the pinned TensorKit/MPSKit structures across
processes, the fallback is an explicit tensor-block schema in JLD2. Julia's
opaque `Serialization` format is not production evidence.

## Observables

The backend must emit the existing table-level contract with backend-neutral
field meanings:

- energy per MPS cell, x period, two-band unit cell, and physical site;
- every site density and every x-ring density;
- Schmidt values, normalized probabilities, entanglement energies, centered
  raw charge, physical excess charge, and charge-sector weights at each
  physical x cut;
- neutral transfer eigenvalues and correlation length in cell and x units;
- transverse-momentum-resolved entanglement levels and validation residuals;
- normalized mixed-transfer fidelity for adjacent flux states.

The MPSKit center tensor `C` is the source of Schmidt data. Transfer and mixed
transfer solves must record convergence counts and residuals; insufficient
eigenpairs produce an invalid diagnostic, not a fabricated correlation length.

For transverse momentum, the six-site ring translation acts within each
Schmidt-charge block. Fermionic permutation signs are included exactly. A
momentum label is valid only when the restricted translation is unitary within
tolerance, its sixth power is the identity within tolerance, and its
commutator/closure residual in the retained Schmidt subspace passes the
predeclared gate. The expected `1,1,2,3,5` counting is evaluated only in the
declared low-entanglement-energy window and never imposed by relabeling levels.

## Flux Pump and Branches

The infinite cylinder has no edges. Pumping is measured from the raw centered
Schmidt polarization at one fixed physical x cut. The workflow stores raw
charge sectors, physical excess charges, sector weights, entanglement levels,
momentum labels, mixed-transfer fidelity, and candidate energy crossings.

No observable is manually shifted, modulo-reduced, or unwrapped. A discrete
relative sector-origin offset may be reported as a separate gauge diagnostic,
but it is never applied to the raw pump curve. The target comparison is one
electron transported over `phi_y = 0 ... 6pi` together with consistent sector
flow and low-lying momentum counting.

## Workflow Adapter and Audit

The benchmark manifest gains an explicit backend identifier. The default
paper-production backend is `mpskit_idmrg_v1`; the legacy backend can only run
when explicitly selected and its results are marked diagnostic.

`MPSKitFig2Adapter.jl` supplies the existing `Fig2BenchmarkOperations`
callbacks for candidate construction, optimization, checkpoint load/audit,
observable extraction, and provenance. Ledger, selection, raw-data, and
acceptance logic remain backend-neutral. Hard-coded assumptions about ITensor
checkpoint internals are moved behind the operation boundary.

Acceptance requires all of:

- manifest and production-checkout provenance;
- exact Hamiltonian/Jordan-Wigner parity gates;
- two-process checkpoint/restart gate;
- converged selected states at every accepted dimension/flux point;
- energy and residual convergence across bond dimensions;
- raw pump, sector flow, momentum counting, and transfer diagnostics;
- auditable PBS resource records and checkpoint hashes;
- no protected finite-DMRG/user file staged or modified.

## W003 and PBS Execution

All Julia commands run on W003. Unit and very small tensor tests may run in an
interactive W003 shell; every solver gate and production run uses PBS queue
`cmt`, one node, 24 requested CPU cores, and an explicit wall time. The Julia
thread count is benchmarked within that allocation; the current measured
default is four Julia threads, not 24, because it was faster for the existing
tensor workload.

Every solver job writes progress atomically and checkpoints often enough to
resume within its wall time. The D=32 gate requests a conservative wall time;
larger dimensions use staged jobs with predecessor checkpoints. Wavefunctions
are never discarded merely because a convergence gate fails.

## Risks and Stop Gates

- **Unregistered API risk:** exact commits are pinned and all adapter calls are
  isolated. A pin mismatch invalidates provenance.
- **Fermionic sign risk:** no solver run is trusted until elementwise dense
  parity passes for nonadjacent and complex hopping.
- **Centered-space risk:** physical spaces, charge conversion, product-state
  closure, and transverse translation are tested independently.
- **Checkpoint risk:** no restart claim is made from a same-process round trip;
  the two-process W003 gate is mandatory.
- **Solver risk:** IDMRG2 diagnostics are not convergence. Final VUMPS Galerkin
  residual and a restart refinement are mandatory.
- **Memory/time risk:** the 36-site MPO bond dimension and D=32 memory are
  measured before higher D. A gate exceeding one 24-core node is redesigned,
  not silently submitted across nodes.
- **Momentum risk:** levels without validated translation closure remain
  charge-resolved only and cannot support the paper counting claim.

## Completion Evidence

The backend is complete only when the repository contains passing W003 test
logs, restart-gate metadata, PBS job provenance, retained checkpoints, final
tables/plots, and `dmrg/idmrg/report/fqahc_fig2_benchmark.md`. That report must
compare numerical energy/residual convergence, one-electron raw sector flow,
and `1,1,2,3,5` momentum counting with the paper, state every discrepancy, and
link the exact result hashes. The `DMRG` branch is then verified clean except
for the protected user files and synchronized to `origin/DMRG`.
