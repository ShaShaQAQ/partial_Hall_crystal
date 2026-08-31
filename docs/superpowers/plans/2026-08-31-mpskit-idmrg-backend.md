# MPSKit iDMRG Backend Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the failed 36-site production solver path with a pinned MPSKit IDMRG2-plus-VUMPS backend while preserving the validated paper Hamiltonian, raw pump policy, momentum entanglement spectrum, checkpoint/restart audit, PBS workflow, and finite-DMRG files.

**Architecture:** Keep the existing ITensor backend as a diagnostic reference and add qualified MPSKit/TensorKit adapters under `dmrg/idmrg`. The adapter consumes `build_infinite_model_terms`, uses one uniform centered U(1)-and-fermion-parity physical space, grows with IDMRG2, gates with an independently recomputed VUMPS Galerkin residual, and plugs into the backend-neutral Fig. 2 ledger through operation callbacks.

**Tech Stack:** Julia 1.12.5 on W003 only; MPSKit commit `811ecf6c06c1f7c1bc656da61abcd679effcd428` (the last numerically equivalent pre-TimerOutputs-1 commit); BlockTensorKit commit `230cec77c9c7527817d2216b9c6de87f6d8bdda8`; TensorKitTensors commit `3755705a1c44a3d5e32086e7d89b2c561b268cb1`; TensorKit 0.17; MatrixAlgebraKit; JLD2; existing pinned ITensors stack; PBS/Torque queue `cmt`, one 24-core node, four Julia threads by measured default.

---

## File Map

- Modify `dmrg/idmrg/Project.toml` and regenerate `dmrg/idmrg/Manifest.toml` on W003: exact backend pins and JLD2.
- Modify `dmrg/idmrg/jobs/bootstrap_w003.sh`: check both backend stacks.
- Create `dmrg/idmrg/src/MPSKitBackend.jl`: dependency constants, centered space, sector helpers, initial product states, and backend provenance.
- Create `dmrg/idmrg/src/MPSKitHamiltonian.jl`: term-list-to-finite/infinite Jordan MPO conversion and dense test helpers.
- Create `dmrg/idmrg/src/MPSKitSolver.jl`: IDMRG2 growth, VUMPS refinement, final Galerkin gate, records, and progress callbacks.
- Create `dmrg/idmrg/src/MPSKitCheckpoint.jl`: versioned atomic JLD2 state payload and compatibility checks.
- Create `dmrg/idmrg/src/MPSKitObservables.jl`: energy, density, Schmidt sectors, transfer spectrum, correlation length, mixed fidelity, and momentum translation data.
- Create `dmrg/idmrg/src/MPSKitFig2Adapter.jl`: `Fig2BenchmarkOperations` callbacks and backend-specific checkpoint audit.
- Modify `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`: qualified imports, includes, and exports only.
- Modify `dmrg/idmrg/src/Fig2Benchmark.jl`: backend-neutral checkpoint path/audit hooks and v5 provenance.
- Modify `dmrg/idmrg/benchmarks/fqahc_fig2.toml`: select `mpskit_idmrg_v1` and pin its gate policy.
- Modify `dmrg/idmrg/bin/run_fig2_benchmark.jl`: select operations from the manifest backend.
- Modify `dmrg/idmrg/jobs/run_tests.pbs`, `run_fig2_stage.pbs`, and `fig2_job_contract.sh`: four-thread default, backend hashes, gate wall times, and restart artifacts.
- Create `dmrg/idmrg/jobs/run_mpskit_gate.pbs`: small tensor, two-process restart, and D=32 solver gates.
- Create tests `test_mpskit_backend.jl`, `test_mpskit_hamiltonian.jl`, `test_mpskit_solver.jl`, `test_mpskit_checkpoint.jl`, `test_mpskit_observables.jl`, and `test_mpskit_fig2_adapter.jl`.
- Modify `dmrg/idmrg/test/runtests.jl`: include the new tests without removing legacy coverage.
- Modify `dmrg/idmrg/README.md`: backend distinction, parameter map, centered charge, W003 commands, and branch diagnosis.
- Create `dmrg/idmrg/report/fqahc_fig2_benchmark.md` and result plots only after audited production data exist.

## Task 1: Pin and Prove the W003 Environment

**Files:**
- Modify: `dmrg/idmrg/Project.toml`
- Modify on W003 through `Pkg`: `dmrg/idmrg/Manifest.toml`
- Modify: `dmrg/idmrg/jobs/bootstrap_w003.sh`
- Test: `dmrg/idmrg/test/test_mpskit_backend.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write the failing dependency/provenance test**

Create a test that imports dependencies with qualified names and asserts the source pins exposed by the adapter:

```julia
using Test
import MPSKit
import TensorKit
import TensorKitTensors
import JLD2

@testset "pinned MPSKit backend" begin
    @test MPSKIT_BACKEND_COMMIT == "811ecf6c06c1f7c1bc656da61abcd679effcd428"
    @test TENSORKITTENSORS_COMMIT == "3755705a1c44a3d5e32086e7d89b2c561b268cb1"
    p = mpskit_backend_provenance()
    @test p.backend == "mpskit_idmrg_v1"
    @test p.mpskit_version == string(Base.pkgversion(MPSKit))
    @test p.tensorkit_version == string(Base.pkgversion(TensorKit))
    @test p.tensorkittensors_version == string(Base.pkgversion(TensorKitTensors))
end
```

Include it at the end of `test/runtests.jl`.

- [ ] **Step 2: Push the RED test and run it on W003**

Commit only the test include and test file, push `DMRG`, fast-forward the clean production checkout, and run:

```bash
cd /home/public/shajy/codex/partial_Hall_crystal-idmrg-production
JULIA_DEPOT_PATH=/home/public/shajy/codex/depots/idmrg-julia-1.12.5 \
  /home/public/shajy/codex/runtime/julia-1.12.5/bin/julia \
  --startup-file=no --project=dmrg/idmrg \
  dmrg/idmrg/test/test_mpskit_backend.jl
```

Expected: fail because MPSKit/JLD2 are not direct dependencies and the adapter constants do not exist. Preserve the log under `/home/public/shajy/codex/results/fqahc-fig2/pbs/manual-red-<commit>/`.

- [ ] **Step 3: Add exact dependency sources and the minimal backend module**

Add direct dependencies for `MPSKit`, `TensorKit`, `TensorKitTensors`, `MatrixAlgebraKit`, and `JLD2`. Pin source commits:

```toml
[sources]
MPSKit = {url = "https://github.com/QuantumKitHub/MPSKit.jl.git", rev = "811ecf6c06c1f7c1bc656da61abcd679effcd428"}
TensorKitTensors = {url = "https://github.com/QuantumKitHub/TensorKitTensors.jl.git", rev = "3755705a1c44a3d5e32086e7d89b2c561b268cb1"}
```

Create `MPSKitBackend.jl` with the two commit constants and a named tuple returned by `mpskit_backend_provenance`. Add only qualified imports and include/export it from the package entry point.

- [ ] **Step 4: Resolve and regenerate only on W003**

On W003 run `Pkg.resolve()`, `Pkg.instantiate()`, and `Pkg.precompile()` in the isolated environment. Update `bootstrap_w003.sh --check` so its import expression includes the five new dependencies. Copy the generated manifest back to the worktree and verify its source revisions with `rg`.

- [ ] **Step 5: Run GREEN and the legacy suite on W003**

Run the targeted test, then submit `jobs/run_tests.pbs` with `JULIA_NUM_THREADS=4`. Expected: new dependency test passes and the complete legacy suite remains green. Commit project, manifest, bootstrap, module, and tests; push `DMRG`.

## Task 2: Uniform Centered Charge and Product Candidates

**Files:**
- Modify: `dmrg/idmrg/src/MPSKitBackend.jl`
- Test: `dmrg/idmrg/test/test_mpskit_backend.jl`

- [ ] **Step 1: Write failing centered-space arithmetic tests**

Test the public API:

```julia
@testset "uniform centered paper space" begin
    c = InfiniteCylinderConfig(; geometry=:paper_straight, Ny=6, x_period=3,
                               filling_num=7, filling_den=9)
    q = mpskit_centered_charges(c)
    @test q == (scale=18, empty=-7, occupied=11)
    pspaces = mpskit_physical_spaces(c)
    @test length(pspaces) == 36
    @test all(==(first(pspaces)), pspaces)
    for candidate in fig2_initial_candidates(c)
        psi = mpskit_product_state(c, candidate.occupied_sites)
        @test mpskit_product_occupations(psi, c) == candidate.occupied_sites
        @test mpskit_total_raw_charge(candidate.occupied_sites, c) == 0
        @test length(psi) == 36
    end
    @test_throws ArgumentError mpskit_product_state(c, 1:13)
end
```

- [ ] **Step 2: Run RED on W003**

Expected: undefined centered-space and product-state functions, not an import or syntax error.

- [ ] **Step 3: Implement the minimal centered representation**

Construct a product-sector `Vect` whose two one-dimensional sectors are `(parity=0, raw=-7)` and `(parity=1, raw=11)`. Build one-dimensional virtual spaces from the cumulative raw charge, with the site tensor convention `Vleft ⊗ P <- Vright`. Reject a path unless the last cumulative charge is zero. Normalize the returned MPS and expose occupation reconstruction only for bond-dimension-one test states.

- [ ] **Step 4: Run GREEN on W003 and commit**

Run the targeted test twice with the same candidate and assert identical tensor-space signatures. Run the configuration, candidate, and transverse-translation legacy tests. Commit and push.

## Task 3: Fermionic Jordan MPO Parity

**Files:**
- Create: `dmrg/idmrg/src/MPSKitHamiltonian.jl`
- Create: `dmrg/idmrg/test/test_mpskit_hamiltonian.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write the independent exact Fock-space oracle**

In the test file, implement test-only bit-basis actions for `c†_i c_j`:

```julia
function exact_hop!(H, target, source, amp, nsites)
    for ket in 0:(2^nsites - 1)
        ((ket >> (source - 1)) & 1) == 1 || continue
        ((ket >> (target - 1)) & 1) == 0 || continue
        between = (min(target, source) + 1):(max(target, source) - 1)
        parity = count(k -> ((ket >> (k - 1)) & 1) == 1, between)
        bra = xor(xor(ket, 1 << (source - 1)), 1 << (target - 1))
        H[bra + 1, ket + 1] += amp * (isodd(parity) ? -1 : 1)
    end
    return H
end
```

Use it to create elementwise expected matrices for two through six sites, including a `1<->4` complex hop, onsite number, and `n_2*n_5` interaction.

- [ ] **Step 2: Write the failing finite-MPO tests**

Specify:

```julia
H = mpskit_finite_hamiltonian(nsites, hops, interactions;
                              charge_mode=:uncentered)
@test mpskit_dense_matrix(H) ≈ H_exact atol=1e-12 rtol=0
@test H_exact ≈ H_exact' atol=1e-12
```

Repeat with `charge_mode=:centered` in a valid fixed-particle cell and compare every matrix element after desymmetrization. Test a nonadjacent hop specifically so a missing Jordan-Wigner sign fails.

- [ ] **Step 3: Run RED on W003**

Expected: missing adapter constructor. Preserve exact test seed and matrix mismatch location in the RED log.

- [ ] **Step 4: Implement ordered Hermitian pair conversion**

Group directed terms by the normalized unordered pair. For `i<j`, emit exactly one local operator

```julia
a * centered_f_plus_f_min(c) - conj(a) * centered_f_min_f_plus(c)
```

and emit `real(a)*centered_f_num(c)` for onsite terms. Emit interactions as `V*(n ⊗ n)`. Construct `MPSKit.FiniteMPOHamiltonian` and `MPSKit.InfiniteMPOHamiltonian` directly from `(indices => operator)` pairs. Desymmetrize only in the small finite test helper.

- [ ] **Step 5: Run GREEN and mutation checks on W003**

Prove the long-range test fails if the JW parity factor in the oracle is deliberately removed, then restore it. Run both centered and uncentered elementwise comparisons, Hermiticity, and the legacy ITensor Hamiltonian tests. Commit and push.

## Task 4: Full Paper Hamiltonian and Flux Gates

**Files:**
- Modify: `dmrg/idmrg/src/MPSKitHamiltonian.jl`
- Modify: `dmrg/idmrg/test/test_mpskit_hamiltonian.jl`

- [ ] **Step 1: Write failing term-list adapter tests**

For `Ny=2,x_period=1` and `Ny=2,x_period=3`, call the existing `build_infinite_model_terms` and compare the MPSKit finite-window matrix with the exact Fock oracle. For the 36-site production cell, avoid densification and assert:

```julia
hops, ints = build_infinite_model_terms(c, p)
H = mpskit_infinite_hamiltonian(c, p)
@test length(H) == sites_per_cell(c) == 36
@test mpskit_local_operators_are_hermitian(c, p; atol=1e-12)
@test mpskit_term_fingerprint(c, p) == canonical_term_fingerprint(hops, ints)
@test mpskit_terms_are_approx(with_flux(c, phi + 2pi),
                              with_flux(c, phi), p; atol=1e-12)
```

- [ ] **Step 2: Run RED on W003**

Expected: missing infinite adapter/fingerprint functions.

- [ ] **Step 3: Implement direct term consumption and fingerprinting**

Do not call the legacy ITensor `infinite_opsum`. Record site pairs and complex amplitudes before MPO construction, plus a SHA-256 over the canonical sorted representation. Expose MPO virtual dimensions for memory provenance.

- [ ] **Step 4: Run GREEN and commit**

Run new tests plus `test_paper_hamiltonian.jl`, `test_hamiltonian.jl`, and finite model regressions on W003. Commit and push.

## Task 5: IDMRG2 Growth and VUMPS Galerkin Gate

**Files:**
- Create: `dmrg/idmrg/src/MPSKitSolver.jl`
- Create: `dmrg/idmrg/test/test_mpskit_solver.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write a failing exactly-solvable solver test**

Use a two- or four-site gapped spinless-fermion chain with known finite-size/free-fermion energy density. Start from a product state and request `D=8`:

```julia
r = run_mpskit_idmrg(H, psi0;
    maxdim_schedule=[4,8], cutoff=1e-10,
    idmrg_maxiter=20, vumps_maxiter=50, galerkin_tol=1e-7)
@test maximum(r.link_dimensions) <= 8
@test r.galerkin_residual <= 1e-7
@test r.recomputed_galerkin_residual <= 1e-7
@test isapprox(r.energy_per_site, exact_energy; atol=1e-6)
@test r.converged
@test all(isfinite, getproperty.(r.records, :energy_per_site))
```

Also test that a low IDMRG center-tensor change with a high recomputed Galerkin residual is rejected.

- [ ] **Step 2: Run RED on W003**

Expected: missing solver API.

- [ ] **Step 3: Implement one schedule stage**

At each target dimension, run `MPSKit.IDMRG2` with `MatrixAlgebraKit.truncrank(target)` and record its returned diagnostic separately. Recreate environments, run `MPSKit.VUMPS`, then call `MPSKit.calc_galerkin` on the final state/environment. Record stage, iteration budget, requested/actual dimensions, energy, IDMRG diagnostic, returned VUMPS residual, recomputed residual, and elapsed seconds.

- [ ] **Step 4: Implement deterministic progress callbacks**

Write records atomically after every completed solver stage. A callback exception must preserve the last valid state and exit nonzero. The convergence predicate must use only the recomputed Galerkin residual plus energy/observable validity.

- [ ] **Step 5: Run GREEN and solver repeatability checks on W003**

Run the exact model from two independent product states. Both must reach the same energy within tolerance, and the second VUMPS call from the returned state must not increase the residual beyond tolerance. Commit and push.

## Task 6: Versioned Two-Process Checkpoint/Restart

**Files:**
- Create: `dmrg/idmrg/src/MPSKitCheckpoint.jl`
- Create: `dmrg/idmrg/test/test_mpskit_checkpoint.jl`
- Create: `dmrg/idmrg/test/mpskit_checkpoint_save.jl`
- Create: `dmrg/idmrg/test/mpskit_checkpoint_resume.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing same-process and negative tests**

Require `save_mpskit_checkpoint(path, state, config, metadata)` and
`load_mpskit_checkpoint(path, config)`. Test wrong config, wrong backend commit,
truncated input, and a legacy ITensor checkpoint. Each failure must be a typed
compatibility/format error.

- [ ] **Step 2: Write the failing two-process scripts**

Process A constructs and refines a nontrivial centered-charge state, writes
`state.h5`, `before.toml`, and exits. Process B loads it without sharing process
memory, recomputes energy/residual, performs one VUMPS refinement, writes
`state-resumed.h5` and `after.toml`. The harness compares configuration,
physical/virtual space fingerprints, energy, density, sector weights, and
checkpoint SHA-256.

- [ ] **Step 3: Run RED on W003**

Expected: missing checkpoint API. Preserve both process logs separately.

- [ ] **Step 4: Implement atomic JLD2 checkpointing**

Write to a sibling temporary file, close and reopen it, verify the completion
marker and metadata, then rename atomically. Store the MPSKit state and an
explicit space/config fingerprint. Never use Julia `Serialization` as the
production format.

- [ ] **Step 5: Run GREEN in two fresh W003 processes**

Run the scripts twice, once interactively and once through PBS. Confirm process
B performs a real solver iteration. Commit and push only after both runs pass.

## Task 7: Energy, Density, Entanglement, and Correlation Length

**Files:**
- Create: `dmrg/idmrg/src/MPSKitObservables.jl`
- Create: `dmrg/idmrg/test/test_mpskit_observables.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing product-state observable tests**

For a known 14-particle product candidate, assert every density, ring sum,
Schmidt probability, raw charge, physical charge `raw/18`, zero entropy, and
energy against direct term evaluation. Require all three physical x cuts.

- [ ] **Step 2: Write a failing transfer-spectrum test**

Use an entangled small state with an analytically or independently known second
transfer eigenvalue. Require two converged neutral eigenpairs, recorded
residuals, and `xi_x = x_period * (-1/log(abs(lambda2/lambda1)))`. A solver that
returns one eigenpair must produce `valid=false`.

- [ ] **Step 3: Run RED on W003**

Expected: missing observable API.

- [ ] **Step 4: Implement backend-neutral result records**

Use `MPSKit.expectation_value` for energy/local number operators; use the center
tensor `C` SVD for Schmidt data; extract the centered product-sector raw U(1)
label with a documented left-bond orientation; use KrylovKit on the one-cell
neutral transfer action and store eigenvalue/residual provenance.

- [ ] **Step 5: Run GREEN, compare with exact small models, and commit**

The MPSKit energy and densities must match the exact Fock-space oracle. Commit
and push.

## Task 8: Momentum ES and Mixed Fidelity

**Files:**
- Modify: `dmrg/idmrg/src/MPSKitObservables.jl`
- Modify: `dmrg/idmrg/test/test_mpskit_observables.jl`

- [ ] **Step 1: Write exact fermionic translation tests**

Reuse `fermionic_ring_translation` as the reference permutation. On small Fock
spaces, prove the MPSKit tensor action includes the wraparound fermion sign and
satisfies `T_y^Ny=I`. Test a deliberately omitted sign to show the assertion is
sensitive.

- [ ] **Step 2: Write failing Schmidt-block momentum tests**

Construct a state with known cyclic-orbit superpositions. Require integer
`ky_index`, phase, unitarity residual, sixth-power residual, commutator/closure
residual, and invalid status when the retained Schmidt subspace is not closed.
Do not accept labels inferred only by rounding phases.

- [ ] **Step 3: Write failing mixed-transfer fidelity tests**

Identical states must give one, a gauge-rotated equivalent state must give one,
and orthogonal product states must give zero within tolerance. Record mixed and
self-transfer residuals.

- [ ] **Step 4: Run RED, implement minimal tensor actions, and run GREEN on W003**

Keep backend-specific tensor contractions in `MPSKitObservables.jl`; share only
the exact site permutation/sign logic. Run all legacy momentum and branch tests
as regressions. Commit and push.

## Task 9: Fig. 2 Workflow Backend Adapter

**Files:**
- Create: `dmrg/idmrg/src/MPSKitFig2Adapter.jl`
- Create: `dmrg/idmrg/test/test_mpskit_fig2_adapter.jl`
- Modify: `dmrg/idmrg/src/Fig2Benchmark.jl`
- Modify: `dmrg/idmrg/bin/run_fig2_benchmark.jl`
- Modify: `dmrg/idmrg/benchmarks/fqahc_fig2.toml`
- Modify: `dmrg/idmrg/test/test_fig2_driver.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing manifest/backend selection tests**

Bump the manifest format to `fqahc_fig2_benchmark_v5` and require:

```toml
[backend]
id = "mpskit_idmrg_v1"
mpskit_commit = "811ecf6c06c1f7c1bc656da61abcd679effcd428"
tensorkittensors_commit = "3755705a1c44a3d5e32086e7d89b2c561b268cb1"
legacy_backend_role = "diagnostic_only"
```

Tests must reject an unknown backend, a source-pin mismatch, or a production
manifest selecting the legacy backend.

- [ ] **Step 2: Write failing operation-adapter tests**

Inject a tiny MPSKit model through `Fig2BenchmarkOperations`. Assert candidate
directory tables, `state.h5`, progress generations, ledger hashes, restart
audit, raw pump, sector-tracked table, and selection replay all use the existing
contract without reading ITensor internals.

- [ ] **Step 3: Run RED on W003**

Expected: v4-only manifest and hard-coded ITensor checkpoint paths/audits.

- [ ] **Step 4: Implement backend-neutral checkpoint hooks**

Move checkpoint filename, loader, max-link-dimension query, provenance, and
restart audit behind `Fig2BenchmarkOperations`. Retain `state.h5` as the
directory filename. Implement `mpskit_fig2_operations(spec)` using the new
solver, checkpoint, observables, momentum, and fidelity APIs.

- [ ] **Step 5: Run GREEN and the complete 900+ assertion driver suite**

Run targeted tests and all `test_fig2_*` files on W003. Confirm legacy injected
operation tests still pass. Commit and push.

## Task 10: PBS Gates and D=32 Paper Model

**Files:**
- Create: `dmrg/idmrg/jobs/run_mpskit_gate.pbs`
- Modify: `dmrg/idmrg/jobs/run_tests.pbs`
- Modify: `dmrg/idmrg/jobs/run_fig2_stage.pbs`
- Modify: `dmrg/idmrg/jobs/fig2_job_contract.sh`
- Modify: `dmrg/idmrg/jobs/submit_fig2_stage.sh`
- Test: `dmrg/idmrg/test/test_fig2_jobs.jl`

- [ ] **Step 1: Write failing PBS contract tests**

Require queue `cmt`, `nodes=1:ppn=24`, explicit wall time, default four Julia
threads, single-thread BLAS, immutable job config/manifest, clean production
checkout equal to `origin/DMRG`, and hashes for every new backend source file.

- [ ] **Step 2: Run RED then implement wrappers**

Add `gate=hamiltonian`, `gate=restart`, `gate=solver-small`, and
`gate=paper-d32` modes. Each writes `provenance.txt`, `job_status.toml`,
`time.txt`, logs, checkpoint hashes, and a gate TOML under a PBS-job-specific
result directory.

- [ ] **Step 3: Run GREEN and submit prerequisite gates**

Submit Hamiltonian, two-process restart, and small-solver gates. Monitor with
`qstat`; do not use fixed sleeps longer than 60 seconds. Audit exit status,
walltime, CPU utilization, maximum RSS, test counts, and result hashes.

- [ ] **Step 4: Submit the full 36-site D=32 gate**

Use the paper-straight parameters, two cold candidates, `D=32`, four Julia
threads, one 24-core `cmt` node, and a conservative wall time. Retain every
wavefunction and progress generation. Acceptance requires finite energy,
achieved dimension, valid observables, and a recomputed Galerkin residual below
the manifest gate without a two-cycle.

- [ ] **Step 5: Audit the architecture decision**

Compare the MPSKit residual trajectory and restart residual with the failed
ITensor event-113 source residual `1.144636e-3`; the MPSKit D=32 gate requires
the predeclared `1e-6` residual. If the MPSKit gate fails, stop
before D=64 and diagnose the exact failed layer. If it passes, record the D=32
gate SHA-256 values and commit the PBS wrapper/provenance references.

## Task 11: Bond-Dimension and Flux Production

**Files:**
- Production outputs under: `/home/public/shajy/codex/results/fqahc-fig2/`
- Modify only if evidence requires: `dmrg/idmrg/benchmarks/fqahc_fig2.toml`

- [ ] **Step 1: Run sparse pilot D=64 and D=128**

Use predecessor checkpoints, the seven predeclared flux points from `0` to
`6pi`, two cold candidates where required, and the existing ground/adiabatic
selection ledger. Audit energy, residual, checkpoint restart, density, transfer,
sector flow, momentum closure, memory, and walltime before continuing.

- [ ] **Step 2: Run D=256 and D=512**

Submit staged PBS jobs with explicit dependencies and immutable manifests.
Preserve nonselected candidates. Confirm convergence trends are not inferred
from energy alone.

- [ ] **Step 3: Run D=1000, D=2000, and D=3000 or the declared early plateau**

Use the manifest's predeclared early-plateau rule only if every lower-D
residual, pump, sector, momentum, restart, and provenance gate passes. Never
claim D=3000 evidence from a lower dimension.

- [ ] **Step 4: Refine flux only from recorded branch evidence**

Add points near fidelity, raw Schmidt-polarization, sector-weight, ES, or energy
crossing events. Keep the original raw grid values unchanged; do not unwrap or
shift any pump value.

- [ ] **Step 5: Freeze the production ledger**

Hash the manifest, ledger, selected checkpoints, tables, PBS provenance, and
all plotting inputs. Make immutable copies under the result root before report
generation.

## Task 12: Figures, Paper Comparison, Documentation, and Sync

**Files:**
- Create: `dmrg/idmrg/report/fqahc_fig2_benchmark.md`
- Create: `dmrg/idmrg/report/figures/fqahc_fig2_energy_convergence.pdf`
- Create: `dmrg/idmrg/report/figures/fqahc_fig2_pump.pdf`
- Create: `dmrg/idmrg/report/figures/fqahc_fig2_momentum_es.pdf`
- Create plotting source and tests under: `dmrg/idmrg/report/`
- Modify: `dmrg/idmrg/README.md`

- [ ] **Step 1: Write plotting tests before plotting code**

Test that plot inputs are selected only from the frozen ledger, raw pump values
are byte-for-byte unchanged, every point has a checkpoint/provenance hash, and
missing/invalid momentum labels cannot enter the counting figure.

- [ ] **Step 2: Generate plots from frozen data on W003**

Generate energy/residual versus D, raw and separately sector-tracked pump versus
`phi/2pi`, low-lying momentum ES with the predeclared cut, and branch diagnostic
panels. Copy only generated report artifacts back to the worktree.

- [ ] **Step 3: Write the evidence report**

Report exact parameters, backend pins, PBS resources, Hamiltonian/JW proof,
checkpoint proof, convergence table, raw endpoint pump, sector shift, momentum
counting, correlation lengths, discrepancies, and all SHA-256 references.
Compare against the digitized Fig. 2(d) curve and `1,1,2,3,5` target without
editing numerical data.

- [ ] **Step 4: Update README parameter mapping and branch diagnosis**

Document finite versus infinite geometry/filling/flux, the 36-site centered
charge rule, why finite edge charge is not the infinite pump, how to interpret
sector-origin jumps, and why ITensorInfiniteMPS results are diagnostic only for
this unit cell.

- [ ] **Step 5: Run completion verification on W003**

Run the full Julia suite, two-process restart gate, manifest/ledger verifier,
plot tests, and report link/hash verifier. Inspect PBS final states and selected
checkpoint hashes. Run finite-DMRG regressions on W003 as a non-regression gate.

- [ ] **Step 6: Protect user files and synchronize `DMRG`**

Before staging, assert that `.codex_diagnostics/`, `dmrg/plot_ground_density.py`,
`dmrg/test_plot_ground_density.py`, `tmp/`, and protected finite-DMRG scripts
are absent from the index. Stage explicit paths only, commit, push
`DMRG:DMRG`, verify local HEAD equals `origin/DMRG`, and leave the worktree
otherwise unchanged.
