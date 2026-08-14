# FQAHC Fig. 2 iDMRG Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce the paper-exact `Ny=6`, `Nx=3`, `nu=7/9` infinite-cylinder FQAHC benchmark in Fig. 2(d-f), with raw adiabatic pump data, validated charge/momentum entanglement spectra, resumable wavefunctions, and a convergence record through the required bond dimensions.

**Architecture:** Extend only `dmrg/idmrg` with a geometry-dispatched paper-straight cylinder, keeping the legacy sheared geometry and finite DMRG behavior unchanged. A benchmark driver owns deterministic product-state candidates, flux continuation, checkpoint ledgers, and acceptance tables; geometry, translation symmetry, observables, and persistence remain separate tested modules. All Julia tests and tensor calculations run in an isolated W003 clone through PBS, while the Mac is limited to edits and Git transport.

**Tech Stack:** Julia 1.12.5, ITensors 0.9.30, ITensorMPS 0.3.45, pinned ITensorInfiniteMPS commit `765f2777703bc1138b009adbed1b97bde1973402`, HDF5, KrylovKit, TOML, PBS Professional `cmt`, 24-core W003 nodes.

---

## File Map

- Create `dmrg/idmrg/jobs/bootstrap_w003.sh`: idempotent W003 Julia/depot/project bootstrap with exact-version checks.
- Create `dmrg/idmrg/jobs/run_tests.pbs`: short W003-only PBS harness for unit and restart tests.
- Modify `dmrg/idmrg/src/Configuration.jl`: tagged legacy/paper geometry, paper coordinates, exact filling arithmetic, and configuration signatures.
- Modify `dmrg/idmrg/src/Hamiltonian.jl`: dispatch legacy and straight-cylinder orbit/winding construction.
- Create `dmrg/idmrg/src/PaperGeometry.jl`: paper two-orbital coordinates, real-space hopping orbits, shell enumeration, and Bloch reconstruction.
- Create `dmrg/idmrg/src/InitialCandidates.jl`: deterministic period-three CDW and doped-particle patterns.
- Create `dmrg/idmrg/src/TransverseTranslation.jl`: exact fermionic ring translation and virtual symmetry validation.
- Create `dmrg/idmrg/src/MomentumEntanglement.jl`: charge- and momentum-resolved ES with residual-bearing labels.
- Modify `dmrg/idmrg/src/Observables.jl`: raw fixed-cut pump/sector observables and momentum-ES integration.
- Modify `dmrg/idmrg/src/BranchTracking.jl`: overlap-supported sector correspondence without changing raw data.
- Modify `dmrg/idmrg/src/Output.jl`: version-2 geometry-aware checkpoints, provenance, atomic candidate records, and scan ledger.
- Modify `dmrg/idmrg/src/CLI.jl`: unambiguous `--geometry`, `--Ny`, benchmark-manifest, and restart options.
- Modify `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`: include and export the new focused modules.
- Create `dmrg/idmrg/bin/run_fig2_benchmark.jl`: manifest-driven candidate/flux/convergence benchmark entry point.
- Create `dmrg/idmrg/benchmarks/fqahc_fig2.toml`: immutable target parameters, gates, grids, and acceptance tolerances.
- Create `dmrg/idmrg/jobs/run_fig2_stage.pbs`: 24-core `cmt` stage runner with timing and provenance.
- Create `dmrg/idmrg/jobs/submit_fig2_stage.sh`: collision-safe PBS submission and dependency helper.
- Create `dmrg/idmrg/test/test_paper_configuration.jl`: exact 36-site/18-cell/14-particle geometry tests.
- Create `dmrg/idmrg/test/test_paper_hamiltonian.jl`: Bloch, seam, shells, Hermiticity, periodicity, and parity tests.
- Create `dmrg/idmrg/test/test_initial_candidates.jl`: exact candidate-set tests.
- Create `dmrg/idmrg/test/test_transverse_translation.jl`: exact small-ring fermionic translation tests.
- Create `dmrg/idmrg/test/test_momentum_entanglement.jl`: virtual-operator validation and invalid-label tests.
- Create `dmrg/idmrg/test/test_checkpoint_restart.jl`: two-process restart regression driver.
- Create `dmrg/idmrg/test/checkpoint_restart_save.jl`: process-one optimize/save fixture.
- Create `dmrg/idmrg/test/checkpoint_restart_resume.jl`: process-two load/iterate/expand fixture.
- Create `dmrg/idmrg/test/test_fig2_driver.jl`: manifest, ledger, raw-data, and acceptance-report tests.
- Create `dmrg/idmrg/report/fqahc_fig2_benchmark.md`: generated human-readable evidence report.
- Modify `dmrg/idmrg/README.md`: finite/iDMRG/paper parameter mapping and branch-jump diagnosis.

The protected user files `dmrg/run_dmrg.jl`, `dmrg/jobs/submit_prl2014_tpl_lx15_ground.pbs`, and `dmrg/jobs/submit_prl2014_tpl_lx15_ground_refine.pbs` are outside this plan and must never be staged.

## Common W003 Test Contract

Every RED/GREEN command below runs from the dedicated W003 checkout stored in `W003_REPO=/home/public/shajy/codex/partial_Hall_crystal-idmrg-benchmark`. The user-owned runtime and depot are explicit:

```bash
export W003_REPO=/home/public/shajy/codex/partial_Hall_crystal-idmrg-benchmark
export JULIA_BIN=/home/public/shajy/codex/runtime/julia-1.12.5/bin/julia
export JULIA_DEPOT_PATH=/home/public/shajy/codex/depots/idmrg-julia-1.12.5
cd "$W003_REPO"
/opt/pbs/bin/qsub -v W003_REPO,JULIA_BIN,JULIA_DEPOT_PATH,TEST_TARGET=test/runtests.jl dmrg/idmrg/jobs/run_tests.pbs
```

`run_tests.pbs` writes a unique log directory under `/home/public/shajy/codex/results/fqahc-fig2/pbs/$PBS_JOBID`, records `hostname`, `git rev-parse HEAD`, `versioninfo()`, `Pkg.status(; mode=PKGMODE_MANIFEST)`, environment thread values, `/usr/bin/time -v`, and `qstat -f`. A test step passes only after `qstat` reports exit status zero and the captured Julia test summary contains no failures or errors.

### Task 1: Establish the isolated W003 runtime and PBS test harness

**Files:**
- Create: `dmrg/idmrg/jobs/bootstrap_w003.sh`
- Create: `dmrg/idmrg/jobs/run_tests.pbs`
- Modify: `dmrg/idmrg/README.md`

- [ ] **Step 1: Write the failing environment check**

Create `dmrg/idmrg/jobs/bootstrap_w003.sh` with a check-only mode whose required assertions are:

```bash
#!/usr/bin/env bash
set -euo pipefail

runtime_root=/home/public/shajy/codex/runtime
julia_bin="$runtime_root/julia-1.12.5/bin/julia"
depot_root=/home/public/shajy/codex/depots/idmrg-julia-1.12.5

if [[ ${1:-} == --check ]]; then
  [[ -x "$julia_bin" ]]
  [[ "$($julia_bin --startup-file=no -e 'print(VERSION)')" == 1.12.5 ]]
  JULIA_DEPOT_PATH="$depot_root" "$julia_bin" --startup-file=no \
    --project=dmrg/idmrg -e '
      using Pkg
      Pkg.instantiate()
      using ITensors, ITensorMPS, ITensorInfiniteMPS
      VERSION == v"1.12.5" || error("wrong Julia version")
    '
  exit 0
fi
```

- [ ] **Step 2: Run RED on W003**

Run:

```bash
ssh W003 'cd /home/public/shajy/codex/partial_Hall_crystal-idmrg-benchmark && bash dmrg/idmrg/jobs/bootstrap_w003.sh --check'
```

Expected: nonzero exit because `/home/public/shajy/codex/runtime/julia-1.12.5/bin/julia` does not yet exist.

- [ ] **Step 3: Implement idempotent bootstrap and the PBS harness**

Complete the non-check branch of `bootstrap_w003.sh` to download the official Linux x86-64 Julia 1.12.5 tarball to a private `mktemp -d` directory, verify the published SHA-256 stored in the script, extract it under `runtime_root`, and instantiate `dmrg/idmrg` with the dedicated depot. Create `run_tests.pbs` with:

```bash
#!/usr/bin/env bash
#PBS -N phc_idmrg_test
#PBS -q cmt
#PBS -l nodes=1:ppn=24
#PBS -l walltime=02:00:00
#PBS -j oe
set -euo pipefail
: "${W003_REPO:?}"
: "${JULIA_BIN:?}"
: "${JULIA_DEPOT_PATH:?}"
: "${TEST_TARGET:=test/runtests.jl}"
result_root=/home/public/shajy/codex/results/fqahc-fig2/pbs/${PBS_JOBID}
mkdir -p "$result_root"
cd "$W003_REPO/dmrg/idmrg"
export JULIA_NUM_THREADS=${NCPUS:-24}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
{
  hostname
  git -C "$W003_REPO" rev-parse HEAD
  "$JULIA_BIN" --startup-file=no -e 'versioninfo()'
  "$JULIA_BIN" --startup-file=no --project=. -e 'using Pkg; Pkg.status(; mode=Pkg.PKGMODE_MANIFEST)'
  env | sort | sed -n '/^JULIA_/p;/^OMP_/p;/^OPENBLAS_/p;/^MKL_/p;/^PBS_/p'
} > "$result_root/provenance.txt"
/usr/bin/time -v "$JULIA_BIN" --startup-file=no --project=. "$TEST_TARGET" \
  > "$result_root/test.log" 2> "$result_root/time.txt"
```

- [ ] **Step 4: Run GREEN on W003**

Run the bootstrap, submit `run_tests.pbs`, and poll with `qstat -f` until completion. Expected: Julia `1.12.5`, manifest instantiation succeeds, and the existing suite exits zero.

- [ ] **Step 5: Commit only the harness files**

```bash
git add dmrg/idmrg/jobs/bootstrap_w003.sh dmrg/idmrg/jobs/run_tests.pbs dmrg/idmrg/README.md
git commit -m "build: add isolated W003 iDMRG test harness"
```

### Task 2: Add unambiguous paper-straight configuration arithmetic

**Files:**
- Modify: `dmrg/idmrg/src/Configuration.jl`
- Modify: `dmrg/idmrg/src/CLI.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Create: `dmrg/idmrg/test/test_paper_configuration.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write the failing configuration tests**

Create `test/test_paper_configuration.jl`:

```julia
using Test
using InfiniteCylinderDMRG

@testset "paper-straight configuration" begin
  c = InfiniteCylinderConfig(; geometry=:paper_straight, Ny=6, x_period=3,
    filling_num=7, filling_den=9)
  @test c.geometry == :paper_straight
  @test circumference_unit_cells(c) == 6
  @test physical_circumference_sites(c) == 12
  @test sites_per_cell(c) == 36
  @test unit_cells_per_cell(c) == 18
  @test particles_per_cell(c) == 14
  @test physical_site_density(c) == 7//18
  @test paper_site(c, 0, 0, :A) == 1
  @test paper_site(c, 0, 0, :B) == 2
  @test paper_site(c, 1, 0, :A) == 13
  @test paper_coordinates(c, 14) == (1, 0, :B)
  @test configuration_signature(c) != configuration_signature(
    InfiniteCylinderConfig(; Ly=12, x_period=3, filling_num=7, filling_den=9))
  @test_throws ArgumentError InfiniteCylinderConfig(; geometry=:paper_straight,
    Ly=12, Ny=6, x_period=3, filling_num=7, filling_den=9)
end

@testset "legacy compatibility" begin
  old = InfiniteCylinderConfig(; Ly=6, x_period=1)
  @test old.geometry == :legacy_sheared
  @test sites_per_cell(old) == 6
  @test circumference_unit_cells(old) == 3
end
```

Include it from `runtests.jl`.

- [ ] **Step 2: Submit RED on W003**

Submit `run_tests.pbs` with `TEST_TARGET=test/test_paper_configuration.jl`. Expected: `MethodError` for the `geometry`/`Ny` constructor.

- [ ] **Step 3: Implement tagged geometry without changing legacy calls**

Introduce:

```julia
const SUPPORTED_GEOMETRIES = (:legacy_sheared, :paper_straight)

struct InfiniteCylinderConfig
  geometry::Symbol
  transverse_cells::Int
  x_period::Int
  filling_num::Int
  filling_den::Int
  phi_y::Float64
end

physical_circumference_sites(c) = c.geometry == :paper_straight ?
  2c.transverse_cells : c.transverse_cells
circumference_unit_cells(c) = c.geometry == :paper_straight ?
  c.transverse_cells : c.transverse_cells ÷ 2
sites_per_cell(c) = physical_circumference_sites(c) * c.x_period
paper_site(c, x, y, orbital) = begin
  c.geometry == :paper_straight || throw(ArgumentError("paper coordinates require paper_straight"))
  orbital_offset = orbital === :A ? 0 : orbital === :B ? 1 :
    throw(ArgumentError("orbital must be :A or :B"))
  x * physical_circumference_sites(c) + 2mod(y, c.transverse_cells) + orbital_offset + 1
end
```

Preserve the old positional constructor as `legacy_sheared`; reject `Ly` together with `Ny`; and include geometry plus circumference convention in the site-index/configuration signature.

- [ ] **Step 4: Submit GREEN and full regression on W003**

Expected: the new paper tests and every pre-existing test pass; existing CLI invocations with only `--Ly` parse identically.

- [ ] **Step 5: Commit**

```bash
git add dmrg/idmrg/src/Configuration.jl dmrg/idmrg/src/CLI.jl \
  dmrg/idmrg/src/InfiniteCylinderDMRG.jl dmrg/idmrg/test/runtests.jl \
  dmrg/idmrg/test/test_paper_configuration.jl
git commit -m "feat: add paper-straight cylinder configuration"
```

### Task 3: Construct and verify the paper Hamiltonian

**Files:**
- Create: `dmrg/idmrg/src/PaperGeometry.jl`
- Modify: `dmrg/idmrg/src/Hamiltonian.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Create: `dmrg/idmrg/test/test_paper_hamiltonian.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing geometry/Hamiltonian tests**

The new test fixes generic momenta before production results are inspected:

```julia
using Test, LinearAlgebra
using InfiniteCylinderDMRG

@testset "paper real-space Hamiltonian" begin
  p = CylinderModelParams(; t1=1, t3=0.2, V1=10, V2=2, V3=2)
  c0 = InfiniteCylinderConfig(; geometry=:paper_straight, Ny=6, x_period=3,
    filling_num=7, filling_den=9, phi_y=0)
  for k in ((0.173, 0.291), (1.137, -0.413), (-2.021, 0.777))
    @test paper_bloch_from_orbits(k, p) ≈ paper_reference_Hk(k, p) atol=1e-11 rtol=0
  end
  hops, ints = build_infinite_model_terms(c0, p)
  @test hamiltonian_is_hermitian(c0, hops)
  @test neighbor_shell_counts(c0, ints) == Dict(1=>108, 2=>108, 3=>108)
  @test straight_seam_is_connected(c0, hops)
  c2pi = with_flux(c0, 2pi)
  @test canonical_term_dict(build_infinite_model_terms(c0, p)...) ≈
    canonical_term_dict(build_infinite_model_terms(c2pi, p)...) atol=1e-11 rtol=0
  @test paper_window_parity(c0, p; x_cells=7, atol=1e-11)
end
```

- [ ] **Step 2: Submit RED on W003**

Expected: missing `paper_bloch_from_orbits` and paper geometry dispatch.

- [ ] **Step 3: Implement paper orbit ownership and winding**

Define focused data:

```julia
struct PaperHoppingOrbit
  dx::Int
  dy::Int
  target_orbital::Int
  source_orbital::Int
  amplitude::ComplexF64
end

function paper_winding_y(source_y::Int, dy::Int, Ny::Int)
  target_unwrapped = source_y + dy
  return fld(target_unwrapped, Ny), mod(target_unwrapped, Ny)
end
```

Generate the finite nonzero Fourier coefficients of the repository two-band `get_Hk` in the paper basis once, assert conjugate orbit closure, then map each orbit to `paper_site`. Apply `cis(winding*phi_y)` during mapping. Enumerate V1/V2/V3 pairs in unwrapped paper coordinates and deduplicate only after translation normalization. Keep the existing `build_infinite_hoppings` and `finite_term_parity` implementation unchanged for `legacy_sheared` dispatch.

- [ ] **Step 4: Submit GREEN and full regression on W003**

Expected: deterministic `H(k)` error below `1e-11`; all shell, seam, flux-periodicity, Hermiticity, and legacy tests pass.

- [ ] **Step 5: Commit**

```bash
git add dmrg/idmrg/src/PaperGeometry.jl dmrg/idmrg/src/Hamiltonian.jl \
  dmrg/idmrg/src/InfiniteCylinderDMRG.jl dmrg/idmrg/test/runtests.jl \
  dmrg/idmrg/test/test_paper_hamiltonian.jl
git commit -m "feat: build paper-straight infinite Hamiltonian"
```

### Task 4: Generate the period-three FQAHC candidate set

**Files:**
- Create: `dmrg/idmrg/src/InitialCandidates.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/src/CLI.jl`
- Create: `dmrg/idmrg/test/test_initial_candidates.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write the failing candidate tests**

```julia
using Test
using InfiniteCylinderDMRG

@testset "Fig. 2 initial candidates" begin
  c = InfiniteCylinderConfig(; geometry=:paper_straight, Ny=6, x_period=3,
    filling_num=7, filling_den=9)
  candidates = fig2_initial_candidates(c)
  @test length(candidates) >= 7
  @test length(unique(x.occupied_sites for x in candidates)) == length(candidates)
  @test all(length(x.occupied_sites) == 14 for x in candidates)
  @test all(length(unique(x.occupied_sites)) == 14 for x in candidates)
  @test Set(x.cdw_translation for x in candidates if x.kind == :cdw_doped) == Set(0:2)
  @test any(x.kind == :uniform_control for x in candidates)
  @test fig2_initial_candidates(c) == candidates
end
```

- [ ] **Step 2: Submit RED on W003**

Expected: `UndefVarError: fig2_initial_candidates not defined`.

- [ ] **Step 3: Implement bounded deterministic candidates**

Create:

```julia
struct InitialCandidate
  id::String
  kind::Symbol
  cdw_translation::Int
  dopant_cells::NTuple{2,Int}
  occupied_sites::Vector{Int}
end
```

Build each period-three background at `nu=2/3` with 12 particles in the 18 two-orbital cells, add two particles at symmetry-inequivalent empty-cell representatives, translate through all three CDW phases, deduplicate by the exact sorted occupied vector, and append the existing deterministic pattern as `uniform_control`. Stable IDs are SHA-free descriptive strings derived only from the pattern metadata.

- [ ] **Step 4: Submit GREEN and full regression on W003**

Expected: every candidate contains exactly 14 unique sites and repeated generation is byte-for-byte deterministic.

- [ ] **Step 5: Commit**

```bash
git add dmrg/idmrg/src/InitialCandidates.jl dmrg/idmrg/src/InfiniteCylinderDMRG.jl \
  dmrg/idmrg/src/CLI.jl dmrg/idmrg/test/runtests.jl \
  dmrg/idmrg/test/test_initial_candidates.jl
git commit -m "feat: add deterministic FQAHC product-state candidates"
```

### Task 5: Version checkpoints and pass the cross-process restart gate

**Files:**
- Modify: `dmrg/idmrg/src/Output.jl`
- Modify: `dmrg/idmrg/src/VUMPSRunner.jl`
- Modify: `dmrg/idmrg/src/CLI.jl`
- Modify: `dmrg/idmrg/test/test_checkpoint.jl`
- Create: `dmrg/idmrg/test/test_checkpoint_restart.jl`
- Create: `dmrg/idmrg/test/checkpoint_restart_save.jl`
- Create: `dmrg/idmrg/test/checkpoint_restart_resume.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Add the two-process failing regression**

`checkpoint_restart_save.jl` constructs a `paper_straight, Ny=2, Nx=3, nu=1/3` state, performs one no-expansion iteration, writes the canonical state plus site-index identity table and observables. `checkpoint_restart_resume.jl` starts in a fresh Julia process, loads it, validates identical site IDs/QNs/directions, performs one no-expansion iteration, expands from `D=4` to `D=8`, advances flux from `0` to `0.1`, and writes the post-restart observables. `test_checkpoint_restart.jl` runs both commands and asserts:

```julia
@test save_cmd.exitcode == 0
@test resume_cmd.exitcode == 0
@test resumed["pre_energy"] ≈ saved["energy"] atol=1e-10 rtol=0
@test resumed["pre_density"] ≈ saved["density"] atol=1e-10 rtol=0
@test resumed["pre_schmidt_qns"] == saved["schmidt_qns"]
@test resumed["post_maxlinkdim"] > resumed["pre_maxlinkdim"]
@test resumed["next_flux_iteration_completed"] === true
```

- [ ] **Step 2: Submit RED on W003**

Submit the restart test alone. Expected: reproduce the pinned backend QN link-direction mismatch during resume/expansion.

- [ ] **Step 3: Diagnose and implement a source-level restart repair**

Record every link index as `(id, tags, direction, qn blocks, dimensions)` before save and after load. Add `repair_loaded_canonical_links(psi)` only if the regression proves HDF5 reconstructs equivalent link spaces with inconsistent arrow orientation. The repair must replace a mismatched link by `sim(old; dir=expected_dir)` and apply `delta(old, dag(new))` to every adjacent tensor so physical tensors, QNs, and amplitudes are unchanged. Validate canonical residuals again after repair; never strip QNs or rebuild a product state.

Upgrade checkpoint metadata to `infinite_cylinder_vumps_v2`, including geometry, `Ny`/legacy circumference, coordinate convention, model-independent site signature, backend commit, Julia version, and package versions. The loader accepts v1 only as `legacy_sheared` after all v1 invariants pass. A paper state and a legacy state are mutually incompatible.

- [ ] **Step 4: Submit GREEN twice and full regression on W003**

Expected: two successive restart jobs pass; pre-restart observables match; no-expansion, expansion, and next-flux iterations complete; v1 legacy checkpoint load tests pass; paper/legacy mismatch tests fail cleanly with `CheckpointCompatibilityError`.

- [ ] **Step 5: Commit**

```bash
git add dmrg/idmrg/src/Output.jl dmrg/idmrg/src/VUMPSRunner.jl \
  dmrg/idmrg/src/CLI.jl dmrg/idmrg/test/test_checkpoint.jl \
  dmrg/idmrg/test/test_checkpoint_restart.jl \
  dmrg/idmrg/test/checkpoint_restart_save.jl \
  dmrg/idmrg/test/checkpoint_restart_resume.jl dmrg/idmrg/test/runtests.jl
git commit -m "fix: make geometry-aware iDMRG checkpoints resumable"
```

### Task 6: Implement exact fermionic transverse translation

**Files:**
- Create: `dmrg/idmrg/src/TransverseTranslation.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Create: `dmrg/idmrg/test/test_transverse_translation.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write exact small-Fock-space tests**

```julia
using Test, LinearAlgebra
using InfiniteCylinderDMRG

@testset "fermionic ring translation" begin
  for Ny in 2:4
    T = fermionic_ring_translation(Ny)
    @test T' * T ≈ I atol=1e-14 rtol=0
    @test T^Ny ≈ I atol=1e-14 rtol=0
    for bits in 0:(2^(2Ny)-1)
      target, sign = translated_fock_basis(bits, Ny)
      @test T[target+1, bits+1] == sign
      @test count(!iszero, T[:, bits+1]) == 1
    end
  end
end
```

The basis orders `(y=0,A),(y=0,B),...`; translating the final two modes to the front contributes the parity of the occupied wrapped modes permuted past the occupied retained modes.

- [ ] **Step 2: Submit RED on W003**

Expected: missing translation functions.

- [ ] **Step 3: Implement exact permutation signs and tensor action**

Implement `translated_fock_basis(bits, Ny)` by explicitly counting occupied inversions in the mode permutation, not by a bosonic reshape. Add a sparse operator constructor and a local-tensor action used by the twisted transfer solve. Reject legacy geometry because a one-cell straight transverse translation is not its symmetry generator.

- [ ] **Step 4: Submit GREEN on W003**

Expected: unitarity and `T^Ny=I` residuals are below `1e-14` for every basis state tested.

- [ ] **Step 5: Commit**

```bash
git add dmrg/idmrg/src/TransverseTranslation.jl \
  dmrg/idmrg/src/InfiniteCylinderDMRG.jl \
  dmrg/idmrg/test/test_transverse_translation.jl dmrg/idmrg/test/runtests.jl
git commit -m "feat: add fermionic transverse translation"
```

### Task 7: Resolve and validate momentum-entanglement labels

**Files:**
- Create: `dmrg/idmrg/src/MomentumEntanglement.jl`
- Modify: `dmrg/idmrg/src/Observables.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Create: `dmrg/idmrg/test/test_momentum_entanglement.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write valid and invalid virtual-symmetry tests**

```julia
using Test, LinearAlgebra
using InfiniteCylinderDMRG

@testset "momentum label validation" begin
  rho = Diagonal([0.4, 0.3, 0.2, 0.1])
  Ty = Diagonal(ComplexF64[1, cis(pi), 1, cis(pi)])
  result = validate_momentum_blocks(rho, Ty; Ny=2, residual_tol=1e-12)
  @test result.valid
  @test all(x -> x.ky_index in 0:1, result.levels)
  @test maximum(x.phase_residual for x in result.levels) < 1e-12

  bad_unitary = validate_momentum_blocks(rho, 1.01Ty; Ny=2, residual_tol=1e-12)
  @test !bad_unitary.valid
  @test all(x -> ismissing(x.ky_index), bad_unitary.levels)

  broken = [1 0.1; 0 1]
  broken_result = validate_momentum_blocks(Diagonal([0.6, 0.4]), broken;
    Ny=2, residual_tol=1e-12)
  @test !broken_result.valid
end
```

- [ ] **Step 2: Submit RED on W003**

Expected: missing virtual translation validation.

- [ ] **Step 3: Implement the twisted transfer and strict label gate**

Solve the mixed transfer problem between the state and its one-cell-translated image in each conserved charge block. Normalize the dominant fixed point into virtual `Ty`. Store these four residuals separately:

```julia
struct MomentumValidation
  unitary_residual::Float64
  power_residual::Float64
  commutator_residual::Float64
  transfer_residual::Float64
  valid::Bool
  reason::String
end
```

Diagonalize `Ty` only inside equal-charge and numerically degenerate Schmidt blocks. Quantize an eigenphase only if its distance to `2pi*m/Ny` is within the configured tolerance. If any global validation fails, return all levels with missing `ky_index`; do not report counting.

- [ ] **Step 4: Submit GREEN plus a Ny=2 tensor smoke on W003**

Expected: algebraic tests pass; the tensor smoke either returns fully validated labels or an explicit invalid reason, never fabricated indices.

- [ ] **Step 5: Commit**

```bash
git add dmrg/idmrg/src/MomentumEntanglement.jl dmrg/idmrg/src/Observables.jl \
  dmrg/idmrg/src/InfiniteCylinderDMRG.jl \
  dmrg/idmrg/test/test_momentum_entanglement.jl dmrg/idmrg/test/runtests.jl
git commit -m "feat: resolve validated transverse momentum entanglement"
```

### Task 8: Add the manifest-driven benchmark, atomic ledger, and acceptance report

**Files:**
- Create: `dmrg/idmrg/benchmarks/fqahc_fig2.toml`
- Create: `dmrg/idmrg/bin/run_fig2_benchmark.jl`
- Modify: `dmrg/idmrg/src/CLI.jl`
- Modify: `dmrg/idmrg/src/BranchTracking.jl`
- Modify: `dmrg/idmrg/src/Output.jl`
- Create: `dmrg/idmrg/test/test_fig2_driver.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`
- Modify: `dmrg/idmrg/README.md`

- [ ] **Step 1: Write failing manifest/ledger/raw-data tests**

The test parses the committed manifest and asserts exact target values:

```julia
@test m["geometry"] == "paper_straight"
@test m["Ny"] == 6
@test m["x_period"] == 3
@test (m["filling_num"], m["filling_den"]) == (7, 9)
@test (m["t1"], m["t3"], m["V1"], m["V2"], m["V3"]) == (1.0, 0.2, 10.0, 2.0, 2.0)
@test m["flux_start"] == 0.0
@test m["flux_stop"] == 6pi
@test m["pump_tolerance"] == 0.05
@test m["counting"] == [1, 1, 2, 3, 5]
```

A fake two-point candidate run asserts that `raw_schmidt_polarization` is byte-identical before/after report generation, every candidate has its own checkpoint, the selected ID is justified by energy at flux zero and mixed fidelity later, and an unvalidated momentum result forces acceptance failure.

- [ ] **Step 2: Submit RED on W003**

Expected: missing manifest and benchmark driver types.

- [ ] **Step 3: Implement the production contract**

For every `(D, flux, candidate)` write into a private directory: `state.h5`, `summary.toml`, convergence/expansion tables, density, raw charge ES, raw Schmidt sectors, momentum ES, mixed fidelities, and candidate metadata. Atomically update `ledger.toml` only after all files validate. Store `selected_candidate` separately from raw candidate rows. Derive a cumulative pump only when overlap-supported sector matching is valid; preserve the raw polarization column without modulo, offset, unwrap, sign flip, or branch translation.

The acceptance table has fixed rows `endpoint_pump`, `sector_shift`, `momentum_counting`, `paper_curve`, `bond_convergence`, `restart`, and `provenance`, each with measured value, tolerance, evidence path, and pass/fail.

- [ ] **Step 4: Submit GREEN and full suite on W003**

Expected: all driver tests pass, the full suite exits zero, and generated fake reports fail rather than hide missing momentum or convergence evidence.

- [ ] **Step 5: Commit**

```bash
git add dmrg/idmrg/benchmarks/fqahc_fig2.toml \
  dmrg/idmrg/bin/run_fig2_benchmark.jl dmrg/idmrg/src/CLI.jl \
  dmrg/idmrg/src/BranchTracking.jl dmrg/idmrg/src/Output.jl \
  dmrg/idmrg/test/test_fig2_driver.jl dmrg/idmrg/test/runtests.jl \
  dmrg/idmrg/README.md
git commit -m "feat: add auditable Fig. 2 benchmark workflow"
```

### Task 9: Run the W003 geometry/restart gates and sparse-flux pilot

**Files:**
- Create: `dmrg/idmrg/jobs/run_fig2_stage.pbs`
- Create: `dmrg/idmrg/jobs/submit_fig2_stage.sh`
- Modify: `dmrg/idmrg/README.md`
- Output only: `/home/public/shajy/codex/results/fqahc-fig2/*`

- [ ] **Step 1: Add and test the production PBS wrapper**

`run_fig2_stage.pbs` must request:

```bash
#PBS -q cmt
#PBS -l nodes=1:ppn=24
#PBS -l walltime=12:00:00
export JULIA_NUM_THREADS=24
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

It rejects output roots already locked by another job, records PBS/node/code/runtime/package provenance, enables ITensor block-sparse threading explicitly, and invokes `run_fig2_benchmark.jl` with an immutable copied manifest and stage name.

- [ ] **Step 2: Submit geometry and restart gates**

Submit the full test suite and two independent restart regressions. Expected: all have PBS exit status zero before any flux job is accepted.

- [ ] **Step 3: Submit D=32/64/128 sparse pilot**

Use the inclusive flux grid `phi_y/(2pi) = 0,0.5,1,1.5,2,2.5,3`, all zero-flux candidates, and warm plus bounded cold candidates later. Walltime is at most 12 hours. Save all candidate wavefunctions and raw observables. Do not treat low-D nonquantization as a failure of topology; do treat invalid state, restart, or branch evidence as a gate failure.

- [ ] **Step 4: Inspect evidence and run the 24-core utilization comparison**

Compare identical D=32 iterations at 4, 12, and 24 Julia threads with BLAS/OMP fixed to one. Record elapsed time, CPU utilization, and peak RSS. Use 24 cores for later stages only if measured throughput is competitive; otherwise document the measured thread count while still keeping the production allocation collision-free.

- [ ] **Step 5: Commit wrappers and pilot provenance links**

```bash
git add dmrg/idmrg/jobs/run_fig2_stage.pbs dmrg/idmrg/jobs/submit_fig2_stage.sh \
  dmrg/idmrg/README.md
git commit -m "build: add staged W003 Fig. 2 jobs"
```

### Task 10: Escalate convergence and produce the paper comparison

**Files:**
- Create: `dmrg/idmrg/report/fqahc_fig2_benchmark.md`
- Modify: `dmrg/idmrg/README.md`
- Output only: `/home/public/shajy/codex/results/fqahc-fig2/*`

- [ ] **Step 1: Declare convergence and counting windows before inspecting high-D output**

Commit the manifest thresholds: per-site energy change `<=1e-6` between final stages, endpoint pump change `<=0.02`, raw/derived sector correspondence residual `<=1e-6`, momentum validation residual `<=1e-6`, and the fixed low-lying ES counting rule (charge sector, momentum sectors, and entanglement-energy window) digitized from Fig. 2(f). Store the Fig. 2(d) digitization with an uncertainty column and immutable source-image checksum.

- [ ] **Step 2: Submit D=256 then D=512 stages**

Use walltimes no larger than 36 and 72 hours respectively. Inspect only after each job and all atomic candidate checkpoints are complete. If a stage fails, apply the design debugging order with one hypothesis/one changed variable and add a failing regression before code changes.

- [ ] **Step 3: Submit D=1000, D=2000, and D=3000 progression**

Use up to 72 hours for D=1000 and 120 hours for D=2000/3000, segmenting with `afterok` only because Task 5 passed. Never allow segments to share a writable candidate directory. A documented early plateau may make D=3000 redundant only if the predeclared energy, pump, sector, momentum-counting, and curve tolerances all pass at both D=1000 and D=2000.

- [ ] **Step 4: Refine flux only where branch evidence requires it**

Add points around fidelity/sector-flow changes while retaining the seven sparse points. Candidate selection uses mixed-MPS fidelity and convergence eligibility, never manual curve alignment. Preserve every raw candidate and produce both `pump_raw.tsv` and `pump_sector_tracked.tsv`.

- [ ] **Step 5: Generate and verify the final acceptance report**

The report must state measured `Delta Q(6pi)`, the mapped one-sector ES shift, actual validated momentum counting, digitized-curve residuals, energy/bond convergence, every PBS job ID, and links to raw evidence. It must mark unmet criteria as failures. Run the complete W003 suite again at the final code commit and verify all referenced checkpoint and table checksums.

- [ ] **Step 6: Commit documentation and synchronize DMRG**

```bash
git add dmrg/idmrg/report/fqahc_fig2_benchmark.md dmrg/idmrg/README.md \
  dmrg/idmrg/benchmarks/fqahc_fig2.toml
git commit -m "docs: report FQAHC Fig. 2 iDMRG benchmark"
git status --short
git push origin DMRG
```

Expected final status: only the three pre-existing protected user paths remain modified/untracked; `origin/DMRG` matches local `DMRG`; every acceptance claim is backed by W003 output at the reported commit.

## Self-Review

- Spec coverage: Tasks 2-3 cover exact geometry, filling, Hamiltonian, twist, and legacy isolation; Tasks 4 and 8 cover candidates and raw branch selection; Tasks 5 and 9 gate checkpoint segmentation; Tasks 6-7 cover true fermionic transverse momentum; Tasks 8-10 cover provenance, staged PBS computation, convergence, digitized comparison, and acceptance.
- Placeholder scan: the plan contains no deferred implementation markers; each test, command, file boundary, tolerance, and failure expectation is explicit.
- Type consistency: `InfiniteCylinderConfig`, `InitialCandidate`, `MomentumValidation`, `paper_site`, `fig2_initial_candidates`, and the geometry dispatch names are used consistently from definition through tests and outputs.
- Safety: every Julia command runs on W003/PBS; W003's dirty production tree and the Mac's compute runtime are excluded; the three protected finite-DMRG paths are never staged.
