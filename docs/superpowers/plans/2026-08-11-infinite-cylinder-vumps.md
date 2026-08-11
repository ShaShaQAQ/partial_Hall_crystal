# Infinite-Cylinder VUMPS Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible VUMPS-based infinite-cylinder workflow that reuses the finite DMRG Hamiltonian, preserves raw flux branches, and outputs energy density, neutral correlation length, charge-resolved entanglement spectra, and Schmidt-sector diagnostics.

**Architecture:** Extract a tensor-network-independent cylinder model core, leave finite DMRG behind its existing API, and add an isolated `dmrg/idmrg` Julia environment. Translation-normalized Hamiltonian terms feed `InfiniteSum{MPO}`; a pinned backend adapter owns VUMPS iteration internals, while separate modules own observables, checkpoints, and branch tracking.

**Tech Stack:** Julia 1.12, ITensors 0.9.30, ITensorMPS 0.3.45 in the isolated environment, pinned ITensorInfiniteMPS 0.2.2 commit `765f2777703bc1138b009adbed1b97bde1973402`, KrylovKit, HDF5, Test, TOML.

---

## File Map

- Create `dmrg/CylinderModelCore.jl`: shared parameter type, neighbor displacements, reciprocal vectors, and two-band hopping transform.
- Modify `dmrg/DMRGFullModel.jl`: import and re-export the shared definitions without changing finite public behavior.
- Create `dmrg/idmrg/Project.toml` and `dmrg/idmrg/Manifest.toml`: isolated pinned backend environment.
- Create `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`: package entry point and exports.
- Create `dmrg/idmrg/src/Configuration.jl`: unit-cell arithmetic, site mapping, centered-charge scale, and initial patterns.
- Create `dmrg/idmrg/src/Hamiltonian.jl`: infinite hopping/interaction orbits, Hermiticity, tiling, parity, and `InfiniteSum{MPO}` construction.
- Create `dmrg/idmrg/src/VUMPSRunner.jl`: staged subspace expansion, per-iteration residual records, convergence gates, and warm starts.
- Create `dmrg/idmrg/src/Observables.jl`: energy normalization, density, entanglement data, Schmidt sectors, transfer spectrum, and neutral correlation length.
- Create `dmrg/idmrg/src/BranchTracking.jl`: mixed-transfer fidelity, sector distances, branch metrics, and candidate selection.
- Create `dmrg/idmrg/src/Output.jl`: TOML/TSV output and HDF5 checkpoints.
- Create `dmrg/idmrg/bin/run_vumps.jl`: single-point CLI.
- Create `dmrg/idmrg/bin/run_flux_scan.jl`: sequential raw flux scan CLI.
- Create `dmrg/idmrg/test/*.jl`: fast, backend, checkpoint, and smoke tests.
- Create `dmrg/idmrg/jobs/run_smoke.pbs`: cluster smoke entry point.
- Create `dmrg/idmrg/README.md`: dependency, parameter, observable, branch, and command documentation.

### Task 1: Extract the shared cylinder model core without changing finite DMRG

**Files:**
- Create: `dmrg/CylinderModelCore.jl`
- Modify: `dmrg/DMRGFullModel.jl`
- Test: `dmrg/test_convergence_logic.jl`
- Test: `dmrg/test_flux_merge.jl`
- Test: `dmrg/test_flux_segment_checkpoint.jl`

- [ ] **Step 1: Record the finite regression baseline**

Run:

```bash
julia --project=. dmrg/test_convergence_logic.jl
julia --project=. dmrg/test_flux_merge.jl
julia --project=. dmrg/test_flux_segment_checkpoint.jl
```

Expected: 4, 6, and 9 passing assertions respectively.

- [ ] **Step 2: Create the shared model module**

Create `dmrg/CylinderModelCore.jl` with the exact public surface:

```julia
module CylinderModelCore

using LinearAlgebra

include("../shared/lattice.jl")
include("../shared/hoppings.jl")

export CylinderModelParams,
  NB_DISPS,
  reciprocal_vectors,
  hopping_matrices_from_get_Hk

struct CylinderModelParams
  t1::Float64
  t3::Float64
  V1::Float64
  V2::Float64
  V3::Float64
end

CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.0, V3=0.0) =
  CylinderModelParams(Float64(t1), Float64(t3), Float64(V1), Float64(V2), Float64(V3))

const NB_DISPS = Dict(
  1 => [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, 1), (1, -1)],
  2 => [(1, 1), (-1, -1), (-1, 2), (1, -2), (-2, 1), (2, -1)],
  3 => [(2, 0), (-2, 0), (0, 2), (0, -2), (-2, 2), (2, -2)],
)

function reciprocal_vectors(a1::NTuple{2,Float64}, a2::NTuple{2,Float64})
  area = a1[1] * a2[2] - a1[2] * a2[1]
  return (
    (2π * a2[2] / area, -2π * a2[1] / area),
    (-2π * a1[2] / area, 2π * a1[1] / area),
  )
end

function hopping_matrices_from_get_Hk(
  t1::Float64,
  t3::Float64;
  nx_aux::Int=16,
  ny_aux::Int=16,
  max_dx::Int=4,
  max_dy::Int=4,
  tol::Float64=1e-10,
)
  a1 = [1.0, 0.0]
  a2 = [0.5, sqrt(3) / 2]
  b1, b2 = reciprocal_vectors(Tuple(a1), Tuple(a2))
  sublat = [[0.0, 0.0], copy(a2)]
  out = Dict{NTuple{2,Int},Matrix{ComplexF64}}()
  Nk = nx_aux * ny_aux
  for n1 in -max_dx:max_dx, n2 in -max_dy:max_dy
    R = n1 .* a1 .+ n2 .* (2 .* a2)
    tmat = zeros(ComplexF64, 2, 2)
    for m1 in 0:(nx_aux - 1), m2 in 0:(ny_aux - 1)
      k = (m1 / nx_aux) .* collect(b1) .+ (m2 / (2 * ny_aux)) .* collect(b2)
      Hk = get_Hk(k, t1, t3)
      for alpha in 1:2, beta in 1:2
        delta = sublat[alpha] .- sublat[beta]
        tmat[alpha, beta] += exp(-1im * dot(k, R .+ delta)) * Hk[alpha, beta]
      end
    end
    tmat ./= Nk
    for alpha in 1:2, beta in 1:2
      abs(real(tmat[alpha, beta])) < tol && (tmat[alpha, beta] -= real(tmat[alpha, beta]))
      abs(imag(tmat[alpha, beta])) < tol && (tmat[alpha, beta] -= 1im * imag(tmat[alpha, beta]))
    end
    maximum(abs, tmat) > tol && (out[(n1, n2)] = tmat)
  end
  return out
end

end
```

- [ ] **Step 3: Replace duplicated finite definitions with shared imports**

At the top of `dmrg/DMRGFullModel.jl`, replace the shared lattice/hopping includes and duplicated `CylinderModelParams`, `NB_DISPS`, `reciprocal_vectors`, and `hopping_matrices_from_get_Hk` definitions with:

```julia
include("CylinderModelCore.jl")
using .CylinderModelCore: CylinderModelParams,
  NB_DISPS,
  reciprocal_vectors,
  hopping_matrices_from_get_Hk

include("../shared/basis.jl")
```

Keep every finite export and call signature unchanged.

- [ ] **Step 4: Run finite regression tests and model-range check**

Run:

```bash
julia --project=. dmrg/test_convergence_logic.jl
julia --project=. dmrg/test_flux_merge.jl
julia --project=. dmrg/test_flux_segment_checkpoint.jl
julia --project=. -e 'include("dmrg/DMRGFullModel.jl"); using .DMRGFullModel; tR=DMRGFullModel.hopping_matrices_from_get_Hk(1.0,0.2); @assert length(tR)==11; @assert maximum(abs(first(k)) for k in keys(tR))==2'
```

Expected: all finite tests pass; the hopping check exits zero.

- [ ] **Step 5: Commit the extraction**

```bash
git add dmrg/CylinderModelCore.jl dmrg/DMRGFullModel.jl
git commit -m "refactor: extract shared cylinder model core"
```

### Task 2: Scaffold the isolated package and validated unit-cell configuration

**Files:**
- Create: `dmrg/idmrg/Project.toml`
- Create: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Create: `dmrg/idmrg/src/Configuration.jl`
- Create: `dmrg/idmrg/test/runtests.jl`
- Create: `dmrg/idmrg/test/test_configuration.jl`

- [ ] **Step 1: Write failing configuration tests**

Create `dmrg/idmrg/test/test_configuration.jl`:

```julia
using Test
using InfiniteCylinderDMRG

@testset "infinite cylinder configuration" begin
  c1 = InfiniteCylinderConfig(; Ly=6, x_period=1, filling_num=1, filling_den=3)
  @test sites_per_cell(c1) == 6
  @test unit_cells_per_cell(c1) == 3
  @test particles_per_cell(c1) == 1
  @test charge_scale(c1) == 6

  c3 = InfiniteCylinderConfig(; Ly=6, x_period=3, filling_num=1, filling_den=3)
  @test sites_per_cell(c3) == 18
  @test particles_per_cell(c3) == 3
  @test charge_scale(c3) == 6

  @test_throws ArgumentError InfiniteCylinderConfig(; Ly=2, x_period=1)
  csmall = InfiniteCylinderConfig(; Ly=2, x_period=3)
  @test particles_per_cell(csmall) == 1

  @test linear_site(c1, 0, 0) == 1
  @test linear_site(c1, 1, 0) == 7
  @test site_coordinates(c1, 7) == (1, 0)
  @test canonical_site(c1, 7) == (1, 1)
end
```

Create `dmrg/idmrg/test/runtests.jl`:

```julia
using Test
include("test_configuration.jl")
```

- [ ] **Step 2: Run the test and observe the missing package failure**

Run:

```bash
julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl
```

Expected: FAIL because `InfiniteCylinderDMRG` and its configuration types do not exist.

- [ ] **Step 3: Add the isolated project metadata**

Create `dmrg/idmrg/Project.toml`:

```toml
name = "InfiniteCylinderDMRG"
uuid = "bb119d70-e6e4-4ad4-a0fd-2dc1be8dca26"
version = "0.1.0"

[deps]
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
ITensorInfiniteMPS = "1dc1fb26-a137-4954-ae60-1bd4106e95ad"
ITensorMPS = "0d1a4710-d33b-49a5-8f18-73bdf49b47e2"
ITensors = "9136182c-28ba-11e9-034c-db9fb085ebd5"
KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"

[compat]
HDF5 = "0.17"
ITensorMPS = "=0.3.45"
ITensors = "=0.9.30"
KrylovKit = "0.10"
julia = "1.10, 1.11, 1.12"

[sources]
ITensorInfiniteMPS = {url = "https://github.com/ITensor/ITensorInfiniteMPS.jl.git", rev = "765f2777703bc1138b009adbed1b97bde1973402"}

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test"]
```

- [ ] **Step 4: Implement configuration arithmetic and mapping**

Create `dmrg/idmrg/src/Configuration.jl` with:

```julia
struct InfiniteCylinderConfig
  Ly::Int
  x_period::Int
  filling_num::Int
  filling_den::Int
  phi_y::Float64
end

function InfiniteCylinderConfig(;
  Ly::Int=6,
  x_period::Int=1,
  filling_num::Int=1,
  filling_den::Int=3,
  phi_y::Real=0.0,
)
  Ly > 0 && iseven(Ly) || throw(ArgumentError("Ly must be positive and even"))
  x_period > 0 || throw(ArgumentError("x_period must be positive"))
  filling_num > 0 || throw(ArgumentError("filling_num must be positive"))
  filling_den > 0 || throw(ArgumentError("filling_den must be positive"))
  c = InfiniteCylinderConfig(Ly, x_period, filling_num, filling_den, Float64(phi_y))
  numerator = filling_num * unit_cells_per_cell(c)
  numerator % filling_den == 0 || throw(ArgumentError("unit cell has noninteger particle number"))
  return c
end

sites_per_cell(c::InfiniteCylinderConfig) = c.Ly * c.x_period
unit_cells_per_cell(c::InfiniteCylinderConfig) = sites_per_cell(c) ÷ 2
particles_per_cell(c::InfiniteCylinderConfig) =
  c.filling_num * unit_cells_per_cell(c) ÷ c.filling_den
physical_site_density(c::InfiniteCylinderConfig) = c.filling_num // (2 * c.filling_den)
charge_scale(c::InfiniteCylinderConfig) =
  sites_per_cell(c) ÷ gcd(sites_per_cell(c), particles_per_cell(c))

linear_site(c::InfiniteCylinderConfig, x::Int, y::Int) = x * c.Ly + mod(y, c.Ly) + 1

function site_coordinates(c::InfiniteCylinderConfig, n::Int)
  x, y = fldmod(n - 1, c.Ly)
  return x, y
end

function canonical_site(c::InfiniteCylinderConfig, n::Int)
  cell, offset = fldmod(n - 1, sites_per_cell(c))
  return offset + 1, cell
end

function default_occupied_sites(c::InfiniteCylinderConfig)
  np = particles_per_cell(c)
  return unique(round.(Int, range(1, sites_per_cell(c); length=np + 2))[2:(end - 1)])
end
```

Create `dmrg/idmrg/src/InfiniteCylinderDMRG.jl` to include the shared core and configuration, and export every function used by the test.

- [ ] **Step 5: Resolve the pinned environment and run configuration tests**

Run:

```bash
julia --project=dmrg/idmrg -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl
```

Expected: the manifest pins ITensorInfiniteMPS commit `765f277...`, ITensorMPS 0.3.45, and ITensors 0.9.30; configuration tests pass.

- [ ] **Step 6: Commit package configuration**

```bash
git add dmrg/idmrg/Project.toml dmrg/idmrg/Manifest.toml dmrg/idmrg/src/InfiniteCylinderDMRG.jl dmrg/idmrg/src/Configuration.jl dmrg/idmrg/test
git commit -m "feat: add infinite cylinder package configuration"
```

### Task 3: Build translation-normalized infinite Hamiltonian terms with finite parity

**Files:**
- Create: `dmrg/idmrg/src/Hamiltonian.jl`
- Create: `dmrg/idmrg/test/test_hamiltonian.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing canonical-term and parity tests**

Create tests covering support, Hermiticity, periodicity, and finite tiling:

```julia
using Test
using InfiniteCylinderDMRG

@testset "infinite Hamiltonian terms" begin
  params = CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.4, V3=0.3)
  for x_period in (1, 3), phi in (0.0, 0.37π, 2π)
    cfg = InfiniteCylinderConfig(; Ly=6, x_period, phi_y=phi)
    hops, ints = build_infinite_model_terms(cfg, params)
    @test !isempty(hops)
    @test !isempty(ints)
    @test all(min(h.target, h.source) in 1:sites_per_cell(cfg) for h in hops)
    @test hamiltonian_is_hermitian(cfg, hops)
  end

  c0 = InfiniteCylinderConfig(; Ly=6, x_period=1, phi_y=0.0)
  c2 = InfiniteCylinderConfig(; Ly=6, x_period=1, phi_y=2π)
  @test canonical_term_dict(build_infinite_model_terms(c0, params)...) ≈
        canonical_term_dict(build_infinite_model_terms(c2, params)...)

  for phi in (0.0, 0.37π, 2π)
    cfg = InfiniteCylinderConfig(; Ly=6, x_period=3, phi_y=phi)
    @test finite_term_parity(cfg, params; Lx=7, atol=1e-10)
  end
end
```

- [ ] **Step 2: Run the Hamiltonian test and verify it fails on missing functions**

Run:

```bash
julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl
```

Expected: FAIL on `build_infinite_model_terms`.

- [ ] **Step 3: Implement canonical term types and translation normalization**

Create `dmrg/idmrg/src/Hamiltonian.jl` with these types and normalization helpers:

```julia
struct InfiniteHoppingTerm
  target::Int
  source::Int
  amp::ComplexF64
end

struct InfiniteInteractionTerm
  i::Int
  j::Int
  V::Float64
  shell::Int
end

function normalize_directed(c::InfiniteCylinderConfig, target::Int, source::Int)
  n = sites_per_cell(c)
  shift = fld(min(target, source) - 1, n) * n
  return target - shift, source - shift
end

function normalize_undirected(c::InfiniteCylinderConfig, i::Int, j::Int)
  i, j = normalize_directed(c, i, j)
  return minmax(i, j)
end
```

Implement hopping construction by enumerating source sites in the reference
MPS cell, applying the shared (t(R)), periodicizing `y`, multiplying
`cis(winding_y * phi_y)`, normalizing by whole MPS cells, and accumulating
translation-equivalent directed pairs.

Implement interaction construction by enumerating `NB_DISPS`, normalizing
unordered pairs, and accumulating one pair per translation orbit and shell.

- [ ] **Step 4: Implement Hermiticity and finite-window tiling**

Use normalized reverse keys for Hermiticity:

```julia
function hamiltonian_is_hermitian(c, hops; atol=1e-10)
  d = Dict((h.target, h.source) => h.amp for h in hops)
  for h in hops
    reverse_key = normalize_directed(c, h.source, h.target)
    haskey(d, reverse_key) || return false
    isapprox(d[reverse_key], conj(h.amp); atol, rtol=0) || return false
  end
  return true
end
```

Implement `tile_to_finite_window` by translating every canonical term by
integer multiples of `sites_per_cell(cfg)`, retaining exactly the terms whose
sites both lie in `1:(Lx*Ly)`. `finite_term_parity` constructs the existing
finite `CylinderLat`, compares accumulated hopping dictionaries, and compares
interaction dictionaries keyed by pair and shell.

- [ ] **Step 5: Build the infinite OpSum adapter**

Add:

```julia
function infinite_opsum(c::InfiniteCylinderConfig, params::CylinderModelParams)
  hops, ints = build_infinite_model_terms(c, params)
  os = OpSum()
  for h in hops
    if h.target == h.source
      os += h.amp, "N", h.source
    else
      os += h.amp, "Cdag", h.target, "C", h.source
    end
  end
  for v in ints
    os += v.V, "N", v.i, "N", v.j
  end
  starts = Set(minimum(ITensors.sites(t)) for t in ITensors.terms(os))
  for j in 1:sites_per_cell(c)
    j in starts || (os += 0.0, "Id", j)
  end
  return os
end

build_infinite_mpo(c, params, sites) = InfiniteSum{MPO}(infinite_opsum(c, params), sites)
```

- [ ] **Step 6: Run parity and finite regressions**

Run:

```bash
julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl
julia --project=. dmrg/test_convergence_logic.jl
julia --project=. dmrg/test_flux_merge.jl
julia --project=. dmrg/test_flux_segment_checkpoint.jl
```

Expected: all Hamiltonian, configuration, and finite regression tests pass.

- [ ] **Step 7: Commit Hamiltonian construction**

```bash
git add dmrg/idmrg/src/Hamiltonian.jl dmrg/idmrg/src/InfiniteCylinderDMRG.jl dmrg/idmrg/test
git commit -m "feat: build infinite cylinder Hamiltonian"
```

### Task 4: Add QN product states and a residual-reporting VUMPS runner

**Files:**
- Create: `dmrg/idmrg/src/VUMPSRunner.jl`
- Create: `dmrg/idmrg/test/test_vumps_runner.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing state and convergence-helper tests**

```julia
using Test
using ITensors
using ITensorMPS
using InfiniteCylinderDMRG

@testset "VUMPS runner helpers" begin
  cfg = InfiniteCylinderConfig(; Ly=2, x_period=3)
  sites, initstate, psi = initial_infinite_mps(cfg)
  @test length(sites) == sites_per_cell(cfg)
  @test flux(psi.AL) == QN()
  @test maximum(link_dimensions(psi)) == 1

  @test vumps_converged(9e-6, 8e-6, 2; vumps_tol=1e-5, energy_tol=1e-4, energy_mismatch_tol=1e-4, stable_iterations=2)
  @test !vumps_converged(2e-5, 8e-6, 2; vumps_tol=1e-5, energy_tol=1e-4, energy_mismatch_tol=1e-4, stable_iterations=2)
end
```

- [ ] **Step 2: Run tests and observe missing runner APIs**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`.

Expected: FAIL on `initial_infinite_mps`.

- [ ] **Step 3: Implement initial QN state construction**

In `VUMPSRunner.jl`, define a deterministic state closure over
`default_occupied_sites(cfg)`, then construct:

```julia
function initial_infinite_mps(c::InfiniteCylinderConfig; occupied_sites=default_occupied_sites(c))
  occupied = Set(occupied_sites)
  length(occupied) == particles_per_cell(c) || throw(ArgumentError("initial pattern particle count mismatch"))
  initstate(n) = mod1(n, sites_per_cell(c)) in occupied ? "Occ" : "Emp"
  sites = infsiteinds("Fermion", sites_per_cell(c); initstate, conserve_qns=true)
  return sites, initstate, InfMPS(sites, initstate)
end
```

Add `link_dimensions`, centered site-charge checks, and the pure
`vumps_converged` helper used by the test.

- [ ] **Step 4: Implement per-iteration records around the pinned backend**

Define:

```julia
struct VUMPSRecord
  stage::Int
  iteration::Int
  maxlinkdim::Int
  energy_left::Float64
  energy_right::Float64
  energy_mismatch::Float64
  delta_energy::Float64
  eps_left::Float64
  eps_right::Float64
  precision_error::Float64
  elapsed_seconds::Float64
  converged::Bool
end

struct VUMPSResult
  psi::InfiniteCanonicalMPS
  records::Vector{VUMPSRecord}
  converged::Bool
  reason::String
end
```

Call the pinned `ITensorInfiniteMPS.tdvp_iteration` with
`ITensorInfiniteMPS.vumps_solver`, explicit left/right epsilon arrays,
`time_step=-Inf`, and the sequential multisite update. Convert scalar or vector
left/right energy returns to unit-cell totals consistently. Require two stable
iterations by default.

- [ ] **Step 5: Implement checked staged subspace expansion**

For each target in the requested schedule:

```julia
before = link_dimensions(psi)
if maximum(before) < target
  expanded = subspace_expansion(psi, H; maxdim=target, cutoff)
  after = link_dimensions(expanded)
  maximum(after) > maximum(before) || error("subspace expansion made no progress toward maxdim=$target")
  psi = expanded
end
```

Run VUMPS iterations after every stage and return nonconverged status instead
of silently accepting a maximum-iteration exit.

- [ ] **Step 6: Run helper tests**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`.

Expected: state, centered-QN, and convergence-helper tests pass. The full smoke
is added after observables and output exist.

- [ ] **Step 7: Commit the runner**

```bash
git add dmrg/idmrg/src/VUMPSRunner.jl dmrg/idmrg/src/InfiniteCylinderDMRG.jl dmrg/idmrg/test
git commit -m "feat: add residual-reporting VUMPS runner"
```

### Task 5: Implement energy, density, entanglement, and transfer observables

**Files:**
- Create: `dmrg/idmrg/src/Observables.jl`
- Create: `dmrg/idmrg/test/test_observables.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing normalization and product-spectrum tests**

```julia
using Test
using InfiniteCylinderDMRG

@testset "observable helpers" begin
  cfg = InfiniteCylinderConfig(; Ly=6, x_period=1)
  e = normalize_energy(cfg, -6.0)
  @test e.per_cell == -6.0
  @test e.per_x == -6.0
  @test e.per_unit_cell == -2.0
  @test e.per_site == -1.0

  sites, initstate, psi = initial_infinite_mps(cfg)
  ent = entanglement_data(psi, cfg; cut_x=1)
  @test isapprox(sum(level.probability for level in ent.levels), 1.0; atol=1e-12)
  @test ent.entropy ≈ 0.0 atol=1e-12
  @test !isempty(ent.sectors)
  @test all(level.physical_charge == level.raw_charge / charge_scale(cfg) for level in ent.levels)
end
```

- [ ] **Step 2: Run and observe missing observable APIs**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`.

Expected: FAIL on `normalize_energy`.

- [ ] **Step 3: Implement energy and density observables**

Define named result types and:

```julia
function normalize_energy(c::InfiniteCylinderConfig, ecell::Real)
  return (
    per_cell=Float64(ecell),
    per_x=Float64(ecell) / c.x_period,
    per_unit_cell=Float64(ecell) / unit_cells_per_cell(c),
    per_site=Float64(ecell) / sites_per_cell(c),
  )
end

energy_data(psi, H, c) = normalize_energy(c, sum(real, expect(psi, H)))
density_data(psi, c) = [
  (site=n, x=site_coordinates(c, n)[1], y=site_coordinates(c, n)[2], density=real(expect(psi, "N", n)))
  for n in 1:sites_per_cell(c)
]
```

- [ ] **Step 4: Implement charge-resolved entanglement data**

At cut `b=cut_x*Ly`, SVD `psi.C[b]`, normalize squared diagonal singular
values, map singular-index block ranges to their QNs, orient the left Schmidt
index consistently, divide raw charge by `charge_scale(cfg)`, and return both
individual levels and sector-weight sums. Assert probability normalization to
`1e-10`.

- [ ] **Step 5: Implement neutral transfer spectrum and correlation length**

Use `TransferMatrix(psi.AL)`, a neutral-flux random input tensor, and
`KrylovKit.eigsolve(..., :LM)`. Return eigenvalues, explicit residual norms,
the Krylov convergence count, and:

```julia
ratio = abs(values[2] / values[1])
xi_cell = -inv(log(ratio))
xi_x = c.x_period * xi_cell
```

Reject fewer than two converged eigenpairs, nonpositive ratios, ratios not
strictly below one within tolerance, and nonfinite correlation lengths.

- [ ] **Step 6: Run observable tests**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`.

Expected: normalization, product entanglement, charge conversion, and transfer
failure-path tests pass.

- [ ] **Step 7: Commit observables**

```bash
git add dmrg/idmrg/src/Observables.jl dmrg/idmrg/src/InfiniteCylinderDMRG.jl dmrg/idmrg/test
git commit -m "feat: add infinite cylinder observables"
```

### Task 6: Add output files and HDF5 checkpoint validation

**Files:**
- Create: `dmrg/idmrg/src/Output.jl`
- Create: `dmrg/idmrg/test/test_checkpoint.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write a failing checkpoint round-trip test**

Use `mktempdir()` to save a product state, reload it with the same config, and
assert that entanglement entropy and density agree to `1e-10`. Also assert that
loading with `Ly=4` throws a configuration-mismatch error.

- [ ] **Step 2: Run and observe the missing output API**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`.

Expected: FAIL on `save_checkpoint`.

- [ ] **Step 3: Implement atomic HDF5 checkpoints**

Write to `state.h5.tmp`, store the infinite state using the backend HDF5
extension, and store these attributes:

```text
format = infinite_cylinder_vumps_v1
backend_commit = 765f2777703bc1138b009adbed1b97bde1973402
Ly, x_period, filling_num, filling_den, phi_y
sites_per_cell, particles_per_cell, charge_scale
```

Close the file and atomically rename it to `state.h5`. On load, compare every
geometry, filling, and backend field before reading the state.

- [ ] **Step 4: Implement TOML and TSV writers**

Add deterministic writers for `summary.toml`, `convergence.tsv`, `density.tsv`,
`entanglement_spectrum.tsv`, `schmidt_sectors.tsv`, and
`transfer_spectrum.tsv`. Invalid or nonconverged values receive explicit
`valid=false` or `converged=false` fields instead of replacement numbers.

- [ ] **Step 5: Run checkpoint and output tests**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`.

Expected: checkpoint round-trip, mismatch rejection, and text-header tests pass.

- [ ] **Step 6: Commit output support**

```bash
git add dmrg/idmrg/src/Output.jl dmrg/idmrg/src/InfiniteCylinderDMRG.jl dmrg/idmrg/test
git commit -m "feat: add infinite state output and checkpoints"
```

### Task 7: Add mixed-transfer branch tracking and raw flux metrics

**Files:**
- Create: `dmrg/idmrg/src/BranchTracking.jl`
- Create: `dmrg/idmrg/test/test_branch_tracking.jl`
- Modify: `dmrg/idmrg/src/InfiniteCylinderDMRG.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing pure-metric tests**

```julia
using Test
using InfiniteCylinderDMRG

@testset "branch metrics" begin
  @test bhattacharyya_distance(Dict(0 => 1.0), Dict(0 => 1.0)) ≈ 0.0
  @test bhattacharyya_distance(Dict(0 => 1.0), Dict(1 => 1.0)) ≈ 1.0
  flags = branch_flags(; fidelity=0.9, charge_step=1 / 3, sector_distance=0.2, spectrum_distance=0.6, energy_crossing=false)
  @test flags.fidelity
  @test flags.charge
  @test flags.sector
  @test flags.spectrum
  @test !flags.energy_crossing
end
```

- [ ] **Step 2: Run and observe missing metric functions**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`.

Expected: FAIL on `bhattacharyya_distance`.

- [ ] **Step 3: Implement sector, spectrum, and flag metrics**

Implement union-of-keys Bhattacharyya distance, discrete integer-QN shift scan,
charge-matched lowest-20 spectrum distance, and the approved default thresholds:

```julia
branch_flags(; fidelity, charge_step, sector_distance, spectrum_distance, energy_crossing,
  fidelity_drop_tol=1e-3, charge_jump_tol=0.1, sector_tol=0.1, spectrum_tol=0.5) = (
    fidelity=(1 - fidelity) > fidelity_drop_tol,
    charge=abs(charge_step) > charge_jump_tol,
    sector=sector_distance > sector_tol,
    spectrum=spectrum_distance > spectrum_tol,
    energy_crossing=energy_crossing,
  )
```

- [ ] **Step 4: Implement mixed-transfer fidelity**

Construct a mixed `ITensorMap` over one full cell from `psi1.AL` and the
primed dagger of `psi2.AL`, solve its dominant eigenvalue, solve the two self
transfer dominant eigenvalues, and return the normalized fidelity and
eigensolver residual. Reuse identical site indices across warm starts. Mark
the result invalid if any dominant solve fails.

- [ ] **Step 5: Implement candidate selection**

For `ground`, select the converged candidate with minimum `energy.per_site`.
For `adiabatic`, require a previous state and select the converged candidate
with maximum valid mixed-transfer fidelity. Return all candidates and the
selected index; do not alter their raw Schmidt data.

- [ ] **Step 6: Run branch tests and commit**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`, expect all pure
metrics and identical-state fidelity tests to pass, then:

```bash
git add dmrg/idmrg/src/BranchTracking.jl dmrg/idmrg/src/InfiniteCylinderDMRG.jl dmrg/idmrg/test
git commit -m "feat: add infinite branch diagnostics"
```

### Task 8: Add single-point and flux-scan command-line workflows

**Files:**
- Create: `dmrg/idmrg/bin/run_vumps.jl`
- Create: `dmrg/idmrg/bin/run_flux_scan.jl`
- Create: `dmrg/idmrg/test/test_cli.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`

- [ ] **Step 1: Write failing CLI parsing tests**

Test `--Ly=6`, `--x_period=1`, `--filling_num=1`, `--filling_den=3`,
`--phi_y=0.0`, `--maxdim=1,4,8`, and invalid non-`key=value` arguments. Assert
that parsed settings reproduce the configuration and schedule exactly.

- [ ] **Step 2: Implement shared CLI parsing in `run_vumps.jl`**

Follow the existing DMRG `--key=value` convention. The single-point command:

1. parses and validates the configuration;
2. configures ITensor/BLAS threading;
3. creates or loads a compatible infinite state;
4. builds the twisted `InfiniteSum{MPO}`;
5. runs staged VUMPS;
6. computes observables only after optimization;
7. writes every output file and checkpoint;
8. exits nonzero when convergence or required observables are invalid unless
   `--allow_nonconverged=true` is explicitly supplied.

- [ ] **Step 3: Implement the raw flux scan**

`run_flux_scan.jl` creates one site-index set, loops through an inclusive flux
grid, warm-starts the previous selected state, optionally adds specified cold
patterns, retains all candidate directories, computes raw Schmidt polarization
relative to the first point, computes branch metrics, and writes
`scan_summary.tsv` and `branch_events.tsv`. It never applies an offset or modulo
unwrap.

- [ ] **Step 4: Run CLI unit tests**

Run `julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl`.

Expected: CLI parsing and invalid-input tests pass without launching a full
smoke run.

- [ ] **Step 5: Commit CLI workflows**

```bash
git add dmrg/idmrg/bin dmrg/idmrg/test
git commit -m "feat: add VUMPS and flux scan commands"
```

### Task 9: Run and lock the minimal and `Ly=6` smoke tests

**Files:**
- Create: `dmrg/idmrg/test/test_smoke.jl`
- Modify: `dmrg/idmrg/test/runtests.jl`
- Create: `dmrg/idmrg/jobs/run_smoke.pbs`

- [ ] **Step 1: Add opt-in smoke test entry points**

Guard expensive tests with `IDMRG_RUN_SMOKE=true`. The minimal smoke command is:

```bash
IDMRG_RUN_SMOKE=true julia --project=dmrg/idmrg dmrg/idmrg/bin/run_vumps.jl \
  --Ly=2 --x_period=3 --filling_num=1 --filling_den=3 \
  --t1=1.0 --t3=0.2 --V1=1.0 --V2=0.0 --V3=0.0 \
  --phi_y=0.0 --maxdim=1,4,8 --cutoff=1e-8 \
  --vumps_tol=1e-5 --maxiter=50 --transfer_tol=1e-8 \
  --output=dmrg/idmrg/smoke_output/Ly2_xp3
```

The target command changes to `--Ly=6 --x_period=1 --vumps_tol=1e-4` and writes
`dmrg/idmrg/smoke_output/Ly6_xp1`.

- [ ] **Step 2: Run the minimal smoke and debug only evidence-backed failures**

Expected: actual bond growth, converged VUMPS residual, finite four-way energy
normalization, two converged neutral transfer eigenpairs, positive finite
correlation length, normalized entanglement probabilities, nonempty sector
data, and a reloadable checkpoint.

- [ ] **Step 3: Run the `Ly=6` target smoke**

Expected: the same data contract and `vumps_tol=1e-4` convergence. A backend
zero-expansion or residual failure remains a failing acceptance criterion.

- [ ] **Step 4: Add a PBS smoke script**

Create a 24-core PBS script that activates `dmrg/idmrg`, sets BLAS threads to
one, runs the minimal smoke first, and runs the `Ly=6` smoke only after the first
command exits zero. It writes logs below `dmrg/idmrg/logs/`.

- [ ] **Step 5: Keep generated smoke outputs out of git and commit tests/jobs**

Add only the output directory pattern to `.gitignore`, then:

```bash
git add .gitignore dmrg/idmrg/test dmrg/idmrg/jobs
git commit -m "test: add infinite cylinder smoke coverage"
```

### Task 10: Document finite/infinite mapping and perform final verification

**Files:**
- Create: `dmrg/idmrg/README.md`
- Modify: `dmrg/README.md`

- [ ] **Step 1: Write the iDMRG/VUMPS README**

Document the pinned experimental dependency, isolated instantiate command,
single-point and flux commands, `Ly/x_period/Np_cell` arithmetic, all four
energy normalizations, neutral-correlation-length limitation, HDF5 compatibility,
and the finite/infinite mapping table from the approved design.

State explicitly:

```text
The finite Lx=15 cumulative charge is raw data. The infinite workflow uses raw
Schmidt-charge polarization and does not add 1/3 or perform modulo unwrapping.
A branch-switch interpretation requires fidelity, Schmidt-sector,
entanglement-spectrum, and candidate-energy evidence.
```

- [ ] **Step 2: Link the new module from `dmrg/README.md`**

Add an `Infinite cylinder` section pointing to `dmrg/idmrg/README.md`, naming
VUMPS as the algorithm and warning that it uses a separate project.

- [ ] **Step 3: Run complete verification**

Run:

```bash
julia --project=. dmrg/test_convergence_logic.jl
julia --project=. dmrg/test_flux_merge.jl
julia --project=. dmrg/test_flux_segment_checkpoint.jl
julia --project=dmrg/idmrg dmrg/idmrg/test/runtests.jl
python3 -m py_compile dmrg/plot_charge_pump.py
git diff --check
git status --short
```

Expected: all finite and fast infinite tests pass, Python compilation succeeds,
no whitespace errors appear, and only intended documentation changes remain.

- [ ] **Step 4: Verify smoke summaries**

Read both `summary.toml` files and assert `converged=true`, required-observable
validity, expected geometry/filling, backend commit, and correct charge scale.

- [ ] **Step 5: Commit documentation**

```bash
git add dmrg/idmrg/README.md dmrg/README.md
git commit -m "docs: explain infinite cylinder VUMPS workflow"
```

- [ ] **Step 6: Request code review and address only verified findings**

Use `superpowers:requesting-code-review`, run the review against the approved
design and this plan, apply technically valid findings with targeted tests, and
rerun the complete verification commands before declaring completion.
