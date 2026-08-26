# Physical Optical Operators on a Finite Torus Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the aliased single-displacement Peierls construction by physically correct finite-torus operators derived from `H(k + Ax*xhat)`, then regenerate the validated `4x6, Np=4, m=2` optical-response benchmark.

**Architecture:** The Bloch Hamiltonian and its analytic first and second Cartesian `kx` derivatives are the source of truth. Each matrix-valued Bloch function is inverse-Fourier transformed on the finite momentum grid with the existing orbital embedding, and only then accumulated into directed finite-torus hopping matrix elements; the existing response-Lanczos and Kubo layers remain unchanged.

**Tech Stack:** Julia 1.10.10, LinearAlgebra, Test, JLD2, Plots, the existing momentum-sector ED code, SSH/rsync to W003.

---

## Execution environment

All Julia execution, ED, benchmark runs, and plotting must run on `W003`. The Mac is limited to source edits, Git, static diff checks, and synchronization.

Use:

```text
Local worktree:     /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/optical-response-lanczos
Remote source:      /home/public/shajy/codex_runs/optical_response_lanczos
Remote Julia env:   /home/public/shajy/codex_envs/optical_response_lanczos
Remote Julia depot: /home/public/shajy/codex_depots/optical_response_lanczos
```

Before every W003 test, synchronize with:

```bash
rsync -az \
  --exclude='.git/' \
  --exclude='.worktrees/' \
  --exclude='case4_4x6_Np4/optical_response_output/' \
  /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/optical-response-lanczos/ \
  W003:/home/public/shajy/codex_runs/optical_response_lanczos/
```

Every Julia command below is executed through:

```bash
ssh W003 'cd /home/public/shajy/codex_runs/optical_response_lanczos && \
  JULIA_DEPOT_PATH=/home/public/shajy/codex_depots/optical_response_lanczos \
  julia --project=/home/public/shajy/codex_envs/optical_response_lanczos <arguments>'
```

## File map

- Modify: `shared/hoppings.jl` — analytic Bloch derivatives, general inverse Fourier helper, finite-torus accumulation, and public `H(Ax)`, `Jx`, `Kxx` constructors.
- Modify: `test/hopping_derivatives_test.jl` — independent Bloch and single-particle finite-torus derivative regression tests.
- Modify: `test/case4_optical_integration_test.jl` — interacting-sector finite-difference and known physical benchmark assertions.
- Modify: `case4_4x6_Np4/benchmark_optical_response.jl` — operator-convention metadata in in-memory and saved results.
- Regenerate on W003: `case4_4x6_Np4/optical_response_output/*` — JLD2, DAT, and PDF artifacts; this directory remains ignored by Git.

### Task 1: Add analytic Bloch derivatives

**Files:**
- Modify: `test/hopping_derivatives_test.jl`
- Modify: `shared/hoppings.jl`

- [ ] **Step 1: Add the failing analytic-derivative test**

Put this test before the finite-torus test in `test/hopping_derivatives_test.jl`:

```julia
@testset "analytic Bloch kx derivatives" begin
    t1, t3 = 1.0, 0.2
    delta1 = 1e-6
    delta2 = 1e-4
    for k in ([0.173, -0.291], [1.137, 0.419], [-0.733, 1.271])
        data = get_Hk_x_derivatives(k, t1, t3)
        hp1 = get_Hk(k .+ [delta1, 0.0], t1, t3)
        hm1 = get_Hk(k .- [delta1, 0.0], t1, t3)
        hp2 = get_Hk(k .+ [delta2, 0.0], t1, t3)
        hm2 = get_Hk(k .- [delta2, 0.0], t1, t3)
        d1_fd = (hp1 - hm1) / (2delta1)
        d2_fd = (hp2 - 2data.Hk + hm2) / delta2^2

        @test data.Hk == get_Hk(k, t1, t3)
        @test norm(data.dHdkx - d1_fd) / max(norm(data.dHdkx), 1.0) < 1e-9
        @test norm(data.d2Hdkx2 - d2_fd) / max(norm(data.d2Hdkx2), 1.0) < 2e-7
    end
end
```

- [ ] **Step 2: Synchronize and verify RED on W003**

Run:

```bash
julia --threads=1 -e 'using Test, LinearAlgebra; include("shared/lattice.jl"); include("shared/hoppings.jl"); include("test/hopping_derivatives_test.jl")'
```

Expected: FAIL with `UndefVarError: get_Hk_x_derivatives not defined`.

- [ ] **Step 3: Implement the minimal analytic derivative API**

In `shared/hoppings.jl`, share the Pauli matrices and paper vectors and define:

```julia
const A_PAPER = [[0.5, -sqrt(3)/2], [0.5, sqrt(3)/2], [-1.0, 0.0]]
const SIGMA_X = ComplexF64[0 1; 1 0]
const SIGMA_Y = ComplexF64[0 -im; im 0]
const SIGMA_Z = ComplexF64[1 0; 0 -1]
const IDENTITY_2 = Matrix{ComplexF64}(I, 2, 2)

function get_Hk_x_derivatives(k::Vector{Float64},
                               t1::Float64, t3::Float64)
    phases = [dot(k, a) for a in A_PAPER]
    g = [2t1 * cos(phase) for phase in phases]
    dg = [-2t1 * a[1] * sin(phase)
          for (a, phase) in zip(A_PAPER, phases)]
    d2g = [-2t1 * a[1]^2 * cos(phase)
           for (a, phase) in zip(A_PAPER, phases)]
    g0 = 2t3 * sum(cos(2phase) for phase in phases)
    dg0 = -4t3 * sum(a[1] * sin(2phase)
                     for (a, phase) in zip(A_PAPER, phases))
    d2g0 = -8t3 * sum(a[1]^2 * cos(2phase)
                      for (a, phase) in zip(A_PAPER, phases))
    assemble(v, v0) = (v[1] * SIGMA_X + v[2] * SIGMA_Y +
                       v[3] * SIGMA_Z + v0 * IDENTITY_2)
    return (Hk=assemble(g, g0),
            dHdkx=assemble(dg, dg0),
            d2Hdkx2=assemble(d2g, d2g0))
end

function get_Hk(k::Vector{Float64}, t1::Float64, t3::Float64)
    return get_Hk_x_derivatives(k, t1, t3).Hk
end
```

- [ ] **Step 4: Synchronize and verify GREEN on W003**

Run the command from Step 2. Expected: the analytic derivative test passes; the still-existing finite-torus test may continue to pass against the old self-consistent but unphysical `build_hops_Ax` and is replaced in Task 2.

- [ ] **Step 5: Commit the Bloch derivative layer**

```bash
git add shared/hoppings.jl test/hopping_derivatives_test.jl
git commit -m "feat: add analytic Bloch current derivatives"
```

### Task 2: Build physical finite-torus operators from Bloch data

**Files:**
- Modify: `test/hopping_derivatives_test.jl`
- Modify: `shared/hoppings.jl`

- [ ] **Step 1: Replace the circular Peierls test by an independent finite-torus oracle**

Define the following test-only helpers and test in `test/hopping_derivatives_test.jl`:

```julia
function single_particle_matrix(lat, hops)
    operator = zeros(ComplexF64, lat.Ns, lat.Ns)
    for (target, source, amplitude) in hops
        operator[target, source] += amplitude
    end
    return operator
end

function torus_matrix_from_bloch(lat, bloch_matrix)
    operator = zeros(ComplexF64, lat.Ns, lat.Ns)
    a1_uc = collect(lat.a1)
    a2_uc = 2 .* collect(lat.a2)
    Nk = length(lat.kpoints)
    for (source, (ix, iy)) in enumerate(lat.sites)
        source_orb = mod(iy, 2) + 1
        source_cell_y = div(iy - (source_orb - 1), 2)
        for (n1_uc, n2_uc) in lat.uc_trans
            Rcart = n1_uc .* a1_uc .+ n2_uc .* a2_uc
            for target_orb in 1:2
                amplitude = sum(lat.kpoints; init=0.0 + 0.0im) do k
                    delta = SUBLAT_POS[target_orb] .- SUBLAT_POS[source_orb]
                    phase = exp(-1im * dot(k, Rcart .+ delta))
                    phase * bloch_matrix(k)[target_orb, source_orb]
                end / Nk
                tix = ix + n1_uc
                tiy = 2(source_cell_y + n2_uc) + target_orb - 1
                target_site = lat.site_idx[canon_prim(lat, tix, tiy)]
                operator[target_site, source] += amplitude
            end
        end
    end
    return operator
end

@testset "physical finite-torus Ax operators" begin
    lat = RectLat4x6()
    t1, t3 = 1.0, 0.2
    hops0, jx_hops, kxx_hops = build_hops_x_derivatives(lat, t1, t3)
    H0 = single_particle_matrix(lat, hops0)
    Jx = single_particle_matrix(lat, jx_hops)
    Kxx = single_particle_matrix(lat, kxx_hops)

    data(k) = get_Hk_x_derivatives(k, t1, t3)
    H0_ref = torus_matrix_from_bloch(lat, k -> data(k).Hk)
    Jx_ref = torus_matrix_from_bloch(lat, k -> data(k).dHdkx)
    Kxx_ref = torus_matrix_from_bloch(lat, k -> data(k).d2Hdkx2)
    relerr(actual, reference) = norm(actual - reference) / max(norm(reference), 1.0)

    @test relerr(H0, H0_ref) < 1e-12
    @test relerr(Jx, Jx_ref) < 1e-12
    @test relerr(Kxx, Kxx_ref) < 1e-12

    Ax = 0.071
    HAx = single_particle_matrix(lat, build_hops_Ax(lat, t1, t3, Ax))
    HAx_ref = torus_matrix_from_bloch(
        lat, k -> get_Hk(k .+ [Ax, 0.0], t1, t3))
    @test relerr(HAx, HAx_ref) < 1e-12

    delta1 = 1e-6
    delta2 = 1e-4
    Hp1 = single_particle_matrix(lat, build_hops_Ax(lat, t1, t3, delta1))
    Hm1 = single_particle_matrix(lat, build_hops_Ax(lat, t1, t3, -delta1))
    Hp2 = single_particle_matrix(lat, build_hops_Ax(lat, t1, t3, delta2))
    Hm2 = single_particle_matrix(lat, build_hops_Ax(lat, t1, t3, -delta2))
    @test relerr((Hp1 - Hm1) / (2delta1), Jx) < 1e-8
    @test relerr((Hp2 - 2H0 + Hm2) / delta2^2, Kxx) < 1e-5

    Hlegacy = single_particle_matrix(lat, build_hops(lat, t1, t3, 0.0))
    @test relerr(H0, Hlegacy) < 1e-12
    for operator in (H0, Jx, Kxx, HAx)
        @test norm(operator - operator') / max(norm(operator), 1.0) < 1e-12
    end
end
```

- [ ] **Step 2: Synchronize and verify the physical regression test is RED on W003**

Run the Task 1 Step 2 command. Expected: FAIL at the nonzero-`Ax`, `Jx`, and `Kxx` comparisons. The old code's previously measured relative errors are approximately `1.10` for `Jx` and `2.16` for `Kxx`.

- [ ] **Step 3: Implement reusable Fourier and finite-torus accumulation helpers**

Refactor `shared/hoppings.jl` so `fourier_to_real` delegates to:

```julia
function fourier_to_real_from_bloch(lat::GenLat, bloch_matrix)
    Nk = length(lat.kpoints)
    a1_uc = collect(lat.a1)
    a2_uc = 2 .* collect(lat.a2)
    tR = Dict{NTuple{2,Int},Matrix{ComplexF64}}()
    for (n1_uc, n2_uc) in lat.uc_trans
        Rcart = n1_uc .* a1_uc .+ n2_uc .* a2_uc
        tmat = zeros(ComplexF64, 2, 2)
        for k in lat.kpoints
            Hk = bloch_matrix(k)
            for alpha in 1:2, beta in 1:2
                delta = SUBLAT_POS[alpha] .- SUBLAT_POS[beta]
                tmat[alpha, beta] +=
                    exp(-1im * dot(k, Rcart .+ delta)) * Hk[alpha, beta]
            end
        end
        tmat ./= Nk
        tR[(n1_uc, n2_uc)] = tmat
    end
    return tR
end

function hopping_list_from_tR(lat::GenLat, tR; threshold::Float64=1e-12)
    amplitudes = Dict{Tuple{Int,Int},ComplexF64}()
    for (source, (ix, iy)) in enumerate(lat.sites)
        source_orb = mod(iy, 2) + 1
        source_cell_y = div(iy - (source_orb - 1), 2)
        for ((n1_uc, n2_uc), tmat) in tR
            for target_orb in 1:2
                tix = ix + n1_uc
                tiy = 2(source_cell_y + n2_uc) + target_orb - 1
                target = lat.site_idx[canon_prim(lat, tix, tiy)]
                key = (target, source)
                amplitudes[key] = get(amplitudes, key, 0.0 + 0.0im) +
                                  tmat[target_orb, source_orb]
            end
        end
    end
    keys_sorted = sort!(collect(keys(amplitudes)))
    return [(target, source, amplitudes[(target, source)])
            for (target, source) in keys_sorted
            if abs(amplitudes[(target, source)]) >= threshold]
end
```

Do not synthesize reverse bonds or attach a single unfolded displacement in this helper.

- [ ] **Step 4: Implement physical `H(Ax)`, `Jx`, and `Kxx` constructors**

Add:

```julia
function fourier_to_real(lat::GenLat, t1::Float64, t3::Float64)
    return fourier_to_real_from_bloch(lat, k -> get_Hk(k, t1, t3))
end

function fourier_to_real_Ax(lat::GenLat, t1::Float64,
                            t3::Float64, Ax::Float64)
    return fourier_to_real_from_bloch(
        lat, k -> get_Hk(k .+ [Ax, 0.0], t1, t3))
end

function fourier_to_real_x_derivatives(lat::GenLat,
                                       t1::Float64, t3::Float64)
    h = fourier_to_real_from_bloch(
        lat, k -> get_Hk_x_derivatives(k, t1, t3).Hk)
    jx = fourier_to_real_from_bloch(
        lat, k -> get_Hk_x_derivatives(k, t1, t3).dHdkx)
    kxx = fourier_to_real_from_bloch(
        lat, k -> get_Hk_x_derivatives(k, t1, t3).d2Hdkx2)
    return h, jx, kxx
end

function build_hops_Ax(lat::GenLat, t1::Float64, t3::Float64,
                       Ax::Float64=0.0)
    return hopping_list_from_tR(lat, fourier_to_real_Ax(lat, t1, t3, Ax))
end

function build_hops_x_derivatives(lat::GenLat, t1::Float64, t3::Float64)
    h, jx, kxx = fourier_to_real_x_derivatives(lat, t1, t3)
    return (hopping_list_from_tR(lat, h),
            hopping_list_from_tR(lat, jx),
            hopping_list_from_tR(lat, kxx))
end
```

Keep the `phi_y` spectral-flow path separate. `build_hops(lat, t1, t3, phi_y)` may retain `RealSpaceBond`, but it is no longer the source of `Ax` optical derivatives.

- [ ] **Step 5: Synchronize and verify GREEN on W003**

Run the Task 1 Step 2 command. Expected: all Bloch and finite-torus tests pass, including `H(0)` legacy equivalence and Hermiticity.

- [ ] **Step 6: Run the full existing suite on W003**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: all testsets pass.

- [ ] **Step 7: Commit the finite-torus operator layer**

```bash
git add shared/hoppings.jl test/hopping_derivatives_test.jl
git commit -m "fix: derive finite-torus optical operators from Bloch flux"
```

### Task 3: Pin the interacting benchmark and operator metadata

**Files:**
- Modify: `test/case4_optical_integration_test.jl`
- Modify: `case4_4x6_Np4/benchmark_optical_response.jl`

- [ ] **Step 1: Add interacting finite-difference and benchmark assertions**

Add a testset that constructs the `Np=4, m=2` sector and materializes `H(+delta)`, `H(-delta)`, `Jx`, and `Kxx` through `materialize_sector_operator`. Use `delta1=1e-6` for `Jx` and `delta2=1e-4` for `Kxx`, and assert:

```julia
relerr(actual, reference) = norm(actual - reference) / max(norm(reference), 1.0)
@test relerr((Hp1 - Hm1) / (2delta1), Jx) < 1e-8
@test relerr((Hp2 - 2H0 + Hm2) / delta2^2, Kxx) < 1e-5
```

In the existing `4x6 optical response integration` test add:

```julia
@test result.Eg ≈ -9.3135758879616 atol=1e-11
@test result.source_norm2 ≈ 10.844308142783673 atol=1e-8
@test result.Kexp ≈ 4.4792283713 atol=1e-7
@test result.operator_convention ==
      "H_A(k)=H(k+A_x*xhat); Jx=dH/dA_x; Kxx=d2H/dA_x2"
```

Also construct the regular conductivity directly from the stored exact
transition energies and weights, independently of `evaluate_curve`:

```julia
function explicit_regular(data, area, omega, eta)
    z = complex(omega, eta)
    Gz = sum(data.weights ./ (z .- data.deltas))
    Gminusz = sum(data.weights ./ (-z .- data.deltas))
    chi0 = -2sum(data.weights ./ data.deltas)
    return (2pi * im / area) * (Gz + Gminusz - chi0) / z
end

explicit_curve = [explicit_regular(result.exact_data, result.area,
                                   omega, result.eta)
                  for omega in result.omegas]
@test scaled_max_error(explicit_curve, result.exact_curve.regular) < 1e-12
```

This is the old sum-over-states content evaluated with the corrected physical
`Jx`; it must agree pointwise with the resolvent/Kubo `regular` curve.

In the output test, load the JLD2 file and assert that the same `operator_convention` string is stored.

- [ ] **Step 2: Synchronize and verify RED on W003**

Run the full suite from Task 2 Step 6. Expected: FAIL because `result.operator_convention` and the saved metadata do not yet exist. The numerical physical assertions should already pass after Task 2.

- [ ] **Step 3: Add operator-convention metadata without changing Kubo code**

In `case4_4x6_Np4/benchmark_optical_response.jl`, define:

```julia
const OPTICAL_OPERATOR_CONVENTION =
    "H_A(k)=H(k+A_x*xhat); Jx=dH/dA_x; Kxx=d2H/dA_x2"
```

Add `operator_convention=OPTICAL_OPERATOR_CONVENTION` to the benchmark result, save it as `operator_convention=result.operator_convention` in JLD2, and prepend this line to both DAT outputs:

```julia
println(io, "# operator_convention: $(result.operator_convention)")
```

Do not modify `shared/optical_response.jl`.

- [ ] **Step 4: Synchronize and verify GREEN on W003**

Run the full suite from Task 2 Step 6. Expected: all testsets pass; the final integration curve error is `< 1e-6` for its reduced frequency grid.

- [ ] **Step 5: Commit benchmark regression coverage and metadata**

```bash
git add test/case4_optical_integration_test.jl \
  case4_4x6_Np4/benchmark_optical_response.jl
git commit -m "test: pin physical optical response benchmark"
```

### Task 4: Run full W003 benchmark and return artifacts

**Files:**
- Regenerate: `case4_4x6_Np4/optical_response_output/optical_response_benchmark.jld2`
- Regenerate: `case4_4x6_Np4/optical_response_output/optical_response_curves.dat`
- Regenerate: `case4_4x6_Np4/optical_response_output/lanczos_convergence.dat`
- Regenerate: `case4_4x6_Np4/optical_response_output/*.pdf`

- [ ] **Step 1: Run the complete frequency-grid benchmark on W003**

Run:

```bash
julia --threads=1 case4_4x6_Np4/benchmark_optical_response.jl
```

Expected printed invariants:

```text
sector=2
E0=-9.313575887962
<Kxx> approximately 4.4792283713
||QJxg||^2 approximately 10.8443081428
final scaled max error < 1e-8
```

- [ ] **Step 2: Synchronize generated artifacts back to the local worktree**

```bash
mkdir -p /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/optical-response-lanczos/case4_4x6_Np4/optical_response_output
rsync -az \
  W003:/home/public/shajy/codex_runs/optical_response_lanczos/case4_4x6_Np4/optical_response_output/ \
  /Users/shajianyu/CMP_manybody/partial_Hall_crystal/.worktrees/optical-response-lanczos/case4_4x6_Np4/optical_response_output/
```

- [ ] **Step 3: Perform fresh final verification**

Synchronize source once more, then rerun on W003:

```bash
julia --threads=1 test/runtests.jl
julia --threads=1 case4_4x6_Np4/benchmark_optical_response.jl --quick --no-plot
```

On the Mac run only:

```bash
git diff --check
git status --short --branch
```

Expected: W003 test and quick benchmark commands exit zero, static diff check is clean, and Git status contains no unintended tracked modifications.

- [ ] **Step 4: Record the final implementation commit**

If plan checkbox updates are the only remaining tracked change:

```bash
git add docs/superpowers/plans/2026-08-26-physical-optical-operators.md
git commit -m "docs: record physical optical operator verification"
```
