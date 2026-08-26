# Optical Conductivity Response-Lanczos Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在现有 `4×6, Np=4, m=2` ED case 上实现包含真实 diamagnetic term 和 Drude weight 的一般纵向光电导，并用完整本征态求和验证 response-Lanczos/resolvent 结果。

**Architecture:** 从统一的均匀矢势 Peierls 耦合生成 `H`、`Jx`、`Kxx`，由 `Q*Jx*g` 生成 response Lanczos 起始态，只保存三对角系数并在后处理中计算 resolvent。完整对角化是小系统专用 reference；生产式 Lanczos 和 exact reference 使用同一组算符，并在 benchmark runner 中逐点比较。

**Tech Stack:** Julia 1.10-compatible Julia、LinearAlgebra、SparseArrays、Test、KrylovKit（现有 ED）、Plots、JLD2；所有数值运行位于 W003。

---

## Execution environment and W003 synchronization

All Julia tests, full diagonalizations, benchmark runs, and figure generation run on SSH host `W003`, not on the Mac. The Mac is limited to source edits, Git operations, `git diff --check`, and file synchronization.

Use these isolated paths on W003:

```text
Remote source copy: /home/public/shajy/codex_runs/optical_response_lanczos
Remote Julia env:   /home/public/shajy/codex_envs/optical_response_lanczos
Existing repo:      /home/public/shajy/partial_Hall_crystal   (read-only for this task)
```

Do not run `rsync --delete` against `/home/public/shajy/partial_Hall_crystal`; it contains unrelated untracked calculations. Before each remote test, synchronize the isolated local worktree with:

```bash
ssh W003 'mkdir -p /home/public/shajy/codex_runs/optical_response_lanczos'
rsync -az \
  --exclude='.git/' \
  --exclude='.worktrees/' \
  --exclude='case4_4x6_Np4/optical_response_output/' \
  "$LOCAL_WORKTREE/" \
  W003:/home/public/shajy/codex_runs/optical_response_lanczos/
```

Every `julia ...` test command in the tasks below means the following W003 wrapper after that sync:

```bash
ssh W003 'cd /home/public/shajy/codex_runs/optical_response_lanczos && \
  JULIA_DEPOT_PATH=/home/public/shajy/codex_depots/optical_response_lanczos \
  julia --project=/home/public/shajy/codex_envs/optical_response_lanczos <arguments>'
```

After a benchmark run, synchronize generated results back into the isolated local worktree:

```bash
mkdir -p "$LOCAL_WORKTREE/case4_4x6_Np4/optical_response_output"
rsync -az \
  W003:/home/public/shajy/codex_runs/optical_response_lanczos/case4_4x6_Np4/optical_response_output/ \
  "$LOCAL_WORKTREE/case4_4x6_Np4/optical_response_output/"
```

The remote runtime uses Julia 1.10.10. All implementation must remain compatible with Julia 1.10; local Julia 1.12-only syntax or APIs are forbidden.

## File map

- Modify: `shared/hoppings.jl` — 保留每条 hopping 的未折叠物理位移，并提供 `A_x` Peierls hopping 及一、二阶导数 hopping。
- Create: `shared/optical_response.jl` — ground projection、response Lanczos、continued-fraction resolvent、exact spectral reference、一般 Kubo 与 Drude/regular 分解。
- Create: `test/runtests.jl` — 统一加载现有 ED 基础文件和新增测试。
- Create: `test/hopping_derivatives_test.jl` — Peierls 导数和 many-body Hermiticity 测试。
- Create: `test/optical_response_test.jl` — resolvent、Lanczos、投影和 Kubo 分解单元测试。
- Create: `test/case4_optical_integration_test.jl` — `4×6, Np=4, m=2` exact/Lanczos 集成测试。
- Create: `case4_4x6_Np4/benchmark_optical_response.jl` — 完整 benchmark、JLD2/文本/PDF 输出和命令行入口。
- Modify: `.gitignore` — 忽略 benchmark 生成的独立输出目录。

## Task 1: Add displacement-aware Peierls hopping derivatives

**Files:**
- Modify: `shared/hoppings.jl`
- Create: `test/runtests.jl`
- Create: `test/hopping_derivatives_test.jl`

- [ ] **Step 1: Create the test runner and failing hopping tests**

Create `test/runtests.jl`:

```julia
using Test
using LinearAlgebra
using SparseArrays

include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")

include("hopping_derivatives_test.jl")
```

Create `test/hopping_derivatives_test.jl`:

```julia
@testset "uniform Ax hopping derivatives" begin
    lat = RectLat4x6()
    t1, t3 = 1.0, 0.2
    hops0, jx_hops, kxx_hops = build_hops_x_derivatives(lat, t1, t3)

    asdict(hops) = Dict((tgt, src) => amp for (tgt, src, amp) in hops)
    @test asdict(hops0) == asdict(build_hops(lat, t1, t3, 0.0))

    delta = 1e-5
    hp = asdict(build_hops_Ax(lat, t1, t3, delta))
    hm = asdict(build_hops_Ax(lat, t1, t3, -delta))
    h0 = asdict(hops0)
    jx = asdict(jx_hops)
    kxx = asdict(kxx_hops)

    for key in keys(h0)
        @test (hp[key] - hm[key]) / (2delta) ≈ jx[key] atol=1e-9 rtol=1e-9
        @test (hp[key] - 2h0[key] + hm[key]) / delta^2 ≈ kxx[key] atol=2e-6 rtol=2e-6
    end

    basis = gen_basis(lat.Ns, 1)
    sec = build_ksector(basis, lat, 2)
    H = build_sparse_H(sec, lat, hops0, 0.0, 0.0, 0.0)
    Jx = build_sparse_H(sec, lat, jx_hops, 0.0, 0.0, 0.0)
    Kxx = build_sparse_H(sec, lat, kxx_hops, 0.0, 0.0, 0.0)
    @test norm(H - H') / max(norm(H), 1.0) < 1e-12
    @test norm(Jx - Jx') / max(norm(Jx), 1.0) < 1e-12
    @test norm(Kxx - Kxx') / max(norm(Kxx), 1.0) < 1e-12
end
```

- [ ] **Step 2: Run the test and verify the new API is missing**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: FAIL with `UndefVarError: build_hops_x_derivatives not defined`.

- [ ] **Step 3: Implement displacement-aware bonds and derivative hopping lists**

Refactor `shared/hoppings.jl` so the existing `build_hops` delegates to one bond enumerator. Preserve its public signature and add the following definitions:

```julia
struct RealSpaceBond
    tgt::Int
    src::Int
    amp::ComplexF64
    winding_T2::Int
    displacement::NTuple{2,Float64}
end

function build_real_space_bonds(lat::GenLat, t1::Float64, t3::Float64)
    tR = fourier_to_real(lat, t1, t3)
    a1_uc = collect(lat.a1)
    a2_uc = 2 .* collect(lat.a2)
    bond_dict = Dict{Tuple{Int,Int},RealSpaceBond}()

    for (si, (ix, iy)) in enumerate(lat.sites)
        source_orb = mod(iy, 2) + 1
        source_cell_y = div(iy - (source_orb - 1), 2)

        for ((n1_uc, n2_uc), tmat) in tR
            Rcart = n1_uc .* a1_uc .+ n2_uc .* a2_uc
            for target_orb in 1:2
                amp = tmat[target_orb, source_orb]
                abs(amp) < 1e-12 && continue

                tix = ix + n1_uc
                tiy = 2 * (source_cell_y + n2_uc) + target_orb - 1
                target = canon_prim(lat, tix, tiy)
                haskey(lat.site_idx, target) || continue
                ti = lat.site_idx[target]
                haskey(bond_dict, (ti, si)) && continue

                delta = SUBLAT_POS[target_orb] .- SUBLAT_POS[source_orb]
                dcart = Rcart .+ delta
                d = (Float64(dcart[1]), Float64(dcart[2]))
                w2 = winding_T2(lat, tix, tiy)

                bond_dict[(ti, si)] = RealSpaceBond(ti, si, amp, w2, d)
                bond_dict[(si, ti)] = RealSpaceBond(si, ti, conj(amp), -w2,
                                                    (-d[1], -d[2]))
            end
        end
    end

    keys_sorted = sort!(collect(keys(bond_dict)))
    return [bond_dict[key] for key in keys_sorted]
end

function build_hops(lat::GenLat, t1::Float64, t3::Float64,
                    phi_y::Float64=0.0)
    return [(b.tgt, b.src, b.amp * cis(b.winding_T2 * phi_y))
            for b in build_real_space_bonds(lat, t1, t3)]
end

function build_hops_Ax(lat::GenLat, t1::Float64, t3::Float64,
                       Ax::Float64=0.0)
    return [(b.tgt, b.src, b.amp * cis(Ax * b.displacement[1]))
            for b in build_real_space_bonds(lat, t1, t3)]
end

function build_hops_x_derivatives(lat::GenLat, t1::Float64, t3::Float64)
    bonds = build_real_space_bonds(lat, t1, t3)
    hops = [(b.tgt, b.src, b.amp) for b in bonds]
    jx = [(b.tgt, b.src, im * b.displacement[1] * b.amp) for b in bonds]
    kxx = [(b.tgt, b.src, -b.displacement[1]^2 * b.amp) for b in bonds]
    return hops, jx, kxx
end
```

Delete the old duplicated body of `build_hops`; retain `get_Hk`, `SUBLAT_POS`, and `fourier_to_real` unchanged.

- [ ] **Step 4: Run the hopping derivative tests**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: PASS for `uniform Ax hopping derivatives`, including the finite-difference and Hermiticity checks.

- [ ] **Step 5: Verify an existing case still builds the same zero-flux spectrum**

Run:

```bash
julia --threads=4 case4_4x6_Np4/main.jl
```

Expected: command exits 0; sector dimensions and the lowest energies agree with the tracked `case4_4x6_Np4/spectrum_Np4.dat` to `1e-8`.

- [ ] **Step 6: Commit the hopping derivative layer**

```bash
git add shared/hoppings.jl test/runtests.jl test/hopping_derivatives_test.jl
git commit -m "feat: derive current operators from Peierls hopping"
```

## Task 2: Implement the projected response-Lanczos kernel

**Files:**
- Create: `shared/optical_response.jl`
- Modify: `test/runtests.jl`
- Create: `test/optical_response_test.jl`

- [ ] **Step 1: Add failing projection and Lanczos tests**

Append to `test/runtests.jl` before the test-file includes:

```julia
include("../shared/optical_response.jl")
```

and append:

```julia
include("optical_response_test.jl")
```

Create `test/optical_response_test.jl`:

```julia
@testset "projected response source" begin
    H = ComplexF64[0 0 0; 0 1 0; 0 0 3]
    J = ComplexF64[2 1-im 0; 1+im 0 2; 0 2 0]
    g = ComplexF64[1, 0, 0]
    f = projected_response_source(J, g)
    @test dot(g, f) ≈ 0.0 atol=1e-14
    @test f ≈ ComplexF64[0, 1+im, 0]
    @test_throws ArgumentError projected_response_source(Matrix{ComplexF64}(I, 3, 3), g)
end

@testset "response Lanczos continued fraction" begin
    H = ComplexF64[0.0 0.2 0.0 0.0;
                   0.2 1.0 0.3 0.0;
                   0.0 0.3 2.0 0.4;
                   0.0 0.0 0.4 4.0]
    E0 = 0.0
    f = ComplexF64[0, 1, 1im, 0]
    kernel = response_lanczos(H, E0, f; mmax=4, breakdown_tol=1e-14)
    @test length(kernel.alpha) == 4
    @test length(kernel.beta) == 3
    @test kernel.source_norm2 ≈ real(dot(f, f))

    for z in (0.3 + 0.2im, 1.4 + 0.1im, -0.7 - 0.15im)
        direct = dot(f, (z * I - (H - E0 * I)) \ f)
        @test lanczos_resolvent(kernel, z) ≈ direct atol=1e-11 rtol=1e-11
    end

    closed = response_lanczos(Diagonal(ComplexF64[0, 1, 2]), 0.0,
                              ComplexF64[0, 1, 0]; mmax=3)
    @test closed.breakdown
    @test length(closed.alpha) == 1
    @test isempty(closed.beta)
end
```

- [ ] **Step 2: Run the test and verify the response API is missing**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: FAIL while including missing `shared/optical_response.jl` or with `projected_response_source not defined`.

- [ ] **Step 3: Implement the response source, three-vector Lanczos, and continued fraction**

Create `shared/optical_response.jl`:

```julia
using LinearAlgebra

struct ResponseLanczosKernel
    alpha::Vector{Float64}
    beta::Vector{Float64}
    source_norm2::Float64
    breakdown::Bool
end

function projected_response_source(J::AbstractMatrix,
                                   ground_state::AbstractVector;
                                   ground_subspace::Union{Nothing,AbstractMatrix}=nothing,
                                   atol::Float64=1e-12)
    g = ComplexF64.(ground_state)
    f = ComplexF64.(J * g)
    states = ground_subspace === nothing ? reshape(g, :, 1) : ground_subspace
    for state in eachcol(states)
        abs(norm(state) - 1) <= 1e-10 || throw(ArgumentError("ground projector states must be normalized"))
        f .-= state .* dot(state, f)
    end
    norm(f) > atol || throw(ArgumentError("projected current source has zero norm"))
    return f
end

function response_lanczos(H, E0::Real, source::AbstractVector;
                          mmax::Int=length(source),
                          breakdown_tol::Float64=1e-13)
    n = length(source)
    1 <= mmax <= n || throw(ArgumentError("mmax must satisfy 1 <= mmax <= length(source)"))
    source_norm = norm(source)
    source_norm > breakdown_tol || throw(ArgumentError("response source has zero norm"))

    q_prev = zeros(ComplexF64, n)
    q = ComplexF64.(source) ./ source_norm
    w = similar(q)
    alpha = Float64[]
    beta = Float64[]
    beta_prev = 0.0
    did_breakdown = false

    for step in 1:mmax
        mul!(w, H, q)
        @. w = w - E0 * q - beta_prev * q_prev
        a = real(dot(q, w))
        push!(alpha, a)
        @. w = w - a * q

        # Repeat local orthogonalization against the two recurrence vectors.
        w .-= q .* dot(q, w)
        step > 1 && (w .-= q_prev .* dot(q_prev, w))
        b = norm(w)

        if b <= breakdown_tol
            did_breakdown = true
            break
        end
        step == mmax && break

        push!(beta, b)
        q_prev, q, w = q, w, q_prev
        q ./= b
        beta_prev = b
    end

    return ResponseLanczosKernel(alpha, beta, source_norm^2, did_breakdown)
end

function lanczos_resolvent(kernel::ResponseLanczosKernel, z::Number)
    isempty(kernel.alpha) && throw(ArgumentError("empty Lanczos kernel"))
    denominator = complex(z) - kernel.alpha[end]
    for j in length(kernel.alpha)-1:-1:1
        denominator = complex(z) - kernel.alpha[j] - kernel.beta[j]^2 / denominator
    end
    return kernel.source_norm2 / denominator
end

function truncate_kernel(kernel::ResponseLanczosKernel, m::Int)
    1 <= m <= length(kernel.alpha) || throw(ArgumentError("invalid truncated dimension"))
    return ResponseLanczosKernel(copy(kernel.alpha[1:m]),
                                 copy(kernel.beta[1:max(0, m-1)]),
                                 kernel.source_norm2,
                                 kernel.breakdown && m == length(kernel.alpha))
end
```

- [ ] **Step 4: Run the response-Lanczos unit tests**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: all hopping and response-Lanczos tests PASS.

- [ ] **Step 5: Commit the reusable Lanczos kernel**

```bash
git add shared/optical_response.jl test/runtests.jl test/optical_response_test.jl
git commit -m "feat: add projected response Lanczos kernel"
```

## Task 3: Implement exact spectral reference and general Kubo decomposition

**Files:**
- Modify: `shared/optical_response.jl`
- Modify: `test/optical_response_test.jl`

- [ ] **Step 1: Add failing tests for exact spectral data and conductivity decomposition**

Append to `test/optical_response_test.jl`:

```julia
@testset "exact spectral reference and Kubo decomposition" begin
    H = ComplexF64[0 0 0; 0 1 0; 0 0 3]
    J = ComplexF64[0 1-im 2; 1+im 0 0; 2 0 0]
    K = ComplexF64[5 0 0; 0 0 0; 0 0 0]
    eig = eigen(Hermitian(H))
    Eg = eig.values[1]
    g = eig.vectors[:, 1]
    f = projected_response_source(J, g)
    data = exact_spectral_data(eig.values, eig.vectors, Eg, f)

    @test data.deltas ≈ [1.0, 3.0]
    @test data.weights ≈ [2.0, 4.0]
    for z in (0.5 + 0.2im, -0.8 - 0.1im)
        direct = dot(f, (z * I - (H - Eg * I)) \ f)
        @test exact_resolvent(data, z) ≈ direct atol=1e-12
    end

    Kexp = real(dot(g, K * g))
    G(z) = exact_resolvent(data, z)
    parts = conductivity_parts(G, Kexp, 7.5, 0.8, 0.065)
    @test parts.total ≈ parts.drude + parts.regular atol=1e-13
    @test parts.drude_weight ≈ Kexp - 2sum(data.weights ./ data.deltas)

    kernel = response_lanczos(H, Eg, f; mmax=3)
    @test lanczos_drude_weight(kernel, Kexp) ≈ parts.drude_weight atol=1e-12
end

@testset "system area" begin
    @test system_area(RectLat4x6()) ≈ 12sqrt(3) atol=1e-12
end
```

- [ ] **Step 2: Run the tests and verify the Kubo helpers are missing**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: FAIL with `UndefVarError: exact_spectral_data not defined`.

- [ ] **Step 3: Implement exact reference, area, Drude weight, and conductivity parts**

Append to `shared/optical_response.jl`:

```julia
struct ExactSpectralData
    deltas::Vector{Float64}
    weights::Vector{Float64}
end

function exact_spectral_data(eigenvalues::AbstractVector,
                             eigenvectors::AbstractMatrix,
                             E0::Real,
                             source::AbstractVector;
                             gap_tol::Float64=1e-10)
    deltas_all = real.(eigenvalues .- E0)
    amplitudes = eigenvectors' * source
    keep = findall(delta -> delta > gap_tol, deltas_all)
    return ExactSpectralData(deltas_all[keep], abs2.(amplitudes[keep]))
end

function exact_resolvent(data::ExactSpectralData, z::Number)
    return sum(data.weights ./ (complex(z) .- data.deltas))
end

function system_area(lat::GenLat)
    a1 = collect(lat.a1)
    a2 = collect(lat.a2)
    T1 = lat.T1[1] .* a1 .+ lat.T1[2] .* a2
    T2 = lat.T2[1] .* a1 .+ lat.T2[2] .* a2
    return abs(T1[1] * T2[2] - T1[2] * T2[1])
end

function lanczos_drude_weight(kernel::ResponseLanczosKernel, Kexp::Real)
    return Float64(Kexp + 2real(lanczos_resolvent(kernel, 0.0 + 0.0im)))
end

function conductivity_parts(G, Kexp::Real, area::Real,
                            omega::Real, eta::Real)
    area > 0 || throw(ArgumentError("system area must be positive"))
    eta > 0 || throw(ArgumentError("eta must be positive"))
    z = complex(omega, eta)
    chi = G(z) + G(-z)
    chi0 = 2real(G(0.0 + 0.0im))
    D = Float64(Kexp + chi0)
    prefactor = 2pi * im / area
    drude = prefactor * D / z
    regular = prefactor * (chi - chi0) / z
    return (total=drude + regular,
            drude=drude,
            regular=regular,
            drude_weight=D)
end

function scaled_max_error(actual::AbstractVector, reference::AbstractVector)
    length(actual) == length(reference) || throw(DimensionMismatch("curve lengths differ"))
    scale = max(1.0, maximum(abs, reference))
    return maximum(abs, actual .- reference) / scale
end
```

- [ ] **Step 4: Run all unit tests**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: all tests PASS, including exact resolvent equality and `total = drude + regular`.

- [ ] **Step 5: Commit the Kubo evaluation layer**

```bash
git add shared/optical_response.jl test/optical_response_test.jl
git commit -m "feat: evaluate general Kubo optical conductivity"
```

## Task 4: Build the `4×6, Np=4, m=2` exact/Lanczos benchmark

**Files:**
- Create: `case4_4x6_Np4/benchmark_optical_response.jl`
- Create: `test/case4_optical_integration_test.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Add a failing integration test for the benchmark function**

Append to `test/runtests.jl`:

```julia
include("case4_optical_integration_test.jl")
```

Create `test/case4_optical_integration_test.jl`:

```julia
include("../case4_4x6_Np4/benchmark_optical_response.jl")

@testset "4x6 optical response integration" begin
    result = run_optical_benchmark(
        omegas=collect(0.05:0.10:2.05),
        eta=0.065,
        requested_m=[50, 100, 200, 400],
        write_outputs=false,
        make_plots=false,
    )
    @test result.sector == 2
    @test result.ground_residual < 1e-10
    @test result.hermiticity.H < 1e-12
    @test result.hermiticity.Jx < 1e-12
    @test result.hermiticity.Kxx < 1e-12
    @test abs(result.source_overlap) < 1e-11
    @test all(diff(result.errors) .<= 1e-10)
    @test result.errors[end] < 1e-6
end
```

- [ ] **Step 2: Run the integration test and verify the benchmark entry is missing**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: FAIL while including missing `case4_4x6_Np4/benchmark_optical_response.jl`.

- [ ] **Step 3: Implement the benchmark computation without output side effects**

Create `case4_4x6_Np4/benchmark_optical_response.jl` with guarded includes and these public entry points:

```julia
if !isdefined(Main, :GenLat)
    include("../shared/lattice.jl")
    include("../shared/hoppings.jl")
    include("../shared/basis.jl")
    include("../shared/ksector.jl")
    include("../shared/hamiltonian.jl")
end
if !isdefined(Main, :ResponseLanczosKernel)
    include("../shared/optical_response.jl")
end

using LinearAlgebra
using SparseArrays
using Printf

function hermiticity_error(operator)
    return norm(operator - operator') / max(norm(operator), 1.0)
end

function evaluate_curve(G, Kexp, area, omegas, eta)
    parts = [conductivity_parts(G, Kexp, area, omega, eta) for omega in omegas]
    return (total=ComplexF64[p.total for p in parts],
            regular=ComplexF64[p.regular for p in parts],
            drude=ComplexF64[p.drude for p in parts],
            drude_weight=parts[1].drude_weight)
end

function run_optical_benchmark(; omegas=collect(0.0:0.001:10.0),
                               eta::Float64=0.065,
                               requested_m::Vector{Int}=[50, 100, 200, 400, 800],
                               write_outputs::Bool=true,
                               make_plots::Bool=true,
                               output_dir::String=joinpath(@__DIR__, "optical_response_output"))
    lat = RectLat4x6()
    basis = gen_basis(lat.Ns, 4)
    sec = build_ksector(basis, lat, 2)
    hops, jx_hops, kxx_hops = build_hops_x_derivatives(lat, 1.0, 0.2)
    H = build_sparse_H(sec, lat, hops, 1.0, 0.0, 0.0)
    Jx = build_sparse_H(sec, lat, jx_hops, 0.0, 0.0, 0.0)
    Kxx = build_sparse_H(sec, lat, kxx_hops, 0.0, 0.0, 0.0)
    dropzeros!(H); dropzeros!(Jx); dropzeros!(Kxx)

    herm = (H=hermiticity_error(H),
            Jx=hermiticity_error(Jx),
            Kxx=hermiticity_error(Kxx))
    maximum(values(herm)) < 1e-10 || error("operator Hermiticity check failed: $herm")

    eig = eigen(Hermitian(Matrix(H)))
    Eg = eig.values[1]
    g = eig.vectors[:, 1]
    residual = norm(H * g - Eg * g)
    f = projected_response_source(Jx, g)
    overlap = dot(g, f)
    Kexp_complex = dot(g, Kxx * g)
    abs(imag(Kexp_complex)) < 1e-10 || error("Kxx expectation is not real")
    Kexp = real(Kexp_complex)
    area = system_area(lat)

    exact_data = exact_spectral_data(eig.values, eig.vectors, Eg, f)
    exact_G(z) = exact_resolvent(exact_data, z)
    exact_curve = evaluate_curve(exact_G, Kexp, area, omegas, eta)

    requested = sort(unique(min.(requested_m, length(f))))
    maximum(requested) >= 1 || throw(ArgumentError("requested_m must contain a positive value"))
    full_kernel = response_lanczos(H, Eg, f; mmax=maximum(requested))
    m_values = sort(unique(min.(requested, length(full_kernel.alpha))))
    kernels = [truncate_kernel(full_kernel, m) for m in m_values]
    curves = map(kernels) do kernel
        G(z) = lanczos_resolvent(kernel, z)
        evaluate_curve(G, Kexp, area, omegas, eta)
    end
    errors = [scaled_max_error(curve.regular, exact_curve.regular) for curve in curves]

    result = (sector=2, Eg=Eg, ground_residual=residual,
              hermiticity=herm, source_overlap=overlap,
              source_norm2=real(dot(f, f)), Kexp=Kexp, area=area,
              exact_data=exact_data, exact_curve=exact_curve,
              m_values=m_values, kernels=kernels, curves=curves,
              errors=errors, omegas=Float64.(omegas), eta=eta)

    write_outputs && write_benchmark_outputs(result, output_dir; make_plots=make_plots)
    return result
end
```

For this task, add a temporary stub at the end so the integration test can keep `write_outputs=false`:

```julia
function write_benchmark_outputs(result, output_dir; make_plots=true)
    error("output writer is implemented in Task 5")
end
```

Do not add a top-level call yet.

- [ ] **Step 4: Run the integration test and inspect convergence**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: all unit tests and the physical integration test PASS. The printed/inspected `errors` must decrease with `M`; if the final coarse-grid error is above `1e-6`, fix recurrence orthogonality before continuing.

- [ ] **Step 5: Commit the side-effect-free benchmark core**

```bash
git add case4_4x6_Np4/benchmark_optical_response.jl test/runtests.jl test/case4_optical_integration_test.jl
git commit -m "test: benchmark optical Lanczos against full ED"
```

## Task 5: Add reproducible data and figure output

**Files:**
- Modify: `case4_4x6_Np4/benchmark_optical_response.jl`
- Modify: `test/case4_optical_integration_test.jl`
- Modify: `.gitignore`

- [ ] **Step 1: Create the isolated Julia environment on W003**

Run:

```bash
ssh W003 'mkdir -p /home/public/shajy/codex_envs/optical_response_lanczos \
    /home/public/shajy/codex_depots/optical_response_lanczos && \
  JULIA_DEPOT_PATH=/home/public/shajy/codex_depots/optical_response_lanczos \
  julia -e '\''using Pkg; \
    Pkg.activate("/home/public/shajy/codex_envs/optical_response_lanczos"); \
    Pkg.add(["KrylovKit", "JLD2", "Plots"])'\'''
ssh W003 'JULIA_DEPOT_PATH=/home/public/shajy/codex_depots/optical_response_lanczos \
  julia --project=/home/public/shajy/codex_envs/optical_response_lanczos \
  -e '\''using KrylovKit, JLD2, Plots; println("output dependencies available")'\'''
```

Expected: the second command prints `output dependencies available`. Only the dedicated environment under `/home/public/shajy/codex_envs/optical_response_lanczos` changes; the existing server repository and its untracked `Project.toml`/`Manifest.toml` remain untouched.

- [ ] **Step 2: Add a failing output smoke test**

Append to `test/case4_optical_integration_test.jl`:

```julia
@testset "optical benchmark output" begin
    mktempdir() do output_dir
        result = run_optical_benchmark(
            omegas=collect(0.05:0.20:1.05),
            eta=0.065,
            requested_m=[50],
            write_outputs=true,
            make_plots=false,
            output_dir=output_dir,
        )
        @test isfile(joinpath(output_dir, "optical_response_benchmark.jld2"))
        @test isfile(joinpath(output_dir, "optical_response_curves.dat"))
        @test isfile(joinpath(output_dir, "lanczos_convergence.dat"))
        @test result.sector == 2
    end
end
```

- [ ] **Step 3: Run the smoke test and verify the stub fails**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: FAIL with `output writer is implemented in Task 5`.

- [ ] **Step 4: Implement JLD2, text, plot, and terminal outputs**

Replace the stub in `case4_4x6_Np4/benchmark_optical_response.jl` with a lazy dependency loader and output writer:

```julia
function write_benchmark_outputs(result, output_dir; make_plots=true)
    @eval using JLD2
    mkpath(output_dir)
    final_index = length(result.m_values)
    final_curve = result.curves[final_index]
    final_kernel = result.kernels[final_index]

    JLD2.jldsave(joinpath(output_dir, "optical_response_benchmark.jld2");
        sector=result.sector, Eg=result.Eg,
        ground_residual=result.ground_residual,
        hermiticity=result.hermiticity,
        source_overlap=result.source_overlap,
        source_norm2=result.source_norm2,
        Kexp=result.Kexp, area=result.area,
        eta=result.eta, omegas=result.omegas,
        m_values=result.m_values, errors=result.errors,
        alpha=final_kernel.alpha, beta=final_kernel.beta,
        exact_deltas=result.exact_data.deltas,
        exact_weights=result.exact_data.weights,
        exact_total=result.exact_curve.total,
        exact_regular=result.exact_curve.regular,
        lanczos_total=final_curve.total,
        lanczos_regular=final_curve.regular,
        exact_drude_weight=result.exact_curve.drude_weight,
        lanczos_drude_weight=final_curve.drude_weight)

    open(joinpath(output_dir, "optical_response_curves.dat"), "w") do io
        println(io, "# omega Re_sigma_exact Im_sigma_exact Re_sigma_L Im_sigma_L Re_reg_exact Im_reg_exact Re_reg_L Im_reg_L")
        for i in eachindex(result.omegas)
            values = (result.omegas[i],
                      real(result.exact_curve.total[i]), imag(result.exact_curve.total[i]),
                      real(final_curve.total[i]), imag(final_curve.total[i]),
                      real(result.exact_curve.regular[i]), imag(result.exact_curve.regular[i]),
                      real(final_curve.regular[i]), imag(final_curve.regular[i]))
            println(io, join(values, ' '))
        end
    end

    open(joinpath(output_dir, "lanczos_convergence.dat"), "w") do io
        println(io, "# M scaled_max_error_regular")
        for (m, err) in zip(result.m_values, result.errors)
            println(io, "$m $err")
        end
    end

    if make_plots
        @eval using Plots
        p1 = Plots.plot(result.omegas, real.(result.exact_curve.regular),
                        label="exact", xlabel="omega", ylabel="Re sigma_xx^reg")
        for (m, curve) in zip(result.m_values, result.curves)
            Plots.plot!(p1, result.omegas, real.(curve.regular), label="M=$m")
        end
        Plots.savefig(p1, joinpath(output_dir, "optical_conductivity_regular.pdf"))

        p2 = Plots.plot(result.omegas, imag.(result.exact_curve.total),
                        label="exact", xlabel="omega", ylabel="Im sigma_xx total")
        Plots.plot!(p2, result.omegas, imag.(final_curve.total),
                    label="Lanczos M=$(result.m_values[end])")
        Plots.savefig(p2, joinpath(output_dir, "optical_conductivity_total.pdf"))

        p3 = Plots.plot(result.m_values, result.errors, marker=:circle,
                        xscale=:log10, yscale=:log10, legend=false,
                        xlabel="Lanczos M", ylabel="scaled max error")
        Plots.savefig(p3, joinpath(output_dir, "lanczos_convergence.pdf"))
    end

    @printf("sector=%d  E0=%.12f  residual=%.3e\n",
            result.sector, result.Eg, result.ground_residual)
    @printf("<Kxx>=%.12g  ||QJxg||^2=%.12g  area=%.12g\n",
            result.Kexp, result.source_norm2, result.area)
    @printf("Dxx exact=%.12g  Lanczos=%.12g\n",
            result.exact_curve.drude_weight, final_curve.drude_weight)
    for (m, err) in zip(result.m_values, result.errors)
        @printf("M=%4d  scaled max error=%.3e\n", m, err)
    end
end
```

Add the guarded CLI entry:

```julia
if abspath(PROGRAM_FILE) == @__FILE__
    quick = "--quick" in ARGS
    no_plot = "--no-plot" in ARGS
    result = run_optical_benchmark(
        omegas=quick ? collect(0.0:0.05:10.0) : collect(0.0:0.001:10.0),
        requested_m=quick ? [50, 100, 200, 400, 600] : [50, 100, 200, 400, 800, 10_000],
        make_plots=!no_plot,
    )
    target_error = quick ? 1e-6 : 1e-8
    result.errors[end] < target_error ||
        error("final Lanczos curve did not reach target accuracy $target_error")
end
```

Append to `.gitignore`:

```gitignore
case4_4x6_Np4/optical_response_output/
```

- [ ] **Step 5: Run the full tests and a quick CLI benchmark**

Run:

```bash
julia --threads=1 test/runtests.jl
julia --threads=4 case4_4x6_Np4/benchmark_optical_response.jl --quick --no-plot
```

Expected: tests PASS; benchmark exits 0, writes the three data files, and reports monotonically decreasing error.

- [ ] **Step 6: Commit output support**

```bash
git add .gitignore case4_4x6_Np4/benchmark_optical_response.jl test/case4_optical_integration_test.jl
git commit -m "feat: save optical response benchmark results"
```

## Task 6: Run the full-resolution benchmark and verify acceptance criteria

**Files:**
- Generated, ignored: `case4_4x6_Np4/optical_response_output/*`
- Modify only if a verified numerical defect is found: `shared/optical_response.jl`, `case4_4x6_Np4/benchmark_optical_response.jl`

- [ ] **Step 1: Run the complete test suite**

Run:

```bash
julia --threads=1 test/runtests.jl
```

Expected: all unit, finite-difference, and physical integration tests PASS.

- [ ] **Step 2: Run the full benchmark with the approved grid**

Run:

```bash
julia --threads=4 case4_4x6_Np4/benchmark_optical_response.jl
```

Expected: command exits 0; output contains JLD2, two conductivity PDFs, one convergence PDF, and two text data files. Final scaled maximum error is below `1e-8`.

- [ ] **Step 3: Synchronize W003 benchmark results back to the Mac worktree**

Run on the Mac:

```bash
mkdir -p "$LOCAL_WORKTREE/case4_4x6_Np4/optical_response_output"
rsync -az \
  W003:/home/public/shajy/codex_runs/optical_response_lanczos/case4_4x6_Np4/optical_response_output/ \
  "$LOCAL_WORKTREE/case4_4x6_Np4/optical_response_output/"
```

Expected: JLD2, DAT, and PDF files exist locally with the same sizes reported by `ssh W003 'ls -l .../optical_response_output'`.

- [ ] **Step 4: Inspect numerical invariants**

Confirm from the terminal summary and JLD2 data:

```text
Hermiticity(H), Hermiticity(Jx), Hermiticity(Kxx) < 1e-12
ground residual < 1e-10
abs(<g|QJx|g>) < 1e-11
abs(Im <Kxx>) < 1e-10
errors decrease as M increases
abs(Dxx_Lanczos - Dxx_exact) / max(1, abs(Dxx_exact)) < 1e-8
```

- [ ] **Step 5: Visually inspect the generated PDFs**

Open:

```text
case4_4x6_Np4/optical_response_output/optical_conductivity_regular.pdf
case4_4x6_Np4/optical_response_output/optical_conductivity_total.pdf
case4_4x6_Np4/optical_response_output/lanczos_convergence.pdf
```

Expected: exact and final Lanczos curves are visually indistinguishable; convergence plot contains finite positive errors and no missing/NaN points.

- [ ] **Step 6: Run repository hygiene checks**

Run:

```bash
git diff --check
git status --short
```

Expected: no whitespace errors; generated benchmark outputs do not appear as untracked files; only intentional source changes are present.

- [ ] **Step 7: Commit any verification-only numerical fix, if one was required**

If Steps 1–5 required a source correction, stage only the files changed for that correction and commit:

```bash
git add shared/optical_response.jl case4_4x6_Np4/benchmark_optical_response.jl
git commit -m "fix: stabilize optical response benchmark convergence"
```

If no correction was required, do not create an empty commit.
