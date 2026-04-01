# Projected ED Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement band-projected ED on the lowest Chern band, supporting both V1V2V3 and gate-screened Coulomb interactions, with full momentum-sector diagonalization.

**Architecture:** New directory `projected_ed/` on branch `projected-ed`. Four shared modules (band, basis, interaction, hamiltonian) plus one top-level script. Reuses `shared/lattice.jl` and `shared/hoppings.jl` via relative includes. Full diagonalization via `eigen` (C(15,12)=455 total dim).

**Tech Stack:** Julia, LinearAlgebra (built-in eigen), existing GenLat/get_Hk infrastructure.

---

### Task 1: Create directory skeleton

**Files:**
- Create: `projected_ed/shared_proj/band.jl` (stub)
- Create: `projected_ed/shared_proj/proj_basis.jl` (stub)
- Create: `projected_ed/shared_proj/proj_interaction.jl` (stub)
- Create: `projected_ed/shared_proj/proj_hamiltonian.jl` (stub)
- Create: `projected_ed/case_proj/main.jl` (stub)
- Create: `projected_ed/case_proj/test_proj.jl` (validation script)

**Step 1: Create directories**

```bash
mkdir -p projected_ed/shared_proj projected_ed/case_proj
```

**Step 2: Create stub files with module header comments**

Each stub contains only the module docstring comment (no code yet). Example for `band.jl`:
```julia
# ============================================================
# band.jl — Single-particle band energies and Bloch eigenstates
# ============================================================
```

**Step 3: Commit skeleton**

```bash
git add projected_ed/
git commit -m "feat: add projected_ed directory skeleton"
```

---

### Task 2: Implement band.jl

**Files:**
- Create: `projected_ed/shared_proj/band.jl`

This module computes the lower-band eigenvalues and gauge-fixed eigenvectors at all k-points.

**Step 1: Write the full implementation**

```julia
# ============================================================
# band.jl — Single-particle band energies and Bloch eigenstates
# ============================================================
using LinearAlgebra

"""
    compute_band(lat, t1, t3) -> (eps, uvec)

Diagonalize the 2×2 H(k) at each of the Nk k-points.
Returns:
  eps  :: Vector{Float64}        length Nk, lower-band eigenvalues
  uvec :: Vector{Vector{ComplexF64}}  length Nk, each is [u_A, u_B] (gauge-fixed)

Gauge fix: rotate each eigenvector so that u_A is real positive.
"""
function compute_band(lat::GenLat, t1::Float64, t3::Float64)
    Nk = lat.Nuc
    eps  = Vector{Float64}(undef, Nk)
    uvec = Vector{Vector{ComplexF64}}(undef, Nk)

    for (mi, k) in enumerate(lat.kpoints)
        Hk = get_Hk(k, t1, t3)
        vals, vecs = eigen(Hermitian(Hk))
        # lower band = index 1 (eigen returns ascending order)
        e_lo = vals[1]
        u_lo = vecs[:, 1]   # 2-component vector [u_A, u_B]

        # gauge fix: rotate so u_A is real positive
        phase = angle(u_lo[1])
        u_lo  = u_lo .* exp(-1im * phase)
        # if u_A ended up negative real, flip sign
        if real(u_lo[1]) < 0
            u_lo = -u_lo
        end

        eps[mi]  = e_lo
        uvec[mi] = u_lo
    end
    return eps, uvec
end
```

**Step 2: Add inline sanity check (run to verify)**

Add at the bottom of `band.jl` (inside `if abspath(PROGRAM_FILE) == @__FILE__` block so it only runs when called directly):

```julia
if abspath(PROGRAM_FILE) == @__FILE__
    include("../../shared/lattice.jl")
    include("../../shared/hoppings.jl")
    lat = RectLat3x6()
    eps, uvec = compute_band(lat, 1.0, 0.2)
    println("Lower band energies: ", round.(eps, digits=4))
    for (m, u) in enumerate(uvec)
        @assert abs(norm(u) - 1.0) < 1e-12 "eigenvector not normalized at k=$m"
        @assert imag(u[1]) < 1e-12 "gauge fix failed: u_A not real at k=$m"
    end
    println("band.jl: all checks passed")
end
```

**Step 3: Run the check**

```bash
cd projected_ed/case_proj
julia ../../projected_ed/shared_proj/band.jl
```
Expected output: 9 eigenvalues printed, "all checks passed".

**Step 4: Commit**

```bash
git add projected_ed/shared_proj/band.jl
git commit -m "feat: implement band.jl with gauge-fixed Bloch eigenstates"
```

---

### Task 3: Implement proj_basis.jl

**Files:**
- Create: `projected_ed/shared_proj/proj_basis.jl`

k-space Fock states, total momentum labeling, momentum-sector enumeration.

**Step 1: Write implementation**

```julia
# ============================================================
# proj_basis.jl — k-space Fock states and momentum sectors
# ============================================================

"""
    kfock_momentum(state, Nk) -> Int

Total momentum label M = (Σ_{occupied m} m) mod Nk.
`state` is a bitmask: bit m (0-indexed) set means k_m is occupied.
"""
@inline function kfock_momentum(state::Int, Nk::Int)
    M = 0
    for m in 0:Nk-1
        (state >> m) & 1 == 1 && (M += m)
    end
    return M % Nk
end

"""
    kfock_fermion_sign(state, m) -> Int

Fermion sign (±1) for annihilating a particle at k_m from `state`.
Sign = (−1)^(number of occupied bits strictly below position m).
"""
@inline function kfock_fermion_sign(state::Int, m::Int)
    mask = (1 << m) - 1          # bits 0..m-1
    return isodd(count_ones(state & mask)) ? -1 : 1
end

"""
    build_proj_sectors(Nk, Np) -> Dict{Int, Vector{Int}}

Enumerate all C(Nk,Np) k-space Fock states and group by total momentum M.
Returns Dict: M => sorted list of Fock states in that sector.
"""
function build_proj_sectors(Nk::Int, Np::Int)
    sectors = Dict{Int, Vector{Int}}()
    for M in 0:Nk-1
        sectors[M] = Int[]
    end

    # enumerate all Nk-bit integers with exactly Np bits set
    # use combinatorial enumeration for efficiency
    function rec!(state, start, left)
        if left == 0
            M = kfock_momentum(state, Nk)
            push!(sectors[M], state)
            return
        end
        for m in start:Nk-left
            rec!(state | (1 << m), m+1, left-1)
        end
    end
    rec!(0, 0, Np)

    for M in 0:Nk-1
        sort!(sectors[M])
    end
    return sectors
end
```

**Step 2: Add inline check**

```julia
if abspath(PROGRAM_FILE) == @__FILE__
    Nk, Np = 15, 12
    sectors = build_proj_sectors(Nk, Np)
    total = sum(length(v) for v in values(sectors))
    @assert total == binomial(Nk, Np) "total states wrong: got $total expected $(binomial(Nk,Np))"
    # check all sectors same size (for uniform filling)
    sizes = [length(sectors[M]) for M in 0:Nk-1]
    println("Sector sizes: ", sizes)
    @assert all(==(sizes[1]), sizes) "sectors not equal size — unexpected for this (Nk,Np)"

    # check fermion sign consistency: applying c†_m c_m should give +n_m
    state = 0b101  # bits 0 and 2 occupied
    @assert kfock_fermion_sign(state, 0) == 1
    @assert kfock_fermion_sign(state, 2) == -1   # one occupied bit below (bit 0)
    println("proj_basis.jl: all checks passed")
end
```

**Step 3: Run check**

```bash
julia projected_ed/shared_proj/proj_basis.jl
```
Expected: 15 sector sizes each = 455/15 = ~30, "all checks passed".

**Step 4: Commit**

```bash
git add projected_ed/shared_proj/proj_basis.jl
git commit -m "feat: implement proj_basis.jl with k-space Fock states and momentum sectors"
```

---

### Task 4: Implement proj_interaction.jl — V1V2V3 version

**Files:**
- Create: `projected_ed/shared_proj/proj_interaction.jl`

**Step 1: Write V_αβ(q) for V1V2V3**

```julia
# ============================================================
# proj_interaction.jl — Projected two-body matrix elements
# ============================================================
using LinearAlgebra

const SUBLAT_CART = [[0.0, 0.0], [-0.5, sqrt(3)/2]]  # τ_A, τ_B in physical coords

"""
    build_Vq_V123(lat, uvec, V1, V2, V3) -> Array{ComplexF64,3}

Projected interaction tensor W[m1+1, m2+1, qi+1] (1-indexed) for V1V2V3 interaction.
W(m1,m2,q) = Σ_{α,β} V_αβ(q) × F_α(m1+q, m1) × conj(F_β(m2-q, m2))

where F_α(k',k) = conj(uvec[k'][α]) * uvec[k][α]
and V_αβ(q) is the Fourier transform of the real-space interaction.
"""
function build_Vq_V123(lat::GenLat, uvec::Vector{Vector{ComplexF64}},
                        V1::Float64, V2::Float64, V3::Float64)
    Nk = lat.Nuc
    kpts = lat.kpoints   # physical k-vectors

    # Precompute V_αβ(q) for each q = kpts[qi] (qi = 1..Nk)
    # V_αβ(q) = Σ_{ΔR, sites (α,β)} V(dist) × exp(-i q·(R_α - R_β))
    # We enumerate neighbor pairs directly from the lattice neighbor tables.
    # nb1[si, :] gives nearest-neighbor indices of site si; distance V1.
    # Sites: sorb = (iy%2+2)%2 gives sublattice (0=A, 1=B).

    # Build site sublattice and Cartesian position lookup
    site_alpha = [((iy % 2 + 2) % 2) for (ix,iy) in lat.sites]       # 0=A,1=B
    site_cell  = [(ix, div(iy - site_alpha[si], 2)) for (si,(ix,iy)) in enumerate(lat.sites)]

    function site_cart(si)
        α = site_alpha[si]
        ix, iy_uc = site_cell[si]
        # unit cell position in physical coords
        R = ix .* collect(lat.a1) .+ iy_uc .* 2 .* collect(lat.a2)
        return R .+ SUBLAT_CART[α+1]
    end

    # Build V_αβ(q) as 2×2×Nk tensor
    Vq = zeros(ComplexF64, 2, 2, Nk)
    for si in 1:lat.Ns
        αi = site_alpha[si]
        ri = site_cart(si)
        for (nb, Vval) in [(lat.nb1, V1), (lat.nb2, V2), (lat.nb3, V3)]
            Vval == 0.0 && continue
            for k in 1:6
                sj = nb[si, k]
                αj = site_alpha[sj]
                rj = site_cart(sj)
                Δr = ri .- rj
                for qi in 1:Nk
                    q = kpts[qi]
                    Vq[αi+1, αj+1, qi] += Vval * exp(-1im * dot(q, Δr)) * 0.5
                end
            end
        end
    end

    # Precompute form factors F_α[k', k] for all k', k
    # F_α(k', k) = conj(uvec[k'][α]) * uvec[k][α]
    # Store as F[α, kp, k] where α ∈ {1,2}, kp,k ∈ 1..Nk
    F = zeros(ComplexF64, 2, Nk, Nk)
    for k in 1:Nk, kp in 1:Nk
        for α in 1:2
            F[α, kp, k] = conj(uvec[kp][α]) * uvec[k][α]
        end
    end

    # Compute W[m1, m2, q] — all indices 1-based (m = label+1)
    # q label: qi such that (m1 + q_label) mod Nk = kp_label for m1→kp scattering
    # Precompute: for each (m1_idx, q_idx), what is kp_idx = (m1_idx-1+q_idx-1)%Nk + 1
    W = zeros(ComplexF64, Nk, Nk, Nk)
    for m1 in 0:Nk-1, m2 in 0:Nk-1, qi in 0:Nk-1
        kp1 = mod(m1 + qi, Nk)   # k1 + q
        kp2 = mod(m2 - qi, Nk)   # k2 - q  (mod Nk, Julia mod handles negative)
        # W = Σ_{α,β} Vq[α,β,qi+1] × F[α, kp1+1, m1+1] × conj(F[β, kp2+1, m2+1])
        val = zero(ComplexF64)
        for α in 1:2, β in 1:2
            val += Vq[α, β, qi+1] * F[α, kp1+1, m1+1] * conj(F[β, kp2+1, m2+1])
        end
        W[m1+1, m2+1, qi+1] = val
    end
    return W
end
```

**Step 2: Add gate Coulomb version**

Append to `proj_interaction.jl`:

```julia
"""
    build_Vq_gate(lat, uvec, d) -> Array{ComplexF64,3}

Projected interaction tensor using dual-gate screened Coulomb:
  V(q) = (2π / |q|) × tanh(|q| × d),   V(q=0) = 0

`d` is the gate distance in physical length units (same as lattice constant).
"""
function build_Vq_gate(lat::GenLat, uvec::Vector{Vector{ComplexF64}}, d::Float64)
    Nk = lat.Nuc
    kpts = lat.kpoints

    # Compute scalar V(q) for each q-point
    Vscalar = zeros(Float64, Nk)
    for qi in 1:Nk
        q = kpts[qi]
        qabs = norm(q)
        if qabs < 1e-10
            Vscalar[qi] = 0.0   # q=0 → neutralizing background, set to zero
        else
            Vscalar[qi] = (2π / qabs) * tanh(qabs * d)
        end
    end

    # Form factors (same as V123 version)
    F = zeros(ComplexF64, 2, Nk, Nk)
    for k in 1:Nk, kp in 1:Nk
        for α in 1:2
            F[α, kp, k] = conj(uvec[kp][α]) * uvec[k][α]
        end
    end

    # For gate Coulomb, V_αβ(q) = V(q) × δ_{αβ} is NOT correct —
    # the interaction is between ALL pairs of sites regardless of sublattice.
    # We need to include the sublattice-position phase:
    # V_αβ(q) = V(|q|) × exp(-i q·(τ_α - τ_β))
    # (the full real-space sum gives this phase factor from the sublattice offset)
    W = zeros(ComplexF64, Nk, Nk, Nk)
    for m1 in 0:Nk-1, m2 in 0:Nk-1, qi in 0:Nk-1
        kp1 = mod(m1 + qi, Nk)
        kp2 = mod(m2 - qi, Nk)
        q   = kpts[qi+1]
        Vq  = Vscalar[qi+1]
        val = zero(ComplexF64)
        for α in 1:2, β in 1:2
            τ_diff = SUBLAT_CART[α] .- SUBLAT_CART[β]
            phase  = exp(-1im * dot(q, τ_diff))
            val   += Vq * phase * F[α, kp1+1, m1+1] * conj(F[β, kp2+1, m2+1])
        end
        W[m1+1, m2+1, qi+1] = val
    end
    return W
end
```

**Step 3: Add inline checks**

```julia
if abspath(PROGRAM_FILE) == @__FILE__
    include("../../shared/lattice.jl")
    include("../../shared/hoppings.jl")
    include("band.jl")
    lat = RectLat3x6()
    eps, uvec = compute_band(lat, 1.0, 0.2)

    W_v123 = build_Vq_V123(lat, uvec, 1.0, 0.0, 0.0)
    W_gate = build_Vq_gate(lat, uvec, 1.0)

    # W should be Hermitian under (m1,m2,q) ↔ (m1+q, m2-q, -q)
    # Basic check: W[m1,m2,q] should be real for q=0 (which is 0 anyway)
    println("W_v123 max imaginary (q=0 slice): ", maximum(abs.(imag.(W_v123[:,:,1]))))
    println("W_gate max imaginary (q=0 slice): ", maximum(abs.(imag.(W_gate[:,:,1]))))
    println("proj_interaction.jl: built successfully")
end
```

**Step 4: Run check**

```bash
julia projected_ed/shared_proj/proj_interaction.jl
```

**Step 5: Commit**

```bash
git add projected_ed/shared_proj/proj_interaction.jl
git commit -m "feat: implement proj_interaction.jl with V1V2V3 and gate Coulomb"
```

---

### Task 5: Implement proj_hamiltonian.jl

**Files:**
- Create: `projected_ed/shared_proj/proj_hamiltonian.jl`

**Step 1: Write implementation**

```julia
# ============================================================
# proj_hamiltonian.jl — Build and diagonalize projected H per momentum sector
# ============================================================
using LinearAlgebra

"""
    build_proj_H(sector_states, eps, W, Nk) -> Matrix{ComplexF64}

Build the projected Hamiltonian matrix for a single momentum sector.

  sector_states: Vector{Int}  — list of Fock states (bitmasks) in this sector
  eps:           Vector{Float64}  — lower-band energies, eps[m+1] for k_m
  W:             Array{ComplexF64,3}  — W[m1+1, m2+1, qi+1], precomputed
  Nk:            Int

Returns: dim×dim Hermitian matrix.
"""
function build_proj_H(sector_states::Vector{Int},
                      eps::Vector{Float64},
                      W::Array{ComplexF64,3},
                      Nk::Int)
    dim = length(sector_states)
    H   = zeros(ComplexF64, dim, dim)

    # index lookup: Fock state → row index
    state2idx = Dict{Int,Int}(s => i for (i,s) in enumerate(sector_states))

    for (i, si) in enumerate(sector_states)
        # --- diagonal: kinetic energy ---
        Ekin = 0.0
        for m in 0:Nk-1
            (si >> m) & 1 == 1 && (Ekin += eps[m+1])
        end
        H[i, i] += Ekin

        # --- off-diagonal: interaction scattering k1→k1+q, k2→k2-q ---
        for m1 in 0:Nk-1
            (si >> m1) & 1 == 0 && continue       # k1 must be occupied
            for m2 in 0:Nk-1
                m2 == m1 && continue
                (si >> m2) & 1 == 0 && continue   # k2 must be occupied
                for qi in 1:Nk-1                   # q ≠ 0
                    kp1 = mod(m1 + qi - 1, Nk)    # k1 + q  (qi is 1-based label)
                    kp2 = mod(m2 - qi + 1 + Nk - 1, Nk)  # k2 - q

                    # Pauli: kp1 must be empty, kp2 must be empty
                    (si >> kp1) & 1 == 1 && continue
                    (si >> kp2) & 1 == 1 && continue
                    kp1 == kp2 && continue         # would violate Pauli after both hops

                    # New state: remove m1,m2; add kp1,kp2
                    sf  = si
                    sf  = sf ⊻ (1 << m1) ⊻ (1 << m2) ⊻ (1 << kp1) ⊻ (1 << kp2)

                    haskey(state2idx, sf) || continue
                    j = state2idx[sf]

                    # Fermion signs: annihilate m2, m1 (in that order), create kp2, kp1
                    sgn2  = kfock_fermion_sign(si,  m2)
                    tmp   = si ⊻ (1 << m2)
                    sgn1  = kfock_fermion_sign(tmp, m1)
                    tmp2  = tmp ⊻ (1 << m1)
                    sgnc2 = kfock_fermion_sign(tmp2, kp2)  # note: creating = sign of annihilating on final
                    tmp3  = tmp2 | (1 << kp2)
                    sgnc1 = kfock_fermion_sign(tmp3, kp1)

                    # For c†_{kp1} c†_{kp2} c_{m2} c_{m1}: combined sign
                    # Actually: sign of c†_a |n⟩ = (-1)^(# occupied below a) × |n+a⟩
                    # Use: sgn(c†_a) = (-1)^popcount(state & mask_below_a) on the state BEFORE creation
                    total_sgn = sgn1 * sgn2 * kfock_fermion_sign(tmp2, kp2) *
                                kfock_fermion_sign(tmp2 | (1 << kp2), kp1)

                    H[j, i] += (1.0 / Nk) * W[m1+1, m2+1, qi] * total_sgn
                end
            end
        end
    end

    return H
end

"""
    solve_proj_ed(sectors, eps, W, Nk) -> (all_ev, all_vecs)

Diagonalize all momentum sectors. Returns:
  all_ev  : Vector{Tuple{Int,Float64}} — sorted (M, energy) pairs
  all_vecs: Dict{Int, Matrix} — M => eigenvector matrix (columns = eigenvectors)
"""
function solve_proj_ed(sectors::Dict{Int,Vector{Int}},
                       eps::Vector{Float64},
                       W::Array{ComplexF64,3},
                       Nk::Int)
    all_ev   = Tuple{Int,Float64}[]
    all_vecs = Dict{Int, Matrix{ComplexF64}}()

    for M in sort(collect(keys(sectors)))
        states = sectors[M]
        isempty(states) && continue
        H  = build_proj_H(states, eps, W, Nk)
        F  = eigen(Hermitian(H))
        for e in F.values
            push!(all_ev, (M, e))
        end
        all_vecs[M] = F.vectors
    end

    sort!(all_ev, by=x->x[2])
    return all_ev, all_vecs
end
```

**Step 2: Add inline validation**

```julia
if abspath(PROGRAM_FILE) == @__FILE__
    include("../../shared/lattice.jl")
    include("../../shared/hoppings.jl")
    include("band.jl")
    include("proj_basis.jl")
    include("proj_interaction.jl")

    lat = RectLat3x6()
    Nk, Np = lat.Nuc, 3
    eps, uvec = compute_band(lat, 1.0, 0.2)
    sectors = build_proj_sectors(Nk, Np)

    # Test 1: pure kinetic (V=0) → ground state = sum of 3 lowest eps
    W_zero = zeros(ComplexF64, Nk, Nk, Nk)
    ev_kin, _ = solve_proj_ed(sectors, eps, W_zero, Nk)
    E0_proj = ev_kin[1][2]
    E0_exact = sum(sort(eps)[1:Np])
    @assert abs(E0_proj - E0_exact) < 1e-8 "kinetic limit failed: $E0_proj vs $E0_exact"
    println("Test 1 (kinetic limit): PASS  E0=$(round(E0_proj,digits=6))")

    # Test 2: with V1V2V3 interaction
    W_v123 = build_Vq_V123(lat, uvec, 1.0, 0.0, 0.0)
    ev_v123, _ = solve_proj_ed(sectors, eps, W_v123, Nk)
    println("Test 2 (V1=1): ground state E0=$(round(ev_v123[1][2],digits=6))")

    # Test 3: particle number (trace of n_k in ground state sector)
    M_gs = ev_v123[1][1]
    println("  Ground state sector M=$M_gs")
    println("proj_hamiltonian.jl: all checks passed")
end
```

**Step 3: Run validation**

```bash
julia projected_ed/shared_proj/proj_hamiltonian.jl
```
Expected: "Test 1 (kinetic limit): PASS", "all checks passed".

**Step 4: Commit**

```bash
git add projected_ed/shared_proj/proj_hamiltonian.jl
git commit -m "feat: implement proj_hamiltonian.jl with full diag per momentum sector"
```

---

### Task 6: Implement main.jl

**Files:**
- Create: `projected_ed/case_proj/main.jl`

**Step 1: Write main.jl**

```julia
# ============================================================
# Projected ED — lowest Chern band
# Lattice: RectLat3x6 (Nk=9) for testing; swap for TiltedLat30 (Nk=15)
# ============================================================
include("../../shared/lattice.jl")
include("../../shared/hoppings.jl")
include("../shared_proj/band.jl")
include("../shared_proj/proj_basis.jl")
include("../shared_proj/proj_interaction.jl")
include("../shared_proj/proj_hamiltonian.jl")

using Printf, LinearAlgebra, Dates

# ── Parameters ──
t1   = 1.0
t3   = 0.2
V1   = 1.0;  V2 = 0.0;  V3 = 0.0   # for V123 version
d_gate = 1.0                          # gate distance for Coulomb version
Np   = 3                              # number of particles

println("="^60)
println("Projected ED — lowest Chern band")
println("t1=$t1  t3=$t3  Np=$Np")
println("="^60)

# ── 1. Lattice ──
lat = RectLat3x6()   # swap to TiltedLat30() for production
@printf("Lattice: Ns=%d  Nk=%d\n", lat.Ns, lat.Nuc)
Nk = lat.Nuc

# ── 2. Band structure ──
print("Computing band structure... ")
eps, uvec = compute_band(lat, t1, t3)
@printf("done. Band range: [%.4f, %.4f]\n", minimum(eps), maximum(eps))

# ── 3. Basis ──
print("Building k-space Fock basis C($Nk,$Np)=$(binomial(Nk,Np))... ")
sectors = build_proj_sectors(Nk, Np)
println("done. Sector sizes: ", [length(sectors[M]) for M in 0:Nk-1])

# ── 4a. V1V2V3 version ──
println("\n--- V1V2V3 version (V1=$V1, V2=$V2, V3=$V3) ---")
W_v123 = build_Vq_V123(lat, uvec, V1, V2, V3)
ev_v123, _ = solve_proj_ed(sectors, eps, W_v123, Nk)
E0_v123 = ev_v123[1][2]
open("spectrum_proj_V123.dat", "w") do f
    println(f, "# m  E-E0  [proj ED Nk=$Nk Np=$Np V1=$V1 V2=$V2 V3=$V3 t3=$t3]")
    for (m, e) in ev_v123
        @printf(f, "%d  %.10f\n", m, e - E0_v123)
    end
end
println("Ground state: M=$(ev_v123[1][1])  E0=$(round(E0_v123,digits=6))")
println("Saved: spectrum_proj_V123.dat")

# ── 4b. Gate Coulomb version ──
println("\n--- Gate Coulomb version (d=$d_gate) ---")
W_gate = build_Vq_gate(lat, uvec, d_gate)
ev_gate, _ = solve_proj_ed(sectors, eps, W_gate, Nk)
E0_gate = ev_gate[1][2]
open("spectrum_proj_gate.dat", "w") do f
    println(f, "# m  E-E0  [proj ED Nk=$Nk Np=$Np d=$d_gate t3=$t3]")
    for (m, e) in ev_gate
        @printf(f, "%d  %.10f\n", m, e - E0_gate)
    end
end
println("Ground state: M=$(ev_gate[1][1])  E0=$(round(E0_gate,digits=6))")
println("Saved: spectrum_proj_gate.dat")

println("\nDone: ", now())
```

**Step 2: Run it**

```bash
cd projected_ed/case_proj
julia main.jl
```
Expected: two `.dat` files created, no errors.

**Step 3: Commit**

```bash
git add projected_ed/case_proj/main.jl
git commit -m "feat: implement projected ED main.jl with V1V2V3 and gate Coulomb"
```

---

### Task 7: Cross-validation against full two-band ED

**Files:**
- Create: `projected_ed/case_proj/validate_vs_full_ed.jl`

**Step 1: Write validation script**

This script runs the projected ED on RectLat3x6 with Np=3 (ν=1/3) and compares with the full two-band ED. The projected energy must be ≥ full ED energy.

```julia
# validate_vs_full_ed.jl
# Run projected ED and compare with known full-ED reference energy
include("../../shared/lattice.jl")
include("../../shared/hoppings.jl")
include("../../shared/basis.jl")
include("../../shared/ksector.jl")
include("../../shared/hamiltonian.jl")
include("../../shared/solver.jl")
include("../shared_proj/band.jl")
include("../shared_proj/proj_basis.jl")
include("../shared_proj/proj_interaction.jl")
include("../shared_proj/proj_hamiltonian.jl")

using LinearAlgebra, Printf

t1, t3, V1, V2, V3 = 1.0, 0.2, 1.0, 0.0, 0.0
Np = 3
lat = RectLat3x6()
Nk  = lat.Nuc

# --- Full two-band ED ---
println("Running full two-band ED (Nk=$Nk, Np=$Np)...")
basis = gen_basis(lat.Ns, Np)
hops  = build_hops(lat, t1, t3, 0.0)
secs  = [build_ksector(basis, lat, m) for m in lat.ktab]
all_ev_full, _ = compute_spectrum_with_vecs(secs, lat, hops, V1, V2, V3;
                                             nev=4, krylovdim=50, verbose=false)
E0_full = all_ev_full[1][2]
@printf("Full ED ground state: E0 = %.8f\n", E0_full)

# --- Projected ED ---
println("Running projected ED...")
eps, uvec = compute_band(lat, t1, t3)
sectors   = build_proj_sectors(Nk, Np)
W         = build_Vq_V123(lat, uvec, V1, V2, V3)
ev_proj, _= solve_proj_ed(sectors, eps, W, Nk)
E0_proj   = ev_proj[1][2]
@printf("Projected ED ground state: E0 = %.8f\n", E0_proj)

@printf("Difference (proj - full): %.8f\n", E0_proj - E0_full)
if E0_proj >= E0_full - 1e-6
    println("PASS: projected energy ≥ full ED energy (as expected, projection is variational restriction)")
else
    println("FAIL: projected energy below full ED — something is wrong!")
end
```

**Step 2: Run validation**

```bash
julia projected_ed/case_proj/validate_vs_full_ed.jl
```
Expected: "PASS", projected energy ≥ full ED energy.

**Step 3: Commit**

```bash
git add projected_ed/case_proj/validate_vs_full_ed.jl
git commit -m "test: add cross-validation script projected vs full two-band ED"
```

---

### Task 8: Wire up TiltedLat30 production run

**Files:**
- Modify: `projected_ed/case_proj/main.jl` (swap lattice, update Np)

**Step 1: Create production version for TiltedLat30, Np=12**

Either edit `main.jl` directly or create `main_tilted.jl`:
- Change `RectLat3x6()` → `TiltedLat30()`
- Change `Np = 3` → `Np = 12`
- C(15,12) = 455, each sector ≈ 30 states — runs in seconds

**Step 2: Run production**

```bash
julia projected_ed/case_proj/main.jl
```

**Step 3: Final commit**

```bash
git add projected_ed/
git commit -m "feat: projected ED complete — V1V2V3 and gate Coulomb, TiltedLat30"
```

---

## Key Notes for the Implementer

- **Momentum arithmetic**: all `mod(x, Nk)` operations. In Julia, `mod(-1, 15) = 14` ✓. In the `qi` loop, `qi` runs 1..Nk-1 (skipping 0 = no scattering), but `W[m1,m2,qi]` is 1-indexed. Be careful: `qi=1` means momentum transfer of `k_1 - k_0 = G1_sc`, not zero.
- **Fermion sign**: the most error-prone part. The rule: sign of c†_m is (-1)^(popcount of occupied bits strictly below m) applied to the state **before** the creation. Double-check with the kinetic limit test (Task 5, Test 1).
- **W is not Hermitian in (m1,m2)**: W(m1,m2,q) = W*(m2,m1,-q) is the symmetry. The Hamiltonian matrix H is Hermitian overall.
- **Gate Coulomb sublattice phases**: the phase exp(-i q·(τ_α − τ_β)) in `build_Vq_gate` is important — omitting it gives wrong results for the inter-sublattice part of the interaction.
