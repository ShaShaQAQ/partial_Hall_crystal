# ============================================================
# proj_interaction.jl — Projected two-body matrix elements
# ============================================================
# Requires: lattice.jl, hoppings.jl, band.jl included before this file.
using LinearAlgebra

if abspath(PROGRAM_FILE) == @__FILE__
    include("../../shared/lattice.jl")
    include("../../shared/hoppings.jl")
    include("band.jl")
end

const SUBLAT_CART_PROJ = [[0.0, 0.0], [-0.5, sqrt(3)/2]]  # τ_A, τ_B

# Triangular-coordinate displacement vectors for nb1/nb2/nb3
# Each (dx, dy) gives Cartesian displacement Δr = dx*a1 + dy*a2
const _D1 = NTuple{2,Int}[(1,0),(-1,0),(1,1),(-1,-1),(0,1),(0,-1)]
const _D2 = NTuple{2,Int}[(2,1),(-2,-1),(1,2),(-1,-2),(1,-1),(-1,1)]
const _D3 = NTuple{2,Int}[(2,0),(-2,0),(0,2),(0,-2),(2,2),(-2,-2)]

"""
    _site_sublat(lat) -> Vector{Int}

Return sublattice index (0=A, 1=B) for each site.
"""
function _site_sublat(lat::GenLat)
    return [((iy % 2 + 2) % 2) for (ix, iy) in lat.sites]
end

"""
    build_k_add(lat) -> Matrix{Int}

Precompute momentum addition table (0-based):
  k_add[m1+1, qi+1] = kp1  such that k_{kp1} ≡ k_{m1} + k_{qi} (mod BZ).
Uses phase_table for a robust match that works for any GenLat geometry.
"""
function build_k_add(lat::GenLat)
    Nk = lat.Nuc
    k_add = zeros(Int, Nk, Nk)
    for m1 in 0:Nk-1, qi in 0:Nk-1
        target = lat.phase_table[m1+1, :] .* lat.phase_table[qi+1, :]
        best = 0
        best_err = Inf
        for m3 in 0:Nk-1
            err = sum(abs2.(lat.phase_table[m3+1, :] .- target))
            if err < best_err
                best_err = err
                best = m3
            end
        end
        @assert best_err < 1e-10 "k_add failed m1=$m1 qi=$qi err=$best_err"
        k_add[m1+1, qi+1] = best
    end
    return k_add
end

function _form_factors(uvec::Vector{Vector{ComplexF64}}, Nk::Int)
    F = zeros(ComplexF64, 2, Nk, Nk)
    for k in 1:Nk, kp in 1:Nk
        for α in 1:2
            F[α, kp, k] = conj(uvec[kp][α]) * uvec[k][α]
        end
    end
    return F
end

"""
    _build_Vq_nb(lat, nb, ds, Vval) -> Array{ComplexF64,3}

Sublattice-resolved Fourier transform for one neighbor shell.
Displacement Δr = dx*a1 + dy*a2 (no periodic-boundary issues).
Vq[αi+1, αj+1, qi] += Vval * exp(-i q·Δr) / 2
"""
function _build_Vq_nb(lat::GenLat, nb::Matrix{Int},
                      ds::Vector{NTuple{2,Int}}, Vval::Float64)
    Nk         = lat.Nuc
    kpts       = lat.kpoints
    site_alpha = _site_sublat(lat)
    a1c        = collect(lat.a1)
    a2c        = collect(lat.a2)
    Vq         = zeros(ComplexF64, 2, 2, Nk)
    for si in 1:lat.Ns
        αi = site_alpha[si]
        for (k, (dx, dy)) in enumerate(ds)
            sj  = nb[si, k]
            αj  = site_alpha[sj]
            Δr  = dx .* a1c .+ dy .* a2c   # exact Cartesian displacement, no wrap
            for qi in 1:Nk
                Vq[αi+1, αj+1, qi] += Vval * exp(-1im * dot(kpts[qi], Δr)) * 0.5
            end
        end
    end
    return Vq
end

"""
    build_Vq_V123(lat, uvec, V1, V2, V3) -> Array{ComplexF64,3}

W[m1, m2, qi] (all 1-based). qi=1 is q=k_0 (zero momentum transfer).
"""
function build_Vq_V123(lat::GenLat, uvec::Vector{Vector{ComplexF64}},
                        V1::Float64, V2::Float64, V3::Float64)
    Nk  = lat.Nuc
    Vq  = zeros(ComplexF64, 2, 2, Nk)
    V1 != 0.0 && (Vq .+= _build_Vq_nb(lat, lat.nb1, collect(_D1), V1))
    V2 != 0.0 && (Vq .+= _build_Vq_nb(lat, lat.nb2, collect(_D2), V2))
    V3 != 0.0 && (Vq .+= _build_Vq_nb(lat, lat.nb3, collect(_D3), V3))

    F     = _form_factors(uvec, Nk)
    kadd  = build_k_add(lat)                       # kadd[m+1, q+1] = (m+q)%BZ  (0-based result)
    kneg  = [kadd[1, qi+1] for qi in 0:Nk-1]      # kneg[qi+1] = (-k_qi)%BZ index (0-based)
    # kneg: -k_{qi} = k_add[k=0, qi_neg], but actually we need: 0 + kneg = 0 → kadd[0,kneg]=0
    # simpler: k_neg[qi] = k such that kadd[qi+1, k+1] == 0
    k_neg_tab = zeros(Int, Nk)
    for qi in 0:Nk-1
        for kn in 0:Nk-1
            kadd[qi+1, kn+1] == 0 && (k_neg_tab[qi+1] = kn; break)
        end
    end

    W = zeros(ComplexF64, Nk, Nk, Nk)
    for m1 in 0:Nk-1, m2 in 0:Nk-1, qi in 0:Nk-1
        kp1 = kadd[m1+1, qi+1]            # k1 + q  (0-based)
        qi_neg = k_neg_tab[qi+1]           # -q  (0-based)
        kp2 = kadd[m2+1, qi_neg+1]        # k2 - q = k2 + (-q)  (0-based)
        val = zero(ComplexF64)
        for α in 1:2, β in 1:2
            val += Vq[α, β, qi+1] * F[α, kp1+1, m1+1] * F[β, kp2+1, m2+1]
        end
        W[m1+1, m2+1, qi+1] = val
    end
    return W
end

"""
    build_Vq_gate(lat, uvec, d) -> Array{ComplexF64,3}

Dual-gate screened Coulomb: V(q) = (2π/|q|) × tanh(|q|×d), V(q=0)=0.
"""
function build_Vq_gate(lat::GenLat, uvec::Vector{Vector{ComplexF64}}, d::Float64)
    Nk   = lat.Nuc
    kpts = lat.kpoints

    Vscalar = [norm(kpts[qi]) < 1e-10 ? 0.0 :
               (2π / norm(kpts[qi])) * tanh(norm(kpts[qi]) * d)
               for qi in 1:Nk]

    F    = _form_factors(uvec, Nk)
    kadd = build_k_add(lat)
    k_neg_tab = zeros(Int, Nk)
    for qi in 0:Nk-1
        for kn in 0:Nk-1
            kadd[qi+1, kn+1] == 0 && (k_neg_tab[qi+1] = kn; break)
        end
    end

    W = zeros(ComplexF64, Nk, Nk, Nk)
    for m1 in 0:Nk-1, m2 in 0:Nk-1, qi in 0:Nk-1
        kp1    = kadd[m1+1, qi+1]
        qi_neg = k_neg_tab[qi+1]
        kp2    = kadd[m2+1, qi_neg+1]
        q      = kpts[qi+1]
        Vq     = Vscalar[qi+1]
        val    = zero(ComplexF64)
        for α in 1:2, β in 1:2
            τ_diff = SUBLAT_CART_PROJ[α] .- SUBLAT_CART_PROJ[β]
            phase  = exp(-1im * dot(q, τ_diff))
            val   += Vq * phase * F[α, kp1+1, m1+1] * F[β, kp2+1, m2+1]
        end
        W[m1+1, m2+1, qi+1] = val
    end
    return W
end

if abspath(PROGRAM_FILE) == @__FILE__
    lat = RectLat3x6()
    eps, uvec = compute_band(lat, 1.0, 0.2)
    Nk = lat.Nuc

    W_v123 = build_Vq_V123(lat, uvec, 1.0, 0.0, 0.0)
    W_gate = build_Vq_gate(lat, uvec, 1.0)

    @assert all(isfinite, W_v123) "W_v123 contains NaN/Inf"
    @assert all(isfinite, W_gate) "W_gate contains NaN/Inf"

    # Hermiticity: W[m1+q, m2-q, -q] = conj(W[m1, m2, q])
    kadd_t = build_k_add(lat)
    k_neg_t = zeros(Int, Nk)
    for qi in 0:Nk-1
        for kn in 0:Nk-1
            kadd_t[qi+1, kn+1] == 0 && (k_neg_t[qi+1] = kn; break)
        end
    end
    local max_err = 0.0
    for m1 in 0:Nk-1, m2 in 0:Nk-1, qi in 0:Nk-1
        m1p    = kadd_t[m1+1, qi+1]
        qi_neg = k_neg_t[qi+1]
        m2p    = kadd_t[m2+1, qi_neg+1]
        err    = abs(W_v123[m1p+1, m2p+1, qi_neg+1] - conj(W_v123[m1+1, m2+1, qi+1]))
        max_err = max(max_err, err)
    end
    println("W_v123 Hermiticity error: ", max_err)
    @assert max_err < 1e-10 "W_v123 not Hermitian: $max_err"

    println("proj_interaction.jl: all checks passed")
end
