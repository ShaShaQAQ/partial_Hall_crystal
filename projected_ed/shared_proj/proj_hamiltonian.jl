# ============================================================
# proj_hamiltonian.jl — Build and diagonalize projected H per momentum sector
# ============================================================
# Requires: proj_basis.jl and proj_interaction.jl included before this file.
using LinearAlgebra

"""
    build_proj_H(sector_states, eps, W, kadd, Nk) -> Matrix{ComplexF64}

Build the projected Hamiltonian matrix for a single momentum sector.

  sector_states: Vector{Int}          — Fock states (bitmasks) in this sector
  eps:           Vector{Float64}      — lower-band energies, eps[m+1] for k_m
  W:             Array{ComplexF64,3}  — W[m1+1, m2+1, qi+1], precomputed
  kadd:          Matrix{Int}          — kadd[m+1, q+1] = (m+q) mod BZ (0-based result)
  Nk:            Int

Returns: dim×dim Hermitian matrix.
"""
function build_proj_H(sector_states::Vector{Int},
                      eps::Vector{Float64},
                      W::Array{ComplexF64,3},
                      kadd::Matrix{Int},
                      Nk::Int)
    dim = length(sector_states)
    H   = zeros(ComplexF64, dim, dim)

    state2idx = Dict{Int,Int}(s => i for (i,s) in enumerate(sector_states))

    # Precompute k_neg table once per call
    k_neg = [_k_neg_qi(kadd, qi, Nk) for qi in 0:Nk-1]  # k_neg[qi+1] = -qi (0-based)

    for (i, si) in enumerate(sector_states)
        # --- diagonal: kinetic energy ---
        Ekin = 0.0
        for m in 0:Nk-1
            (si >> m) & 1 == 1 && (Ekin += eps[m+1])
        end
        H[i, i] += Ekin

        # --- diagonal: q=0 Hartree interaction ---
        # H += (1/Nk) * Σ_{m1<m2, both occupied} W[m1,m2,q=0]
        # (the (m1,m2,0)+(m2,m1,0) pair together give 2×, cancelled by summing m1<m2 only)
        for m1 in 0:Nk-1
            (si >> m1) & 1 == 0 && continue
            for m2 in m1+1:Nk-1
                (si >> m2) & 1 == 0 && continue
                H[i, i] += (1.0/Nk) * real(W[m1+1, m2+1, 1])
            end
        end

        # --- off-diagonal: c†_{m1+q} c†_{m2-q} c_{m2} c_{m1} for q≠0 ---
        # Factor 0.5/Nk: loop over ALL (m1,m2,q) double-counts each scattering event
        # (the pair (m1,m2,q) and (m2,m1,-q) both map si→sf with equal contribution).
        for m1 in 0:Nk-1
            (si >> m1) & 1 == 0 && continue
            for m2 in 0:Nk-1
                m2 == m1 && continue
                (si >> m2) & 1 == 0 && continue
                for qi in 1:Nk-1   # 0-based q-label, skip q=0
                    kp1 = kadd[m1+1, qi+1]             # k1 + q  (0-based)
                    kp2 = kadd[m2+1, k_neg[qi+1]+1]   # k2 - q  (0-based); k_neg is 1-indexed

                    # Pauli blocking
                    (si >> kp1) & 1 == 1 && continue
                    (si >> kp2) & 1 == 1 && continue
                    kp1 == kp2 && continue

                    sf = si ⊻ (1 << m1) ⊻ (1 << m2) ⊻ (1 << kp1) ⊻ (1 << kp2)
                    haskey(state2idx, sf) || continue
                    j = state2idx[sf]

                    # Fermion signs for c†_{kp1} c†_{kp2} c_{m2} c_{m1}
                    sgn_m1  = kfock_fermion_sign(si,   m1)
                    tmp1    = si ⊻ (1 << m1)
                    sgn_m2  = kfock_fermion_sign(tmp1, m2)
                    tmp2    = tmp1 ⊻ (1 << m2)
                    sgn_kp2 = kfock_fermion_sign(tmp2, kp2)
                    tmp3    = tmp2 | (1 << kp2)
                    sgn_kp1 = kfock_fermion_sign(tmp3, kp1)
                    total_sgn = sgn_m1 * sgn_m2 * sgn_kp2 * sgn_kp1

                    H[j, i] += (0.5/Nk) * W[m1+1, m2+1, qi+1] * total_sgn
                end
            end
        end
    end

    return H
end

# Helper: return 0-based index of -q given 0-based q index qi
@inline function _k_neg_qi(kadd::Matrix{Int}, qi::Int, Nk::Int)
    for kn in 0:Nk-1
        kadd[qi+1, kn+1] == 0 && return kn
    end
    error("k_neg not found for qi=$qi")
end

"""
    solve_proj_ed(sectors, eps, W, kadd, Nk) -> (all_ev, all_vecs)

Diagonalize all momentum sectors. Returns:
  all_ev  : Vector{Tuple{Int,Float64}} — (M, energy) pairs, sorted by energy
  all_vecs: Dict{Int, Matrix} — M => eigenvector matrix
"""
function solve_proj_ed(sectors::Dict{Int,Vector{Int}},
                       eps::Vector{Float64},
                       W::Array{ComplexF64,3},
                       kadd::Matrix{Int},
                       Nk::Int)
    all_ev   = Tuple{Int,Float64}[]
    all_vecs = Dict{Int, Matrix{ComplexF64}}()

    for M in sort(collect(keys(sectors)))
        states = sectors[M]
        isempty(states) && continue
        H = build_proj_H(states, eps, W, kadd, Nk)
        F = eigen(Hermitian(H))
        for e in F.values
            push!(all_ev, (M, e))
        end
        all_vecs[M] = F.vectors
    end

    sort!(all_ev, by=x->x[2])
    return all_ev, all_vecs
end

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
    kadd = build_k_add(lat)

    # Test 1: pure kinetic (V=0) → ground state = sum of 3 lowest eps
    W_zero = zeros(ComplexF64, Nk, Nk, Nk)
    ev_kin, _ = solve_proj_ed(sectors, eps, W_zero, kadd, Nk)
    E0_proj  = ev_kin[1][2]
    E0_exact = sum(sort(eps)[1:Np])
    @assert abs(E0_proj - E0_exact) < 1e-8 "kinetic limit failed: $E0_proj vs $E0_exact"
    println("Test 1 (kinetic limit): PASS  E0=$(round(E0_proj, digits=6))")

    # Test 2: with V1 interaction
    W_v123 = build_Vq_V123(lat, uvec, 1.0, 0.0, 0.0)
    ev_v123, _ = solve_proj_ed(sectors, eps, W_v123, kadd, Nk)
    println("Test 2 (V1=1): E0=$(round(ev_v123[1][2], digits=6))  sector M=$(ev_v123[1][1])")

    # Test 3: Hermitian check on H matrix in one sector
    states_0 = sectors[0]
    H_test = build_proj_H(states_0, eps, W_v123, kadd, Nk)
    herm_err = maximum(abs.(H_test .- H_test'))
    @assert herm_err < 1e-10 "H not Hermitian: $herm_err"
    println("Test 3 (Hermitian H): PASS  err=$herm_err")

    println("proj_hamiltonian.jl: all checks passed")
end
