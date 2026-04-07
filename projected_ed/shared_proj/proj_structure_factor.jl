# ============================================================
# proj_structure_factor.jl — Structure factor for projected ED
# ============================================================
# Requires: proj_basis.jl (kfock_fermion_sign), proj_interaction.jl
#           (build_k_add) included before this file.
#
# Computes the projected density structure factor:
#
#   N_proj(q) = (1/Nk) ||ρ_q^proj |ψ_GS>||²
#
# where the projected density operator is
#
#   ρ_q^proj = Σ_k F(k+q, k) c†_{k+q} c_k
#
# with band form factor F(k', k) = Σ_α conj(u_{k'α}) u_{kα}.
#
# This is the lowest-band projection of the real-space density
# Fourier component, so N_proj(q) is directly comparable to the
# full-ED structure factor N_full(q).
# ============================================================

"""
    build_form_factors(uvec, Nk) -> Matrix{ComplexF64}

F[kp+1, k+1] = Σ_α conj(u_{kp,α}) u_{k,α}   (band form factor)
"""
function build_form_factors(uvec::Vector{Vector{ComplexF64}}, Nk::Int)
    nb = length(uvec[1])   # number of orbitals (2 for honeycomb)
    F  = zeros(ComplexF64, Nk, Nk)
    for k in 0:Nk-1, kp in 0:Nk-1
        for α in 1:nb
            F[kp+1, k+1] += conj(uvec[kp+1][α]) * uvec[k+1][α]
        end
    end
    return F
end

"""
    compute_proj_sq(gs_vec, sector_states, kadd, Nk, uvec) -> Vector{Float64}

Compute the projected structure factor N_proj(q) for all Nk BZ momenta.

  gs_vec:         length-dim eigenvector (coefficients over sector_states)
  sector_states:  Fock states for the GS momentum sector
  kadd:           kadd[m+1, qi+1] = (m+qi) mod BZ  (0-based)
  Nk:             number of k-points
  uvec:           uvec[k+1] = lower-band spinor at k_m

Returns sq[qi+1] = N_proj(k_{qi}) for qi = 0,...,Nk-1.
sq[1] (q=0) is set to 0.0 by convention (charge conservation).
"""
function compute_proj_sq(gs_vec::Vector{ComplexF64},
                          sector_states::Vector{Int},
                          kadd::Matrix{Int},
                          Nk::Int,
                          uvec::Vector{Vector{ComplexF64}})
    F_ff = build_form_factors(uvec, Nk)
    sq   = zeros(Float64, Nk)

    for qi in 1:Nk-1     # skip qi=0
        # Apply ρ_q^proj = Σ_k F[k+q,k] c†_{k+q} c_k to |ψ>
        # Accumulate output in a sparse dict (key = Fock bitmask)
        out = Dict{Int, ComplexF64}()

        for (idx, s) in enumerate(sector_states)
            α_c = gs_vec[idx]
            abs2(α_c) < 1e-30 && continue

            for k in 0:Nk-1
                (s >> k) & 1 == 0 && continue    # k must be occupied
                kpq = kadd[k+1, qi+1]             # k + q  (0-based)
                (s >> kpq) & 1 == 1 && continue   # k+q must be unoccupied

                # Apply c†_{kpq} c_k: annihilate k, then create kpq
                sgn_k   = kfock_fermion_sign(s,               k)
                tmp     = s ⊻ (1 << k)
                sgn_kpq = kfock_fermion_sign(tmp,             kpq)
                t       = tmp | (1 << kpq)

                val     = α_c * F_ff[kpq+1, k+1] * (sgn_k * sgn_kpq)
                out[t]  = get(out, t, zero(ComplexF64)) + val
            end
        end

        sq[qi+1] = isempty(out) ? 0.0 : sum(abs2(v) for v in values(out)) / Nk
    end
    # sq[1] remains 0.0  (N_proj(q=0) = 0 by charge conservation)
    return sq
end

"""
    compute_proj_sq_allsectors(all_vecs, sectors, kadd, Nk, uvec)
        -> Matrix{Float64}

Compute N_proj(q) for ALL momentum sectors' lowest eigenstates.
Returns sq_all[M+1, qi+1] for M in 0:Nk-1, qi in 0:Nk-1.
"""
function compute_proj_sq_allsectors(all_vecs::Dict{Int, Matrix{ComplexF64}},
                                     sectors::Dict{Int, Vector{Int}},
                                     kadd::Matrix{Int},
                                     Nk::Int,
                                     uvec::Vector{Vector{ComplexF64}})
    sq_all = zeros(Float64, Nk, Nk)
    for M in 0:Nk-1
        haskey(all_vecs, M) || continue
        isempty(sectors[M]) && continue
        vec_M = all_vecs[M][:, 1]    # lowest eigenvector in sector M
        sq_all[M+1, :] = compute_proj_sq(vec_M, sectors[M], kadd, Nk, uvec)
    end
    return sq_all
end
