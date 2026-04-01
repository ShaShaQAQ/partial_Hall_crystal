# ============================================================
# band.jl — Single-particle band energies and Bloch eigenstates
# ============================================================
using LinearAlgebra

if abspath(PROGRAM_FILE) == @__FILE__
    include("../../shared/lattice.jl")
    include("../../shared/hoppings.jl")
end

"""
    compute_band(lat, t1, t3) -> (eps, uvec)

Diagonalize the 2×2 H(k) at each of the Nk k-points.
Returns:
  eps  :: Vector{Float64}             length Nk, lower-band eigenvalues
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

if abspath(PROGRAM_FILE) == @__FILE__
    lat = RectLat3x6()
    eps, uvec = compute_band(lat, 1.0, 0.2)
    println("Lower band energies: ", round.(eps, digits=4))
    for (m, u) in enumerate(uvec)
        @assert abs(norm(u) - 1.0) < 1e-12 "eigenvector not normalized at k=$m"
        @assert abs(imag(u[1])) < 1e-12 "gauge fix failed: u_A not real at k=$m"
    end
    println("band.jl: all checks passed")
end
