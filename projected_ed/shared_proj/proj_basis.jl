# ============================================================
# proj_basis.jl — k-space Fock states and momentum sectors
# ============================================================
# Requires: no external dependencies

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

if abspath(PROGRAM_FILE) == @__FILE__
    Nk, Np = 15, 12
    sectors = build_proj_sectors(Nk, Np)
    total = sum(length(v) for v in values(sectors))
    @assert total == binomial(Nk, Np) "total states wrong: got $total expected $(binomial(Nk,Np))"
    sizes = [length(sectors[M]) for M in 0:Nk-1]
    println("Sector sizes (Nk=$Nk, Np=$Np): ", sizes)

    # fermion sign check: state = 0b101 (bits 0 and 2 occupied)
    state = Int(0b101)
    @assert kfock_fermion_sign(state, 0) == 1   "sign at bit 0 wrong"
    @assert kfock_fermion_sign(state, 2) == -1  "sign at bit 2 wrong (one occupied below)"
    println("proj_basis.jl: all checks passed")
end
