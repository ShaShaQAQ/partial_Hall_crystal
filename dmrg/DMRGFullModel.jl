module DMRGFullModel

using LinearAlgebra
using Printf
using Random
using Statistics

using KrylovKit
using ITensors
using ITensorMPS
import ITensorMPS: checkdone!, measure!, truncerror

include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")

export CylinderLat,
    CylinderModelParams,
    HoppingPath,
    HoppingTerm,
    InteractionTerm,
    ConvergenceRecord,
    DMRGRunResult,
    cylinder_lattice,
    filling_np,
    site_index,
    site_xy,
    site_phys_pos,
    neighbor_bonds,
    interaction_terms,
    build_cylinder_hopping_paths,
    build_cylinder_hoppings,
    build_model_terms,
    full_ed_ground_energy,
    build_mpo,
    convergence_satisfied,
    run_dmrg,
    density_profile,
    y_averaged_density,
    green_function,
    connected_density_correlation,
    cumulative_charge,
    write_density,
    write_convergence,
    write_vector,
    write_matrix,
    write_complex_matrix

struct CylinderLat
    Lx::Int
    Ly::Int
    Ly_uc::Int
    Ns::Int
    Nuc::Int
    a1::NTuple{2,Float64}
    a2::NTuple{2,Float64}
    b1::NTuple{2,Float64}
    b2::NTuple{2,Float64}
    sites::Vector{NTuple{2,Int}}
    site_idx::Dict{NTuple{2,Int},Int}
end

struct CylinderModelParams
    t1::Float64
    t3::Float64
    V1::Float64
    V2::Float64
    V3::Float64
end

CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.0, V3=0.0) =
    CylinderModelParams(Float64(t1), Float64(t3), Float64(V1), Float64(V2), Float64(V3))

struct HoppingPath
    target::Int
    source::Int
    amp0::ComplexF64
    winding_y::Int
    dx_uc::Int
    dy_uc::Int
end

struct HoppingTerm
    target::Int
    source::Int
    amp::ComplexF64
end

struct InteractionTerm
    i::Int
    j::Int
    V::Float64
    shell::Int
end

struct ConvergenceRecord
    sweep::Int
    energy::Float64
    delta_energy::Float64
    max_density_delta::Float64
    maxlinkdim::Int
    maxerr::Float64
    converged::Bool
    reason::String
end

struct DMRGRunResult
    energy::Float64
    psi::MPS
    sites::Vector{Index}
    density::Vector{Float64}
    convergence::Vector{ConvergenceRecord}
end

const NB_DISPS = Dict(
    1 => [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, 1), (1, -1)],
    2 => [(1, 1), (-1, -1), (-1, 2), (1, -2), (-2, 1), (2, -1)],
    3 => [(2, 0), (-2, 0), (0, 2), (0, -2), (-2, 2), (2, -2)],
)

function reciprocal_vectors(a1::NTuple{2,Float64}, a2::NTuple{2,Float64})
    area = a1[1] * a2[2] - a1[2] * a2[1]
    b1 = (2π * a2[2] / area, -2π * a2[1] / area)
    b2 = (-2π * a1[2] / area, 2π * a1[1] / area)
    return b1, b2
end

function cylinder_lattice(Lx::Int; Ly::Int=6)
    Lx > 0 || throw(ArgumentError("Lx must be positive"))
    Ly > 0 || throw(ArgumentError("Ly must be positive"))
    iseven(Ly) || throw(ArgumentError("Ly must be even because each unit cell has two sites"))
    a1 = (1.0, 0.0)
    a2 = (0.5, sqrt(3) / 2)
    b1, b2 = reciprocal_vectors(a1, a2)
    sites = [(x, y) for x in 0:(Lx - 1) for y in 0:(Ly - 1)]
    site_idx = Dict(s => i for (i, s) in enumerate(sites))
    Ns = Lx * Ly
    return CylinderLat(Lx, Ly, Ly ÷ 2, Ns, Ns ÷ 2, a1, a2, b1, b2, sites, site_idx)
end

filling_np(lat::CylinderLat; numerator::Int=1, denominator::Int=3) =
    numerator * lat.Nuc ÷ denominator

function site_index(lat::CylinderLat, x::Int, y::Int)
    0 <= x < lat.Lx || throw(BoundsError(lat.sites, (x, y)))
    return lat.site_idx[(x, mod(y, lat.Ly))]
end

site_xy(lat::CylinderLat, i::Int) = lat.sites[i]

function site_phys_pos(lat::CylinderLat, x::Int, y::Int)
    rx = Float64(x) * lat.a1[1] + Float64(y) * lat.a2[1]
    ry = Float64(x) * lat.a1[2] + Float64(y) * lat.a2[2]
    return rx, ry
end

function neighbor_bonds(lat::CylinderLat, shell::Int)
    haskey(NB_DISPS, shell) || throw(ArgumentError("unknown neighbor shell $shell"))
    bonds = Set{Tuple{Int,Int}}()
    for (i, (x, y)) in enumerate(lat.sites)
        for (dx, dy) in NB_DISPS[shell]
            tx = x + dx
            0 <= tx < lat.Lx || continue
            ty = mod(y + dy, lat.Ly)
            j = lat.site_idx[(tx, ty)]
            i == j && continue
            push!(bonds, minmax(i, j))
        end
    end
    return sort!(collect(bonds))
end

function interaction_terms(lat::CylinderLat, params::CylinderModelParams)
    terms = InteractionTerm[]
    for (shell, V) in ((1, params.V1), (2, params.V2), (3, params.V3))
        abs(V) < 1e-14 && continue
        for (i, j) in neighbor_bonds(lat, shell)
            push!(terms, InteractionTerm(i, j, V, shell))
        end
    end
    return terms
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
    b1, b2 = reciprocal_vectors((a1[1], a1[2]), (a2[1], a2[2]))
    b1v = collect(b1)
    b2v = collect(b2)
    sublat = [[0.0, 0.0], copy(a2)]
    out = Dict{NTuple{2,Int},Matrix{ComplexF64}}()
    Nk = nx_aux * ny_aux
    for n1 in -max_dx:max_dx, n2 in -max_dy:max_dy
        R = n1 .* a1 .+ n2 .* (2 .* a2)
        tmat = zeros(ComplexF64, 2, 2)
        for m1 in 0:(nx_aux - 1), m2 in 0:(ny_aux - 1)
            k = (m1 / nx_aux) .* b1v .+ (m2 / (2 * ny_aux)) .* b2v
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

function build_cylinder_hopping_paths(lat::CylinderLat, t1::Float64, t3::Float64; kwargs...)
    tR = hopping_matrices_from_get_Hk(t1, t3; kwargs...)
    paths = HoppingPath[]
    for (source, (x, y)) in enumerate(lat.sites)
        source_orb0 = mod(y, 2)
        source_cell_y = (y - source_orb0) ÷ 2
        for ((dx_uc, dy_uc), tmat) in tR
            tx = x + dx_uc
            0 <= tx < lat.Lx || continue
            target_cell_y_unwrapped = source_cell_y + dy_uc
            winding_y = fld(target_cell_y_unwrapped, lat.Ly_uc)
            target_cell_y = mod(target_cell_y_unwrapped, lat.Ly_uc)
            for target_orb in 1:2
                amp0 = tmat[target_orb, source_orb0 + 1]
                abs(amp0) < 1e-12 && continue
                ty = 2 * target_cell_y + (target_orb - 1)
                target = site_index(lat, tx, ty)
                push!(paths, HoppingPath(target, source, amp0, winding_y, dx_uc, dy_uc))
            end
        end
    end
    return paths
end

function build_cylinder_hoppings(lat::CylinderLat, t1::Float64, t3::Float64, phi_y::Float64=0.0; kwargs...)
    paths = build_cylinder_hopping_paths(lat, t1, t3; kwargs...)
    accum = Dict{Tuple{Int,Int},ComplexF64}()
    for p in paths
        key = (p.target, p.source)
        accum[key] = get(accum, key, 0.0 + 0.0im) + p.amp0 * cis(p.winding_y * phi_y)
    end
    terms = HoppingTerm[]
    for ((target, source), amp) in accum
        abs(amp) < 1e-12 && continue
        push!(terms, HoppingTerm(target, source, amp))
    end
    sort!(terms, by=t -> (t.source, t.target, real(t.amp), imag(t.amp)))
    return terms
end

function build_model_terms(lat::CylinderLat, params::CylinderModelParams; phi_y::Float64=0.0)
    hops = build_cylinder_hoppings(lat, params.t1, params.t3, phi_y)
    ints = interaction_terms(lat, params)
    return hops, ints
end

@inline occupied(state::Int64, site::Int) = ((state >> (site - 1)) & 1) == 1

function interaction_energy(state::Int64, interactions::Vector{InteractionTerm})
    E = 0.0
    for term in interactions
        occupied(state, term.i) && occupied(state, term.j) && (E += term.V)
    end
    return E
end

function apply_full_hamiltonian!(
    w::Vector{ComplexF64},
    v::Vector{ComplexF64},
    basis::Vector{Int64},
    index::Dict{Int64,Int},
    hops::Vector{HoppingTerm},
    interactions::Vector{InteractionTerm},
)
    fill!(w, 0.0 + 0.0im)
    for (bi, F) in enumerate(basis)
        vi = v[bi]
        abs2(vi) < 1e-28 && continue
        w[bi] += interaction_energy(F, interactions) * vi
        for h in hops
            if h.target == h.source
                occupied(F, h.source) && (w[bi] += h.amp * vi)
                continue
            end
            occupied(F, h.source) || continue
            occupied(F, h.target) && continue
            newF, sgn = hop(F, h.target, h.source)
            j = index[newF]
            w[j] += h.amp * sgn * vi
        end
    end
    return w
end

function full_ed_ground_energy(
    lat::CylinderLat,
    params::CylinderModelParams,
    Np::Int;
    phi_y::Float64=0.0,
    krylovdim::Int=60,
    tol::Float64=1e-10,
    maxiter::Int=400,
    seed::Int=17,
)
    hops, ints = build_model_terms(lat, params; phi_y)
    basis = gen_basis(lat.Ns, Np)
    index = Dict(F => i for (i, F) in enumerate(basis))
    rng = MersenneTwister(seed)
    v0 = randn(rng, ComplexF64, length(basis))
    v0 ./= norm(v0)
    function Hf(v)
        w = similar(v)
        apply_full_hamiltonian!(w, v, basis, index, hops, ints)
        return w
    end
    vals, vecs, info = eigsolve(Hf, v0, 1, :SR; ishermitian=true, krylovdim, tol, maxiter)
    info.converged < 1 && @warn "ED did not converge" info
    return real(vals[1]), vecs[1]
end

function build_mpo(
    lat::CylinderLat,
    params::CylinderModelParams;
    phi_y::Float64=0.0,
    sites=nothing,
)
    sites === nothing && (sites = siteinds("Fermion", lat.Ns; conserve_qns=true))
    hops, ints = build_model_terms(lat, params; phi_y)
    os = OpSum()
    for h in hops
        if h.target == h.source
            os += h.amp, "N", h.source
        else
            os += h.amp, "Cdag", h.target, "C", h.source
        end
    end
    for term in ints
        os += term.V, "N", term.i, "N", term.j
    end
    return MPO(ComplexF64, os, sites), sites
end

function convergence_satisfied(
    sweep::Int,
    delta_energy::Real,
    max_density_delta::Real,
    maxerr::Real;
    min_sweeps::Int,
    target_sweep::Int=1,
    stable_sweeps::Int=1,
    satisfied_streak::Int=1,
    energy_tol::Real,
    density_tol::Real,
    truncerr_tol::Real,
)
    return sweep >= min_sweeps &&
        sweep >= target_sweep &&
        isfinite(delta_energy) &&
        isfinite(max_density_delta) &&
        delta_energy <= energy_tol &&
        max_density_delta <= density_tol &&
        maxerr <= truncerr_tol &&
        satisfied_streak >= stable_sweeps
end

mutable struct DMRGConvergenceObserver <: AbstractObserver
    min_sweeps::Int
    target_sweep::Int
    stable_sweeps::Int
    energy_tol::Float64
    density_tol::Float64
    truncerr_tol::Float64
    previous_energy::Union{Nothing,Float64}
    previous_density::Union{Nothing,Vector{Float64}}
    current_sweep::Int
    current_maxerr::Float64
    satisfied_streak::Int
    records::Vector{ConvergenceRecord}
end

function DMRGConvergenceObserver(;
    min_sweeps::Int,
    target_sweep::Int,
    stable_sweeps::Int,
    energy_tol::Real,
    density_tol::Real,
    truncerr_tol::Real,
)
    min_sweeps > 0 || throw(ArgumentError("min_sweeps must be positive"))
    target_sweep > 0 || throw(ArgumentError("target_sweep must be positive"))
    stable_sweeps > 0 || throw(ArgumentError("stable_sweeps must be positive"))
    energy_tol >= 0 || throw(ArgumentError("energy_tol must be non-negative"))
    density_tol >= 0 || throw(ArgumentError("density_tol must be non-negative"))
    truncerr_tol >= 0 || throw(ArgumentError("truncerr_tol must be non-negative"))
    return DMRGConvergenceObserver(
        min_sweeps,
        target_sweep,
        stable_sweeps,
        Float64(energy_tol),
        Float64(density_tol),
        Float64(truncerr_tol),
        nothing,
        nothing,
        0,
        0.0,
        0,
        ConvergenceRecord[],
    )
end

function measure!(obs::DMRGConvergenceObserver; kwargs...)
    sweep = kwargs[:sweep]
    if sweep != obs.current_sweep
        obs.current_sweep = sweep
        obs.current_maxerr = 0.0
    end
    if haskey(kwargs, :spec)
        obs.current_maxerr = max(obs.current_maxerr, truncerror(kwargs[:spec]))
    end
    return nothing
end

function convergence_reason(
    sweep::Int,
    delta_energy::Real,
    max_density_delta::Real,
    maxerr::Real,
    next_streak::Int,
    obs::DMRGConvergenceObserver,
)
    sweep < obs.min_sweeps && return "min_sweeps"
    sweep < obs.target_sweep && return "target_sweep"
    !isfinite(delta_energy) && return "first_sweep"
    delta_energy > obs.energy_tol && return "energy"
    max_density_delta > obs.density_tol && return "density"
    maxerr > obs.truncerr_tol && return "truncerr"
    next_streak < obs.stable_sweeps && return "stable_sweeps"
    return "energy_density_truncerr"
end

function checkdone!(obs::DMRGConvergenceObserver; energy, psi, sweep, outputlevel=false)
    e = real(energy)
    dens = density_profile(psi)
    delta_energy = obs.previous_energy === nothing ? Inf : abs(e - obs.previous_energy)
    max_density_delta = obs.previous_density === nothing ? Inf : maximum(abs.(dens .- obs.previous_density))
    maxerr = obs.current_maxerr
    thresholds_met = sweep >= obs.min_sweeps &&
        sweep >= obs.target_sweep &&
        isfinite(delta_energy) &&
        isfinite(max_density_delta) &&
        delta_energy <= obs.energy_tol &&
        max_density_delta <= obs.density_tol &&
        maxerr <= obs.truncerr_tol
    next_streak = thresholds_met ? obs.satisfied_streak + 1 : 0
    converged = convergence_satisfied(
        sweep,
        delta_energy,
        max_density_delta,
        maxerr;
        min_sweeps=obs.min_sweeps,
        target_sweep=obs.target_sweep,
        stable_sweeps=obs.stable_sweeps,
        satisfied_streak=next_streak,
        energy_tol=obs.energy_tol,
        density_tol=obs.density_tol,
        truncerr_tol=obs.truncerr_tol,
    )
    reason = convergence_reason(sweep, delta_energy, max_density_delta, maxerr, next_streak, obs)
    push!(
        obs.records,
        ConvergenceRecord(sweep, e, delta_energy, max_density_delta, maxlinkdim(psi), maxerr, converged, reason),
    )
    obs.previous_energy = e
    obs.previous_density = copy(dens)
    obs.satisfied_streak = next_streak
    if outputlevel > 0
        @printf(
            "Convergence check sweep %d: |dE|=%.3e max|dn|=%.3e maxerr=%.3e streak=%d/%d reason=%s\n",
            sweep,
            delta_energy,
            max_density_delta,
            maxerr,
            next_streak,
            obs.stable_sweeps,
            reason,
        )
        flush(stdout)
    end
    if converged && outputlevel > 0
        @printf(
            "Adaptive convergence reached at sweep %d: |dE|=%.3e max|dn|=%.3e maxerr=%.3e\n",
            sweep,
            delta_energy,
            max_density_delta,
            maxerr,
        )
        flush(stdout)
    end
    return converged
end

function initial_state(lat::CylinderLat, Np::Int)
    0 <= Np <= lat.Ns || throw(ArgumentError("Np must be between 0 and Ns"))
    state = fill("Emp", lat.Ns)
    Np == 0 && return state
    picks = unique(round.(Int, range(1, lat.Ns; length=Np + 2))[2:(end - 1)])
    candidate = 1
    while length(picks) < Np
        candidate in picks || push!(picks, candidate)
        candidate += 1
    end
    for i in sort(picks[1:Np])
        state[i] = "Occ"
    end
    return state
end

function run_dmrg(
    lat::CylinderLat,
    params::CylinderModelParams,
    Np::Int;
    phi_y::Float64=0.0,
    nsweeps::Int=8,
    maxdim::Vector{Int}=[20, 50, 100, 200],
    cutoff::Float64=1e-8,
    outputlevel::Int=1,
    seed::Int=1234,
    sites=nothing,
    psi0=nothing,
    init_linkdim::Int=8,
    adaptive::Bool=false,
    max_sweeps=nothing,
    min_sweeps::Int=4,
    target_sweep::Int=1,
    stable_sweeps::Int=1,
    energy_tol::Float64=1e-7,
    density_tol::Float64=1e-5,
    truncerr_tol::Float64=Inf,
)
    H, sites = build_mpo(lat, params; phi_y, sites)
    if psi0 === nothing
        state = initial_state(lat, Np)
        rng = MersenneTwister(seed)
        psi0 = random_mps(rng, ComplexF64, sites, state; linkdims=init_linkdim)
    end
    records = ConvergenceRecord[]
    if adaptive
        sweeps_to_run = max_sweeps === nothing ? nsweeps : Int(max_sweeps)
        observer = DMRGConvergenceObserver(;
            min_sweeps,
            target_sweep,
            stable_sweeps,
            energy_tol,
            density_tol,
            truncerr_tol,
        )
        energy, psi = dmrg(H, psi0; nsweeps=sweeps_to_run, maxdim, cutoff, outputlevel, observer)
        records = observer.records
    else
        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, outputlevel)
    end
    dens = density_profile(psi)
    return DMRGRunResult(real(energy), psi, sites, dens, records)
end

density_profile(psi::MPS) = real.(expect(psi, "N"))

function y_averaged_density(lat::CylinderLat, dens::AbstractVector{<:Real})
    out = zeros(Float64, lat.Lx)
    for x in 0:(lat.Lx - 1)
        inds = [site_index(lat, x, y) for y in 0:(lat.Ly - 1)]
        out[x + 1] = mean(dens[inds])
    end
    return out
end

green_function(psi::MPS) = correlation_matrix(psi, "Cdag", "C")

function connected_density_correlation(psi::MPS)
    dens = density_profile(psi)
    nn = real.(correlation_matrix(psi, "N", "N"))
    return nn .- dens * dens'
end

function cumulative_charge(lat::CylinderLat, density_delta::AbstractVector{<:Real})
    out = zeros(Float64, lat.Lx)
    running = 0.0
    for x in 0:(lat.Lx - 1)
        for y in 0:(lat.Ly - 1)
            running += density_delta[site_index(lat, x, y)]
        end
        out[x + 1] = running
    end
    return out
end

function write_density(path::AbstractString, lat::CylinderLat, dens::AbstractVector{<:Real})
    open(path, "w") do io
        println(io, "# site x y rx ry density")
        for (i, (x, y)) in enumerate(lat.sites)
            rx, ry = site_phys_pos(lat, x, y)
            @printf(io, "%d %d %d %.12g %.12g %.16g\n", i, x, y, rx, ry, dens[i])
        end
    end
end

function write_convergence(path::AbstractString, records::AbstractVector{ConvergenceRecord})
    open(path, "w") do io
        println(io, "# sweep energy delta_energy max_density_delta maxlinkdim maxerr converged reason")
        for r in records
            @printf(
                io,
                "%d %.16g %.16g %.16g %d %.16g %s %s\n",
                r.sweep,
                r.energy,
                r.delta_energy,
                r.max_density_delta,
                r.maxlinkdim,
                r.maxerr,
                string(r.converged),
                r.reason,
            )
        end
    end
end

function write_vector(path::AbstractString, xs::AbstractVector, ys::AbstractVector; header::AbstractString="# x value")
    open(path, "w") do io
        println(io, header)
        for (x, y) in zip(xs, ys)
            @printf(io, "%s %.16g\n", string(x), y)
        end
    end
end

function write_matrix(path::AbstractString, M::AbstractMatrix{<:Real}; header::AbstractString="# i j value")
    open(path, "w") do io
        println(io, header)
        for i in axes(M, 1), j in axes(M, 2)
            @printf(io, "%d %d %.16g\n", i, j, M[i, j])
        end
    end
end

function write_complex_matrix(path::AbstractString, M::AbstractMatrix; header::AbstractString="# i j real imag")
    open(path, "w") do io
        println(io, header)
        for i in axes(M, 1), j in axes(M, 2)
            @printf(io, "%d %d %.16g %.16g\n", i, j, real(M[i, j]), imag(M[i, j]))
        end
    end
end

end
