#!/usr/bin/env julia

include("../../shared/lattice.jl")
include("../../shared/hoppings.jl")
include("../../shared/basis.jl")
include("../../shared/ksector.jl")
include("../../shared/hamiltonian.jl")
include("../../shared/optical_response.jl")
include("../../shared/solver.jl")

using Dates
using JLD2
using LinearAlgebra
using Printf
using Random
using SparseArrays

function parse_spectrum_argument(name, default)
    index = findfirst(==(name), ARGS)
    index === nothing && return default
    index < length(ARGS) || error("missing value after $name")
    return ARGS[index + 1]
end

const SPECTRUM_SECTOR_START = parse(
    Int, parse_spectrum_argument("--sector-start", "0"))
const SPECTRUM_SECTOR_END = parse(
    Int, parse_spectrum_argument(
        "--sector-end", string(SPECTRUM_SECTOR_START)))
const SPECTRUM_OUTPUT_DIR = parse_spectrum_argument(
    "--output-dir", joinpath(@__DIR__, "output"))
const SPECTRUM_LATTICE = parse_spectrum_argument(
    "--lattice", "corrected")
const SPECTRUM_T1 = parse(
    Float64, parse_spectrum_argument("--t1", "1.0"))
const SPECTRUM_T3 = parse(
    Float64, parse_spectrum_argument("--t3", "0.2"))
const SPECTRUM_V1 = parse(
    Float64, parse_spectrum_argument("--V1", "10.0"))
const SPECTRUM_V2 = parse(
    Float64, parse_spectrum_argument("--V2", "2.0"))
const SPECTRUM_V3 = parse(
    Float64, parse_spectrum_argument("--V3", "2.0"))
const SPECTRUM_NP = parse(
    Int, parse_spectrum_argument("--Np", "12"))
const SPECTRUM_NEV = parse(
    Int, parse_spectrum_argument("--nev", "8"))
const SPECTRUM_KRYLOVDIM = parse(
    Int, parse_spectrum_argument("--krylovdim", "60"))
const SPECTRUM_SEED = parse(
    Int, parse_spectrum_argument(
        "--seed", string(20260919 + SPECTRUM_SECTOR_START)))

0 <= SPECTRUM_SECTOR_START <= SPECTRUM_SECTOR_END <= 14 ||
    error("sector range must satisfy 0 <= start <= end <= 14")
SPECTRUM_SECTOR_END == SPECTRUM_SECTOR_START ||
    error("30-site corrected spectrum jobs must contain one sector")
SPECTRUM_LATTICE in ("corrected", "legacy") ||
    error("lattice must be corrected or legacy")
1 <= SPECTRUM_NP < 30 || error("Np must satisfy 1 <= Np < 30")
SPECTRUM_NEV > 0 || error("nev must be positive")
SPECTRUM_KRYLOVDIM > SPECTRUM_NEV ||
    error("krylovdim must be greater than nev")

function spectrum_lattice(convention::AbstractString=SPECTRUM_LATTICE)
    convention == "corrected" && return TiltedLat30()
    convention == "legacy" && return LegacyTiltedLat30()
    throw(ArgumentError("unknown lattice convention: $convention"))
end

function spectrum_partial_path(
        first_sector::Int=SPECTRUM_SECTOR_START)
    return joinpath(
        SPECTRUM_OUTPUT_DIR, "partial_$(first_sector).jld2")
end

spectrum_lanczos_matrix(matrix::SparseMatrixCSC{ComplexF64,I}) where {I<:Integer} =
    threaded_csr(matrix)

function run_large_spectrum()
    n_local = SPECTRUM_SECTOR_END - SPECTRUM_SECTOR_START + 1
    blas_threads = max(1, Threads.nthreads() ÷ n_local)
    BLAS.set_num_threads(blas_threads)
    Random.seed!(SPECTRUM_SEED)

    println("="^72)
    println("30-site ED spectrum, sectors " *
            "$SPECTRUM_SECTOR_START:$SPECTRUM_SECTOR_END")
    println("start: ", now())
    println("lattice=$SPECTRUM_LATTICE Np=$SPECTRUM_NP " *
            "t1=$SPECTRUM_T1 t3=$SPECTRUM_T3 " *
            "V1=$SPECTRUM_V1 V2=$SPECTRUM_V2 V3=$SPECTRUM_V3")
    println("nev=$SPECTRUM_NEV krylovdim=$SPECTRUM_KRYLOVDIM " *
            "Julia_threads=$(Threads.nthreads()) " *
            "BLAS_threads=$(BLAS.get_num_threads()) seed=$SPECTRUM_SEED")
    flush(stdout)

    lattice = spectrum_lattice()
    println("[1] generating C($(lattice.Ns),$SPECTRUM_NP) basis")
    basis_time = @elapsed basis = gen_basis(lattice.Ns, SPECTRUM_NP)
    @printf("    states=%d time=%.2f s RSS=%.2f GB\n",
            length(basis), basis_time, mem_rss_gb())

    println("[2] building corrected real-space hopping")
    hopping_time = @elapsed hops = build_hops(
        lattice, SPECTRUM_T1, SPECTRUM_T3, 0.0)
    character_error = translation_character_error(lattice)
    @printf("    hoppings=%d character_error=%.3e time=%.2f s\n",
            length(hops), character_error, hopping_time)
    character_error < 1e-12 || error(
        "momentum mesh is inconsistent with translation characters")
    length(hops) == 300 || error(
        "unexpected tilted-30 hopping count: $(length(hops))")

    translation_maps = translation_site_maps(lattice)
    sectors = KSector[]
    hamiltonians = SparseMatrixCSC{ComplexF64,Int32}[]
    println("[3] building momentum sectors and sparse Hamiltonians")
    pipeline_time = @elapsed begin
        for momentum in SPECTRUM_SECTOR_START:SPECTRUM_SECTOR_END
            sector_time = @elapsed sector = build_ksector(
                basis, lattice, momentum; trans_maps=translation_maps)
            matrix_time = @elapsed matrix = build_sparse_H(
                sector, lattice, hops,
                SPECTRUM_V1, SPECTRUM_V2, SPECTRUM_V3)
            @printf("    k=%2d dim=%d nnz=%d sector=%.2f s H=%.2f s RSS=%.2f GB\n",
                    momentum, length(sector.reps), nnz(matrix),
                    sector_time, matrix_time, mem_rss_gb())
            flush(stdout)
            empty!(sector.fock2rep)
            empty!(sector.orbit_data)
            push!(sectors, sector)
            push!(hamiltonians, matrix)
            GC.gc()
        end
    end
    empty!(basis)
    GC.gc()
    @printf("    pipeline=%.2f s RSS=%.2f GB\n",
            pipeline_time, mem_rss_gb())

    println("[4] converting sparse Hamiltonian to threaded CSR")
    conversion_time = @elapsed lanczos_matrices = [
        spectrum_lanczos_matrix(matrix) for matrix in hamiltonians
    ]
    empty!(hamiltonians)
    GC.gc()
    @printf("    conversion=%.2f s RSS=%.2f GB\n",
            conversion_time, mem_rss_gb())

    println("[5] solving lowest eigenpairs")
    solve_time = @elapsed eigenvalues, ground_states =
        lanczos_sparse_sectors(
            sectors, lanczos_matrices;
            nev=SPECTRUM_NEV,
            krylovdim=SPECTRUM_KRYLOVDIM,
            verbose=true)

    residuals = Dict{Int,Float64}()
    sector_ground_energies = Dict{Int,Float64}()
    for (index, sector) in enumerate(sectors)
        momentum = sector.m
        energies = [energy for (label, energy) in eigenvalues
                    if label == momentum]
        isempty(energies) && error("sector $momentum has no eigenvalue")
        energy = minimum(energies)
        state = ground_states[momentum]
        product = similar(state)
        mul!(product, lanczos_matrices[index], state)
        residuals[momentum] = norm(product .- energy .* state)
        sector_ground_energies[momentum] = energy
        @printf("    k=%2d E0=%.12f residual=%.3e norm=%.12f\n",
                momentum, energy, residuals[momentum], norm(state))
        residuals[momentum] < 1e-7 ||
            error("sector $momentum ground-state residual is too large")
    end
    @printf("    solve=%.2f s\n", solve_time)

    println("[6] saving states")
    mkpath(SPECTRUM_OUTPUT_DIR)
    output_path = spectrum_partial_path()
    jldsave(
        output_path;
        ev_pairs=eigenvalues,
        gs_vecs=ground_states,
        residuals,
        sector_ground_energies,
        Np=SPECTRUM_NP,
        V1=SPECTRUM_V1,
        V2=SPECTRUM_V2,
        V3=SPECTRUM_V3,
        t1=SPECTRUM_T1,
        t3=SPECTRUM_T3,
        lattice_convention=SPECTRUM_LATTICE,
        momentum_step=copy(lattice.kpoints[2]),
        translation_character_error=character_error,
        hopping_count=length(hops),
        sector_start=SPECTRUM_SECTOR_START,
        sector_end=SPECTRUM_SECTOR_END,
        nev=SPECTRUM_NEV,
        krylovdim=SPECTRUM_KRYLOVDIM,
        random_seed=SPECTRUM_SEED,
    )
    println("saved: $output_path")
    println("done: ", now())
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_large_spectrum()
end
