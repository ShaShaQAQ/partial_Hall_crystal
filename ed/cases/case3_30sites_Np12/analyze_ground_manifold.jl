#!/usr/bin/env julia

include("../../shared/lattice.jl")
include("../../shared/basis.jl")
include("../../shared/ksector.jl")
include("../../shared/structure_factor.jl")

using Dates
using JLD2
using LinearAlgebra
using Printf
using SHA

function parse_manifold_argument(name, default)
    index = findfirst(==(name), ARGS)
    index === nothing && return default
    index < length(ARGS) || error("missing value after $name")
    return ARGS[index + 1]
end

const MANIFOLD_DATA_DIR = parse_manifold_argument(
    "--data-dir", joinpath(@__DIR__, "output"))
const MANIFOLD_RESULT_ID = parse_manifold_argument(
    "--result-id", "phc30_np12_v1_10_v2_2_v3_2_fqahc")
const MANIFOLD_RESULT_DIR = parse_manifold_argument(
    "--result-dir",
    normpath(joinpath(
        @__DIR__, "..", "..", "..", "results",
        MANIFOLD_RESULT_ID)))
const MANIFOLD_PHASE_LABEL = parse_manifold_argument(
    "--phase-label", "FQAHC candidate")
const MANIFOLD_SIZE = parse(
    Int, parse_manifold_argument("--manifold-size", "15"))
const MANIFOLD_LATTICE = parse_manifold_argument(
    "--lattice", "corrected")
const MANIFOLD_NP = parse(
    Int, parse_manifold_argument("--Np", "12"))
const MANIFOLD_T1 = parse(
    Float64, parse_manifold_argument("--t1", "1.0"))
const MANIFOLD_T3 = parse(
    Float64, parse_manifold_argument("--t3", "0.2"))
const MANIFOLD_V1 = parse(
    Float64, parse_manifold_argument("--V1", "10.0"))
const MANIFOLD_V2 = parse(
    Float64, parse_manifold_argument("--V2", "2.0"))
const MANIFOLD_V3 = parse(
    Float64, parse_manifold_argument("--V3", "2.0"))

MANIFOLD_SIZE > 0 || error("manifold size must be positive")
1 <= MANIFOLD_NP < 30 || error("Np must satisfy 1 <= Np < 30")
MANIFOLD_LATTICE == "corrected" ||
    error("new 30-site manifold results must use corrected lattice")

function manifold_expected_parameters()
    return (
        Np=MANIFOLD_NP,
        t1=MANIFOLD_T1,
        t3=MANIFOLD_T3,
        V1=MANIFOLD_V1,
        V2=MANIFOLD_V2,
        V3=MANIFOLD_V3,
        lattice_convention=MANIFOLD_LATTICE,
    )
end

function manifold_partial_paths()
    return [joinpath(MANIFOLD_DATA_DIR, "partial_$first.jld2")
            for first in 0:14]
end

function classify_candidate_manifold(
        levels::Vector{Tuple{Int,Float64}}, manifold_size::Int)
    manifold_size > 0 || throw(ArgumentError(
        "manifold size must be positive"))
    length(levels) > manifold_size || throw(ArgumentError(
        "need at least manifold_size + 1 levels"))
    ordered = sort(levels; by=level -> level[2])
    ground_energy = ordered[1][2]
    top_energy = ordered[manifold_size][2]
    next_absolute = ordered[manifold_size + 1][2]
    candidate_sectors = [level[1]
                         for level in ordered[1:manifold_size]]
    repeated_sectors = sort([
        sector for sector in unique(candidate_sectors)
        if count(==(sector), candidate_sectors) > 1
    ])
    isempty(repeated_sectors) || throw(ArgumentError(
        "candidate manifold contains repeated sector " *
        "$(first(repeated_sectors)); partial files save only one ground " *
        "state per sector, so an excited-state wavefunction is " *
        "unavailable"))
    sectors = sort(candidate_sectors)
    return (
        sectors=sectors,
        ground_energy=ground_energy,
        width=top_energy - ground_energy,
        next_energy=next_absolute - ground_energy,
        separation=next_absolute - top_energy,
        isolated=(next_absolute - top_energy) >
                 (top_energy - ground_energy),
    )
end

function file_sha256(path::AbstractString)
    return open(path, "r") do input
        bytes2hex(sha256(input))
    end
end

function validate_manifold_partial(
        data, path, expected, expected_sector::Integer)
    observed = (
        Np=Int(data["Np"]),
        t1=Float64(data["t1"]),
        t3=Float64(data["t3"]),
        V1=Float64(data["V1"]),
        V2=Float64(data["V2"]),
        V3=Float64(data["V3"]),
        lattice_convention=String(data["lattice_convention"]),
    )
    observed == expected || error(
        "parameter mismatch in $path: expected $expected, got $observed")

    sector_start = Int(data["sector_start"])
    sector_end = Int(data["sector_end"])
    sector_start == expected_sector == sector_end || error(
        "sector range mismatch in $path: expected " *
        "$expected_sector:$expected_sector, got $sector_start:$sector_end")

    ev_pairs = data["ev_pairs"]
    ev_sectors = Set(Int(momentum)
                     for (momentum, _) in ev_pairs)
    ev_sectors == Set([Int(expected_sector)]) || error(
        "ev_pairs sectors mismatch in $path: expected " *
        "[$expected_sector], got $(sort!(collect(ev_sectors)))")

    ground_states = data["gs_vecs"]
    state_sectors = Set(Int(momentum) for momentum in keys(ground_states))
    state_sectors == Set([Int(expected_sector)]) || error(
        "gs_vecs sectors mismatch in $path: expected " *
        "[$expected_sector], got $(sort!(collect(state_sectors)))")

    residuals = data["residuals"]
    residual_sectors = Set(Int(momentum) for momentum in keys(residuals))
    residual_sectors == Set([Int(expected_sector)]) || error(
        "residuals sectors mismatch in $path: expected " *
        "[$expected_sector], got $(sort!(collect(residual_sectors)))")

    state_norm = norm(ground_states[expected_sector])
    isfinite(state_norm) &&
        isapprox(state_norm, 1.0; atol=1e-10, rtol=0) || error(
        "sector $expected_sector ground state in $path is not normalized: " *
        "norm=$state_norm")

    residual = Float64(residuals[expected_sector])
    isfinite(residual) && 0 <= residual < 1e-7 || error(
        "sector $expected_sector residual in $path must be finite and " *
        "< 1e-7, got $residual")

    sector_ground_energies = data["sector_ground_energies"]
    energy_sectors = Set(
        Int(momentum) for momentum in keys(sector_ground_energies))
    energy_sectors == Set([Int(expected_sector)]) || error(
        "sector_ground_energies sectors mismatch in $path: expected " *
        "[$expected_sector], got $(sort!(collect(energy_sectors)))")
    saved_ground_energy = Float64(
        sector_ground_energies[expected_sector])
    sector_energies = Float64[
        energy for (momentum, energy) in ev_pairs
        if Int(momentum) == expected_sector
    ]
    ev_ground_energy = minimum(sector_energies)
    isfinite(saved_ground_energy) && isfinite(ev_ground_energy) &&
        isapprox(saved_ground_energy, ev_ground_energy;
                 atol=1e-10, rtol=1e-12) || error(
        "sector $expected_sector ground energy mismatch in $path: " *
        "saved=$saved_ground_energy ev_pairs=$ev_ground_energy")

    for key in ("momentum_step", "translation_character_error", "hopping_count")
        haskey(data, key) || error("missing lattice fingerprint $key in $path")
    end
    momentum_step = Float64.(data["momentum_step"])
    expected_step = TiltedLat30().kpoints[2]
    isapprox(momentum_step, expected_step; atol=1e-12, rtol=0) || error(
        "momentum-step mismatch in $path")
    character_error = Float64(data["translation_character_error"])
    character_error < 1e-12 || error(
        "translation-character mismatch in $path: $character_error")
    hopping_count = Int(data["hopping_count"])
    hopping_count == 300 || error(
        "unexpected hopping count in $path: $hopping_count")
    return merge(observed, (
        momentum_step=momentum_step,
        translation_character_error=character_error,
        hopping_count=hopping_count,
    ))
end

function load_manifold_partials()
    all_levels = Tuple{Int,Float64}[]
    states = Dict{Int,Vector{ComplexF64}}()
    residuals = Dict{Int,Float64}()
    metadata = NamedTuple[]
    seen_sectors = Set{Int}()

    expected = manifold_expected_parameters()
    for (expected_sector, path) in zip(0:14, manifold_partial_paths())
        isfile(path) || error("missing partial spectrum file: $path")
        data = load(path)
        parameters = validate_manifold_partial(
            data, path, expected, expected_sector)
        append!(all_levels, [
            (Int(momentum), Float64(energy))
            for (momentum, energy) in data["ev_pairs"]
        ])
        for (momentum, state) in data["gs_vecs"]
            label = Int(momentum)
            label in seen_sectors && error(
                "duplicate sector $label in partial files")
            push!(seen_sectors, label)
            states[label] = ComplexF64.(state)
        end
        for (momentum, residual) in data["residuals"]
            residuals[Int(momentum)] = Float64(residual)
        end
        push!(metadata, merge(parameters, (
            path=path,
            size_bytes=filesize(path),
            sha256=file_sha256(path),
            sector_start=Int(data["sector_start"]),
            sector_end=Int(data["sector_end"]),
            nev=Int(data["nev"]),
            krylovdim=Int(data["krylovdim"]),
        )))
        data = nothing
        GC.gc()
    end

    seen_sectors == Set(0:14) || error(
        "partial files do not contain exactly sectors 0:14")
    return all_levels, states, residuals, metadata
end

function write_manifold_parameters(path)
    open(path, "w") do output
        println(output, "result_id = \"$MANIFOLD_RESULT_ID\"")
        println(output, "phase_label = \"$MANIFOLD_PHASE_LABEL\"")
        println(output, "Ns = 30")
        println(output, "Nuc = 15")
        println(output, "Np = $MANIFOLD_NP")
        println(output, "t1 = $MANIFOLD_T1")
        println(output, "t3 = $MANIFOLD_T3")
        println(output, "V1 = $MANIFOLD_V1")
        println(output, "V2 = $MANIFOLD_V2")
        println(output, "V3 = $MANIFOLD_V3")
        println(output, "lattice_convention = \"$MANIFOLD_LATTICE\"")
        println(output, "candidate_manifold_size = $MANIFOLD_SIZE")
    end
end

function write_manifold_manifest(path, metadata)
    open(path, "w") do output
        println(output, "schema_version = 1")
        println(output, "result_id = \"$MANIFOLD_RESULT_ID\"")
        println(output, "storage_host = \"W003\"")
        println(output, "data_directory = \"$MANIFOLD_DATA_DIR\"")
        println(output, "generated_at = \"$(now())\"")
        for item in metadata
            println(output, "\n[[partial_files]]")
            println(output, "sector_start = $(item.sector_start)")
            println(output, "sector_end = $(item.sector_end)")
            println(output, "path = \"$(item.path)\"")
            println(output, "size_bytes = $(item.size_bytes)")
            println(output, "sha256 = \"$(item.sha256)\"")
            println(output, "nev = $(item.nev)")
            println(output, "krylovdim = $(item.krylovdim)")
            println(output, "momentum_step = [$(join(item.momentum_step, ", "))]")
            println(output,
                    "translation_character_error = $(item.translation_character_error)")
            println(output, "hopping_count = $(item.hopping_count)")
        end
    end
end

function run_ground_manifold_analysis()
    println("="^72)
    println("30-site candidate-manifold analysis")
    println("start: ", now())
    println("data: $MANIFOLD_DATA_DIR")
    println("result: $MANIFOLD_RESULT_ID phase: $MANIFOLD_PHASE_LABEL")
    println("Np=$MANIFOLD_NP t1=$MANIFOLD_T1 t3=$MANIFOLD_T3 " *
            "V1=$MANIFOLD_V1 V2=$MANIFOLD_V2 V3=$MANIFOLD_V3")
    flush(stdout)

    levels, states, residuals, metadata = load_manifold_partials()
    summary = classify_candidate_manifold(levels, MANIFOLD_SIZE)
    ordered_levels = sort(levels; by=level -> level[2])
    @printf("E0=%.12f width=%.12g next=%.12g separation=%.12g\n",
            summary.ground_energy, summary.width,
            summary.next_energy, summary.separation)
    println("candidate sectors: ", join(summary.sectors, ","))
    println("isolated by separation > width: ", summary.isolated)

    data_dir = joinpath(MANIFOLD_RESULT_DIR, "data")
    mkpath(data_dir)
    spectrum_path = joinpath(data_dir, "spectrum.dat")
    open(spectrum_path, "w") do output
        println(output,
                "# k E-E0; $MANIFOLD_LATTICE TiltedLat30, " *
                "Np=$MANIFOLD_NP, V1=$MANIFOLD_V1, " *
                "V2=$MANIFOLD_V2, V3=$MANIFOLD_V3")
        for (momentum, energy) in ordered_levels
            @printf(output, "%d %.12f\n",
                    momentum, energy - summary.ground_energy)
        end
    end

    lattice = TiltedLat30()
    basis = gen_basis(lattice.Ns, MANIFOLD_NP)
    translation_maps = translation_site_maps(lattice)
    curves = Dict{Int,Vector{Tuple{Int,Float64,Float64,Float64}}}()
    println("computing structure factors for $MANIFOLD_SIZE candidate states")
    for momentum in summary.sectors
        state = states[momentum]
        sector_time = @elapsed sector = build_ksector(
            basis, lattice, momentum; trans_maps=translation_maps)
        length(state) == length(sector.reps) || error(
            "state dimension does not match sector $momentum")
        survey_time = @elapsed ordered, _, norm2, particle_number =
            structure_factor_survey_representatives(
                sector, state, lattice, MANIFOLD_NP)
        curves[momentum] = ordered
        @printf("  k=%2d dim=%d residual=%.3e norm=%.12f Np=%.9f sector=%.2f s survey=%.2f s\n",
                momentum, length(state), residuals[momentum], norm2,
                particle_number, sector_time, survey_time)
        sector = nothing
        GC.gc()
    end

    structure_path = joinpath(
        data_dir, "structure_factor_ground_manifold.dat")
    open(structure_path, "w") do output
        println(output,
                "# state_sector q N(q) qx qy; unit-cell normalization 1/Nuc")
        println(output, "# corrected TiltedLat30 Fourier mesh")
        for momentum in summary.sectors
            for (q, value, qx, qy) in curves[momentum]
                @printf(output, "%d %d %.12f %.12f %.12f\n",
                        momentum, q, value, qx, qy)
            end
        end
        for q_index in eachindex(lattice.ktab)
            q = lattice.ktab[q_index]
            qx, qy = lattice.kpoints[q_index]
            average = sum(curves[m][q_index][2]
                          for m in summary.sectors) / length(summary.sectors)
            @printf(output, "average %d %.12f %.12f %.12f\n",
                    q, average, qx, qy)
        end
    end

    diagnostics_path = joinpath(data_dir, "manifold_diagnostics.txt")
    open(diagnostics_path, "w") do output
        diagnostic_label = replace(MANIFOLD_PHASE_LABEL, ' ' => '_')
        println(output, "phase_label=$diagnostic_label")
        println(output, "topology_confirmed_by_this_run=false")
        println(output, "lattice_convention=$MANIFOLD_LATTICE")
        println(output, "Np=$MANIFOLD_NP")
        println(output, "t1=$MANIFOLD_T1")
        println(output, "t3=$MANIFOLD_T3")
        println(output, "V1=$MANIFOLD_V1")
        println(output, "V2=$MANIFOLD_V2")
        println(output, "V3=$MANIFOLD_V3")
        println(output, "ground_energy=$(summary.ground_energy)")
        println(output, "candidate_manifold_size=$MANIFOLD_SIZE")
        println(output, "candidate_sectors=$(join(summary.sectors, ','))")
        println(output, "manifold_width=$(summary.width)")
        println(output, "next_energy=$(summary.next_energy)")
        println(output, "separation=$(summary.separation)")
        println(output, "isolated=$(summary.isolated)")
        for momentum in summary.sectors
            println(output, "residual_k$(momentum)=$(residuals[momentum])")
        end
    end

    write_manifold_parameters(joinpath(MANIFOLD_RESULT_DIR, "parameters.toml"))
    write_manifold_manifest(joinpath(MANIFOLD_RESULT_DIR, "manifest.toml"), metadata)
    println("saved: $spectrum_path")
    println("saved: $structure_path")
    println("saved: $diagnostics_path")
    println("done: ", now())
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_ground_manifold_analysis()
end
