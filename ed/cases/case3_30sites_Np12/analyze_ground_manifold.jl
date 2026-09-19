#!/usr/bin/env julia

include("../../shared/lattice.jl")
include("../../shared/basis.jl")
include("../../shared/ksector.jl")
include("../../shared/structure_factor.jl")

using Dates
using JLD2
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
const MANIFOLD_RESULT_DIR = parse_manifold_argument(
    "--result-dir",
    normpath(joinpath(
        @__DIR__, "..", "..", "..", "results",
        "phc30_np12_v1_10_v2_2_v3_2_fqahc")))
const MANIFOLD_SIZE = parse(
    Int, parse_manifold_argument("--manifold-size", "15"))
const MANIFOLD_LATTICE = parse_manifold_argument(
    "--lattice", "corrected")

MANIFOLD_SIZE > 0 || error("manifold size must be positive")
MANIFOLD_LATTICE == "corrected" ||
    error("new moderate-coupling results must use corrected lattice")

function manifold_partial_paths()
    return [joinpath(MANIFOLD_DATA_DIR, "partial_$first.jld2")
            for first in (0, 5, 10)]
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
    sectors = sort([level[1] for level in ordered[1:manifold_size]])
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

function validate_manifold_partial(data, path)
    expected = (
        Np=12,
        t1=1.0,
        t3=0.2,
        V1=10.0,
        V2=2.0,
        V3=2.0,
        lattice_convention="corrected",
    )
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
    return observed
end

function load_manifold_partials()
    all_levels = Tuple{Int,Float64}[]
    states = Dict{Int,Vector{ComplexF64}}()
    residuals = Dict{Int,Float64}()
    metadata = NamedTuple[]
    seen_sectors = Set{Int}()

    for path in manifold_partial_paths()
        isfile(path) || error("missing partial spectrum file: $path")
        data = load(path)
        parameters = validate_manifold_partial(data, path)
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
        println(output, "result_id = \"phc30_np12_v1_10_v2_2_v3_2_fqahc\"")
        println(output, "phase_label = \"FQAHC candidate\"")
        println(output, "Ns = 30")
        println(output, "Nuc = 15")
        println(output, "Np = 12")
        println(output, "t1 = 1.0")
        println(output, "t3 = 0.2")
        println(output, "V1 = 10.0")
        println(output, "V2 = 2.0")
        println(output, "V3 = 2.0")
        println(output, "lattice_convention = \"corrected\"")
        println(output, "candidate_manifold_size = $MANIFOLD_SIZE")
    end
end

function write_manifold_manifest(path, metadata)
    open(path, "w") do output
        println(output, "schema_version = 1")
        println(output, "result_id = \"phc30_np12_v1_10_v2_2_v3_2_fqahc\"")
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
        end
    end
end

function run_ground_manifold_analysis()
    println("="^72)
    println("30-site moderate-coupling candidate-manifold analysis")
    println("start: ", now())
    println("data: $MANIFOLD_DATA_DIR")
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
                "# k E-E0; corrected TiltedLat30, Np=12, " *
                "V1=10, V2=V3=2")
        for (momentum, energy) in ordered_levels
            @printf(output, "%d %.12f\n",
                    momentum, energy - summary.ground_energy)
        end
    end

    lattice = TiltedLat30()
    basis = gen_basis(lattice.Ns, 12)
    translation_maps = translation_site_maps(lattice)
    curves = Dict{Int,Vector{Tuple{Int,Float64,Float64,Float64}}}()
    println("computing structure factors for 15 candidate states")
    for momentum in summary.sectors
        state = states[momentum]
        sector_time = @elapsed sector = build_ksector(
            basis, lattice, momentum; trans_maps=translation_maps)
        length(state) == length(sector.reps) || error(
            "state dimension does not match sector $momentum")
        survey_time = @elapsed ordered, _, norm2, particle_number =
            structure_factor_survey_representatives(
                sector, state, lattice, 12)
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
                          for m in summary.sectors) / MANIFOLD_SIZE
            @printf(output, "average %d %.12f %.12f %.12f\n",
                    q, average, qx, qy)
        end
    end

    diagnostics_path = joinpath(data_dir, "manifold_diagnostics.txt")
    open(diagnostics_path, "w") do output
        println(output, "phase_label=FQAHC_candidate")
        println(output, "topology_confirmed_by_this_run=false")
        println(output, "lattice_convention=corrected")
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
