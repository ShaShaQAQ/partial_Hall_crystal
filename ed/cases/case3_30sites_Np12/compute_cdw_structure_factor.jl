#!/usr/bin/env julia

include("../../shared/lattice.jl")
include("../../shared/basis.jl")
include("../../shared/ksector.jl")
include("../../shared/structure_factor.jl")

using Dates
using JLD2
using Printf

function parse_structure_factor_argument(name, default)
    index = findfirst(==(name), ARGS)
    index === nothing && return default
    index < length(ARGS) || error("missing value after $name")
    return ARGS[index + 1]
end

const STRUCTURE_FACTOR_DATA_DIR = parse_structure_factor_argument(
    "--data-dir", joinpath(@__DIR__, "output"))
const STRUCTURE_FACTOR_RESULT_DIR = parse_structure_factor_argument(
    "--result-dir",
    normpath(joinpath(
        @__DIR__, "..", "..", "..", "results",
        "phc30_np12_v1_100_cdw")))
const CDW_GROUND_SECTORS = (0, 5, 10)

function structure_factor_partial_path(sector::Int)
    first_sector = 5 * div(sector, 5)
    return joinpath(
        STRUCTURE_FACTOR_DATA_DIR, "partial_$(first_sector).jld2")
end

function load_cdw_state(sector::Int)
    path = structure_factor_partial_path(sector)
    isfile(path) || error("missing saved data: $path")
    data = load(path)
    parameters = (
        Np=Int(data["Np"]),
        t1=Float64(data["t1"]),
        t3=Float64(data["t3"]),
        V1=Float64(data["V1"]),
        V2=Float64(data["V2"]),
        V3=Float64(data["V3"]),
    )
    expected = (Np=12, t1=1.0, t3=0.2, V1=100.0, V2=0.0, V3=0.0)
    parameters == expected ||
        error("saved parameters do not match CDW result package: $parameters")
    states = data["gs_vecs"]
    haskey(states, sector) || error("sector $sector is absent from $path")
    energies = [Float64(energy)
                for (momentum, energy) in data["ev_pairs"]
                if Int(momentum) == sector]
    isempty(energies) && error("sector $sector energy is absent from $path")
    return (state=ComplexF64.(states[sector]),
            energy=minimum(energies), path=path)
end

function run_cdw_structure_factor()
    println("="^72)
    println("30-site strong-coupling CDW ground-state structure factors")
    println("start: ", now())
    println("data: $STRUCTURE_FACTOR_DATA_DIR")
    flush(stdout)

    lattice = LegacyTiltedLat30()
    basis = gen_basis(lattice.Ns, 12)
    translation_maps = translation_site_maps(lattice)
    curves = Dict{Int,Vector{Tuple{Int,Float64,Float64,Float64}}}()
    diagnostics = Dict{Int,NamedTuple}()

    for sector_label in CDW_GROUND_SECTORS
        saved = load_cdw_state(sector_label)
        @printf("sector k=%d source=%s\n", sector_label, saved.path)
        sector_time = @elapsed sector = build_ksector(
            basis, lattice, sector_label; trans_maps=translation_maps)
        length(saved.state) == length(sector.reps) ||
            error("saved vector dimension does not match sector $sector_label")
        survey_time = @elapsed begin
            ordered, sorted, norm2, particle_number =
                structure_factor_survey_representatives(
                    sector, saved.state, lattice, 12)
        end
        curves[sector_label] = ordered
        diagnostics[sector_label] = (
            energy=saved.energy,
            norm2=norm2,
            particle_number=particle_number,
            peak_momentum_1=sorted[1][1],
            peak_value_1=sorted[1][2],
            peak_momentum_2=sorted[2][1],
            peak_value_2=sorted[2][2],
        )
        @printf("  dim=%d sector=%.2f s survey=%.2f s N(q=%d)=%.6f N(q=%d)=%.6f\n",
                length(saved.state), sector_time, survey_time,
                sorted[1][1], sorted[1][2], sorted[2][1], sorted[2][2])
        sector = nothing
        saved = nothing
        GC.gc()
    end

    result_data_dir = joinpath(STRUCTURE_FACTOR_RESULT_DIR, "data")
    mkpath(result_data_dir)
    output_path = joinpath(
        result_data_dir, "structure_factor_ground_states.dat")
    open(output_path, "w") do output
        println(output,
                "# state_sector q N(q) qx qy; unit-cell normalization 1/Nuc")
        println(output,
                "# legacy TiltedLat30 Fourier mesh matching saved partial_*.jld2")
        for sector_label in CDW_GROUND_SECTORS
            for (momentum, value, qx, qy) in curves[sector_label]
                @printf(output, "%d %d %.12f %.12f %.12f\n",
                        sector_label, momentum, value, qx, qy)
            end
        end
        for momentum_index in eachindex(lattice.ktab)
            momentum = lattice.ktab[momentum_index]
            qx = lattice.kpoints[momentum_index][1]
            qy = lattice.kpoints[momentum_index][2]
            average = sum(
                curves[sector_label][momentum_index][2]
                for sector_label in CDW_GROUND_SECTORS
            ) / length(CDW_GROUND_SECTORS)
            @printf(output, "average %d %.12f %.12f %.12f\n",
                    momentum, average, qx, qy)
        end
    end

    diagnostic_path = joinpath(
        result_data_dir, "structure_factor_diagnostics.txt")
    open(diagnostic_path, "w") do output
        println(output, "normalization=unit_cell")
        println(output, "normalization_denominator=$(lattice.Nuc)")
        println(output, "fourier_mesh=legacy_tilted30_april_2026")
        for sector_label in CDW_GROUND_SECTORS
            values = diagnostics[sector_label]
            println(output, "sector=$sector_label")
            println(output, "energy=$(values.energy)")
            println(output, "norm2=$(values.norm2)")
            println(output, "particle_number=$(values.particle_number)")
            println(output,
                    "peak_1=$(values.peak_momentum_1),$(values.peak_value_1)")
            println(output,
                    "peak_2=$(values.peak_momentum_2),$(values.peak_value_2)")
        end
    end
    println("saved: $output_path")
    println("saved: $diagnostic_path")
    println("done: ", now())
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cdw_structure_factor()
end
