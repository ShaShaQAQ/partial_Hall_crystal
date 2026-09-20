using Test

include(joinpath(
    @__DIR__, "..", "cases", "case3_30sites_Np12",
    "analyze_ground_manifold.jl"))

@testset "candidate manifold classification" begin
    levels = Tuple{Int,Float64}[
        (momentum, 2.0 + 0.01 * momentum) for momentum in 0:14
    ]
    push!(levels, (5, 3.0))
    push!(levels, (10, 3.2))

    summary = classify_candidate_manifold(levels, 15)
    @test summary.sectors == collect(0:14)
    @test summary.ground_energy == 2.0
    @test summary.width ≈ 0.14
    @test summary.next_energy ≈ 1.0
    @test summary.separation ≈ 0.86
end

@testset "non-fifteen candidate manifold classification" begin
    levels = Tuple{Int,Float64}[
        (2, 1.0),
        (7, 1.05),
        (12, 1.1),
        (4, 2.0),
    ]

    summary = classify_candidate_manifold(levels, 3)
    @test summary.sectors == [2, 7, 12]
    @test summary.width ≈ 0.1
    @test summary.next_energy ≈ 1.0
    @test summary.separation ≈ 0.9
end

@testset "candidate manifold input validation" begin
    @test_throws ArgumentError classify_candidate_manifold(
        Tuple{Int,Float64}[(0, 0.0)], 1)
    @test_throws ArgumentError classify_candidate_manifold(
        Tuple{Int,Float64}[(0, 0.0), (1, 1.0)], 0)
end

function valid_partial_metadata(expected, sector)
    physical = TiltedLat30()
    return Dict{String,Any}(
        "Np" => expected.Np,
        "t1" => expected.t1,
        "t3" => expected.t3,
        "V1" => expected.V1,
        "V2" => expected.V2,
        "V3" => expected.V3,
        "lattice_convention" => "corrected",
        "momentum_step" => copy(physical.kpoints[2]),
        "translation_character_error" => 1e-14,
        "hopping_count" => 300,
        "sector_start" => sector,
        "sector_end" => sector,
        "ev_pairs" => Tuple{Int,Float64}[(sector, 1.0), (sector, 1.5)],
        "gs_vecs" => Dict(sector => ComplexF64[1.0, 0.0]),
        "residuals" => Dict(sector => 1e-12),
        "sector_ground_energies" => Dict(sector => 1.0),
    )
end

function validation_error(metadata, expected, sector)
    try
        validate_manifold_partial(
            metadata, "synthetic.jld2", expected, sector)
        return ""
    catch error
        return sprint(showerror, error)
    end
end

@testset "saved partial validation" begin
    expected = (
        Np=13,
        t1=1.0,
        t3=0.2,
        V1=100.0,
        V2=0.0,
        V3=0.0,
        lattice_convention="corrected",
    )
    sector = 3
    metadata = valid_partial_metadata(expected, sector)
    validated = try
        validate_manifold_partial(
            metadata, "synthetic.jld2", expected, sector)
    catch error
        error
    end
    @test validated isa NamedTuple
    if validated isa NamedTuple
        @test validated.Np == 13
        @test validated.V1 == 100.0
        @test validated.hopping_count == 300
        @test validated.translation_character_error < 1e-12
    end

    mismatch = valid_partial_metadata(expected, sector)
    mismatch["Np"] = 12
    @test occursin(
        "parameter mismatch", validation_error(mismatch, expected, sector))

    for key in ("sector_start", "sector_end")
        mismatch = valid_partial_metadata(expected, sector)
        mismatch[key] = sector - 1
        @test occursin(
            "sector range", validation_error(mismatch, expected, sector))
    end

    for key in ("ev_pairs", "gs_vecs", "residuals")
        mismatch = valid_partial_metadata(expected, sector)
        if key == "ev_pairs"
            mismatch[key] = Tuple{Int,Float64}[(sector + 1, 1.0)]
        elseif key == "gs_vecs"
            mismatch[key] = Dict(
                sector + 1 => ComplexF64[1.0, 0.0])
        else
            mismatch[key] = Dict(sector + 1 => 1e-12)
        end
        @test occursin(
            "$key sectors", validation_error(mismatch, expected, sector))
    end

    mismatch = valid_partial_metadata(expected, sector)
    mismatch["gs_vecs"] = Dict(sector => ComplexF64[2.0, 0.0])
    @test occursin(
        "not normalized", validation_error(mismatch, expected, sector))

    for residual in (1e-7, 1e-6, Inf, NaN)
        mismatch = valid_partial_metadata(expected, sector)
        mismatch["residuals"] = Dict(sector => residual)
        message = validation_error(mismatch, expected, sector)
        @test occursin("residual", message)
        @test occursin("finite and < 1e-7", message)
    end

    mismatch = valid_partial_metadata(expected, sector)
    mismatch["sector_ground_energies"] = Dict(sector => 1.1)
    @test occursin(
        "ground energy", validation_error(mismatch, expected, sector))
end

@testset "saved lattice fingerprint validation" begin
    physical = TiltedLat30()
    expected = (
        Np=13,
        t1=1.0,
        t3=0.2,
        V1=100.0,
        V2=0.0,
        V3=0.0,
        lattice_convention="corrected",
    )
    sector = 3
    metadata = valid_partial_metadata(expected, sector)
    delete!(metadata, "momentum_step")
    @test_throws ErrorException validate_manifold_partial(
        metadata, "synthetic.jld2", expected, sector)

    metadata["momentum_step"] = physical.kpoints[2]
    validated = validate_manifold_partial(
        metadata, "synthetic.jld2", expected, sector)
    @test validated.momentum_step ≈ physical.kpoints[2]

    metadata["translation_character_error"] = 1.0
    @test_throws ErrorException validate_manifold_partial(
        metadata, "synthetic.jld2", expected, sector)
end

@testset "candidate manifold rejects repeated sectors" begin
    repeated = try
        classify_candidate_manifold(
            Tuple{Int,Float64}[(0, 0.0), (0, 0.1), (1, 1.0)], 2)
        nothing
    catch error
        error
    end
    @test repeated isa ArgumentError
    if repeated isa Exception
        message = sprint(showerror, repeated)
        @test occursin("repeated sector 0", message)
        @test occursin("only one ground state per sector", message)
    end
end
