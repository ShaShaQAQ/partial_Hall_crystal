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

@testset "candidate manifold input validation" begin
    @test_throws ArgumentError classify_candidate_manifold(
        Tuple{Int,Float64}[(0, 0.0)], 1)
    @test_throws ArgumentError classify_candidate_manifold(
        Tuple{Int,Float64}[(0, 0.0), (1, 1.0)], 0)
end

@testset "saved lattice fingerprint validation" begin
    physical = TiltedLat30()
    metadata = Dict{String,Any}(
        "Np" => 12,
        "t1" => 1.0,
        "t3" => 0.2,
        "V1" => 10.0,
        "V2" => 2.0,
        "V3" => 2.0,
        "lattice_convention" => "corrected",
    )
    @test_throws ErrorException validate_manifold_partial(
        metadata, "synthetic.jld2")

    metadata["momentum_step"] = physical.kpoints[2]
    metadata["translation_character_error"] = 1e-14
    metadata["hopping_count"] = 300
    validated = validate_manifold_partial(metadata, "synthetic.jld2")
    @test validated.hopping_count == 300
    @test validated.translation_character_error < 1e-12

    metadata["translation_character_error"] = 1.0
    @test_throws ErrorException validate_manifold_partial(
        metadata, "synthetic.jld2")
end
