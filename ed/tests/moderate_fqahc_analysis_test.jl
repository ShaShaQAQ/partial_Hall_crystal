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
