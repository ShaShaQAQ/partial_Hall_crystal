using Test

@testset "ground-manifold analyzer loads" begin
    driver = joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "analyze_ground_manifold.jl")
    @test isfile(driver)
    if isfile(driver)
        include(driver)
        @test isdefined(Main, :run_ground_manifold_analysis)
        @test MANIFOLD_SIZE == 15
        @test MANIFOLD_LATTICE == "corrected"
        @test basename(MANIFOLD_RESULT_DIR) ==
              "phc30_np12_v1_10_v2_2_v3_2_fqahc"
        @test manifold_partial_paths() == [
            joinpath(MANIFOLD_DATA_DIR, "partial_$sector.jld2")
            for sector in 0:14
        ]
    end
end
