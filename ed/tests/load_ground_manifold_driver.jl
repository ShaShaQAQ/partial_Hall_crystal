using Test

@testset "ground-manifold analyzer loads" begin
    driver = joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "analyze_ground_manifold.jl")
    @test isfile(driver)
    if isfile(driver)
        include(driver)
        @test isdefined(Main, :run_ground_manifold_analysis)
        @test isdefined(Main, :MANIFOLD_NP)
        @test isdefined(Main, :MANIFOLD_RESULT_ID)
        @test isdefined(Main, :MANIFOLD_PHASE_LABEL)
        @test MANIFOLD_SIZE == 15
        @test MANIFOLD_LATTICE == "corrected"
        if isdefined(Main, :MANIFOLD_NP)
            @test MANIFOLD_NP == parse(
                Int, parse_manifold_argument("--Np", "12"))
        end
        if isdefined(Main, :MANIFOLD_RESULT_ID)
            @test MANIFOLD_RESULT_ID == parse_manifold_argument(
                "--result-id", "phc30_np12_v1_10_v2_2_v3_2_fqahc")
        end
        if isdefined(Main, :MANIFOLD_PHASE_LABEL)
            @test MANIFOLD_PHASE_LABEL == parse_manifold_argument(
                "--phase-label", "FQAHC candidate")
        end
        @test manifold_partial_paths() == [
            joinpath(MANIFOLD_DATA_DIR, "partial_$sector.jld2")
            for sector in 0:14
        ]
    end
end
