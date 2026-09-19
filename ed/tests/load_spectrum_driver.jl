using Test

@testset "30-site spectrum driver loads" begin
    driver = joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "run_spectrum.jl")
    @test isfile(driver)
    if isfile(driver)
        include(driver)
        @test isdefined(Main, :run_large_spectrum)
        @test SPECTRUM_NP == 12
        @test SPECTRUM_V1 == 10.0
        @test SPECTRUM_V2 == 2.0
        @test SPECTRUM_V3 == 2.0
        @test SPECTRUM_LATTICE == "corrected"
        @test endswith(spectrum_partial_path(5), "partial_5.jld2")
    end
end
