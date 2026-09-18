using Test

@testset "CDW structure-factor driver loads" begin
    include(joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "compute_cdw_structure_factor.jl"))
    @test isdefined(Main, :run_cdw_structure_factor)
    @test CDW_GROUND_SECTORS == (0, 5, 10)
end
