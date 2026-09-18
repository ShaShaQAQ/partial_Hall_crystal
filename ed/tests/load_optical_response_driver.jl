using Test

ENV["PHC_PREFLIGHT"] = "1"

@testset "large optical-response driver loads" begin
    include(joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "run_optical_response.jl"))
    @test isdefined(Main, :run_large_optical_response)
    @test RESPONSE_PREFLIGHT
    @test endswith(response_partial_path(5), "partial_5.jld2")
end
