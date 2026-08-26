include("../case4_4x6_Np4/benchmark_optical_response.jl")

@testset "4x6 optical response integration" begin
    result = run_optical_benchmark(
        omegas=collect(0.05:0.10:2.05),
        eta=0.065,
        requested_m=[50, 100, 200, 400],
        write_outputs=false,
        make_plots=false,
    )
    @test result.sector == 2
    @test result.ground_residual < 1e-10
    @test result.hermiticity.H < 1e-12
    @test result.hermiticity.Jx < 1e-12
    @test result.hermiticity.Kxx < 1e-12
    @test abs(result.source_overlap) < 1e-11
    @test all(diff(result.errors) .<= 1e-10)
    @test result.errors[end] < 1e-6
end
