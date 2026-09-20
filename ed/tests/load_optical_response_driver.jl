using Test

ENV["PHC_PREFLIGHT"] = "1"

@testset "large optical-response driver loads" begin
    include(joinpath(
        @__DIR__, "..", "cases", "case3_30sites_Np12",
        "run_optical_response.jl"))
    @test isdefined(Main, :run_large_optical_response)
    @test isdefined(Main, :RESPONSE_NP)
    if isdefined(Main, :RESPONSE_NP)
        expected_np = parse(
            Int, parse_response_argument("--Np", "12"))
        @test RESPONSE_NP == expected_np
    end
    @test RESPONSE_PREFLIGHT
    @test RESPONSE_LATTICE == "legacy"
    @test response_lattice("legacy") isa GenLat
    @test response_lattice("corrected") isa GenLat
    @test response_lattice("legacy").kpoints ≈
          response_lattice("corrected").kpoints atol=1e-14
    @test_throws ArgumentError response_lattice("unknown")
    @test endswith(
        response_partial_path(7, "legacy"), "partial_5.jld2")
    @test endswith(
        response_partial_path(7, "corrected"), "partial_7.jld2")
end

@testset "large optical-response output records model parameters" begin
    mktempdir() do directory
        output_path = joinpath(directory, "response.jld2")
        frequencies = [0.0, 0.1]
        curve = (
            drude_weight=0.25,
            total=ComplexF64[1 + 2im, 3 + 4im],
            regular=ComplexF64[0.5 + 0.2im, 0.4 + 0.1im],
            drude=ComplexF64[0.5 + 1.8im, 2.6 + 3.9im],
        )
        kernel = (alpha=[1.0, 2.0], beta=[0.3], breakdown=false)
        expected_np = parse(
            Int, parse_response_argument("--Np", "12"))
        model = (Np=expected_np, t1=1.0, t3=0.2,
                 V1=10.0, V2=2.0, V3=2.0)

        save_response_curve(
            output_path, 5, 1.25, 1e-12, 0.75, 0.5, 12sqrt(3),
            0.065, frequencies, curve, kernel, "corrected";
            model_parameters=model, requested_mmax=600)

        saved = load(output_path)
        @test saved["Np"] == expected_np
        @test saved["t1"] == 1.0
        @test saved["t3"] == 0.2
        @test saved["V1"] == 10.0
        @test saved["V2"] == 2.0
        @test saved["V3"] == 2.0
        @test saved["requested_mmax"] == 600
        @test saved["lanczos_steps"] == 2
        @test !saved["lanczos_breakdown"]
        header = first(readlines(replace(output_path, ".jld2" => ".dat")))
        @test occursin("drude_weight=0.25", header)
        @test occursin("lanczos_breakdown=false", header)
    end
end
