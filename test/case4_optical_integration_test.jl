include("../case4_4x6_Np4/benchmark_optical_response.jl")

@testset "4x6 interacting Ax derivatives" begin
    lat = RectLat4x6()
    basis = gen_basis(lat.Ns, 4)
    sec = build_ksector(basis, lat, 2)
    hops, jx_hops, kxx_hops = build_hops_x_derivatives(lat, 1.0, 0.2)
    H0 = materialize_sector_operator(sec, lat, hops, 1.0, 0.0, 0.0)
    Jx = materialize_sector_operator(sec, lat, jx_hops, 0.0, 0.0, 0.0)
    Kxx = materialize_sector_operator(sec, lat, kxx_hops, 0.0, 0.0, 0.0)

    delta1 = 1e-6
    delta2 = 1e-4
    Hp1 = materialize_sector_operator(
        sec, lat, build_hops_Ax(lat, 1.0, 0.2, delta1),
        1.0, 0.0, 0.0)
    Hm1 = materialize_sector_operator(
        sec, lat, build_hops_Ax(lat, 1.0, 0.2, -delta1),
        1.0, 0.0, 0.0)
    Hp2 = materialize_sector_operator(
        sec, lat, build_hops_Ax(lat, 1.0, 0.2, delta2),
        1.0, 0.0, 0.0)
    Hm2 = materialize_sector_operator(
        sec, lat, build_hops_Ax(lat, 1.0, 0.2, -delta2),
        1.0, 0.0, 0.0)

    relerr(actual, reference) =
        norm(actual - reference) / max(norm(reference), 1.0)
    @test relerr((Hp1 - Hm1) / (2delta1), Jx) < 1e-8
    @test relerr((Hp2 - 2H0 + Hm2) / delta2^2, Kxx) < 1e-5
end

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
    @test result.Eg ≈ -9.3135758879616 atol=1e-11
    @test result.source_norm2 ≈ 10.844308142783673 atol=1e-8
    @test result.Kexp ≈ 4.479235242914201 atol=1e-10
    @test result.operator_convention ==
          "H_A(k)=H(k+A_x*xhat); Jx=dH/dA_x; Kxx=d2H/dA_x2"

    function explicit_regular(data, area, omega, eta)
        z = complex(omega, eta)
        Gz = sum(data.weights ./ (z .- data.deltas))
        Gminusz = sum(data.weights ./ (-z .- data.deltas))
        chi0 = -2sum(data.weights ./ data.deltas)
        return (2pi * im / area) * (Gz + Gminusz - chi0) / z
    end
    explicit_curve = [explicit_regular(
                          result.exact_data, result.area, omega, result.eta)
                      for omega in result.omegas]
    @test scaled_max_error(
        explicit_curve, result.exact_curve.regular) < 1e-12
end

@testset "optical benchmark output" begin
    mktempdir() do output_dir
        result = run_optical_benchmark(
            omegas=collect(0.05:0.20:1.05),
            eta=0.065,
            requested_m=[50],
            write_outputs=true,
            make_plots=false,
            output_dir=output_dir,
        )
        @test isfile(joinpath(output_dir, "optical_response_benchmark.jld2"))
        @test isfile(joinpath(output_dir, "optical_response_curves.dat"))
        @test isfile(joinpath(output_dir, "lanczos_convergence.dat"))
        @test result.sector == 2
        saved = Base.invokelatest(
            JLD2.load,
            joinpath(output_dir, "optical_response_benchmark.jld2"))
        @test saved["operator_convention"] ==
              "H_A(k)=H(k+A_x*xhat); Jx=dH/dA_x; Kxx=d2H/dA_x2"
    end
end
