using Test

const ANALYSIS_PATH = joinpath(@__DIR__, "analyze_ground_state.jl")

@testset "ground-state diagnostics entry point exists" begin
    @test isfile(ANALYSIS_PATH)
end

if isfile(ANALYSIS_PATH)
    include(ANALYSIS_PATH)
    using .GroundStateDiagnostics

    function write_fixture(path::AbstractString, body::AbstractString)
        mkpath(dirname(path))
        open(path, "w") do io
            write(io, body)
        end
    end

    function fixture_output()
        root = mktempdir()
        write_fixture(
            joinpath(root, "summary.dat"),
            """
            energy 1.25
            Lx 2
            Ly 2
            Ns 4
            Np 2
            converged true
            final_delta_energy 1e-8
            final_max_density_delta 1e-6
            final_maxerr 1e-10
            checkpoint_out checkpoint_final.jls
            """,
        )
        write_fixture(
            joinpath(root, "density.dat"),
            """
            # site x y rx ry density
            1 0 0 0.0 0.0 0.8
            2 0 1 0.5 0.866025403784 0.2
            3 1 0 1.0 0.0 0.8
            4 1 1 1.5 0.866025403784 0.2
            """,
        )
        write_fixture(
            joinpath(root, "density_x.dat"),
            """
            # x y_averaged_density
            0 0.5
            1 0.5
            """,
        )
        write_fixture(
            joinpath(root, "connected_density.dat"),
            """
            # i j <n_i n_j>-<n_i><n_j>
            1 1 0.16
            1 2 -0.16
            1 3 0.16
            1 4 -0.16
            2 1 -0.16
            2 2 0.16
            2 3 -0.16
            2 4 0.16
            3 1 0.16
            3 2 -0.16
            3 3 0.16
            3 4 -0.16
            4 1 -0.16
            4 2 0.16
            4 3 -0.16
            4 4 0.16
            """,
        )
        write_fixture(
            joinpath(root, "green.dat"),
            """
            # i j Re(<c_i^dag c_j>) Im(<c_i^dag c_j>)
            1 1 1.0 0.0
            1 2 0.0 0.0
            1 3 0.25 0.0
            1 4 0.0 0.0
            2 1 0.0 0.0
            2 2 1.0 0.0
            2 3 0.0 0.0
            2 4 0.5 0.0
            3 1 0.25 0.0
            3 2 0.0 0.0
            3 3 1.0 0.0
            3 4 0.0 0.0
            4 1 0.0 0.0
            4 2 0.5 0.0
            4 3 0.0 0.0
            4 4 1.0 0.0
            """,
        )
        write_fixture(
            joinpath(root, "convergence.dat"),
            """
            # sweep energy delta_energy max_density_delta maxlinkdim maxerr converged reason
            1 1.30 Inf Inf 4 1e-8 false min_sweeps
            2 1.25 0.05 1e-6 4 1e-10 true energy_density_truncerr
            """,
        )
        return root
    end

    @testset "diagnostics parse finite-DMRG text outputs" begin
        root = fixture_output()
        report = analyze_output(root)

        @test report.summary["converged"] == "true"
        @test report.density.total_density ≈ 2.0 atol=1e-12
        @test report.density.mean_density ≈ 0.5 atol=1e-12
        @test report.density.cdw_contrast ≈ 0.6 atol=1e-12
        @test report.density.sublattice[0].mean ≈ 0.8 atol=1e-12
        @test report.density.sublattice[1].mean ≈ 0.2 atol=1e-12
        @test report.convergence[end].converged
        @test report.convergence[end].reason == "energy_density_truncerr"
    end

    @testset "diagnostics compute structure-factor and Green-decay summaries" begin
        root = fixture_output()
        report = analyze_output(root)

        peak = first(report.structure_factor_peaks)
        @test peak.my == 1
        @test peak.value ≈ 0.64 atol=1e-12

        dx1 = only(filter(row -> row.dx == 1, report.green_decay))
        @test dx1.mean_abs ≈ 0.375 atol=1e-12
        @test dx1.count == 4
    end

    @testset "diagnostics CLI writes portable artifacts" begin
        root = fixture_output()
        outdir = joinpath(root, "diagnostics")
        GroundStateDiagnostics.main(["--input=$root", "--outdir=$outdir", "--top_peaks=3"])

        @test isfile(joinpath(outdir, "diagnostics_summary.dat"))
        @test isfile(joinpath(outdir, "density_by_sublattice.csv"))
        @test isfile(joinpath(outdir, "structure_factor.csv"))
        @test isfile(joinpath(outdir, "structure_factor_peaks.csv"))
        @test isfile(joinpath(outdir, "green_decay_x.csv"))
        @test isfile(joinpath(outdir, "ground_state_diagnostics.svg"))

        summary = read(joinpath(outdir, "diagnostics_summary.dat"), String)
        @test occursin("density_total 2", summary)
        @test occursin("structure_peak_1_my 1", summary)
    end

    @testset "diagnostics can read a saved MPS checkpoint for entanglement" begin
        using ITensors
        using ITensorMPS
        using Random
        using Serialization

        root = fixture_output()
        sites = siteinds("Fermion", 4; conserve_qns=true)
        state = ["Occ", "Emp", "Occ", "Emp"]
        psi = random_mps(MersenneTwister(11), ComplexF64, sites, state; linkdims=2)
        checkpoint = joinpath(root, "checkpoint_final.jls")
        serialize(checkpoint, (psi=psi,))

        report = analyze_output(root; checkpoint, entanglement_bonds="center", spectrum_levels=4)

        @test isempty(report.entanglement_error)
        @test length(report.entanglement) == 1
        @test report.entanglement[1].bond == 2
        @test isfinite(report.entanglement[1].entropy)
        @test report.entanglement[1].schmidt_rank >= 1
        @test !isempty(report.entanglement_spectrum)
    end
end
