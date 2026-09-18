#!/usr/bin/env julia

using Test

include("flux_segment_utils.jl")

@testset "segmented flux helpers" begin
    @test point_dir("out", 7) == joinpath("out", "phi_007")
    @test checkpoint_path("out", 7) == joinpath("out", "checkpoints", "state_007.jls")

    @test default_phi_for_step(0; cycles=3, steps_per_2pi=24, direction="forward") == 0.0
    @test default_phi_for_step(18; cycles=3, steps_per_2pi=24, direction="forward") ≈ 1.5π
    @test default_phi_for_step(0; cycles=3, steps_per_2pi=24, direction="backward") ≈ 6π
    @test default_phi_for_step(18; cycles=3, steps_per_2pi=24, direction="backward") ≈ 4.5π

    tmp = mktempdir()
    write(
        joinpath(tmp, "pumping_segment_003_005.dat"),
        """
        # step phi energy cumulative_left_edge cumulative_mid cumulative_right_edge actual_sweeps attempts converged reason
        3 0.3 -3.0 0.01 0.03 0.0 10 1 true energy_density_truncerr
        5 0.5 -5.0 0.01 0.05 0.0 12 1 true energy_density_truncerr
        """,
    )
    write(
        joinpath(tmp, "pumping_segment_004_004.dat"),
        """
        # step phi energy cumulative_left_edge cumulative_mid cumulative_right_edge actual_sweeps attempts converged reason
        4 0.4 -4.0 0.01 0.04 0.0 11 1 true energy_density_truncerr
        """,
    )
    rows = collect_pumping_rows(tmp)
    @test [row.step for row in rows] == [3, 4, 5]
    refresh_pumping_index(tmp)
    merged = readlines(joinpath(tmp, "pumping.dat"))
    @test length(filter(line -> !startswith(line, "#"), merged)) == 3
    @test occursin("4 0.4 -4.0", join(merged, "\n"))
end
