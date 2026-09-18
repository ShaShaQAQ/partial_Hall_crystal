include("DMRGFullModel.jl")
using .DMRGFullModel
using Test

@testset "adaptive convergence gating" begin
    @test !convergence_satisfied(
        5,
        1e-9,
        1e-9,
        1e-9;
        min_sweeps=1,
        target_sweep=6,
        stable_sweeps=1,
        satisfied_streak=1,
        energy_tol=1e-8,
        density_tol=1e-8,
        truncerr_tol=1e-8,
    )

    @test !convergence_satisfied(
        6,
        1e-9,
        1e-9,
        1e-9;
        min_sweeps=1,
        target_sweep=6,
        stable_sweeps=3,
        satisfied_streak=2,
        energy_tol=1e-8,
        density_tol=1e-8,
        truncerr_tol=1e-8,
    )

    @test convergence_satisfied(
        6,
        1e-9,
        1e-9,
        1e-9;
        min_sweeps=1,
        target_sweep=6,
        stable_sweeps=3,
        satisfied_streak=3,
        energy_tol=1e-8,
        density_tol=1e-8,
        truncerr_tol=1e-8,
    )

    @test !convergence_satisfied(
        6,
        1e-6,
        1e-9,
        1e-9;
        min_sweeps=1,
        target_sweep=6,
        stable_sweeps=1,
        satisfied_streak=1,
        energy_tol=1e-8,
        density_tol=1e-8,
        truncerr_tol=1e-8,
    )
end
