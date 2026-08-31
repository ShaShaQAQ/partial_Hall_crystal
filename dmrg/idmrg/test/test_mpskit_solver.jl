using Test
using InfiniteCylinderDMRG

import MPSKit

function exact_dimer_solver_fixture(occupied_sites)
    config = InfiniteCylinderConfig(;
        Ly=4,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    physical_spaces = mpskit_physical_spaces(config)
    operators = InfiniteCylinderDMRG._mpskit_fermion_operators(
        first(physical_spaces);
        centered=true,
    )
    amplitude = -1.0 + 0.0im
    hopping = amplitude * operators.plus_min -
        conj(amplitude) * operators.min_plus
    hamiltonian = MPSKit.InfiniteMPOHamiltonian(
        physical_spaces,
        (1, 2) => hopping,
        (3, 4) => hopping,
    )
    state = mpskit_product_state(config, occupied_sites)
    return config, hamiltonian, state
end

@testset "MPSKit solver convergence uses the recomputed residual" begin
    @test mpskit_solver_converged(;
        idmrg_diagnostic=1e-2,
        recomputed_galerkin_residual=1e-9,
        energy_per_site=-0.5,
        energy_imaginary=1e-14,
        final_stage_reached=true,
        galerkin_tol=1e-7,
        energy_imag_tol=1e-12,
    )
    @test !mpskit_solver_converged(;
        idmrg_diagnostic=1e-12,
        recomputed_galerkin_residual=1e-3,
        energy_per_site=-0.5,
        energy_imaginary=0.0,
        final_stage_reached=true,
        galerkin_tol=1e-7,
        energy_imag_tol=1e-12,
    )
    @test !mpskit_solver_converged(;
        idmrg_diagnostic=1e-12,
        recomputed_galerkin_residual=1e-9,
        energy_per_site=NaN,
        energy_imaginary=0.0,
        final_stage_reached=true,
        galerkin_tol=1e-7,
        energy_imag_tol=1e-12,
    )
end

@testset "IDMRG2 growth and VUMPS Galerkin gate" begin
    _, hamiltonian, initial_state = exact_dimer_solver_fixture([1, 3])
    callbacks = Any[]
    result = run_mpskit_idmrg(
        hamiltonian,
        initial_state;
        maxdim_schedule=[4, 8],
        cutoff=1e-10,
        idmrg_maxiter=20,
        vumps_maxiter=50,
        galerkin_tol=1e-7,
        progress_callback=(state, environments, records) -> push!(
            callbacks,
            (
                state=state,
                environments=environments,
                records=copy(records),
            ),
        ),
    )

    @test maximum(result.link_dimensions) <= 8
    @test result.galerkin_residual <= 1e-7
    @test result.recomputed_galerkin_residual <= 1e-7
    @test isapprox(result.energy_per_site, -0.5; atol=1e-6, rtol=0)
    @test result.converged
    @test length(result.records) == 2
    @test getproperty.(result.records, :requested_maxdim) == [4, 8]
    @test all(
        record -> record.actual_maxdim <= record.requested_maxdim,
        result.records,
    )
    @test all(isfinite, getproperty.(result.records, :energy_per_site))
    @test all(isfinite, getproperty.(result.records, :elapsed_seconds))
    @test length(callbacks) == 2
    @test length(last(callbacks).records) == 2
end

@testset "independent dimer product states reach the same energy" begin
    _, hamiltonian, first_state = exact_dimer_solver_fixture([1, 3])
    _, _, second_state = exact_dimer_solver_fixture([2, 4])
    first_result = run_mpskit_idmrg(
        hamiltonian,
        first_state;
        maxdim_schedule=[4],
        cutoff=1e-10,
        idmrg_maxiter=20,
        vumps_maxiter=50,
        galerkin_tol=1e-7,
    )
    second_result = run_mpskit_idmrg(
        hamiltonian,
        second_state;
        maxdim_schedule=[4],
        cutoff=1e-10,
        idmrg_maxiter=20,
        vumps_maxiter=50,
        galerkin_tol=1e-7,
    )
    @test first_result.converged
    @test second_result.converged
    @test isapprox(
        first_result.energy_per_site,
        second_result.energy_per_site;
        atol=1e-8,
        rtol=0,
    )
    @test isapprox(first_result.energy_per_site, -0.5; atol=1e-6, rtol=0)
end

@testset "progress callback receives the last valid stage before failure" begin
    _, hamiltonian, initial_state = exact_dimer_solver_fixture([1, 3])
    captured = Ref{Any}(nothing)
    @test_throws ErrorException run_mpskit_idmrg(
        hamiltonian,
        initial_state;
        maxdim_schedule=[2],
        cutoff=1e-8,
        idmrg_maxiter=2,
        vumps_maxiter=2,
        galerkin_tol=1e-7,
        progress_callback=(state, environments, records) -> begin
            captured[] = (
                state=state,
                environments=environments,
                records=copy(records),
            )
            error("intentional callback failure")
        end,
    )
    @test captured[] !== nothing
    @test length(captured[].records) == 1
    @test isfinite(last(captured[].records).energy_per_site)
end
