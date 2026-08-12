using Test
using ITensors
using ITensorMPS
using ITensorInfiniteMPS
using InfiniteCylinderDMRG

@testset "VUMPS runner helpers" begin
    cfg = InfiniteCylinderConfig(; Ly=2, x_period=3)
    sites, initstate, psi = initial_infinite_mps(cfg)
    @test length(sites) == sites_per_cell(cfg)
    @test flux(psi.AL) == QN()
    @test maximum(link_dimensions(psi)) == 1

    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=Int[])
    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=[1, 1])
    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=[0])
    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=[7])
    @test_throws ArgumentError initial_infinite_mps(cfg; occupied_sites=[true])

    @test centered_site_charges(cfg) == (empty=-1, occupied=5)
    @test isnothing(validate_centered_site_charges(cfg, sites))
    uncentered = siteinds("Fermion", sites_per_cell(cfg); conserve_qns=true)
    @test_throws ArgumentError validate_centered_site_charges(cfg, uncentered)
    dense_sites = siteinds("Fermion", sites_per_cell(cfg))
    @test_throws ArgumentError validate_centered_site_charges(cfg, dense_sites)
    @test_throws ArgumentError validate_centered_site_charges(cfg, sites[1:(end - 1)])

    @test vumps_converged(
        9e-6,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        2e-5,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        9e-6,
        2e-4,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        9e-6,
        8e-6,
        2e-4,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        9e-6,
        8e-6,
        8e-6,
        1;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        NaN,
        8e-6,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test !vumps_converged(
        -1e-6,
        8e-6,
        8e-6,
        2;
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
    @test_throws ArgumentError vumps_converged(
        1e-6,
        1e-6,
        1e-6,
        2;
        vumps_tol=0.0,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )
end

@testset "pinned VUMPS iteration boundary" begin
    initspin(_) = "↑"
    sites = infsiteinds("S=1/2", 1; initstate=initspin)
    psi = InfMPS(sites, initspin)
    H = InfiniteSum{MPO}(Model("ising"), sites; J=1.0, h=1.2)

    step = vumps_iteration(H, psi; vumps_tol=1e-5)
    @test step.psi isa InfiniteCanonicalMPS
    @test length(step.eps_left) == 1
    @test length(step.eps_right) == 1
    @test all(isfinite, step.eps_left)
    @test all(isfinite, step.eps_right)
    @test step.energy_left isa Float64
    @test step.energy_right isa Float64
    @test isfinite(step.elapsed_seconds)
    @test step.elapsed_seconds >= 0
    @test_throws ArgumentError vumps_iteration(
        H, psi; vumps_tol=1e-5, imaginary_tol=-1e-12
    )
end

@testset "actual cylinder complex-energy boundary" begin
    cfg = InfiniteCylinderConfig(; Ly=2, x_period=3)
    params = CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.0, V3=0.0)
    sites, _, psi = initial_infinite_mps(cfg)
    H = build_infinite_mpo(cfg, params, sites)

    result = run_vumps(
        H,
        psi;
        maxdim_schedule=[1],
        cutoff=1e-8,
        max_iterations=1,
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
        imaginary_tol=1e-12,
    )
    @test !result.converged
    @test length(result.records) == 1
    @test result.records[1].stage == 1
    @test result.records[1].iteration == 1
    @test isfinite(result.records[1].energy_left)
    @test isfinite(result.records[1].energy_right)
    @test result.records[1].energy_left ≈ 0.4 atol=1e-12
    @test result.records[1].energy_right ≈ 0.4 atol=1e-12
    @test result.records[1].energy_mismatch <= 1e-12
    @test_throws ArgumentError run_vumps(
        H,
        psi;
        maxdim_schedule=[1],
        cutoff=1e-8,
        max_iterations=1,
        vumps_tol=1e-5,
        imaginary_tol=Inf,
    )
end

@testset "checked subspace expansion" begin
    initspin(_) = "↑"
    sites = infsiteinds("S=1/2", 1; initstate=initspin)
    psi = InfMPS(sites, initspin)
    H = InfiniteSum{MPO}(Model("ising"), sites; J=1.0, h=1.2)
    expanded = expand_subspace(psi, H, 2; cutoff=1e-8)
    @test maximum(link_dimensions(expanded)) > maximum(link_dimensions(psi))

    stalled_H = InfiniteSum{MPO}(Model("ising"), sites; J=0.0, h=0.0)
    error = try
        expand_subspace(psi, stalled_H, 2; cutoff=1e-8)
        nothing
    catch caught
        caught
    end
    @test error isa ErrorException
    @test occursin("target maxdim=2", sprint(showerror, error))
    @test expand_subspace(expanded, H, 2; cutoff=1e-8) === expanded
    @test_throws ArgumentError expand_subspace(psi, H, 0; cutoff=1e-8)
    @test_throws ArgumentError expand_subspace(psi, H, 2; cutoff=0.0)
end

@testset "staged VUMPS runner" begin
    initspin(_) = "↑"
    sites = infsiteinds("S=1/2", 1; initstate=initspin)
    psi = InfMPS(sites, initspin)
    H = InfiniteSum{MPO}(Model("ising"), sites; J=1.0, h=1.2)
    options = (
        maxdim_schedule=[1],
        cutoff=1e-8,
        vumps_tol=1e-5,
        energy_tol=1e-4,
        energy_mismatch_tol=1e-4,
        stable_iterations=2,
    )

    stopped = run_vumps(H, psi; options..., max_iterations=1)
    @test !stopped.converged
    @test stopped.reason == "stage 1 reached maximum iterations (1) without convergence"
    @test length(stopped.records) == 1
    @test stopped.records[1].stage == 1
    @test stopped.records[1].iteration == 1
    @test !stopped.records[1].converged

    expanded_stop = run_vumps(
        H, psi; options..., maxdim_schedule=[2], max_iterations=1
    )
    @test !expanded_stop.converged
    @test expanded_stop.records[1].stage == 1
    @test expanded_stop.records[1].maxlinkdim == 2

    converged = run_vumps(H, psi; options..., max_iterations=3)
    @test converged.converged
    @test converged.reason == "converged after 1 stage"
    @test length(converged.records) == 3
    @test converged.records[end].converged
    @test converged.records[end].precision_error < 1e-5
    @test converged.records[end].delta_energy < 1e-4
    @test converged.records[end].energy_mismatch < 1e-4

    @test_throws ArgumentError run_vumps(
        H, psi; options..., maxdim_schedule=Int[], max_iterations=1
    )
    @test_throws ArgumentError run_vumps(
        H, psi; options..., maxdim_schedule=[2, 1], max_iterations=1
    )
    @test_throws ArgumentError run_vumps(
        H, psi; options..., maxdim_schedule=[true], max_iterations=1
    )
    @test_throws ArgumentError run_vumps(H, psi; options..., max_iterations=0)
    @test_throws ArgumentError run_vumps(H, psi; options..., max_iterations=true)
    @test_throws ArgumentError run_vumps(
        H, psi; options..., stable_iterations=0, max_iterations=1
    )
    @test_throws ArgumentError run_vumps(
        H, psi; options..., stable_iterations=true, max_iterations=1
    )
end

@testset "VUMPS records and energy normalization" begin
    cfg = InfiniteCylinderConfig(; Ly=2, x_period=3)
    _, _, psi = initial_infinite_mps(cfg)
    record = VUMPSRecord(
        2, 3, 4, -1.0, -0.9, 0.1, 0.01, 1e-6, 2e-6, 2e-6, 0.25, false
    )
    @test record.stage == 2
    @test record.iteration == 3
    @test record.maxlinkdim == 4
    @test record.energy_mismatch == 0.1
    result = VUMPSResult(psi, [record], false, "maximum iterations reached")
    @test result.psi === psi
    @test result.records == [record]
    @test !result.converged
    @test result.reason == "maximum iterations reached"

    @test unit_cell_energy(2.5, 6) == 15.0
    @test unit_cell_energy([1.0, 2.0, 3.0], 3) == 6.0
    @test unit_cell_energy(2.5 + 1e-13im, 6; imaginary_tol=1e-12) == 15.0
    @test unit_cell_energy(
        ComplexF64[1.0 + 1e-13im, 2.0 - 2e-13im, 3.0],
        3;
        imaginary_tol=1e-12,
    ) == 6.0
    @test_throws ArgumentError unit_cell_energy([1.0, 2.0], 3)
    @test_throws ArgumentError unit_cell_energy(Float64[], 3)
    @test_throws ArgumentError unit_cell_energy(Inf, 3)
    @test_throws ArgumentError unit_cell_energy(floatmax(Float64), 2)
    @test_throws ArgumentError unit_cell_energy(
        [floatmax(Float64), floatmax(Float64)], 2
    )
    @test_throws ArgumentError unit_cell_energy(
        BigFloat(floatmax(Float64)) * 2, 1
    )
    @test_throws ArgumentError unit_cell_energy(
        2.5 + 2e-12im, 6; imaginary_tol=1e-12
    )
    @test_throws ArgumentError unit_cell_energy(
        ComplexF64[1.0, 2.0 + 2e-12im, 3.0], 3; imaginary_tol=1e-12
    )
    for invalid_energy in (ComplexF64(Inf, 0.0), ComplexF64(0.0, Inf), ComplexF64(NaN, 0.0))
        @test_throws ArgumentError unit_cell_energy(invalid_energy, 3)
        @test_throws ArgumentError unit_cell_energy([invalid_energy], 1)
    end
    oversized_tol = BigFloat(floatmax(Float64)) * 2
    for invalid_tol in (-1e-12, NaN, Inf, oversized_tol)
        @test_throws ArgumentError unit_cell_energy(
            1.0 + 0.0im, 1; imaginary_tol=invalid_tol
        )
        @test_throws ArgumentError unit_cell_energy(
            ComplexF64[1.0], 1; imaginary_tol=invalid_tol
        )
    end
end
