using Test
using LinearAlgebra
using InfiniteCylinderDMRG

@testset "momentum label validation" begin
    rho = Diagonal([0.4, 0.3, 0.2, 0.1])
    Ty = Diagonal(ComplexF64[1, cis(pi), 1, cis(pi)])
    result = validate_momentum_blocks(rho, Ty; Ny=2, residual_tol=1e-12)

    @test result.valid
    @test result.validation.valid
    @test all(level -> level.ky_index in 0:1, result.levels)
    @test maximum(level.phase_residual for level in result.levels) < 1e-12
    @test result.validation.unitary_residual < 1e-12
    @test result.validation.power_residual < 1e-12
    @test result.validation.commutator_residual < 1e-12

    bad_unitary = validate_momentum_blocks(rho, 1.01Ty; Ny=2, residual_tol=1e-12)
    @test !bad_unitary.valid
    @test all(level -> ismissing(level.ky_index), bad_unitary.levels)

    broken = ComplexF64[1 0.1; 0 1]
    broken_result = validate_momentum_blocks(
        Diagonal([0.6, 0.4]),
        broken;
        Ny=2,
        residual_tol=1e-12,
    )
    @test !broken_result.valid
    @test all(level -> ismissing(level.ky_index), broken_result.levels)
end

@testset "degenerate charge blocks and strict invalidation" begin
    rho = Diagonal([0.4, 0.4, 0.1, 0.1])
    Ty = ComplexF64[
        0 1 0 0
        1 0 0 0
        0 0 1 0
        0 0 0 -1
    ]
    result = validate_momentum_blocks(
        rho,
        Ty;
        Ny=2,
        raw_charges=[0, 0, 1, 1],
        residual_tol=1e-12,
    )
    @test result.valid
    @test count(level -> level.ky_index == 0, result.levels) == 2
    @test count(level -> level.ky_index == 1, result.levels) == 2
    @test [level.raw_charge for level in result.levels] == [0, 0, 1, 1]

    charge_mixing = validate_momentum_blocks(
        Diagonal([0.5, 0.5]),
        ComplexF64[0 1; 1 0];
        Ny=2,
        raw_charges=[0, 1],
        residual_tol=1e-12,
    )
    @test !charge_mixing.valid
    @test occursin("charge", lowercase(charge_mixing.reason))
    @test all(level -> ismissing(level.ky_index), charge_mixing.levels)

    transfer_failure = validate_momentum_blocks(
        Diagonal([0.6, 0.4]),
        Diagonal(ComplexF64[1, -1]);
        Ny=2,
        transfer_residual=1e-4,
        residual_tol=1e-12,
    )
    @test !transfer_failure.valid
    @test occursin("transfer", lowercase(transfer_failure.reason))
    @test all(level -> ismissing(level.ky_index), transfer_failure.levels)
end

@testset "momentum validation input checks" begin
    rho = Diagonal([0.6, 0.4])
    Ty = Diagonal(ComplexF64[1, -1])
    @test_throws ArgumentError validate_momentum_blocks(rho, Ty; Ny=0)
    @test_throws ArgumentError validate_momentum_blocks(rho, Ty; Ny=true)
    @test_throws ArgumentError validate_momentum_blocks(rho, Ty; Ny=2, residual_tol=0)
    @test_throws DimensionMismatch validate_momentum_blocks(rho, Matrix{ComplexF64}(I, 3, 3); Ny=2)
    @test_throws DimensionMismatch validate_momentum_blocks(
        rho,
        Ty;
        Ny=2,
        raw_charges=[0],
    )
end

@testset "Ny=2 tensor momentum smoke" begin
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=2,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    _, _, psi = initial_infinite_mps(config; occupied_sites=[1, 3])
    result = momentum_entanglement_data(
        psi,
        config;
        cut_x=1,
        residual_tol=1e-8,
        transfer_tol=1e-8,
    )

    @test result isa MomentumEntanglementData
    @test !isempty(result.levels)
    @test if result.valid
        all(level -> !ismissing(level.ky_index), result.levels)
    else
        !isempty(strip(result.reason)) &&
            all(level -> ismissing(level.ky_index), result.levels)
    end

    legacy = InfiniteCylinderConfig(; Ly=4, x_period=1, filling_num=1, filling_den=1)
    _, _, legacy_psi = initial_infinite_mps(legacy; occupied_sites=[1, 3])
    @test_throws ArgumentError momentum_entanglement_data(
        legacy_psi,
        legacy;
        cut_x=1,
    )
end
