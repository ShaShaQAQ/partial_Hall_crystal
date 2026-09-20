using Test
using LinearAlgebra
using InfiniteCylinderDMRG

import MPSKit

function small_paper_momentum_fixture()
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=2,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    first_state = mpskit_product_state(config, [1, 3])
    orthogonal_state = mpskit_product_state(config, [2, 4])
    return config, first_state, orthogonal_state
end

function signless_ring_translation(Ny::Int)
    modes = 2Ny
    dimension = 1 << modes
    result = zeros(ComplexF64, dimension, dimension)
    for bits in 0:(dimension - 1)
        target = 0
        for mode in 0:(modes - 1)
            iszero(bits & (1 << mode)) && continue
            target |= 1 << mod(mode + 2, modes)
        end
        result[target + 1, bits + 1] = 1
    end
    return result
end

function phase_redistributed_state(state)
    phases = cis.([0.17, -0.41, 0.29, -0.05])
    @assert length(state) == length(phases)
    @assert isapprox(prod(phases), 1; atol=1e-14, rtol=0)
    tensors = [phases[site] * state.AL[site] for site in 1:length(state)]
    return MPSKit.InfiniteMPS(tensors)
end

@testset "MPSKit fermionic transverse-translation MPO" begin
    config, _, _ = small_paper_momentum_fixture()
    finite_translation =
        InfiniteCylinderDMRG._mpskit_finite_transverse_translation_mpo(config)
    dense_translation = mpskit_dense_matrix(finite_translation)
    exact_translation = Matrix(fermionic_ring_translation(config))

    @test dense_translation ≈ exact_translation atol=1e-12 rtol=0
    @test dense_translation^config.Ny ≈
        Matrix{ComplexF64}(I, size(dense_translation)...) atol=1e-12 rtol=0

    # Occupying modes 1 and 3 wraps one fermion past the other under +a2.
    wrapped_bits = 0b0101
    translated_bits, wrapped_sign = translated_fock_basis(wrapped_bits, config.Ny)
    @test translated_bits == wrapped_bits
    @test wrapped_sign == -1
    @test dense_translation[wrapped_bits + 1, wrapped_bits + 1] ≈ -1

    # Mutation sensitivity: a bare bit permutation without the fermion sign is wrong.
    @test norm(dense_translation - signless_ring_translation(config.Ny)) > 1

    direction_config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=3,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    direction_translation = mpskit_dense_matrix(
        InfiniteCylinderDMRG._mpskit_finite_transverse_translation_mpo(
            direction_config,
        ),
    )
    @test direction_translation ≈
        Matrix(fermionic_ring_translation(direction_config)) atol=1e-12 rtol=0
end

@testset "MPSKit Schmidt-block transverse momentum" begin
    Ny = 3
    cyclic_translation = ComplexF64[
        0 0 1
        1 0 0
        0 1 0
    ]
    density = Matrix{ComplexF64}(I, 3, 3) / 3
    result = mpskit_schmidt_momentum_data(
        [density],
        [cyclic_translation];
        Ny,
        raw_charges=[0],
        residual_tol=1e-12,
    )

    @test result isa MPSKitSchmidtMomentumData
    @test result.valid
    @test result.spectrum.valid
    @test result.closure_residual < 1e-12
    @test result.spectrum.validation.unitary_residual < 1e-12
    @test result.spectrum.validation.power_residual < 1e-12
    @test result.spectrum.validation.commutator_residual < 1e-12
    @test result.spectrum.validation.transfer_residual == 0.0
    @test sort(getproperty.(result.spectrum.levels, :ky_index)) == [0, 1, 2]
    @test all(level -> level.phase_residual < 1e-12, result.spectrum.levels)

    open_density = Diagonal(ComplexF64[0.5, 0.5, 0.0])
    open_result = mpskit_schmidt_momentum_data(
        [open_density],
        [cyclic_translation];
        Ny,
        raw_charges=[0],
        residual_tol=1e-12,
        probability_tol=1e-13,
    )
    @test !open_result.valid
    @test !open_result.spectrum.valid
    @test open_result.closure_residual > 0.1
    @test occursin("closure", lowercase(open_result.reason))
    @test all(ismissing, getproperty.(open_result.spectrum.levels, :ky_index))

    rounded_only = copy(cyclic_translation)
    rounded_only[1, 3] *= cis(1e-3)
    rounded_result = mpskit_schmidt_momentum_data(
        [density],
        [rounded_only];
        Ny,
        raw_charges=[0],
        residual_tol=1e-12,
    )
    @test !rounded_result.valid
    @test rounded_result.spectrum.validation.power_residual > 1e-5
    @test all(ismissing, getproperty.(rounded_result.spectrum.levels, :ky_index))
end

@testset "MPSKit product-state momentum contraction" begin
    config, state, _ = small_paper_momentum_fixture()
    result = mpskit_momentum_entanglement_data(
        state,
        config;
        cut_x=1,
        tol=1e-10,
        maxiter=100,
        krylovdim=8,
        seed=0x4d4f4d,
        residual_tol=1e-9,
    )

    @test result isa MPSKitMomentumData
    @test result.valid
    @test result.spectrum.valid
    @test result.cut_x == 1
    @test result.bond == sites_per_cell(config)
    @test result.transfer_converged
    @test result.transfer_residual < 1e-9
    @test abs(abs(result.transfer_eigenvalue) - 1) < 1e-9
    @test result.closure_residual < 1e-9
    @test length(result.spectrum.levels) == 1
    @test only(result.spectrum.levels).ky_index == 0
end

@testset "MPSKit normalized mixed-transfer fidelity" begin
    config, state, orthogonal_state = small_paper_momentum_fixture()
    same = mpskit_mixed_transfer_fidelity(
        state,
        state,
        config;
        tol=1e-10,
        maxiter=100,
        krylovdim=8,
        seed=0x4d4958,
    )
    @test same isa MixedTransferFidelity
    @test same.valid
    @test same.converged
    @test same.fidelity_cell ≈ 1.0 atol=1e-10 rtol=0
    @test maximum((same.mixed_residual, same.self_residual1, same.self_residual2)) < 1e-9

    gauge_state = phase_redistributed_state(state)
    gauge = mpskit_mixed_transfer_fidelity(
        state,
        gauge_state,
        config;
        tol=1e-10,
        maxiter=100,
        krylovdim=8,
        seed=0x4d4958,
    )
    @test gauge.valid
    @test gauge.converged
    @test gauge.fidelity_cell ≈ 1.0 atol=1e-10 rtol=0
    @test maximum((gauge.mixed_residual, gauge.self_residual1, gauge.self_residual2)) < 1e-9

    orthogonal = mpskit_mixed_transfer_fidelity(
        state,
        orthogonal_state,
        config;
        tol=1e-10,
        maxiter=100,
        krylovdim=8,
        seed=0x4d4958,
    )
    @test orthogonal.valid
    @test orthogonal.converged
    @test orthogonal.fidelity_cell ≈ 0.0 atol=1e-12 rtol=0
    @test maximum((
        orthogonal.mixed_residual,
        orthogonal.self_residual1,
        orthogonal.self_residual2,
    )) < 1e-9
end
