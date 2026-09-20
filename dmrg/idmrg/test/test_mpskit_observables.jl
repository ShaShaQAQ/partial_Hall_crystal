using Test
using LinearAlgebra
using InfiniteCylinderDMRG

import MPSKit
import TensorKit

function paper_product_observable_fixture()
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=6,
        x_period=3,
        filling_num=7,
        filling_den=9,
    )
    parameters = CylinderModelParams(;
        t1=1.0,
        t3=0.2,
        V1=1.0,
        V2=0.0,
        V3=0.0,
    )
    candidate = first(fig2_initial_candidates(config))
    state = mpskit_product_state(config, candidate.occupied_sites)
    hamiltonian = mpskit_infinite_hamiltonian(config, parameters)
    return config, parameters, candidate, state, hamiltonian
end

function direct_product_energy_per_cell(config, parameters, occupied_sites)
    occupied = Set(occupied_sites)
    occupation(site) = mod1(site, sites_per_cell(config)) in occupied ? 1.0 : 0.0
    hoppings, interactions = build_infinite_model_terms(config, parameters)
    onsite = sum((
        real(term.amp) * occupation(term.source) for term in hoppings if
        term.target == term.source
    ); init=0.0)
    interaction = sum((
        term.V * occupation(term.i) * occupation(term.j) for
        term in interactions
    ); init=0.0)
    return onsite + interaction
end

function expected_product_raw_charge(config, occupied_sites, bond)
    occupied_left = count(<=(bond), occupied_sites)
    charges = mpskit_centered_charges(config)
    return occupied_left * charges.occupied +
        (bond - occupied_left) * charges.empty
end

function analytic_channel_state()
    a = 0.8
    b = 0.6
    first_kraus = ComplexF64[
        sqrt(a) 0
        0 sqrt(b)
    ]
    second_kraus = ComplexF64[
        0 sqrt(1 - b)
        sqrt(1 - a) 0
    ]
    virtual = TensorKit.ComplexSpace(2)
    physical = TensorKit.ComplexSpace(2)
    tensor = TensorKit.TensorMap(
        vcat(first_kraus, second_kraus),
        TensorKit.:(⊗)(virtual, physical),
        virtual,
    )
    state = MPSKit.InfiniteMPS([tensor, tensor])
    one_site_transfer = kron(conj(first_kraus), first_kraus) +
        kron(conj(second_kraus), second_kraus)
    exact_cell_eigenvalues = eigvals(one_site_transfer * one_site_transfer)
    sort!(exact_cell_eigenvalues; by=abs, rev=true)
    return state, exact_cell_eigenvalues
end

@testset "MPSKit paper product observables" begin
    config, parameters, candidate, state, hamiltonian =
        paper_product_observable_fixture()
    environments = MPSKit.environments(state, hamiltonian, state)

    energy = mpskit_energy_data(
        state,
        hamiltonian,
        config;
        environments,
    )
    exact_energy = direct_product_energy_per_cell(
        config,
        parameters,
        candidate.occupied_sites,
    )
    @test energy isa MPSKitEnergyData
    @test energy.per_cell ≈ exact_energy atol=1e-11 rtol=0
    @test energy.per_x ≈ exact_energy / config.x_period atol=1e-11 rtol=0
    @test energy.per_unit_cell ≈
        exact_energy / unit_cells_per_cell(config) atol=1e-11 rtol=0
    @test energy.per_site ≈
        exact_energy / sites_per_cell(config) atol=1e-11 rtol=0
    @test abs(energy.imaginary_per_cell) <= 1e-12

    density = mpskit_density_data(state, config)
    expected_densities = [
        site in candidate.occupied_sites ? 1.0 : 0.0 for
        site in 1:sites_per_cell(config)
    ]
    @test density isa MPSKitDensityData
    @test density.densities ≈ expected_densities atol=1e-12 rtol=0
    @test density.ring_sums ≈ [
        sum(expected_densities[(x * config.Ly + 1):((x + 1) * config.Ly)]) for
        x in 0:(config.x_period - 1)
    ] atol=1e-12 rtol=0
    @test density.total_particles ≈ particles_per_cell(config) atol=1e-12
    @test density.expected_particles == particles_per_cell(config)
    @test abs(density.particle_error) <= 1e-12
    @test length(density.rows) == sites_per_cell(config)
    @test getproperty.(density.rows, :density) == density.densities

    entanglement = mpskit_entanglement_data(state, config)
    @test length(entanglement) == config.x_period == 3
    for cut in entanglement
        expected_bond = cut.cut_x * config.Ly
        expected_raw_charge = expected_product_raw_charge(
            config,
            candidate.occupied_sites,
            expected_bond,
        )
        @test cut isa MPSKitEntanglementData
        @test cut.bond == expected_bond
        @test cut.entropy ≈ 0.0 atol=1e-12
        @test length(cut.levels) == 1
        @test only(cut.levels).probability ≈ 1.0 atol=1e-12
        @test only(cut.levels).raw_charge == expected_raw_charge
        @test only(cut.levels).physical_charge ==
            expected_raw_charge / charge_scale(config)
        @test length(cut.sectors) == 1
        @test only(cut.sectors).raw_charge == expected_raw_charge
        @test only(cut.sectors).weight ≈ 1.0 atol=1e-12
        @test cut.mean_raw_charge ≈ expected_raw_charge atol=1e-12
        @test cut.mean_physical_charge ≈
            expected_raw_charge / charge_scale(config) atol=1e-12
    end
end

@testset "MPSKit neutral transfer spectrum and correlation length" begin
    config = InfiniteCylinderConfig(;
        Ly=2,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    state, exact_eigenvalues = analytic_channel_state()
    transfer = mpskit_neutral_transfer_data(
        state,
        config;
        neigs=2,
        tol=1e-10,
        maxiter=500,
        krylovdim=12,
        seed=0x4d50534b,
    )
    exact_ratio = abs(exact_eigenvalues[2] / exact_eigenvalues[1])
    @test transfer isa MPSKitNeutralTransferData
    @test transfer.valid
    @test transfer.converged >= 2
    @test length(transfer.eigenvalues) == 2
    @test length(transfer.residual_norms) == 2
    @test all(<=(1e-9), transfer.residual_norms)
    @test transfer.ratio ≈ exact_ratio atol=1e-8 rtol=0
    @test transfer.xi_cell ≈ -1 / log(exact_ratio) atol=1e-8 rtol=0
    @test transfer.xi_x ≈
        config.x_period * (-1 / log(exact_ratio)) atol=1e-8 rtol=0

    product_config = InfiniteCylinderConfig(;
        Ly=4,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    product_state = mpskit_product_state(product_config, [1, 3])
    product_transfer = mpskit_neutral_transfer_data(
        product_state,
        product_config;
        neigs=2,
        tol=1e-10,
        maxiter=100,
        krylovdim=8,
        seed=0x4d50534b,
    )
    @test !product_transfer.valid
    @test product_transfer.converged == 1
    @test length(product_transfer.eigenvalues) == 1
    @test occursin("fewer than two", product_transfer.reason)
    @test isnan(product_transfer.xi_cell)
    @test isnan(product_transfer.xi_x)
end
