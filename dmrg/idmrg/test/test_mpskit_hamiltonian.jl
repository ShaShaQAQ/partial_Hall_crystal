using Test
using LinearAlgebra
using InfiniteCylinderDMRG

function exact_hop!(
    matrix::AbstractMatrix,
    target::Int,
    source::Int,
    amplitude::Number,
    nsites::Int;
    jordan_wigner::Bool=true,
)
    target == source && return exact_number!(
        matrix,
        target,
        amplitude,
        nsites,
    )
    for ket in 0:(2^nsites - 1)
        ((ket >> (source - 1)) & 1) == 1 || continue
        ((ket >> (target - 1)) & 1) == 0 || continue
        between = (min(target, source) + 1):(max(target, source) - 1)
        parity = count(
            site -> ((ket >> (site - 1)) & 1) == 1,
            between,
        )
        sign = jordan_wigner && isodd(parity) ? -1 : 1
        bra = xor(
            xor(ket, 1 << (source - 1)),
            1 << (target - 1),
        )
        matrix[bra + 1, ket + 1] += amplitude * sign
    end
    return matrix
end

function exact_number!(
    matrix::AbstractMatrix,
    site::Int,
    amplitude::Number,
    nsites::Int,
)
    for ket in 0:(2^nsites - 1)
        ((ket >> (site - 1)) & 1) == 1 || continue
        matrix[ket + 1, ket + 1] += amplitude
    end
    return matrix
end

function exact_interaction!(
    matrix::AbstractMatrix,
    first_site::Int,
    second_site::Int,
    interaction::Number,
    nsites::Int,
)
    for ket in 0:(2^nsites - 1)
        ((ket >> (first_site - 1)) & 1) == 1 || continue
        ((ket >> (second_site - 1)) & 1) == 1 || continue
        matrix[ket + 1, ket + 1] += interaction
    end
    return matrix
end

function exact_fock_hamiltonian(
    nsites::Int,
    hoppings::AbstractVector{<:InfiniteHoppingTerm},
    interactions::AbstractVector{<:InfiniteInteractionTerm};
    jordan_wigner::Bool=true,
)
    matrix = zeros(ComplexF64, 2^nsites, 2^nsites)
    for hopping in hoppings
        exact_hop!(
            matrix,
            hopping.target,
            hopping.source,
            hopping.amp,
            nsites;
            jordan_wigner,
        )
    end
    for interaction in interactions
        exact_interaction!(
            matrix,
            interaction.i,
            interaction.j,
            interaction.V,
            nsites,
        )
    end
    return matrix
end

function exact_test_terms(nsites::Int)
    amplitude = 0.37 + 0.19im
    hoppings = InfiniteHoppingTerm[
        InfiniteHoppingTerm(1, nsites, amplitude),
        InfiniteHoppingTerm(nsites, 1, conj(amplitude)),
        InfiniteHoppingTerm(min(2, nsites), min(2, nsites), -0.41),
    ]
    if nsites >= 3
        push!(hoppings, InfiniteHoppingTerm(2, 3, -0.73))
        push!(hoppings, InfiniteHoppingTerm(3, 2, -0.73))
    end
    interactions = InfiniteInteractionTerm[
        InfiniteInteractionTerm(1, nsites, 0.61, 1),
    ]
    return hoppings, interactions
end

@testset "independent exact fermionic Fock oracle" begin
    nsites = 4
    hoppings, interactions = exact_test_terms(nsites)
    exact = exact_fock_hamiltonian(nsites, hoppings, interactions)
    no_jordan_wigner = exact_fock_hamiltonian(
        nsites,
        hoppings,
        interactions;
        jordan_wigner=false,
    )
    @test exact ≈ exact' atol=1e-12 rtol=0
    @test norm(exact - no_jordan_wigner) > 1e-6
    @test exact[0b1011 + 1, 0b1011 + 1] ≈ 0.61 - 0.41
end

@testset "finite MPSKit Jordan-Wigner Hamiltonian" begin
    for nsites in 2:6
        hoppings, interactions = exact_test_terms(nsites)
        expected = exact_fock_hamiltonian(nsites, hoppings, interactions)
        hamiltonian = mpskit_finite_hamiltonian(
            nsites,
            hoppings,
            interactions;
            charge_mode=:uncentered,
        )
        actual = mpskit_dense_matrix(hamiltonian)
        @test actual ≈ expected atol=1e-12 rtol=0
        @test actual ≈ actual' atol=1e-12 rtol=0
    end
end

@testset "centered and uncentered finite Hamiltonians agree" begin
    config = InfiniteCylinderConfig(;
        Ly=4,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    hoppings, interactions = exact_test_terms(sites_per_cell(config))
    expected = exact_fock_hamiltonian(
        sites_per_cell(config),
        hoppings,
        interactions,
    )
    centered = mpskit_finite_hamiltonian(
        config,
        hoppings,
        interactions;
        charge_mode=:centered,
    )
    uncentered = mpskit_finite_hamiltonian(
        sites_per_cell(config),
        hoppings,
        interactions;
        charge_mode=:uncentered,
    )
    @test mpskit_dense_matrix(centered) ≈ expected atol=1e-12 rtol=0
    @test mpskit_dense_matrix(uncentered) ≈ expected atol=1e-12 rtol=0
end

@testset "paper term lists map directly to MPSKit" begin
    parameters = CylinderModelParams(;
        t1=1.0,
        t3=0.2,
        V1=1.0,
        V2=0.0,
        V3=0.0,
    )
    for x_period in (1, 3)
        config = InfiniteCylinderConfig(;
            geometry=:paper_straight,
            Ny=2,
            x_period,
            filling_num=1,
            filling_den=1,
            phi_y=0.37π,
        )
        hoppings, interactions = build_infinite_model_terms(config, parameters)
        finite_hoppings, finite_interactions = tile_to_finite_window(
            config,
            hoppings,
            interactions;
            Lx=x_period,
        )
        expected = exact_fock_hamiltonian(
            sites_per_cell(config),
            finite_hoppings,
            finite_interactions,
        )
        hamiltonian = mpskit_finite_hamiltonian(
            config,
            finite_hoppings,
            finite_interactions;
            charge_mode=:centered,
        )
        @test mpskit_dense_matrix(hamiltonian) ≈ expected atol=1e-12 rtol=0
    end
end

@testset "paper infinite MPSKit Hamiltonian provenance" begin
    parameters = CylinderModelParams(;
        t1=1.0,
        t3=0.2,
        V1=1.0,
        V2=0.0,
        V3=0.0,
    )
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=6,
        x_period=3,
        filling_num=7,
        filling_den=9,
        phi_y=0.413,
    )
    hoppings, interactions = build_infinite_model_terms(config, parameters)
    hamiltonian = mpskit_infinite_hamiltonian(config, parameters)

    @test length(hamiltonian) == sites_per_cell(config) == 36
    @test mpskit_local_operators_are_hermitian(
        config,
        parameters;
        atol=1e-12,
    )
    @test mpskit_term_fingerprint(config, parameters) ==
        canonical_term_fingerprint(hoppings, interactions)
    @test occursin(r"^[0-9a-f]{64}$", mpskit_term_fingerprint(config, parameters))
    @test mpskit_terms_are_approx(
        with_flux(config, config.phi_y + 2π),
        config,
        parameters;
        atol=1e-12,
    )
    virtual_dimensions = mpskit_mpo_virtual_dimensions(hamiltonian)
    @test length(virtual_dimensions) == 36
    @test all(>(1), virtual_dimensions)
end
