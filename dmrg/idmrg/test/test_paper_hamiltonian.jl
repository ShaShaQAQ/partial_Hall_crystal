using Test
using LinearAlgebra
using InfiniteCylinderDMRG

function analytic_paper_Hk(k, t1, t3)
    sigma_x = ComplexF64[0 1; 1 0]
    sigma_y = ComplexF64[0 -im; im 0]
    sigma_z = ComplexF64[1 0; 0 -1]
    identity_2 = Matrix{ComplexF64}(I, 2, 2)
    directions = (
        (0.5, -sqrt(3) / 2),
        (0.5, sqrt(3) / 2),
        (-1.0, 0.0),
    )
    dot2(a, b) = a[1] * b[1] + a[2] * b[2]
    g = [2t1 * cos(dot2(k, direction)) for direction in directions]
    g0 = 2t3 * sum(cos(2dot2(k, direction)) for direction in directions)
    return g[1] * sigma_x + g[2] * sigma_y + g[3] * sigma_z + g0 * identity_2
end

@testset "paper hopping orbits reconstruct the two-band model" begin
    params = CylinderModelParams(; t1=1.0, t3=0.2, V1=10.0, V2=2.0, V3=2.0)
    orbits = paper_hopping_orbits(params)
    @test !isempty(orbits)
    @test paper_orbits_are_hermitian(orbits; atol=1e-12)
    for k in ((0.173, 0.291), (1.137, -0.413), (-2.021, 0.777))
        @test paper_bloch_from_orbits(k, params) ≈
            analytic_paper_Hk(k, params.t1, params.t3) atol = 1e-11 rtol = 0
    end
end

@testset "paper straight-cylinder Hamiltonian" begin
    params = CylinderModelParams(; t1=1.0, t3=0.2, V1=10.0, V2=2.0, V3=2.0)
    c0 = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=6,
        x_period=3,
        filling_num=7,
        filling_den=9,
        phi_y=0.0,
    )
    hops, interactions = build_infinite_model_terms(c0, params)
    @test !isempty(hops)
    @test hamiltonian_is_hermitian(c0, hops)
    @test neighbor_shell_counts(interactions) == Dict(1 => 108, 2 => 108, 3 => 108)
    expected_distance_squared = Dict(1 => 1.0, 2 => 3.0, 3 => 4.0)
    @test all(
        isapprox(
            paper_periodic_distance_squared(c0, interaction.i, interaction.j),
            expected_distance_squared[interaction.shell];
            atol=1e-12,
            rtol=0,
        ) for interaction in interactions
    )
    @test straight_seam_is_connected(c0, hops)
    @test paper_window_parity(c0, params; x_cells=7, atol=1e-11)

    twisted = with_flux(c0, 0.37pi)
    twisted_hops, _ = build_infinite_model_terms(twisted, params)
    @test hamiltonian_is_hermitian(twisted, twisted_hops)

    c2pi = with_flux(c0, 2pi)
    @test canonical_term_dict(build_infinite_model_terms(c0, params)...) ≈
        canonical_term_dict(build_infinite_model_terms(c2pi, params)...) atol = 1e-11 rtol = 0
end

@testset "legacy and paper Hamiltonian dispatch remain distinct" begin
    params = CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.4, V3=0.3)
    legacy = InfiniteCylinderConfig(; Ly=12, x_period=1, filling_num=1, filling_den=3)
    paper = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=6,
        x_period=1,
        filling_num=1,
        filling_den=3,
    )
    legacy_terms = canonical_term_dict(build_infinite_model_terms(legacy, params)...)
    paper_terms = canonical_term_dict(build_infinite_model_terms(paper, params)...)
    @test !isapprox(legacy_terms, paper_terms; atol=1e-12, rtol=0)
end
