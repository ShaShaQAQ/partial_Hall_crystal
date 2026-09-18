using Test
using ITensors
using ITensorMPS
using ITensorInfiniteMPS
using InfiniteCylinderDMRG

@testset "infinite Hamiltonian terms" begin
    params = CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.4, V3=0.3)
    for x_period in (1, 3), phi in (0.0, 0.37π, 2π)
        cfg = InfiniteCylinderConfig(; Ly=6, x_period, phi_y=phi)
        hops, ints = build_infinite_model_terms(cfg, params)
        @test !isempty(hops)
        @test !isempty(ints)
        @test all(min(h.target, h.source) in 1:sites_per_cell(cfg) for h in hops)
        @test hamiltonian_is_hermitian(cfg, hops)
    end

    c0 = InfiniteCylinderConfig(; Ly=6, x_period=1, phi_y=0.0)
    c2 = InfiniteCylinderConfig(; Ly=6, x_period=1, phi_y=2π)
    @test canonical_term_dict(build_infinite_model_terms(c0, params)...) ≈
        canonical_term_dict(build_infinite_model_terms(c2, params)...)

    for x_period in (1, 3), phi in (0.0, 0.37π, 2π)
        cfg = InfiniteCylinderConfig(; Ly=6, x_period, phi_y=phi)
        @test finite_term_parity(cfg, params; Lx=7, atol=1e-10)
    end
end

@testset "infinite Hamiltonian backend adapter" begin
    params = CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.4, V3=0.3)
    cfg = InfiniteCylinderConfig(; Ly=6, x_period=1)
    os = infinite_opsum(cfg, params)
    starts = Set(minimum(ITensors.sites(term)) for term in ITensors.terms(os))
    @test starts == Set(1:sites_per_cell(cfg))

    occupied = Set(default_occupied_sites(cfg))
    initstate(n) = mod1(n, sites_per_cell(cfg)) in occupied ? "Occ" : "Emp"
    sites = infsiteinds(
        "Fermion",
        sites_per_cell(cfg);
        initstate,
        conserve_qns=true,
    )
    H = build_infinite_mpo(cfg, params, sites)
    @test H isa InfiniteSum{MPO}
    @test nsites(H) == sites_per_cell(cfg)

    zero_params = CylinderModelParams(; t1=0.0, t3=0.0, V1=0.0, V2=0.0, V3=0.0)
    zero_os = infinite_opsum(cfg, zero_params)
    zero_starts =
        Set(minimum(ITensors.sites(term)) for term in ITensors.terms(zero_os))
    @test zero_starts == Set(1:sites_per_cell(cfg))
    zero_H = build_infinite_mpo(cfg, zero_params, sites)
    @test zero_H isa InfiniteSum{MPO}
    @test nsites(zero_H) == sites_per_cell(cfg)
end
