using Test
using InfiniteCylinderDMRG

@testset "deterministic Fig. 2 product-state candidates" begin
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=6,
        x_period=3,
        filling_num=7,
        filling_den=9,
    )
    candidates = fig2_initial_candidates(config)
    @test length(candidates) == 10
    @test length(unique(candidate.id for candidate in candidates)) == 10
    @test length(unique(candidate.occupied_sites for candidate in candidates)) == 10
    @test all(length(candidate.occupied_sites) == 14 for candidate in candidates)
    @test all(
        length(unique(candidate.occupied_sites)) == 14 for candidate in candidates
    )
    @test all(
        issorted(candidate.occupied_sites) &&
            all(in(1:sites_per_cell(config)), candidate.occupied_sites) for
        candidate in candidates
    )

    cdw_candidates = filter(candidate -> candidate.kind == :cdw_doped, candidates)
    @test length(cdw_candidates) == 9
    @test Set(candidate.cdw_translation for candidate in cdw_candidates) == Set(0:2)
    for candidate in cdw_candidates
        background = filter(
            site -> paper_cdw_color(config, site) == candidate.cdw_translation,
            candidate.occupied_sites,
        )
        @test length(background) == 12
        @test length(setdiff(candidate.occupied_sites, background)) == 2
    end

    controls = filter(candidate -> candidate.kind == :uniform_control, candidates)
    @test length(controls) == 1
    @test only(controls).occupied_sites == default_occupied_sites(config)
    @test fig2_initial_candidates(config) == candidates
end

@testset "Fig. 2 candidates reject a different cell" begin
    @test_throws ArgumentError fig2_initial_candidates(
        InfiniteCylinderConfig(;
            geometry=:paper_straight,
            Ny=3,
            x_period=3,
            filling_num=7,
            filling_den=9,
        ),
    )
    @test_throws ArgumentError fig2_initial_candidates(
        InfiniteCylinderConfig(; Ly=12, x_period=3, filling_num=7, filling_den=9),
    )
end
