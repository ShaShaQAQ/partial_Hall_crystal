using Test
using InfiniteCylinderDMRG

@testset "infinite cylinder configuration" begin
    c1 = InfiniteCylinderConfig(; Ly=6, x_period=1, filling_num=1, filling_den=3)
    @test sites_per_cell(c1) == 6
    @test unit_cells_per_cell(c1) == 3
    @test particles_per_cell(c1) == 1
    @test charge_scale(c1) == 6

    c3 = InfiniteCylinderConfig(; Ly=6, x_period=3, filling_num=1, filling_den=3)
    @test sites_per_cell(c3) == 18
    @test particles_per_cell(c3) == 3
    @test charge_scale(c3) == 6

    @test_throws ArgumentError InfiniteCylinderConfig(; Ly=2, x_period=1)
    csmall = InfiniteCylinderConfig(; Ly=2, x_period=3)
    @test particles_per_cell(csmall) == 1

    @test linear_site(c1, 0, 0) == 1
    @test linear_site(c1, 1, 0) == 7
    @test site_coordinates(c1, 7) == (1, 0)
    @test canonical_site(c1, 7) == (1, 1)
end

@testset "constructor invariants" begin
    @test_throws ArgumentError InfiniteCylinderConfig(3, 1, 1, 3, 0.0)
    @test_throws ArgumentError InfiniteCylinderConfig(6, 0, 1, 3, 0.0)
    @test_throws ArgumentError InfiniteCylinderConfig(6, 1, 0, 3, 0.0)
    @test_throws ArgumentError InfiniteCylinderConfig(6, 1, 1, 0, 0.0)
    @test_throws ArgumentError InfiniteCylinderConfig(2, 1, 1, 3, 0.0)

    for phi_y in (NaN, Inf, -Inf)
        @test_throws ArgumentError InfiniteCylinderConfig(6, 1, 1, 3, phi_y)
        @test_throws ArgumentError InfiniteCylinderConfig(; phi_y)
    end

    @test_throws ArgumentError InfiniteCylinderConfig(6, 1, 7, 3, 0.0)
    @test_throws ArgumentError InfiniteCylinderConfig(; Ly=6, filling_num=7, filling_den=3)

    positional = InfiniteCylinderConfig(6, 1, 1, 3, 0.25)
    @test positional.phi_y == 0.25
    @test physical_site_density(positional) == 1 // 6
end

@testset "overflow-safe geometry" begin
    @test_throws ArgumentError InfiniteCylinderConfig(
        typemax(Int) - 1,
        2,
        2,
        1,
        0.0,
    )

    large_ratio_term = typemax(Int) ÷ 2 + 1
    large_ratio = InfiniteCylinderConfig(2, 2, large_ratio_term, large_ratio_term, 0.0)
    @test particles_per_cell(large_ratio) == 2
    @test physical_site_density(large_ratio) == 1 // 2

    c = InfiniteCylinderConfig()
    @test_throws OverflowError linear_site(c, typemax(Int), 0)
    @test_throws OverflowError site_coordinates(c, typemin(Int))
    @test_throws OverflowError canonical_site(c, typemin(Int))
end

@testset "exact occupied-site patterns" begin
    @test default_occupied_sites(InfiniteCylinderConfig(; Ly=6, x_period=1)) == [1]
    @test default_occupied_sites(InfiniteCylinderConfig(; Ly=2, x_period=3)) == [1]

    full = InfiniteCylinderConfig(; Ly=6, filling_num=2, filling_den=1)
    near_full = InfiniteCylinderConfig(; Ly=6, filling_num=5, filling_den=3)

    for c in (full, near_full)
        occupied = default_occupied_sites(c)
        @test length(occupied) == particles_per_cell(c)
        @test length(unique(occupied)) == particles_per_cell(c)
        @test all(in(1:sites_per_cell(c)), occupied)
        @test occupied == default_occupied_sites(c)
    end
end
