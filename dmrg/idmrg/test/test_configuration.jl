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
