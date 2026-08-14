using Test
using InfiniteCylinderDMRG

@testset "paper-straight configuration arithmetic" begin
    c = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=6,
        x_period=3,
        filling_num=7,
        filling_den=9,
    )
    @test c.geometry == :paper_straight
    @test c.Ny == 6
    @test c.Ly == 12
    @test circumference_unit_cells(c) == 6
    @test physical_circumference_sites(c) == 12
    @test sites_per_cell(c) == 36
    @test unit_cells_per_cell(c) == 18
    @test particles_per_cell(c) == 14
    @test physical_site_density(c) == 7 // 18
    @test charge_scale(c) == 18

    @test paper_site(c, 0, 0, :A) == 1
    @test paper_site(c, 0, 0, :B) == 2
    @test paper_site(c, 0, 6, :A) == 1
    @test paper_site(c, 1, 0, :A) == 13
    @test paper_coordinates(c, 14) == (1, 0, :B)
    @test paper_coordinates(c, 37) == (3, 0, :A)

    next_flux = with_flux(c, 0.25)
    @test next_flux.geometry == c.geometry
    @test next_flux.Ny == c.Ny
    @test next_flux.phi_y == 0.25
    @test configuration_signature(next_flux) == configuration_signature(c)
end

@testset "paper and legacy geometry are unambiguous" begin
    legacy = InfiniteCylinderConfig(; Ly=6, x_period=1)
    @test legacy.geometry == :legacy_sheared
    @test legacy.Ny == 3
    @test legacy.Ly == 6
    @test circumference_unit_cells(legacy) == 3
    @test physical_circumference_sites(legacy) == 6

    paper = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=6,
        x_period=3,
        filling_num=7,
        filling_den=9,
    )
    legacy_same_sites = InfiniteCylinderConfig(;
        Ly=12,
        x_period=3,
        filling_num=7,
        filling_den=9,
    )
    @test configuration_signature(paper) != configuration_signature(legacy_same_sites)
    @test_throws ArgumentError InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ly=12,
        Ny=6,
        x_period=3,
        filling_num=7,
        filling_den=9,
    )
    @test_throws ArgumentError InfiniteCylinderConfig(;
        geometry=:paper_straight,
        x_period=3,
        filling_num=7,
        filling_den=9,
    )
    @test_throws ArgumentError InfiniteCylinderConfig(; geometry=:legacy_sheared, Ny=3)
    @test_throws ArgumentError InfiniteCylinderConfig(; geometry=:unknown, Ly=6)
    @test_throws ArgumentError paper_site(legacy, 0, 0, :A)
    @test_throws ArgumentError paper_coordinates(legacy, 1)
end

@testset "paper CLI requires geometry and Ny" begin
    settings = parse_single_point_args([
        "--geometry=paper_straight",
        "--Ny=6",
        "--x_period=3",
        "--filling_num=7",
        "--filling_den=9",
        "--phi_y=0.0",
        "--maxdim=8",
        "--output=paper-test",
    ])
    @test settings.config.geometry == :paper_straight
    @test sites_per_cell(settings.config) == 36
    @test particles_per_cell(settings.config) == 14

    @test_throws ArgumentError parse_single_point_args([
        "--Ny=6",
        "--x_period=3",
        "--filling_num=7",
        "--filling_den=9",
        "--phi_y=0.0",
        "--maxdim=8",
        "--output=paper-test",
    ])
    @test_throws ArgumentError parse_single_point_args([
        "--geometry=paper_straight",
        "--Ly=12",
        "--Ny=6",
        "--x_period=3",
        "--filling_num=7",
        "--filling_den=9",
        "--phi_y=0.0",
        "--maxdim=8",
        "--output=paper-test",
    ])
end
