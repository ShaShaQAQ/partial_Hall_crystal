using Test
using InfiniteCylinderDMRG

import JLD2
import BlockTensorKit
import MPSKit
import TensorKit
import TensorKitTensors

@testset "pinned MPSKit backend" begin
    @test BLOCKTENSORKIT_COMMIT ==
        "230cec77c9c7527817d2216b9c6de87f6d8bdda8"
    @test MPSKIT_BACKEND_COMMIT ==
        "811ecf6c06c1f7c1bc656da61abcd679effcd428"
    @test TENSORKITTENSORS_COMMIT ==
        "3755705a1c44a3d5e32086e7d89b2c561b268cb1"

    provenance = mpskit_backend_provenance()
    @test provenance.backend == "mpskit_idmrg_v1"
    @test provenance.blocktensorkit_version ==
        string(Base.pkgversion(BlockTensorKit))
    @test provenance.mpskit_version == string(Base.pkgversion(MPSKit))
    @test provenance.tensorkit_version == string(Base.pkgversion(TensorKit))
    @test provenance.tensorkittensors_version ==
        string(Base.pkgversion(TensorKitTensors))
    @test provenance.jld2_version == string(Base.pkgversion(JLD2))
end

@testset "uniform centered paper charge" begin
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=6,
        x_period=3,
        filling_num=7,
        filling_den=9,
    )

    @test mpskit_centered_charges(config) ==
        (scale=18, empty=-7, occupied=11)
    @test mpskit_total_raw_charge(default_occupied_sites(config), config) == 0

    physical_spaces = mpskit_physical_spaces(config)
    @test length(physical_spaces) == 36
    @test all(==(first(physical_spaces)), physical_spaces)
    @test TensorKit.dim(first(physical_spaces)) == 2
    @test mpskit_physical_raw_charges(first(physical_spaces)) == [-7, 11]

    for candidate in fig2_initial_candidates(config)
        psi = mpskit_product_state(config, candidate.occupied_sites)
        @test length(psi) == sites_per_cell(config)
        @test mpskit_product_occupations(psi, config) ==
            candidate.occupied_sites
        @test mpskit_total_raw_charge(candidate.occupied_sites, config) == 0
        @test all(
            TensorKit.dim(MPSKit.right_virtualspace(psi, site)) == 1 for
            site in 1:length(psi)
        )
        @test all(
            MPSKit.physicalspace(psi, site) == first(physical_spaces) for
            site in 1:length(psi)
        )
    end

    @test_throws ArgumentError mpskit_product_state(config, collect(1:13))
    @test_throws ArgumentError mpskit_product_state(
        config,
        vcat(collect(1:13), 13),
    )
    @test_throws ArgumentError mpskit_product_state(
        config,
        vcat(collect(1:13), 37),
    )
end
