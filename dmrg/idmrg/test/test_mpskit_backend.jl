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
