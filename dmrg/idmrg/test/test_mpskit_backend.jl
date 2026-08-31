using Test
using InfiniteCylinderDMRG

import JLD2
import MPSKit
import TensorKit
import TensorKitTensors

@testset "pinned MPSKit backend" begin
    @test MPSKIT_BACKEND_COMMIT ==
        "4b579944066e9b4bf814205a99f9c24014f76885"
    @test TENSORKITTENSORS_COMMIT ==
        "3755705a1c44a3d5e32086e7d89b2c561b268cb1"

    provenance = mpskit_backend_provenance()
    @test provenance.backend == "mpskit_idmrg_v1"
    @test provenance.mpskit_version == string(Base.pkgversion(MPSKit))
    @test provenance.tensorkit_version == string(Base.pkgversion(TensorKit))
    @test provenance.tensorkittensors_version ==
        string(Base.pkgversion(TensorKitTensors))
    @test provenance.jld2_version == string(Base.pkgversion(JLD2))
end
