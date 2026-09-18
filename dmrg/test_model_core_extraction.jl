using Test

include("DMRGFullModel.jl")
using .DMRGFullModel

@testset "shared cylinder model core" begin
    @test isdefined(DMRGFullModel, :CylinderModelCore)
    if isdefined(DMRGFullModel, :CylinderModelCore)
        @test DMRGFullModel.CylinderModelParams ===
            DMRGFullModel.CylinderModelCore.CylinderModelParams
        tR = DMRGFullModel.CylinderModelCore.hopping_matrices_from_get_Hk(1.0, 0.2)
        @test length(tR) == 11
        @test maximum(abs(first(k)) for k in keys(tR)) == 2
    end
end
