const BLOCKTENSORKIT_COMMIT =
    "230cec77c9c7527817d2216b9c6de87f6d8bdda8"
const MPSKIT_BACKEND_COMMIT =
    "4b579944066e9b4bf814205a99f9c24014f76885"
const TENSORKITTENSORS_COMMIT =
    "3755705a1c44a3d5e32086e7d89b2c561b268cb1"

function mpskit_backend_provenance()
    return (
        backend="mpskit_idmrg_v1",
        blocktensorkit_commit=BLOCKTENSORKIT_COMMIT,
        blocktensorkit_version=string(Base.pkgversion(BlockTensorKit)),
        mpskit_commit=MPSKIT_BACKEND_COMMIT,
        mpskit_version=string(Base.pkgversion(MPSKit)),
        tensorkit_version=string(Base.pkgversion(TensorKit)),
        tensorkittensors_commit=TENSORKITTENSORS_COMMIT,
        tensorkittensors_version=string(Base.pkgversion(TensorKitTensors)),
        jld2_version=string(Base.pkgversion(JLD2)),
    )
end
