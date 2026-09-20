using Test
using TOML
using InfiniteCylinderDMRG

const MPSKIT_FIG2_MANIFEST_PATH = normpath(joinpath(
    @__DIR__, "..", "benchmarks", "fqahc_fig2.toml"
))
const EXPECTED_MPSKIT_FIG2_BACKEND = Dict{String,Any}(
    "id" => "mpskit_idmrg_v1",
    "mpskit_commit" =>
        "811ecf6c06c1f7c1bc656da61abcd679effcd428",
    "tensorkittensors_commit" =>
        "3755705a1c44a3d5e32086e7d89b2c561b268cb1",
    "legacy_backend_role" => "diagnostic_only",
)

function write_mpskit_fig2_manifest_fixture!(directory, data)
    path = joinpath(directory, "manifest.toml")
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function capture_mpskit_fig2_manifest_error(data)
    return mktempdir() do directory
        path = write_mpskit_fig2_manifest_fixture!(directory, data)
        try
            load_fig2_benchmark(path)
            nothing
        catch error
            error
        end
    end
end

@testset "MPSKit Fig. 2 workflow adapter" begin
@testset "production backend manifest contract" begin
    manifest = TOML.parsefile(MPSKIT_FIG2_MANIFEST_PATH)
    @test get(manifest, "format", "") == "fqahc_fig2_benchmark_v5"
    @test get(manifest, "backend", nothing) == EXPECTED_MPSKIT_FIG2_BACKEND

    spec = try
        load_fig2_benchmark(MPSKIT_FIG2_MANIFEST_PATH)
    catch
        nothing
    end
    @test !isnothing(spec)
    if !isnothing(spec)
        @test get(spec.data, "backend", nothing) ==
            EXPECTED_MPSKIT_FIG2_BACKEND
    end
end

@testset "backend selection fails closed" begin
    baseline = TOML.parsefile(MPSKIT_FIG2_MANIFEST_PATH)
    baseline["format"] = "fqahc_fig2_benchmark_v5"
    baseline["backend"] = deepcopy(EXPECTED_MPSKIT_FIG2_BACKEND)

    cases = [
        (
            "missing backend",
            data -> delete!(data, "backend"),
            "missing backend",
        ),
        (
            "unknown backend",
            data -> data["backend"]["id"] = "unknown_idmrg_backend",
            "backend id",
        ),
        (
            "MPSKit pin drift",
            data -> data["backend"]["mpskit_commit"] = repeat("0", 40),
            "mpskit commit",
        ),
        (
            "TensorKitTensors pin drift",
            data -> data["backend"]["tensorkittensors_commit"] = repeat("0", 40),
            "tensorkittensors commit",
        ),
        (
            "legacy role drift",
            data -> data["backend"]["legacy_backend_role"] = "production",
            "legacy backend role",
        ),
        (
            "legacy backend selected for production",
            data -> data["backend"]["id"] = "itensor_infinite_mps_v1",
            "production manifest",
        ),
    ]

    for (label, tamper, reason_fragment) in cases
        @testset "$label" begin
            candidate = deepcopy(baseline)
            tamper(candidate)
            error = capture_mpskit_fig2_manifest_error(candidate)
            @test error isa ArgumentError
            if error isa ArgumentError
                @test occursin(
                    lowercase(reason_fragment),
                    lowercase(sprint(showerror, error)),
                )
            end
        end
    end
end
end
