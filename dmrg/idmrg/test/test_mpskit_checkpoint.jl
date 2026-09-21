using Test
using TOML
using InfiniteCylinderDMRG

import JLD2

include(joinpath(@__DIR__, "mpskit_checkpoint_helpers.jl"))

const MPSKIT_CHECKPOINT_FIG2_MANIFEST = normpath(joinpath(
    @__DIR__, "..", "benchmarks", "fqahc_fig2.toml"
))

function captured_mpskit_checkpoint_exception(f)
    try
        f()
        return nothing
    catch error
        return error
    end
end

function run_mpskit_checkpoint_process(
    script::AbstractString,
    output_directory::AbstractString,
)
    project = normpath(joinpath(@__DIR__, ".."))
    stdout_path = joinpath(output_directory, basename(script) * ".stdout.log")
    stderr_path = joinpath(output_directory, basename(script) * ".stderr.log")
    command = `$(Base.julia_cmd()) --startup-file=no --project=$project $script $output_directory`
    return open(stdout_path, "w") do stdout_io
        open(stderr_path, "w") do stderr_io
            run(
                pipeline(
                    ignorestatus(command);
                    stdout=stdout_io,
                    stderr=stderr_io,
                ),
            )
        end
    end
end

function with_mpskit_checkpoint_directory(f::Function)
    artifact_root = get(ENV, "IDMRG_TEST_ARTIFACTS", "")
    if isempty(artifact_root)
        return mktempdir(f)
    end
    directory = joinpath(artifact_root, "mpskit_checkpoint_restart")
    mkpath(directory)
    return f(directory)
end

@testset "versioned atomic MPSKit checkpoint" begin
    mktempdir() do directory
        config, _, state, _ = mpskit_checkpoint_test_fixture()
        path = joinpath(directory, "state.h5")
        metadata = (
            candidate_id="product_fixture",
            generation=2,
            requested_maxdim=4,
        )
        @test save_mpskit_checkpoint(path, state, config, metadata) == path
        @test readdir(directory) == ["state.h5"]

        loaded = load_mpskit_checkpoint(path, config)
        @test loaded isa MPSKitCheckpointData
        @test loaded.metadata == metadata
        @test loaded.provenance.format == "mpskit_idmrg_checkpoint_v1"
        @test loaded.provenance.backend == "mpskit_idmrg_v1"
        @test loaded.provenance.mpskit_commit == MPSKIT_BACKEND_COMMIT
        @test loaded.provenance.configuration_signature ==
            configuration_signature(config)
        @test mpskit_product_occupations(loaded.state, config) == [1, 3]
        @test mpskit_state_space_fingerprints(loaded.state, config) ==
            mpskit_state_space_fingerprints(state, config)

        wrong_config = with_flux(config, 0.25)
        config_error = captured_mpskit_checkpoint_exception(
            () -> load_mpskit_checkpoint(path, wrong_config),
        )
        @test config_error isa CheckpointCompatibilityError
        @test occursin("phi_y", sprint(showerror, config_error))

        wrong_backend = joinpath(directory, "wrong-backend.h5")
        cp(path, wrong_backend)
        JLD2.jldopen(wrong_backend, "r+") do file
            delete!(file, "mpskit_commit")
            file["mpskit_commit"] = "not-the-pinned-commit"
        end
        backend_error = captured_mpskit_checkpoint_exception(
            () -> load_mpskit_checkpoint(wrong_backend, config),
        )
        @test backend_error isa CheckpointCompatibilityError
        @test occursin("MPSKit commit", sprint(showerror, backend_error))

        incomplete = joinpath(directory, "incomplete.h5")
        cp(path, incomplete)
        JLD2.jldopen(incomplete, "r+") do file
            delete!(file, "complete")
            file["complete"] = false
        end
        incomplete_error = captured_mpskit_checkpoint_exception(
            () -> load_mpskit_checkpoint(incomplete, config),
        )
        @test incomplete_error isa CheckpointFormatError
        @test occursin("complete", lowercase(sprint(showerror, incomplete_error)))

        truncated = joinpath(directory, "truncated.h5")
        bytes = read(path)
        open(truncated, "w") do io
            write(io, @view(bytes[1:min(64, length(bytes))]))
        end
        truncated_error = captured_mpskit_checkpoint_exception(
            () -> load_mpskit_checkpoint(truncated, config),
        )
        @test truncated_error isa CheckpointFormatError

        legacy = joinpath(directory, "legacy-itensor.h5")
        _, _, legacy_state = initial_infinite_mps(config)
        save_checkpoint(legacy, legacy_state, config)
        legacy_error = @test_logs (
            :warn,
            r"File likely not written by JLD2",
        ) captured_mpskit_checkpoint_exception(
            () -> load_mpskit_checkpoint(legacy, config),
        )
        @test legacy_error isa CheckpointFormatError
        @test occursin("MPSKit", sprint(showerror, legacy_error))
    end
end

@testset "two-process MPSKit checkpoint restart" begin
    with_mpskit_checkpoint_directory() do directory
        save_script = joinpath(@__DIR__, "mpskit_checkpoint_save.jl")
        resume_script = joinpath(@__DIR__, "mpskit_checkpoint_resume.jl")

        save_process = run_mpskit_checkpoint_process(save_script, directory)
        @test save_process.exitcode == 0
        if save_process.exitcode != 0
            println(
                stderr,
                read(
                    joinpath(directory, basename(save_script) * ".stderr.log"),
                    String,
                ),
            )
            return
        end

        resume_process = run_mpskit_checkpoint_process(resume_script, directory)
        @test resume_process.exitcode == 0
        if resume_process.exitcode != 0
            println(
                stderr,
                read(
                    joinpath(directory, basename(resume_script) * ".stderr.log"),
                    String,
                ),
            )
            return
        end

        before = TOML.parsefile(joinpath(directory, "before.toml"))
        after = TOML.parsefile(joinpath(directory, "after.toml"))
        @test before["process_id"] != after["process_id"]
        @test after["load_completed"] === true
        @test after["no_expansion_iteration_completed"] === true
        @test after["controlled_expansion_completed"] === true
        @test after["next_flux_iteration_completed"] === true
        @test after["configuration_signature"] ==
            before["configuration_signature"]
        @test after["physical_space_fingerprint"] ==
            before["physical_space_fingerprint"]
        @test after["virtual_space_fingerprint"] ==
            before["virtual_space_fingerprint"]
        @test after["pre_energy_per_site"] ≈
            before["energy_per_site"] atol=1e-10 rtol=0
        @test after["pre_galerkin_residual"] ≈
            before["galerkin_residual"] atol=1e-10 rtol=0
        @test after["pre_densities"] ≈ before["densities"] atol=1e-10 rtol=0
        @test after["pre_sector_weights"] == before["sector_weights"]
        @test after["pre_entanglement_spectrum"] ==
            before["entanglement_spectrum"]
        @test after["original_checkpoint_sha256"] ==
            before["checkpoint_sha256"]
        @test after["refinement_iterations"] == 1
        @test before["maxlinkdim"] == 4
        @test after["pre_maxlinkdim"] == 4
        @test after["post_maxlinkdim"] == 8
        @test after["next_flux_phi_y"] ≈ 0.1 atol=0 rtol=0
        @test isfile(joinpath(directory, "state-resumed.h5"))
        @test isfile(joinpath(directory, "state-next-flux.h5"))
        @test !isempty(after["resumed_checkpoint_sha256"])
        @test !isempty(after["next_flux_checkpoint_sha256"])

        spec = load_fig2_benchmark(MPSKIT_CHECKPOINT_FIG2_MANIFEST)
        gate = write_fig2_restart_gate!(
            spec,
            directory;
            before_path=joinpath(directory, "before.toml"),
            after_path=joinpath(directory, "after.toml"),
            save_process_exit_code=save_process.exitcode,
            resume_process_exit_code=resume_process.exitcode,
        )
        @test gate["pass"] === true
        @test gate["pre_maxlinkdim"] == 4
        @test gate["post_maxlinkdim"] == 8
        @test gate["energy_absolute_error"] <= 1e-10
        @test gate["density_max_absolute_error"] <= 1e-10
        @test TOML.parsefile(joinpath(directory, "restart_gate.toml")) == gate
    end
end
