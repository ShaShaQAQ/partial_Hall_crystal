using Test
using TOML

function run_restart_process(script::AbstractString, output_directory::AbstractString)
    project = normpath(joinpath(@__DIR__, ".."))
    stdout_path = joinpath(output_directory, basename(script) * ".stdout.log")
    stderr_path = joinpath(output_directory, basename(script) * ".stderr.log")
    command = `$(Base.julia_cmd()) --startup-file=no --project=$project $script $output_directory`
    process = open(stdout_path, "w") do stdout_io
        open(stderr_path, "w") do stderr_io
            run(pipeline(ignorestatus(command); stdout=stdout_io, stderr=stderr_io))
        end
    end
    return process
end

function with_restart_directory(f::Function)
    artifact_root = get(ENV, "IDMRG_TEST_ARTIFACTS", "")
    if isempty(artifact_root)
        return mktempdir(f)
    end
    directory = joinpath(artifact_root, "checkpoint_restart")
    mkpath(directory)
    return f(directory)
end

@testset "cross-process checkpoint restart, expansion, and next flux" begin
    with_restart_directory() do directory
        save_script = joinpath(@__DIR__, "checkpoint_restart_save.jl")
        resume_script = joinpath(@__DIR__, "checkpoint_restart_resume.jl")

        save_cmd = run_restart_process(save_script, directory)
        @test save_cmd.exitcode == 0
        if save_cmd.exitcode != 0
            println(
                stderr,
                "checkpoint save stderr:\n",
                read(joinpath(directory, basename(save_script) * ".stderr.log"), String),
            )
            return
        end

        resume_cmd = run_restart_process(resume_script, directory)
        @test resume_cmd.exitcode == 0
        if resume_cmd.exitcode != 0
            println(
                stderr,
                "checkpoint resume stderr:\n",
                read(joinpath(directory, basename(resume_script) * ".stderr.log"), String),
            )
            return
        end

        saved = TOML.parsefile(joinpath(directory, "saved.toml"))
        resumed = TOML.parsefile(joinpath(directory, "resumed.toml"))
        @test resumed["pre_energy"] ≈ saved["energy"] atol=1e-10 rtol=0
        @test resumed["pre_density"] ≈ saved["density"] atol=1e-10 rtol=0
        @test resumed["pre_schmidt_qns"] == saved["schmidt_qns"]
        @test resumed["pre_site_indices"] == saved["site_indices"]
        @test resumed["post_maxlinkdim"] > resumed["pre_maxlinkdim"]
        @test resumed["next_flux_iteration_completed"] === true
    end
end
