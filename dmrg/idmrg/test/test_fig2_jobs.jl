using Test

function run_fig2_submitter_dry_run(
    submitter::AbstractString,
    repository::AbstractString,
    threads::Integer,
)
    return mktempdir() do directory
        stdout_path = joinpath(directory, "stdout.log")
        stderr_path = joinpath(directory, "stderr.log")
        command = addenv(
            `bash $submitter`,
            "W003_REPO" => repository,
            "DRY_RUN" => "1",
            "FIG2_MANIFEST" => joinpath(
                repository,
                "dmrg",
                "idmrg",
                "benchmarks",
                "fqahc_fig2.toml",
            ),
            "FIG2_STAGE" => "thread_allowlist_test",
            "FIG2_OUTPUT" =>
                "/home/public/shajy/codex/results/fqahc-fig2/thread_allowlist_test",
            "FIG2_DIMENSIONS" => "32",
            "FIG2_FLUX_UNITS" => "0",
            "FIG2_THREADS" => string(threads),
            "FIG2_DEPENDENCY" => "",
            "FIG2_WALLTIME" => "12:00:00",
        )
        process = open(stdout_path, "w") do stdout_io
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
        return (
            exitcode=process.exitcode,
            stdout=read(stdout_path, String),
            stderr=read(stderr_path, String),
        )
    end
end

@testset "Fig. 2 W003 production job contract" begin
    job_directory = joinpath(@__DIR__, "..", "jobs")
    runner = joinpath(job_directory, "run_fig2_stage.pbs")
    submitter = joinpath(job_directory, "submit_fig2_stage.sh")

    @test isfile(runner)
    @test isfile(submitter)

    if isfile(runner)
        source = read(runner, String)
        @test success(`bash -n $runner`)
        @test occursin("#PBS -q cmt", source)
        @test occursin("#PBS -l nodes=1:ppn=24", source)
        @test occursin("#PBS -l walltime=12:00:00", source)
        @test occursin("JULIA_NUM_THREADS=24", source)
        @test occursin("OMP_NUM_THREADS=1", source)
        @test occursin("OPENBLAS_NUM_THREADS=1", source)
        @test occursin("MKL_NUM_THREADS=1", source)
        thread_allowlist_contract =
            "case \"\$FIG2_THREADS\" in\n" *
            "  4|12|24) ;;\n" *
            "  *)\n" *
            "    echo \"FIG2_THREADS must be one of 4, 12, or 24\" >&2\n" *
            "    exit 2\n" *
            "    ;;\n" *
            "esac\n"
        @test occursin(thread_allowlist_contract, source)
        @test occursin("/usr/bin/flock -n", source)
        @test occursin(
            "exec 9>> \"\$output/.fig2_stage.lock\"\n" *
            "if ! /usr/bin/flock -n 9; then\n" *
            "  echo \"another PBS job holds the Fig. 2 output lock: \$output\" >&2\n" *
            "  exit 3\n" *
            "fi\n" *
            "printf '%s\\n' \"\$PBS_JOBID\" > \"\$output/.fig2_stage.lock\"",
            source,
        )
        @test occursin("git diff --quiet", source)
        @test occursin("git diff --cached --quiet", source)
        @test occursin("origin/DMRG", source)
        @test occursin("Pkg.status", source)
        @test occursin("sha256sum", source)
        @test occursin("qstat -f", source)
        @test occursin("/usr/bin/time -v", source)
        @test occursin("run_fig2_benchmark.jl", source)
        for argument in (
            "--manifest=",
            "--stage=",
            "--output=",
            "--dimensions=",
            "--flux_units_2pi=",
            "--threads=",
        )
            @test occursin(argument, source)
        end
        @test occursin("job_status.toml", source)
    end

    if isfile(submitter)
        source = read(submitter, String)
        @test success(`bash -n $submitter`)
        @test occursin("DRY_RUN", source)
        @test occursin("afterok", source)
        @test occursin("FIG2_STAGE", source)
        @test occursin("FIG2_OUTPUT", source)
        @test occursin("FIG2_DIMENSIONS", source)
        @test occursin("FIG2_FLUX_UNITS", source)
        @test occursin("FIG2_WALLTIME", source)
        @test occursin("qsub", source)
        @test occursin("run_fig2_stage.pbs", source)
        @test occursin(
            "/home/public/shajy/codex/results/fqahc-fig2/", source
        )

        repository = normpath(joinpath(@__DIR__, "..", "..", ".."))
        @testset "submitter accepts the declared thread counts" begin
            for threads in (4, 12, 24)
                result = run_fig2_submitter_dry_run(
                    submitter,
                    repository,
                    threads,
                )
                @test result.exitcode == 0
                @test occursin("FIG2_THREADS=$threads\n", result.stdout)
                @test isempty(result.stderr)
            end
        end

        @testset "submitter rejects undeclared thread counts" begin
            rejection_message =
                "FIG2_THREADS must be one of 4, 12, or 24\n"
            for threads in (1, 5, 0, 25)
                result = run_fig2_submitter_dry_run(
                    submitter,
                    repository,
                    threads,
                )
                @test result.exitcode == 2
                @test result.stderr == rejection_message
            end
        end
    end
end
