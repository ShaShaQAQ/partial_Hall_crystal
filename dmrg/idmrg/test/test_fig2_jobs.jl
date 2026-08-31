using Test
using LibGit2
using SHA

function run_fig2_submitter_dry_run(
    submitter::AbstractString,
    repository::AbstractString,
    threads::Union{Nothing,Integer},
    ;
    dimensions::AbstractString="32",
    walltime::AbstractString="12:00:00",
)
    return mktempdir() do directory
        stdout_path = joinpath(directory, "stdout.log")
        stderr_path = joinpath(directory, "stderr.log")
        environment = Pair{String,Union{Nothing,String}}[
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
            "FIG2_DIMENSIONS" => dimensions,
            "FIG2_FLUX_UNITS" => "0",
            "FIG2_DEPENDENCY" => "",
            "FIG2_WALLTIME" => walltime,
        ]
        push!(
            environment,
            "FIG2_THREADS" =>
                (isnothing(threads) ? nothing : string(threads)),
        )
        command = addenv(`bash $submitter`, environment...)
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

function run_fig2_submitter_generation(
    submitter::AbstractString,
    repository::AbstractString;
    qsub_exit_code::Integer=0,
)
    return mktempdir() do directory
        qsub_capture = joinpath(directory, "qsub_args.txt")
        launcher_snapshot = joinpath(directory, "launcher.pbs")
        mock_qsub = joinpath(directory, "qsub")
        write(
            mock_qsub,
            """#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "\$@" > "\$QSUB_CAPTURE"
last_argument=\${!#}
cp "\$last_argument" "\$QSUB_LAUNCHER_SNAPSHOT"
if [[ "\$QSUB_EXIT_CODE" != 0 ]]; then
  exit "\$QSUB_EXIT_CODE"
fi
printf '999999.w003\\n'
""",
        )
        chmod(mock_qsub, 0o700)

        stdout_path = joinpath(directory, "stdout.log")
        stderr_path = joinpath(directory, "stderr.log")
        stage = "job_contract_$(getpid())_$(time_ns())"
        command = addenv(
            `bash $submitter`,
            "W003_REPO" => repository,
            "QSUB_BIN" => mock_qsub,
            "QSUB_CAPTURE" => qsub_capture,
            "QSUB_LAUNCHER_SNAPSHOT" => launcher_snapshot,
            "QSUB_EXIT_CODE" => string(qsub_exit_code),
            "DRY_RUN" => "0",
            "FIG2_MANIFEST" => joinpath(
                repository,
                "dmrg",
                "idmrg",
                "benchmarks",
                "fqahc_fig2.toml",
            ),
            "FIG2_STAGE" => stage,
            "FIG2_OUTPUT" =>
                "/home/public/shajy/codex/results/fqahc-fig2/$stage",
            "FIG2_DIMENSIONS" => "32",
            "FIG2_FLUX_UNITS" => "0",
            "FIG2_THREADS" => "4",
            "FIG2_DEPENDENCY" => "",
            "FIG2_WALLTIME" => "00:20:00",
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

        stdout = read(stdout_path, String)
        stderr = read(stderr_path, String)
        qsub_arguments =
            isfile(qsub_capture) ? readlines(qsub_capture) : String[]
        snapshot =
            isfile(launcher_snapshot) ? read(launcher_snapshot, String) : ""
        config_match = match(
            r"^export FIG2_JOB_CONFIG=(.+)$"m,
            snapshot,
        )
        config = isnothing(config_match) ? "" : only(config_match.captures)
        launcher = isempty(qsub_arguments) ? "" : last(qsub_arguments)
        config_source = isfile(config) ? read(config, String) : ""
        manifest_match = match(r"^FIG2_MANIFEST=(.+)$"m, config_source)
        job_manifest = isnothing(manifest_match) ? "" :
            only(manifest_match.captures)
        config_mode = isfile(config) ? stat(config).mode & 0o777 : nothing
        launcher_mode =
            isfile(launcher) ? stat(launcher).mode & 0o777 : nothing
        manifest_mode =
            isfile(job_manifest) ? stat(job_manifest).mode & 0o777 : nothing
        manifest_bytes =
            isfile(job_manifest) ? read(job_manifest) : UInt8[]
        config_exists_after = isfile(config)
        launcher_exists_after = isfile(launcher)
        manifest_exists_after = isfile(job_manifest)

        job_config_root =
            "/home/public/shajy/codex/results/fqahc-fig2/job_configs/"
        submission_files_after = if isdir(job_config_root)
            sort(filter(
                path -> startswith(basename(path), "$stage."),
                readdir(job_config_root; join=true),
            ))
        else
            String[]
        end
        for path in (config, launcher, job_manifest)
            if startswith(path, job_config_root) &&
                    startswith(basename(path), "$stage.")
                rm(path; force=true)
            end
        end

        return (
            exitcode=process.exitcode,
            stdout,
            stderr,
            qsub_arguments,
            snapshot,
            config,
            config_source,
            launcher,
            job_manifest,
            config_mode,
            launcher_mode,
            manifest_mode,
            manifest_bytes,
            config_exists_after,
            launcher_exists_after,
            manifest_exists_after,
            submission_files_after,
            stage,
        )
    end
end

function invoke_checkout_auditor(
    auditor::AbstractString,
    repository::AbstractString,
)
    return mktempdir() do directory
        stdout_path = joinpath(directory, "stdout.log")
        stderr_path = joinpath(directory, "stderr.log")
        command = `$(Base.julia_cmd()) --startup-file=no $auditor $repository`
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

function exercise_checkout_auditor(
    auditor::AbstractString,
    source_repository::AbstractString,
)
    return mktempdir() do directory
        checkout = joinpath(directory, "checkout")
        cloned_repository = LibGit2.clone(
            source_repository,
            checkout;
            branch="DMRG",
        )
        close(cloned_repository)

        clean = invoke_checkout_auditor(auditor, checkout)
        tracked_path = joinpath(checkout, "dmrg", "idmrg", "Project.toml")
        original = read(tracked_path, String)
        open(tracked_path, "a") do io
            write(io, "\n# tracked dirty audit fixture\n")
        end
        tracked_dirty = invoke_checkout_auditor(auditor, checkout)

        open(tracked_path, "w") do io
            write(io, original)
        end
        open(tracked_path, "a") do io
            write(io, "\n# staged dirty audit fixture\n")
        end
        staged_repository = LibGit2.GitRepo(checkout)
        try
            LibGit2.add!(
                staged_repository,
                joinpath("dmrg", "idmrg", "Project.toml"),
            )
        finally
            close(staged_repository)
        end
        staged_dirty = invoke_checkout_auditor(auditor, checkout)

        diverged_checkout = joinpath(directory, "diverged_checkout")
        diverged_repository = LibGit2.clone(
            source_repository,
            diverged_checkout;
            branch="DMRG",
        )
        try
            parent_commit = LibGit2.revparseid(
                diverged_repository,
                "HEAD^",
            )
            LibGit2.reset!(
                diverged_repository,
                parent_commit,
                LibGit2.Consts.RESET_HARD,
            )
        finally
            close(diverged_repository)
        end
        diverged = invoke_checkout_auditor(auditor, diverged_checkout)

        return (; clean, tracked_dirty, staged_dirty, diverged)
    end
end

@testset "Fig. 2 W003 production job contract" begin
    production_repository =
        "/home/public/shajy/codex/partial_Hall_crystal-idmrg-production"
    retired_benchmark_repository =
        "/home/public/shajy/codex/partial_Hall_crystal-idmrg-benchmark"
    job_directory = joinpath(@__DIR__, "..", "jobs")
    runner = joinpath(job_directory, "run_fig2_stage.pbs")
    test_runner = joinpath(job_directory, "run_tests.pbs")
    submitter = joinpath(job_directory, "submit_fig2_stage.sh")
    contract = joinpath(job_directory, "fig2_job_contract.sh")
    auditor = joinpath(
        @__DIR__,
        "..",
        "bin",
        "audit_production_checkout.jl",
    )
    benchmark_source = joinpath(
        @__DIR__,
        "..",
        "src",
        "Fig2Benchmark.jl",
    )

    @test isfile(runner)
    @test isfile(test_runner)
    @test isfile(submitter)
    @test isfile(contract)
    @test isfile(auditor)
    @test isfile(benchmark_source)

    if isfile(runner)
        source = read(runner, String)
        @test success(`bash -n $runner`)
        @test occursin("#PBS -q cmt", source)
        @test occursin("#PBS -l nodes=1:ppn=24", source)
        @test occursin("#PBS -l walltime=12:00:00", source)
        @test occursin("JULIA_NUM_THREADS=24", source)
        @test occursin("export JULIA_DEPOT_PATH\n", source)
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
        @test occursin("origin/DMRG", source)
        @test occursin("Pkg.status", source)
        @test occursin("sha256sum", source)
        @test occursin("qstat -f", source)
        @test occursin("/usr/bin/time -v", source)
        @test occursin("run_fig2_benchmark.jl", source)
        @test occursin("audit_production_checkout.jl", source)
        @test occursin("checkout_audit", source)
        @test !occursin("git diff", source)
        @test !occursin("git status", source)
        @test !occursin("git rev-parse", source)
        @test occursin(
            ": \"\${FIG2_JOB_LAUNCHER:?FIG2_JOB_LAUNCHER is required}\"",
            source,
        )
        @test occursin("job_launcher=%s", source)
        result_root_position = findfirst("result_root=", source)
        config_requirement_position = findfirst(
            ": \"\${FIG2_JOB_CONFIG:?FIG2_JOB_CONFIG is required}\"",
            source,
        )
        @test !isnothing(result_root_position)
        @test !isnothing(config_requirement_position)
        if !isnothing(result_root_position) &&
                !isnothing(config_requirement_position)
            @test first(result_root_position) < first(config_requirement_position)
        end
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
        @test occursin("FIG2_MANIFEST_SHA256", source)
        @test occursin("fig2_validate_walltime", source)
        @test occursin("pbs_walltime", source)
        @test occursin(production_repository, source)
        @test !occursin(retired_benchmark_repository, source)
    end

    if isfile(contract)
        @test success(`bash -n $contract`)
    end

    if isfile(test_runner)
        source = read(test_runner, String)
        @test success(`bash -n $test_runner`)
        @test occursin(production_repository, source)
        @test !occursin(retired_benchmark_repository, source)
        @test occursin(": \"\${IDMRG_TEST_THREADS:=4}\"", source)
        @test occursin(
            "export JULIA_NUM_THREADS=\"\$IDMRG_TEST_THREADS\"",
            source,
        )
        @test !occursin("JULIA_NUM_THREADS=\${NCPUS:-24}", source)
    end

    if isfile(auditor)
        source = read(auditor, String)
        @test occursin("using LibGit2", source)
        @test occursin("LibGit2.isdirty", source)
        @test occursin("refs/remotes/origin/DMRG", source)

        source_repository = get(
            ENV,
            "W003_REPO",
            normpath(joinpath(@__DIR__, "..", "..", "..")),
        )
        @testset "compute-node checkout audit is fail closed" begin
            result = exercise_checkout_auditor(auditor, source_repository)
            @test result.clean.exitcode == 0
            @test occursin("tracked_dirty=false\n", result.clean.stdout)
            @test occursin("staged_dirty=false\n", result.clean.stdout)
            clean_stderr = replace(
                result.clean.stderr,
                "WARNING: failed to select UTF-8 encoding, using ASCII\n" => "",
            )
            @test isempty(clean_stderr)

            @test result.tracked_dirty.exitcode == 4
            @test occursin(
                "tracked_dirty=true\n",
                result.tracked_dirty.stdout,
            )
            @test occursin(
                "tracked worktree changes",
                result.tracked_dirty.stderr,
            )

            @test result.staged_dirty.exitcode == 4
            @test occursin(
                "staged_dirty=true\n",
                result.staged_dirty.stdout,
            )
            @test occursin(
                "staged index changes",
                result.staged_dirty.stderr,
            )

            @test result.diverged.exitcode == 4
            @test occursin(
                "W003 HEAD must equal origin/DMRG",
                result.diverged.stderr,
            )
        end
    end

    if isfile(benchmark_source)
        source = read(benchmark_source, String)
        @test !occursin("Cmd([\"git\"", source)
        @test occursin("audit_production_checkout.jl", source)
        @test occursin("_fig2_checkout_audit", source)
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
        @test occursin(production_repository, source)
        @test !occursin(retired_benchmark_repository, source)
        @test occursin("launcher=\$(mktemp", source)
        @test occursin("export FIG2_JOB_CONFIG=%q", source)
        @test occursin("export FIG2_JOB_LAUNCHER=%q", source)
        @test occursin(
            "chmod 0444 \"\$job_manifest\" \"\$config\" \"\$launcher\"",
            source,
        )
        @test !occursin("-v \"FIG2_JOB_CONFIG=", source)
        qsub_override_supported = occursin("QSUB_BIN", source)
        @test qsub_override_supported

        repository = normpath(joinpath(@__DIR__, "..", "..", ".."))
        if qsub_override_supported
            @testset "submitter generates the launcher passed to qsub" begin
                result = run_fig2_submitter_generation(
                    submitter,
                    repository,
                )
                @test result.exitcode == 0
                @test result.config_mode == 0o444
                @test result.launcher_mode == 0o444
                @test result.manifest_mode == 0o444
                @test result.config_exists_after
                @test result.launcher_exists_after
                @test result.manifest_exists_after
                @test last(result.qsub_arguments) == result.launcher
                @test !any(==("-v"), result.qsub_arguments)
                source_manifest = joinpath(
                    repository,
                    "dmrg",
                    "idmrg",
                    "benchmarks",
                    "fqahc_fig2.toml",
                )
                @test result.job_manifest != source_manifest
                @test startswith(
                    result.job_manifest,
                    "/home/public/shajy/codex/results/fqahc-fig2/job_configs/",
                )
                @test endswith(result.job_manifest, ".manifest.toml")
                @test result.manifest_bytes == read(source_manifest)
                @test occursin(
                    "FIG2_MANIFEST=$(result.job_manifest)\n",
                    result.config_source,
                )
                @test occursin(
                    "FIG2_MANIFEST_SHA256=$(bytes2hex(sha256(result.manifest_bytes)))\n",
                    result.config_source,
                )
                @test occursin(
                    "job_manifest=$(result.job_manifest)\n",
                    result.stdout,
                )
                @test occursin("#PBS -q cmt\n", result.snapshot)
                @test occursin(
                    "#PBS -l nodes=1:ppn=24\n",
                    result.snapshot,
                )
                @test occursin(
                    "#PBS -l walltime=12:00:00\n",
                    result.snapshot,
                )
                @test occursin(
                    "export FIG2_JOB_CONFIG=$(result.config)\n",
                    result.snapshot,
                )
                @test occursin(
                    "export FIG2_JOB_LAUNCHER=$(result.launcher)\n",
                    result.snapshot,
                )
                exec_match = match(r"^exec (.+)$"m, result.snapshot)
                @test !isnothing(exec_match)
                if !isnothing(exec_match)
                    @test normpath(only(exec_match.captures)) ==
                        normpath(runner)
                end
            end

            @testset "failed qsub removes orphan submission files" begin
                result = run_fig2_submitter_generation(
                    submitter,
                    repository;
                    qsub_exit_code=17,
                )
                @test result.exitcode == 17
                @test !result.config_exists_after
                @test !result.launcher_exists_after
                @test !result.manifest_exists_after
                @test isempty(result.submission_files_after)
            end
        end

        @testset "submitter accepts the declared thread counts" begin
            for threads in (4, 12, 24)
                result = run_fig2_submitter_dry_run(
                    submitter,
                    repository,
                    threads,
                )
                @test result.exitcode == 0
                @test occursin("FIG2_THREADS=$threads\n", result.stdout)
                @test occursin(
                    "FIG2_JOB_LAUNCHER=<generated-on-submit>\n",
                    result.stdout,
                )
                @test !occursin("-v FIG2_JOB_CONFIG", result.stdout)
                @test isempty(result.stderr)
            end
        end

        @testset "submitter defaults to the measured 4-thread D=32 setting" begin
            result = run_fig2_submitter_dry_run(
                submitter,
                repository,
                nothing,
            )
            @test result.exitcode == 0
            @test occursin("FIG2_THREADS=4\n", result.stdout)
            @test isempty(result.stderr)
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

        @testset "submitter enforces dimension-specific walltime caps" begin
            for (dimensions, cap, over_cap) in (
                ("32,64,128", "12:00:00", "12:00:01"),
                ("256", "36:00:00", "36:00:01"),
                ("512", "72:00:00", "72:00:01"),
                ("1000", "72:00:00", "72:00:01"),
                ("2000", "120:00:00", "120:00:01"),
                ("3000", "120:00:00", "120:00:01"),
            )
                accepted = run_fig2_submitter_dry_run(
                    submitter,
                    repository,
                    24;
                    dimensions,
                    walltime=cap,
                )
                @test accepted.exitcode == 0
                @test occursin("FIG2_WALLTIME=$cap\n", accepted.stdout)
                @test isempty(accepted.stderr)

                rejected = run_fig2_submitter_dry_run(
                    submitter,
                    repository,
                    24;
                    dimensions,
                    walltime=over_cap,
                )
                @test rejected.exitcode == 2
                @test occursin("walltime", lowercase(rejected.stderr))
                @test occursin("cap", lowercase(rejected.stderr))
            end
        end
    end
end
