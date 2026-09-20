using Test
using TOML
using InfiniteCylinderDMRG

import MPSKit

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

function mpskit_fig2_fixture_solver(
    hamiltonian,
    initial_state;
    maxdim_schedule,
    cutoff,
    idmrg_maxiter,
    vumps_maxiter,
    galerkin_tol,
    energy_imag_tol=1e-10,
    verbosity=0,
    progress_callback=(args...) -> nothing,
)
    _ = (galerkin_tol, energy_imag_tol, verbosity)
    length(maxdim_schedule) == 1 || throw(
        ArgumentError("fixture solver requires exactly one schedule stage"),
    )
    environments = MPSKit.environments(
        initial_state,
        hamiltonian,
        initial_state,
    )
    energy = ComplexF64(MPSKit.expectation_value(
        initial_state,
        hamiltonian,
        environments,
    )) / length(initial_state)
    dimensions = InfiniteCylinderDMRG._mpskit_solver_link_dimensions(
        initial_state
    )
    record = MPSKitSolverStageRecord(
        1,
        only(maxdim_schedule),
        maximum(dimensions),
        Float64(cutoff),
        idmrg_maxiter,
        vumps_maxiter,
        real(energy),
        imag(energy),
        0.0,
        0.0,
        0.0,
        0.01,
    )
    progress_callback(initial_state, environments, [record])
    return MPSKitSolverResult(
        initial_state,
        environments,
        [record],
        real(energy),
        imag(energy),
        0.0,
        0.0,
        dimensions,
        true,
    )
end

function mpskit_fig2_fixture_provenance(spec, output, runtime_seconds)
    active_project = abspath(Base.active_project())
    project_manifest = joinpath(dirname(active_project), "Manifest.toml")
    benchmark_source = normpath(joinpath(
        @__DIR__, "..", "src", "Fig2Benchmark.jl"
    ))
    base = Dict{String,Any}(
        "format" => "fqahc_fig2_provenance_v2",
        "manifest_sha256" => spec.sha256,
        "git_commit" => InfiniteCylinderDMRG._fig2_repository_commit(),
        "git_tree_clean" => true,
        "julia_version" => string(VERSION),
        "pbs_job_id" => get(ENV, "PBS_JOBID", "adapter-test.w003"),
        "threads" => Threads.nthreads(),
        "blas_threads" => 1,
        "strided_threads" => 1,
        "blocksparse_threaded" => Threads.nthreads() > 1,
        "runtime_seconds" => runtime_seconds,
        "active_project" => active_project,
        "project_manifest" => project_manifest,
        "project_manifest_sha256" =>
            InfiniteCylinderDMRG._fig2_file_sha256(project_manifest),
        "benchmark_source" => benchmark_source,
        "benchmark_source_sha256" =>
            InfiniteCylinderDMRG._fig2_file_sha256(benchmark_source),
    )
    return InfiniteCylinderDMRG._mpskit_fig2_provenance(
        spec,
        output,
        runtime_seconds;
        base_provenance=(args...) -> deepcopy(base),
    )
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

@testset "MPSKit checkpoint operations are backend-neutral" begin
    required_fields = Set((
        :backend_id,
        :checkpoint_filename,
        :checkpoint_save,
        :checkpoint_load,
        :state_maxlinkdim,
    ))
    @test required_fields ⊆ Set(fieldnames(Fig2BenchmarkOperations))
    @test isdefined(InfiniteCylinderDMRG, :mpskit_fig2_operations)
    @test isdefined(InfiniteCylinderDMRG, :fig2_operations_for_backend)

    if all(
        name -> isdefined(InfiniteCylinderDMRG, name),
        (:mpskit_fig2_operations, :fig2_operations_for_backend),
    ) && required_fields ⊆ Set(fieldnames(Fig2BenchmarkOperations))
        spec = load_fig2_benchmark(MPSKIT_FIG2_MANIFEST_PATH)
        direct = mpskit_fig2_operations(spec)
        selected = fig2_operations_for_backend(spec)
        @test direct.backend_id == "mpskit_idmrg_v1"
        @test selected.backend_id == direct.backend_id
        @test direct.checkpoint_filename == "state.h5"

        occupied_sites = first(fig2_initial_candidates(spec.config)).occupied_sites
        state = mpskit_product_state(spec.config, occupied_sites)
        metadata = (
            candidate_id="adapter_fixture",
            requested_maxdim=1,
            completed_stage=0,
        )
        mktempdir() do directory
            checkpoint = joinpath(directory, direct.checkpoint_filename)
            @test direct.checkpoint_save(
                checkpoint,
                state,
                spec.config,
                metadata,
            ) == checkpoint
            reloaded = direct.checkpoint_load(checkpoint, spec.config)
            @test !(reloaded isa MPSKitCheckpointData)
            @test direct.state_maxlinkdim(reloaded) == 1
            @test mpskit_state_space_fingerprints(reloaded, spec.config) ==
                mpskit_state_space_fingerprints(state, spec.config)

            audit = direct.checkpoint_audit(
                spec,
                checkpoint,
                spec.config.phi_y,
            )
            @test audit == (
                restart_valid=true,
                checkpoint_maxlinkdim=1,
            )
        end
    end
end

@testset "MPSKit progress generations reload through backend hooks" begin
    required = (
        :MPSKitFig2ProgressEvent,
        :_mpskit_fig2_persist_progress_event!,
        :_mpskit_fig2_finalize_progress!,
    )
    @test all(name -> isdefined(InfiniteCylinderDMRG, name), required)
    if all(name -> isdefined(InfiniteCylinderDMRG, name), required)
        spec = load_fig2_benchmark(MPSKIT_FIG2_MANIFEST_PATH)
        operations = mpskit_fig2_operations(spec)
        occupied_sites = first(fig2_initial_candidates(spec.config)).occupied_sites
        state = mpskit_product_state(spec.config, occupied_sites)
        record = MPSKitSolverStageRecord(
            1,
            1,
            1,
            1.0e-9,
            1,
            1,
            -0.25,
            0.0,
            0.0,
            0.0,
            0.0,
            0.01,
        )
        event = InfiniteCylinderDMRG.MPSKitFig2ProgressEvent(
            1,
            state,
            record,
        )
        candidate_id = "adapter_progress"
        mktempdir() do directory
            persisted =
                InfiniteCylinderDMRG._mpskit_fig2_persist_progress_event!(
                    spec,
                    directory,
                    1,
                    1,
                    spec.config.phi_y,
                    candidate_id,
                    [1],
                    event;
                    resume_count=0,
                )
            @test persisted.event_sequence == 1
            @test persisted.maxlinkdim == 1
            @test isfile(joinpath(directory, persisted.state_path))

            final_checkpoint = joinpath(
                directory,
                operations.checkpoint_filename,
            )
            metadata = (
                candidate_id,
                requested_maxdim=1,
                completed_stage=1,
            )
            operations.checkpoint_save(
                final_checkpoint,
                state,
                spec.config,
                metadata,
            )
            finalized =
                InfiniteCylinderDMRG._mpskit_fig2_finalize_progress!(
                    spec,
                    directory,
                    1,
                    1,
                    spec.config.phi_y,
                    candidate_id,
                    final_checkpoint,
                )
            audited = operations.progress_audit(
                spec,
                directory,
                1,
                1,
                spec.config.phi_y,
                candidate_id,
                final_checkpoint,
            )
            @test finalized == audited
            @test audited.complete
            @test audited.event_count == 1
            @test audited.resume_count == 0
            @test audited.latest_maxlinkdim == 1
            @test audited.final_state_sha256 ==
                InfiniteCylinderDMRG._fig2_file_sha256(final_checkpoint)
            @test occursin(r"^[0-9a-f]{64}$", audited.progress_sha256)
        end
    end
end

@testset "MPSKit provenance and restart pin audit" begin
    spec = load_fig2_benchmark(MPSKIT_FIG2_MANIFEST_PATH)
    @test spec.data["restart_gate"]["backend_commit"] ==
        MPSKIT_BACKEND_COMMIT
    @test isdefined(InfiniteCylinderDMRG, :_mpskit_fig2_provenance)

    active_project = abspath(Base.active_project())
    project_manifest = joinpath(dirname(active_project), "Manifest.toml")
    backend = spec.data["backend"]
    pin_auditor = InfiniteCylinderDMRG._fig2_pinned_backend_revision
    @test applicable(pin_auditor, active_project, project_manifest, backend)
    if applicable(pin_auditor, active_project, project_manifest, backend)
        @test pin_auditor(active_project, project_manifest, backend) ==
            MPSKIT_BACKEND_COMMIT
    end

    if isdefined(InfiniteCylinderDMRG, :_mpskit_fig2_provenance)
        base = Dict{String,Any}(
            "format" => "fqahc_fig2_provenance_v2",
            "manifest_sha256" => spec.sha256,
            "git_commit" => repeat("a", 40),
        )
        provenance = InfiniteCylinderDMRG._mpskit_fig2_provenance(
            spec,
            "unused",
            1.25;
            base_provenance=(args...) -> deepcopy(base),
        )
        @test provenance["backend_id"] == "mpskit_idmrg_v1"
        @test provenance["mpskit_commit"] == MPSKIT_BACKEND_COMMIT
        @test provenance["tensorkittensors_commit"] ==
            TENSORKITTENSORS_COMMIT
        @test provenance["blocktensorkit_commit"] == BLOCKTENSORKIT_COMMIT
        adapter_source = normpath(joinpath(
            @__DIR__, "..", "src", "MPSKitFig2Adapter.jl"
        ))
        @test provenance["backend_adapter_source"] == adapter_source
        @test provenance["backend_adapter_source_sha256"] ==
            InfiniteCylinderDMRG._fig2_file_sha256(adapter_source)
    end
end

@testset "MPSKit adapter runs the audited candidate workflow" begin
    required = (
        :MPSKIT_FIG2_ALGORITHM,
        :_mpskit_fig2_write_solver_outputs!,
    )
    @test all(name -> isdefined(InfiniteCylinderDMRG, name), required)
    if all(name -> isdefined(InfiniteCylinderDMRG, name), required)
        spec = load_fig2_benchmark(MPSKIT_FIG2_MANIFEST_PATH)
        candidate_id = first(fig2_initial_candidates(spec.config)).id
        operations = mpskit_fig2_operations(
            spec;
            solver=mpskit_fig2_fixture_solver,
            provenance=mpskit_fig2_fixture_provenance,
            candidate_ids=(args...) -> [candidate_id],
        )
        mktempdir() do directory
            run = run_fig2_benchmark(
                spec,
                directory;
                stage="adapter_fixture",
                dimensions=[1],
                fluxes=[0.0],
                operations,
            )
            @test length(run.selections) == 1
            selection = only(run.selections)
            @test selection.candidate_id == candidate_id
            @test selection.restart_valid

            candidate_directory = joinpath(directory, selection.directory)
            for filename in InfiniteCylinderDMRG.FIG2_REQUIRED_CANDIDATE_FILES
                @test isfile(joinpath(candidate_directory, filename))
            end
            summary = TOML.parsefile(joinpath(
                candidate_directory, "summary.toml"
            ))
            @test summary["algorithm"] ==
                InfiniteCylinderDMRG.MPSKIT_FIG2_ALGORITHM
            @test summary["backend"]["id"] == "mpskit_idmrg_v1"
            @test summary["optimization"]["maxdim_schedule"] == [1]
            @test summary["optimization"]["recomputed_galerkin_gate"] ==
                1.0e-6

            convergence = readlines(joinpath(
                candidate_directory, "convergence.tsv"
            ))
            @test length(convergence) == 2
            @test first(convergence) ==
                InfiniteCylinderDMRG.FIG2_ARTIFACT_HEADERS["convergence.tsv"]
            @test length(readlines(joinpath(
                candidate_directory, "density.tsv"
            ))) == sites_per_cell(spec.config) + 1

            ledger = TOML.parsefile(joinpath(directory, "ledger.toml"))
            @test length(ledger["candidate"]) == 1
            @test length(ledger["selection"]) == 1
            @test only(ledger["candidate"])["restart_valid"]
            @test only(ledger["candidate"])["progress_event_count"] == 1
            @test isfile(joinpath(directory, "pump_raw.tsv"))
            raw_rows = readlines(joinpath(directory, "pump_raw.tsv"))
            @test length(raw_rows) == 2
        end
    end
end
end
