using Test
using TOML
using InfiniteCylinderDMRG
using ITensors
using ITensorMPS
using Random

const FIG2_MANIFEST_PATH = normpath(joinpath(
    @__DIR__, "..", "benchmarks", "fqahc_fig2.toml"
))
const FIG2_INTERRUPT_PROBES = Dict{Symbol,Any}()
const FIG2_SYNTHETIC_MANIFEST = TOML.parsefile(FIG2_MANIFEST_PATH)
const FIG2_SYNTHETIC_LY = 2 * Int(FIG2_SYNTHETIC_MANIFEST["Ny"])
const FIG2_SYNTHETIC_SITES_PER_CELL =
    FIG2_SYNTHETIC_LY * Int(FIG2_SYNTHETIC_MANIFEST["x_period"])

struct Fig2WarmScheduleProbeState{T}
    AL::T
    maxlinkdim::Int
end

InfiniteCylinderDMRG.link_dimensions(state::Fig2WarmScheduleProbeState) =
    [state.maxlinkdim]

struct Fig2WarmScheduleCaptured <: Exception end

@testset "immutable Fig. 2 benchmark manifest" begin
    manifest = TOML.parsefile(FIG2_MANIFEST_PATH)
    @test manifest["format"] == "fqahc_fig2_benchmark_v4"
    @test manifest["geometry"] == "paper_straight"
    @test manifest["Ny"] == 6
    @test manifest["x_period"] == 3
    @test (manifest["filling_num"], manifest["filling_den"]) == (7, 9)
    @test (
        manifest["t1"],
        manifest["t3"],
        manifest["V1"],
        manifest["V2"],
        manifest["V3"],
    ) == (1.0, 0.2, 10.0, 2.0, 2.0)
    @test manifest["flux_start"] == 0.0
    @test manifest["flux_stop"] == 6pi
    @test manifest["expected_endpoint_pump"] == 1.0
    @test manifest["expected_sector_shift"] == 1
    @test manifest["pump_tolerance"] == 0.05
    @test manifest["counting"] == [1, 1, 2, 3, 5]
    @test manifest["momentum_counting_rule"] == Dict(
        "format" => "fqahc_fig2_momentum_counting_v1",
        "reference_sector" => "largest_eigenvalue_raw_charge_sector",
        "delta_ky" => [0, -1, -2, -3, -4],
        "cut_intercept" => 2.15,
        "cut_slope" => -1.58,
        "digitization_uncertainty" => 0.20,
        "source_pdf_sha256" =>
            "23f0691570b906ae7a9394cef2e3dd6eebd4e36d79e882a4bbd2f17824cc89a0",
        "source_page" => 3,
        "source_panel" => "f",
    )
    @test manifest["paper_curve_rule"] == Dict(
        "format" => "fqahc_fig2_paper_curve_comparison_v3",
        "source_pdf_sha256" =>
            "23f0691570b906ae7a9394cef2e3dd6eebd4e36d79e882a4bbd2f17824cc89a0",
        "source_page" => 3,
        "source_panel" => "d",
        "source_vector_path_sha256" =>
            "54e5c57f88a735e5efe0c2a24de455bbaa987bbceb5a9affa44f85fab58a892c",
        "target" => "nu_star=1/3,V1=10,V2=V3=2",
        "extraction" => "pdftocairo_svg_vector_path",
        "flux_units_2pi" => collect(0.0:0.125:3.0),
        "delta_q" => [
            0.000000, 0.039286, 0.080329, 0.122619, 0.165646,
            0.208617, 0.250794, 0.292290, 0.332370, 0.371882,
            0.410998, 0.450680, 0.490760, 0.532426, 0.575907,
            0.620918, 0.667290, 0.713152, 0.760601, 0.807029,
            0.851927, 0.889739, 0.925283, 0.961735, 1.000000,
        ],
        "uncertainty" => fill(0.02, 25),
    )
    @test manifest["convergence"]["maxdims"] == [256, 512, 1000, 2000, 3000]
    @test manifest["momentum_convention"] ==
        "sector_relative_canonical_cyclic_orbit"
    @test manifest["optimization"] == Dict(
        "multisite_update_alg" => "sequential",
        "progress_generations_to_keep" => 2,
        "cutoff" => 1.0e-9,
        "vumps_tol" => 1.0e-6,
        "energy_tol" => 1.0e-6,
        "energy_mismatch_tol" => 1.0e-6,
        "transfer_tol" => 1.0e-8,
        "momentum_residual_tol" => 1.0e-6,
        "max_iterations" => 50,
        "stable_iterations" => 2,
    )
    @test manifest["pilot"] == Dict(
        "maxdims" => [32, 64, 128],
        "flux_units_2pi" => [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        "cold_candidates" => 2,
    )
    @test manifest["convergence"] == Dict(
        "maxdims" => [256, 512, 1000, 2000, 3000],
        "energy_per_site_change_tol" => 1.0e-6,
        "endpoint_pump_change_tol" => 0.02,
        "sector_residual_tol" => 1.0e-6,
        "momentum_validation_residual_tol" => 1.0e-6,
    )
    @test get(manifest, "early_plateau", nothing) == Dict(
        "policy" => "allowed_if_all_gates_pass",
        "lower_dimension" => 1000,
        "upper_dimension" => 2000,
        "omitted_dimension" => 3000,
        "evidence_filename" => "early_plateau.toml",
    )
    @test get(manifest, "restart_gate", nothing) == Dict(
        "format" => "fqahc_fig2_restart_gate_v1",
        "evidence_filename" => "restart_gate.toml",
        "backend_commit" =>
            "765f2777703bc1138b009adbed1b97bde1973402",
        "energy_tolerance" => 1.0e-10,
        "density_tolerance" => 1.0e-10,
        "required_checks" => [
            "load_completed",
            "no_expansion_iteration_completed",
            "controlled_expansion_completed",
            "energy_reproduced",
            "density_reproduced",
            "schmidt_sectors_reproduced",
            "entanglement_spectrum_reproduced",
            "site_indices_reproduced",
            "next_flux_iteration_completed",
        ],
    )
    @test manifest["raw_pump_policy"] ==
        "no_offset_no_unwrap_no_sign_flip_no_branch_translation"
    @test manifest["acceptance"]["rows"] == [
        "endpoint_pump",
        "sector_shift",
        "momentum_counting",
        "paper_curve",
        "bond_convergence",
        "restart",
        "provenance",
    ]
    @test InfiniteCylinderDMRG.FIG2_LEDGER_FORMAT ==
        "fqahc_fig2_ledger_v3"
end

@testset "Fig. 2 uses sequential long-range QN updates" begin
    spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
    @test spec.data["optimization"]["multisite_update_alg"] == "sequential"
    @test InfiniteCylinderDMRG._fig2_runner_optimization(
        spec
    ).multisite_update_alg == :sequential
end

@testset "Fig. 2 manifest rejects scientific and numerical goalpost drift" begin
    baseline = TOML.parsefile(FIG2_MANIFEST_PATH)
    function rejects_manifest(tamper::Function)
        mktempdir() do directory
            candidate = deepcopy(baseline)
            tamper(candidate)
            path = joinpath(directory, "tampered.toml")
            open(path, "w") do io
                TOML.print(io, candidate; sorted=true)
            end
            @test_throws ArgumentError load_fig2_benchmark(path)
        end
    end

    scientific_cases = [
        ("endpoint pump", data -> data["expected_endpoint_pump"] = 1.1),
        ("sector shift", data -> data["expected_sector_shift"] = 2),
        ("pump tolerance", data -> data["pump_tolerance"] = 0.06),
        ("momentum convention", data -> data["momentum_convention"] = "absolute_raw"),
    ]
    for (label, tamper) in scientific_cases
        @testset "$label" begin
            rejects_manifest(tamper)
        end
    end

    optimization_alternatives = Dict{String,Any}(
        "cutoff" => 1.0e-8,
        "vumps_tol" => 2.0e-6,
        "energy_tol" => 2.0e-6,
        "energy_mismatch_tol" => 2.0e-6,
        "transfer_tol" => 2.0e-8,
        "momentum_residual_tol" => 2.0e-6,
        "max_iterations" => 51,
        "stable_iterations" => 3,
        "progress_generations_to_keep" => 3,
    )
    for (key, value) in optimization_alternatives
        @testset "optimization.$key" begin
            rejects_manifest(data -> data["optimization"][key] = value)
        end
    end

    pilot_cases = [
        ("pilot maxdims", data -> data["pilot"]["maxdims"] = [32, 64, 129]),
        ("pilot flux grid", data -> data["pilot"]["flux_units_2pi"] =
            [0.0, 0.5, 1.0, 1.5, 2.0, 2.6, 3.0]),
        ("pilot cold candidates", data -> data["pilot"]["cold_candidates"] = 3),
    ]
    for (label, tamper) in pilot_cases
        @testset "$label" begin
            rejects_manifest(tamper)
        end
    end

    convergence_alternatives = Dict{String,Any}(
        "maxdims" => [256, 512, 1000, 2000, 3001],
        "energy_per_site_change_tol" => 2.0e-6,
        "endpoint_pump_change_tol" => 0.03,
        "sector_residual_tol" => 2.0e-6,
        "momentum_validation_residual_tol" => 2.0e-6,
    )
    for (key, value) in convergence_alternatives
        @testset "convergence.$key" begin
            rejects_manifest(data -> data["convergence"][key] = value)
        end
    end

    early_plateau_contract = Dict{String,Any}(
        "policy" => "allowed_if_all_gates_pass",
        "lower_dimension" => 1000,
        "upper_dimension" => 2000,
        "omitted_dimension" => 3000,
        "evidence_filename" => "early_plateau.toml",
    )
    early_plateau_alternatives = Dict{String,Any}(
        "policy" => "allowed_after_run",
        "lower_dimension" => 999,
        "upper_dimension" => 2001,
        "omitted_dimension" => 3001,
        "evidence_filename" => "posthoc_plateau.toml",
    )
    for (key, value) in early_plateau_alternatives
        @testset "early_plateau.$key" begin
            rejects_manifest() do data
                data["early_plateau"] = merge(
                    copy(early_plateau_contract), Dict(key => value)
                )
            end
        end
    end
    @testset "early_plateau dimension type" begin
        rejects_manifest() do data
            data["early_plateau"] = merge(
                copy(early_plateau_contract),
                Dict("lower_dimension" => 1000.0),
            )
        end
    end

    restart_gate_alternatives = Dict{String,Any}(
        "format" => "fqahc_fig2_restart_gate_v2",
        "evidence_filename" => "restart.txt",
        "backend_commit" => repeat("0", 40),
        "required_checks" => ["load_completed"],
    )
    for (key, value) in restart_gate_alternatives
        @testset "restart_gate.$key" begin
            rejects_manifest(data -> data["restart_gate"][key] = value)
        end
    end
    strict_cases = [
        ("integer target type", data -> data["expected_sector_shift"] = 1.0),
        ("boolean numerical value", data -> data["optimization"]["cutoff"] = true),
        ("nonfinite numerical value", data -> data["optimization"]["cutoff"] = Inf),
        ("nonpositive numerical value", data -> data["optimization"]["cutoff"] = 0.0),
        ("integer schedule type", data -> data["pilot"]["maxdims"] =
            [32.0, 64.0, 128.0]),
        ("pilot maxdim ordering", data -> data["pilot"]["maxdims"] =
            [32, 128, 64]),
        ("pilot flux ordering", data -> data["pilot"]["flux_units_2pi"] =
            [0.0, 0.5, 1.5, 1.0, 2.0, 2.5, 3.0]),
        ("convergence maxdim ordering", data -> data["convergence"]["maxdims"] =
            [256, 1000, 512, 2000, 3000]),
        ("nonpositive convergence tolerance", data ->
            data["convergence"]["sector_residual_tol"] = 0.0),
        ("counting integer type", data -> data["counting"] =
            [1.0, 1.0, 2.0, 3.0, 5.0]),
        ("momentum sector integer type", data ->
            data["momentum_counting_rule"]["delta_ky"] =
                [0.0, -1.0, -2.0, -3.0, -4.0]),
    ]
    for (label, tamper) in strict_cases
        @testset "$label" begin
            rejects_manifest(tamper)
        end
    end
end

@testset "Fig. 2 production uses bounded powers-of-two maxdim staging" begin
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(2) == [2]
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(4) == [4]
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(32) ==
        [4, 8, 16, 32]
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(1000) ==
        [4, 8, 16, 32, 64, 128, 256, 512, 1000]
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(2000) ==
        [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2000]
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(8, 8) == [8]
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(8, 16) == [8, 16]
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(6, 16) == [8, 16]
    @test InfiniteCylinderDMRG._fig2_maxdim_schedule(1, 16) ==
        [4, 8, 16]
    @test_throws ArgumentError InfiniteCylinderDMRG._fig2_maxdim_schedule(0)
    @test_throws ArgumentError InfiniteCylinderDMRG._fig2_maxdim_schedule(true)
    @test_throws ArgumentError InfiniteCylinderDMRG._fig2_maxdim_schedule(17, 16)
    @test_throws ArgumentError InfiniteCylinderDMRG._fig2_maxdim_schedule(true, 16)
end


@testset "default Fig. 2 runner wires warm-compatible maxdim staging" begin
    spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
    occupied_sites = first(fig2_initial_candidates(spec.config)).occupied_sites
    _, _, product_state = initial_infinite_mps(
        spec.config; occupied_sites
    )
    warm_state = Fig2WarmScheduleProbeState(product_state.AL, 8)
    captured_schedules = Vector{Int}[]
    captured_algorithms = Symbol[]

    function capture_settings(settings, H, psi, operations; transfer_rng)
        push!(captured_schedules, copy(settings.maxdim_schedule))
        push!(captured_algorithms, settings.multisite_update_alg)
        throw(Fig2WarmScheduleCaptured())
    end
    capture_error(f::Function) = try
        f()
        nothing
    catch error
        error
    end

    mktempdir() do directory
        for (point, dimension, expected) in (
            (2, 8, [8]),
            (3, 16, [8, 16]),
        )
            error = capture_error() do
                InfiniteCylinderDMRG._default_fig2_run_candidate(
                    spec,
                    dimension,
                    point,
                    0.0,
                    "warm",
                    warm_state,
                    joinpath(directory, "warm_$point");
                    build_hamiltonian=(args...) -> :synthetic_hamiltonian,
                    run_prepared_point=capture_settings,
                )
            end
            @test error isa Fig2WarmScheduleCaptured
            @test last(captured_schedules) == expected
            @test last(captured_algorithms) == :sequential
        end

        cold_id = "cold_$(last(fig2_initial_candidates(spec.config)).id)"
        cold_error = capture_error() do
            InfiniteCylinderDMRG._default_fig2_run_candidate(
                spec,
                16,
                3,
                0.0,
                cold_id,
                warm_state,
                joinpath(directory, "cold");
                build_hamiltonian=(args...) -> :synthetic_hamiltonian,
                run_prepared_point=capture_settings,
            )
        end
        @test cold_error isa Fig2WarmScheduleCaptured
        @test last(captured_schedules) == [4, 8, 16]
        @test last(captured_algorithms) == :sequential

        oversized_error = capture_error() do
            InfiniteCylinderDMRG._default_fig2_run_candidate(
                spec,
                16,
                3,
                0.0,
                "warm",
                Fig2WarmScheduleProbeState(product_state.AL, 17),
                joinpath(directory, "oversized");
                build_hamiltonian=(args...) -> :synthetic_hamiltonian,
                run_prepared_point=capture_settings,
            )
        end
        @test oversized_error isa ArgumentError
    end
end

if get(ENV, "IDMRG_FIG2_REAL_SMOKE", "0") == "1"
    @testset "real Fig. 2 default runner proves checkpoint bond dimension" begin
        haskey(ENV, "IDMRG_TEST_ARTIFACTS") || error(
            "IDMRG_TEST_ARTIFACTS is required for the real Fig. 2 smoke"
        )
        output = joinpath(
            ENV["IDMRG_TEST_ARTIFACTS"], "fig2_default_runner_D2"
        )
        mkpath(output)
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        runner_override = InfiniteCylinderDMRG.Fig2RunnerOverride(;
            cutoff=1.0e-8,
            vumps_tol=1.0e2,
            energy_tol=1.0e6,
            energy_mismatch_tol=1.0e6,
            transfer_tol=1.0e-4,
            momentum_residual_tol=1.0e-4,
            max_iterations=3,
            stable_iterations=2,
        )
        dimension = 2
        candidate_id = first(InfiniteCylinderDMRG._default_fig2_candidate_ids(
            spec, dimension, 1, nothing
        ))
        evidence = InfiniteCylinderDMRG._default_fig2_run_candidate(
            spec,
            dimension,
            1,
            0.0,
            candidate_id,
            nothing,
            output,
            runner_override=runner_override,
        )

        checkpoint = joinpath(output, "state.h5")
        @test evidence.requested_maxdim == dimension
        @test evidence.achieved_maxlinkdim == dimension
        @test evidence.converged
        @test evidence.restart_valid
        @test isfile(checkpoint) && filesize(checkpoint) > 0
        reloaded = load_checkpoint(checkpoint, spec.config)
        checkpoint_maxlinkdim = maximum(link_dimensions(reloaded))
        @test evidence.checkpoint_maxlinkdim == checkpoint_maxlinkdim
        @test checkpoint_maxlinkdim == evidence.achieved_maxlinkdim

        summary = TOML.parsefile(joinpath(output, "summary.toml"))
        @test summary["optimization"]["multisite_update_alg"] == "sequential"
        @test summary["optimization"]["maxdim_schedule"] ==
            InfiniteCylinderDMRG._fig2_maxdim_schedule(dimension)
        smoke_optimization_data = summary["optimization"]
        smoke_optimization = (
            vumps_tol=Float64(smoke_optimization_data["vumps_tol"]),
            energy_tol=Float64(smoke_optimization_data["energy_tol"]),
            energy_mismatch_tol=Float64(
                smoke_optimization_data["energy_mismatch_tol"]
            ),
            stable_iterations=Int(
                smoke_optimization_data["stable_iterations"]
            ),
        )
        final_convergence = InfiniteCylinderDMRG._validate_fig2_convergence_tsv(
            joinpath(output, "convergence.tsv"),
            smoke_optimization_data["maxdim_schedule"];
            optimization=smoke_optimization,
            energy_normalization_sites=sites_per_cell(spec.config),
        )
        @test final_convergence.converged
        @test final_convergence.maxlinkdim == evidence.achieved_maxlinkdim

        progress_path = joinpath(output, "progress.toml")
        @test isfile(progress_path)
        progress = TOML.parsefile(progress_path)
        @test progress["format"] == "fqahc_fig2_progress_summary_v1"
        @test progress["complete"] === true
        @test progress["event_count"] >= 2
        @test progress["resume_count"] == 0
        progress_audit = InfiniteCylinderDMRG._fig2_validate_progress_artifact(
            spec,
            output,
            dimension,
            1,
            0.0,
            candidate_id,
            checkpoint,
        )
        @test progress_audit.complete
        @test progress_audit.latest_maxlinkdim == evidence.achieved_maxlinkdim

        report = Dict(
            "format" => "fqahc_fig2_real_claimed_dimension_smoke_v1",
            "candidate_id" => candidate_id,
            "requested_maxdim" => evidence.requested_maxdim,
            "achieved_maxlinkdim" => evidence.achieved_maxlinkdim,
            "checkpoint_maxlinkdim" => checkpoint_maxlinkdim,
            "evidence_checkpoint_maxlinkdim" =>
                evidence.checkpoint_maxlinkdim,
            "summary_maxdim_schedule" =>
                summary["optimization"]["maxdim_schedule"],
            "final_convergence_maxlinkdim" => final_convergence.maxlinkdim,
            "final_convergence_converged" => final_convergence.converged,
            "restart_valid" => evidence.restart_valid,
            "state_sha256" => InfiniteCylinderDMRG._fig2_file_sha256(checkpoint),
            "progress_sha256" =>
                InfiniteCylinderDMRG._fig2_file_sha256(progress_path),
        )
        open(joinpath(output, "claimed_dimension_smoke.toml"), "w") do io
            TOML.print(io, report; sorted=true)
        end
    end
end

@testset "Fig. 2 file SHA-256 is streaming" begin
    mktempdir() do directory
        small_path = joinpath(directory, "abc.bin")
        open(small_path, "w") do io
            write(io, "abc")
        end
        @test InfiniteCylinderDMRG._fig2_file_sha256(small_path) ==
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

        large_path = joinpath(directory, "large.bin")
        file_bytes = 32 * 1024 * 1024
        chunk = fill(UInt8(0xa5), 4096)
        open(large_path, "w") do io
            for _ in 1:(file_bytes ÷ length(chunk))
                write(io, chunk)
            end
        end
        @test filesize(large_path) == file_bytes
        InfiniteCylinderDMRG._fig2_file_sha256(large_path)
        GC.gc()
        allocated_bytes = @allocated InfiniteCylinderDMRG._fig2_file_sha256(
            large_path
        )
        println("fig2_file_sha256_allocated_bytes=$allocated_bytes")
        @test allocated_bytes <= 8 * 1024 * 1024
    end
end

@testset "Fig. 2 driver API exists" begin
    for name in (
        :Fig2BenchmarkSpec,
        :Fig2CandidateEvidence,
        :Fig2BenchmarkOperations,
        :Fig2BenchmarkRun,
        :load_fig2_benchmark,
        :run_fig2_benchmark,
        :write_fig2_acceptance_report!,
    )
        @test isdefined(InfiniteCylinderDMRG, name)
    end
end

@testset "Fig. 2 candidate evidence retains complete audit objects" begin
    fields = fieldnames(Fig2CandidateEvidence)
    @test :requested_maxdim in fields
    @test :achieved_maxlinkdim in fields
    @test :checkpoint_maxlinkdim in fields
    @test :mixed_fidelity in fields
    @test :momentum_counting_evidence in fields
    @test :provenance_valid ∉ fields
    @test :progress_audit in fieldnames(Fig2BenchmarkOperations)
    @test "progress.toml" in InfiniteCylinderDMRG.FIG2_REQUIRED_CANDIDATE_FILES
end

@testset "Fig. 2 restart validation executes the checkpoint loader" begin
    @test isdefined(InfiniteCylinderDMRG, :_fig2_checkpoint_restart_valid)
    if isdefined(InfiniteCylinderDMRG, :_fig2_checkpoint_restart_valid)
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        calls = NamedTuple[]
        loader = (path, config) -> begin
            push!(calls, (; path, config))
            return :restarted_state
        end
        @test InfiniteCylinderDMRG._fig2_checkpoint_restart_valid(
            "state.h5", spec.config; loader
        )
        @test calls == [(; path="state.h5", config=spec.config)]
        @test !InfiniteCylinderDMRG._fig2_checkpoint_restart_valid(
            "corrupt.h5",
            spec.config;
            loader=(args...) -> error("corrupt checkpoint"),
        )
        FIG2_INTERRUPT_PROBES[:checkpoint_restart] = try
            InfiniteCylinderDMRG._fig2_checkpoint_restart_valid(
                "interrupted.h5",
                spec.config;
                loader=(args...) -> throw(InterruptException()),
            )
        catch error
            error
        end
    end
end

@testset "Fig. 2 checkpoint audit measures the independently loaded state" begin
    @test isdefined(InfiniteCylinderDMRG, :_fig2_checkpoint_audit)
    if isdefined(InfiniteCylinderDMRG, :_fig2_checkpoint_audit)
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        calls = NamedTuple[]
        audit = InfiniteCylinderDMRG._fig2_checkpoint_audit(
            "state.h5",
            spec.config;
            loader=(path, config) -> begin
                push!(calls, (; operation=:load, path, config))
                return :independently_reloaded_state
            end,
            maxlinkdim_of=state -> begin
                push!(calls, (; operation=:measure, state))
                return 37
            end,
        )
        @test audit == (restart_valid=true, checkpoint_maxlinkdim=37)
        @test calls == [
            (; operation=:load, path="state.h5", config=spec.config),
            (; operation=:measure, state=:independently_reloaded_state),
        ]

        load_failure = InfiniteCylinderDMRG._fig2_checkpoint_audit(
            "corrupt.h5",
            spec.config;
            loader=(args...) -> error("corrupt checkpoint"),
            maxlinkdim_of=state -> error("must not measure"),
        )
        @test load_failure == (restart_valid=false, checkpoint_maxlinkdim=0)
        measurement_failure = InfiniteCylinderDMRG._fig2_checkpoint_audit(
            "state.h5",
            spec.config;
            loader=(args...) -> :state,
            maxlinkdim_of=state -> error("cannot measure checkpoint"),
        )
        @test measurement_failure ==
            (restart_valid=false, checkpoint_maxlinkdim=0)
        FIG2_INTERRUPT_PROBES[:checkpoint_audit_load] = try
            InfiniteCylinderDMRG._fig2_checkpoint_audit(
                "interrupted-load.h5",
                spec.config;
                loader=(args...) -> throw(InterruptException()),
                maxlinkdim_of=state -> error("must not measure"),
            )
        catch error
            error
        end
        FIG2_INTERRUPT_PROBES[:checkpoint_audit_measure] = try
            InfiniteCylinderDMRG._fig2_checkpoint_audit(
                "interrupted-measure.h5",
                spec.config;
                loader=(args...) -> :state,
                maxlinkdim_of=state -> throw(InterruptException()),
            )
        catch error
            error
        end
    end
end

@testset "Fig. 2 integrity digest propagates interruption" begin
    mktempdir() do directory
        path = joinpath(directory, "evidence.toml")
        open(path, "w") do io
            write(io, "format = \"synthetic\"\n")
        end
        @test isfile(path)
        FIG2_INTERRUPT_PROBES[:integrity_digest] = try
            InfiniteCylinderDMRG._fig2_integrity_digest(
                path;
                sha256=(args...) -> throw(InterruptException()),
            )
            nothing
        catch error
            error
        end
        FIG2_INTERRUPT_PROBES[:integrity_digest_ordinary] = try
            InfiniteCylinderDMRG._fig2_integrity_digest(
                path;
                sha256=(args...) -> error("synthetic hash failure"),
            )
        catch error
            error
        end
    end
end

@testset "Fig. 2 acceptance helpers propagate interruption fail closed" begin
    spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
    capture(f::Function) = try
        f()
    catch error
        error
    end

    snapshot_parse_interrupt = capture() do
        InfiniteCylinderDMRG._fig2_validated_snapshot(
            spec;
            parse_snapshot=(args...) -> throw(InterruptException()),
        )
    end
    @test snapshot_parse_interrupt isa InterruptException
    snapshot_validation_interrupt = capture() do
        InfiniteCylinderDMRG._fig2_validated_snapshot(
            spec;
            validate_snapshot=(args...) -> throw(InterruptException()),
        )
    end
    @test snapshot_validation_interrupt isa InterruptException
    @test capture() do
        InfiniteCylinderDMRG._fig2_validated_snapshot(
            spec;
            parse_snapshot=(args...) -> error("synthetic parse failure"),
        )
    end isa ArgumentError
    @test capture() do
        InfiniteCylinderDMRG._fig2_validated_snapshot(
            spec;
            validate_snapshot=(args...) -> error("synthetic validation failure"),
        )
    end isa ArgumentError

    mktempdir() do directory
        for filename in (
            "pump_raw.tsv",
            "pump_sector_tracked.tsv",
            "provenance.toml",
            "early_plateau.toml",
        )
            open(joinpath(directory, filename), "w") do io
                write(io, "synthetic\n")
            end
        end

        paper_interrupt = capture() do
            InfiniteCylinderDMRG._fig2_paper_curve_result(
                spec,
                directory,
                0,
                NamedTuple[];
                raw_reader=(args...) -> throw(InterruptException()),
            )
        end
        @test paper_interrupt isa InterruptException
        paper_ordinary = capture() do
            InfiniteCylinderDMRG._fig2_paper_curve_result(
                spec,
                directory,
                0,
                NamedTuple[];
                raw_reader=(args...) -> error("synthetic raw read failure"),
            )
        end
        @test paper_ordinary isa NamedTuple && !paper_ordinary.pass

        provenance_interrupt = capture() do
            InfiniteCylinderDMRG._fig2_provenance_result(
                spec,
                directory;
                toml_parsefile=(args...) -> throw(InterruptException()),
            )
        end
        @test provenance_interrupt isa InterruptException
        provenance_ordinary = capture() do
            InfiniteCylinderDMRG._fig2_provenance_result(
                spec,
                directory;
                toml_parsefile=(args...) -> error("synthetic provenance parse failure"),
            )
        end
        @test provenance_ordinary isa NamedTuple && !provenance_ordinary.valid

        plateau_interrupt = capture() do
            InfiniteCylinderDMRG._fig2_early_plateau_evidence(
                spec,
                directory;
                toml_parsefile=(args...) -> throw(InterruptException()),
            )
        end
        @test plateau_interrupt isa InterruptException
        plateau_ordinary = capture() do
            InfiniteCylinderDMRG._fig2_early_plateau_evidence(
                spec,
                directory;
                toml_parsefile=(args...) -> error("synthetic plateau parse failure"),
            )
        end
        @test plateau_ordinary isa NamedTuple && !plateau_ordinary.valid
    end
end

@testset "cold Fig. 2 candidates reuse previous QN site identities" begin
    spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
    candidates = fig2_initial_candidates(spec.config)
    _, _, previous_state = initial_infinite_mps(
        spec.config; occupied_sites=first(candidates).occupied_sites
    )
    previous_sites = siteinds(only, previous_state.AL)
    cold_candidate = candidates[2]
    prepared = InfiniteCylinderDMRG._prepare_fig2_candidate_state(
        spec,
        spec.config,
        "cold_$(cold_candidate.id)",
        previous_state,
    )

    @test prepared.occupied_sites == cold_candidate.occupied_sites
    @test prepared.sites == previous_sites
    @test siteinds(only, prepared.psi.AL) == previous_sites
    fidelity = mixed_transfer_fidelity(
        previous_state,
        prepared.psi,
        spec.config;
        tol=1e-8,
        rng=Random.Xoshiro(0x51e1d),
        maxiter=8,
    )
    @test fidelity isa MixedTransferFidelity
    @test !occursin("exact same site indices", fidelity.reason)
end

if all(
    name -> isdefined(InfiniteCylinderDMRG, name),
    (
        :Fig2BenchmarkSpec,
        :Fig2CandidateEvidence,
        :Fig2BenchmarkOperations,
        :Fig2BenchmarkRun,
        :load_fig2_benchmark,
        :run_fig2_benchmark,
        :write_fig2_acceptance_report!,
    ),
)
    function synthetic_momentum(; valid=true)
        validation = MomentumValidation(
            valid ? 0.0 : Inf,
            valid ? 0.0 : Inf,
            valid ? 0.0 : Inf,
            valid ? 0.0 : Inf,
            valid,
            valid ? "synthetic valid momentum" : "synthetic invalid momentum",
        )
        levels = MomentumEntanglementLevel[
            MomentumEntanglementLevel(
                level,
                probability,
                -log(probability),
                0,
                valid ? momentum : missing,
                valid ? 2pi * momentum / 6 : missing,
                valid ? 2pi * momentum / 6 : missing,
                valid ? 0.0 : NaN,
            ) for (level, (probability, momentum)) in enumerate((
                (0.40, 0),
                (0.25, 1),
                (0.15, 2),
                (0.12, 3),
                (0.08, 4),
            ))
        ]
        return MomentumEntanglementData(
            levels,
            validation,
            valid,
            validation.reason,
            valid ? :sector_relative_canonical_cyclic_orbit : :invalid,
            valid ? Dict(0 => 0) : Dict{Int,Int}(),
        )
    end

    function synthetic_counting_momentum(;
        offset=0,
        ambiguous=false,
        transfer_residual=0.0,
    )
        energies = Dict(
            0 => [0.0, 2.7],
            -1 => [1.2, 4.3],
            -2 => [1.0, 3.0, ambiguous ? 5.4 : 5.8],
            -3 => [1.0, 2.5, 4.5, 7.4],
            -4 => [1.0, 2.0, 3.0, 4.0, 5.0, 9.0],
        )
        reference_ky = mod(2 + offset, 6)
        raw_levels = NamedTuple[]
        for delta_ky in (0, -1, -2, -3, -4)
            ky_index = mod(reference_ky + delta_ky, 6)
            for relative_energy in energies[delta_ky]
                push!(raw_levels, (;
                    relative_energy,
                    raw_charge=4,
                    ky_index,
                ))
            end
        end
        push!(raw_levels, (;
            relative_energy=0.5,
            raw_charge=3,
            ky_index=mod(reference_ky - 1, 6),
        ))
        normalization = sum(exp(-row.relative_energy) for row in raw_levels)
        levels = MomentumEntanglementLevel[
            let probability = exp(-row.relative_energy) / normalization
                MomentumEntanglementLevel(
                    level,
                    probability,
                    -log(probability),
                    row.raw_charge,
                    row.ky_index,
                    2pi * row.ky_index / 6,
                    2pi * row.ky_index / 6,
                    0.0,
                )
            end for (level, row) in enumerate(raw_levels)
        ]
        validation = MomentumValidation(
            0.0, 0.0, 0.0, transfer_residual, true, "valid"
        )
        return MomentumEntanglementData(
            levels,
            validation,
            true,
            "valid",
            :sector_relative_canonical_cyclic_orbit,
            Dict(3 => offset, 4 => offset),
        )
    end

    function synthetic_mixed_fidelity(value)
        return MixedTransferFidelity(
            value,
            ComplexF64(value),
            1.0 + 0.0im,
            1.0 + 0.0im,
            1.0e-10,
            2.0e-10,
            3.0e-10,
            true,
            true,
            "valid normalized mixed-transfer fidelity",
        )
    end

    function synthetic_repository_commit()
        repository = normpath(joinpath(@__DIR__, "..", "..", ".."))
        dotgit = joinpath(repository, ".git")
        git_directory = if isdir(dotgit)
            dotgit
        else
            declaration = strip(read(dotgit, String))
            @assert startswith(declaration, "gitdir: ")
            target = strip(declaration[9:end])
            isabspath(target) ? target : normpath(joinpath(repository, target))
        end
        common_directory = if isfile(joinpath(git_directory, "commondir"))
            target = strip(read(joinpath(git_directory, "commondir"), String))
            isabspath(target) ? target :
                normpath(joinpath(git_directory, target))
        else
            git_directory
        end
        head = strip(read(joinpath(git_directory, "HEAD"), String))
        commit = if startswith(head, "ref: ")
            reference = strip(head[6:end])
            loose_paths = (
                joinpath(git_directory, split(reference, '/')...),
                joinpath(common_directory, split(reference, '/')...),
            )
            loose = findfirst(isfile, loose_paths)
            if !isnothing(loose)
                strip(read(loose_paths[loose], String))
            else
                packed = only(filter(
                    line -> endswith(line, " $reference"),
                    readlines(joinpath(common_directory, "packed-refs")),
                ))
                first(split(packed; limit=2))
            end
        else
            head
        end
        @assert occursin(r"^[0-9a-f]{40}$", commit)
        return commit
    end

    function synthetic_fig2_provenance(
        spec,
        output,
        runtime_seconds;
        git_commit=synthetic_repository_commit(),
        pbs_job_id=get(ENV, "PBS_JOBID", "synthetic.w003"),
    )
        active_project = abspath(Base.active_project())
        project_manifest = joinpath(dirname(active_project), "Manifest.toml")
        benchmark_source = joinpath(
            dirname(pathof(InfiniteCylinderDMRG)), "Fig2Benchmark.jl"
        )
        return Dict{String,Any}(
            "format" => "fqahc_fig2_provenance_v2",
            "manifest_sha256" => spec.sha256,
            "git_commit" => git_commit,
            "git_tree_clean" => true,
            "julia_version" => string(VERSION),
            "pbs_job_id" => pbs_job_id,
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
    end

    function synthetic_restart_gate(
        spec,
        directory;
        overrides=Dict{String,Any}(),
    )
        provenance = synthetic_fig2_provenance(spec, directory, 0.0)
        gate = Dict{String,Any}(
            "format" => "fqahc_fig2_restart_gate_v1",
            "pass" => true,
            "manifest_sha256" => spec.sha256,
            "git_commit" => provenance["git_commit"],
            "julia_version" => provenance["julia_version"],
            "project_manifest_sha256" =>
                provenance["project_manifest_sha256"],
            "backend_commit" =>
                "765f2777703bc1138b009adbed1b97bde1973402",
            "save_process_exit_code" => 0,
            "resume_process_exit_code" => 0,
            "save_process_id" => "save:101",
            "resume_process_id" => "resume:202",
            "load_completed" => true,
            "no_expansion_iteration_completed" => true,
            "controlled_expansion_completed" => true,
            "energy_reproduced" => true,
            "density_reproduced" => true,
            "schmidt_sectors_reproduced" => true,
            "entanglement_spectrum_reproduced" => true,
            "site_indices_reproduced" => true,
            "next_flux_iteration_completed" => true,
            "pre_maxlinkdim" => 4,
            "post_maxlinkdim" => 8,
            "energy_absolute_error" => 1.0e-12,
            "energy_tolerance" => 1.0e-10,
            "density_max_absolute_error" => 2.0e-12,
            "density_tolerance" => 1.0e-10,
        )
        merge!(gate, deepcopy(overrides))
        return gate
    end

    function write_synthetic_restart_gate!(
        spec,
        directory;
        overrides=Dict{String,Any}(),
    )
        path = joinpath(directory, "restart_gate.toml")
        InfiniteCylinderDMRG._write_fig2_toml(
            path, synthetic_restart_gate(spec, directory; overrides)
        )
        return path
    end

    function synthetic_checkpoint_audit(spec, checkpoint, phi_y)
        _ = (spec, phi_y)
        matched = match(
            r"^synthetic-checkpoint:maxlinkdim=([0-9]+):",
            read(checkpoint, String),
        )
        isnothing(matched) && return (
            restart_valid=false,
            checkpoint_maxlinkdim=0,
        )
        checkpoint_maxlinkdim = parse(Int, only(matched.captures))
        checkpoint_maxlinkdim > 0 || return (
            restart_valid=false,
            checkpoint_maxlinkdim=0,
        )
        return (
            restart_valid=true,
            checkpoint_maxlinkdim,
        )
    end

    function synthetic_progress_audit(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        final_checkpoint,
    )
        _ = (spec, dimension, point, phi_y, candidate_id, final_checkpoint)
        path = joinpath(candidate_directory, "progress.toml")
        data = TOML.parsefile(path)
        data["format"] == "fqahc_fig2_synthetic_progress_v1" || error(
            "invalid synthetic progress format"
        )
        return (;
            complete=data["complete"],
            event_count=Int(data["event_count"]),
            resume_count=Int(data["resume_count"]),
            latest_maxlinkdim=Int(data["latest_maxlinkdim"]),
            final_state_sha256=String(data["final_state_sha256"]),
            progress_sha256=InfiniteCylinderDMRG._fig2_file_sha256(path),
        )
    end

    function write_synthetic_fig2_acceptance_report!(
        spec,
        directory;
        candidate_ids_provider,
    )
        return write_fig2_acceptance_report!(
            spec,
            directory;
            checkpoint_audit=synthetic_checkpoint_audit,
            progress_audit=synthetic_progress_audit,
            candidate_ids_provider,
        )
    end

    function synthetic_fig2_acceptance_integrity_inputs(
        spec,
        directory;
        candidate_ids_provider,
    )
        return InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
            spec,
            directory;
            checkpoint_audit=synthetic_checkpoint_audit,
            progress_audit=synthetic_progress_audit,
            candidate_ids_provider,
        )
    end

    @testset "synthetic candidate contracts are fixed before ledger replay" begin
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        fixed_candidate_ids = (args...) -> ["warm", "cold_a", "cold_b"]
        @test InfiniteCylinderDMRG._validate_fig2_candidate_id_set(
            fixed_candidate_ids,
            spec,
            4,
            2,
            ["warm"];
            require_complete=false,
        ) == ["warm", "cold_a", "cold_b"]
        @test_throws ArgumentError InfiniteCylinderDMRG._validate_fig2_candidate_id_set(
            fixed_candidate_ids,
            spec,
            4,
            2,
            ["warm", "extra"];
            require_complete=false,
        )
        @test_throws ArgumentError InfiniteCylinderDMRG._validate_fig2_candidate_id_set(
            fixed_candidate_ids,
            spec,
            4,
            2,
            ["warm"];
            require_complete=true,
        )
    end

    @testset "Fig. 2 manifest snapshot rejects post-load runtime drift" begin
        @test TOML.parsefile(FIG2_MANIFEST_PATH)["format"] ==
            "fqahc_fig2_benchmark_v4"

        function drift_error(tamper::Function, action::Function)
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            tamper(spec.data)
            exception = try
                action(spec)
                nothing
            catch error
                error
            end
            @test exception isa ArgumentError
            message = isnothing(exception) ? "" : lowercase(sprint(showerror, exception))
            @test occursin("manifest", message)
            @test occursin("in-memory drift", message)
        end

        drift_error(
            data -> data["expected_endpoint_pump"] = 1.1,
            spec -> mktempdir() do directory
                operations = Fig2BenchmarkOperations(
                    progress_audit=synthetic_progress_audit,
                    provenance=synthetic_fig2_provenance,
                )
                run_fig2_benchmark(
                    spec,
                    directory;
                    stage="runtime_drift",
                    dimensions=Int[],
                    fluxes=Float64[],
                    operations,
                )
            end,
        )
        drift_error(
            data -> data["counting"][1] = 2,
            spec -> mktempdir() do directory
                write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=(args...) -> ["fixed"],
                )
            end,
        )
        drift_error(
            data -> data["paper_curve_rule"]["uncertainty"][1] = 0.2,
            spec -> InfiniteCylinderDMRG._default_fig2_candidate_ids(
                spec, 4, 1, nothing
            ),
        )
        drift_error(
            data -> data["optimization"]["max_iterations"] = 3,
            spec -> InfiniteCylinderDMRG._default_fig2_run_candidate(
                spec,
                2,
                1,
                0.0,
                "not-a-candidate",
                nothing,
                "unused",
            ),
        )
        drift_error(
            data -> data["pump_tolerance"] = 0.5,
            spec -> InfiniteCylinderDMRG._fig2_momentum_counting(
                spec, synthetic_counting_momentum()
            ),
        )
        drift_error(
            data -> begin
                plateau = get!(
                    data,
                    "early_plateau",
                    Dict{String,Any}(
                        "policy" => "allowed_if_all_gates_pass",
                        "lower_dimension" => 1000,
                        "upper_dimension" => 2000,
                        "omitted_dimension" => 3000,
                        "evidence_filename" => "early_plateau.toml",
                    ),
                )
                plateau["policy"] = "allowed_after_run"
            end,
            spec -> mktempdir() do directory
                InfiniteCylinderDMRG._fig2_early_plateau_evidence(
                    spec, directory
                )
            end,
        )

        @test isdefined(InfiniteCylinderDMRG, :Fig2RunnerOverride)
        @test isdefined(InfiniteCylinderDMRG, :_fig2_runner_optimization)
        if isdefined(InfiniteCylinderDMRG, :Fig2RunnerOverride) &&
           isdefined(InfiniteCylinderDMRG, :_fig2_runner_optimization)
            override = InfiniteCylinderDMRG.Fig2RunnerOverride(;
                cutoff=1.0e-8,
                vumps_tol=1.0e2,
                energy_tol=1.0e6,
                energy_mismatch_tol=1.0e6,
                transfer_tol=1.0e-4,
                momentum_residual_tol=1.0e-4,
                max_iterations=3,
                stable_iterations=1,
            )
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            withenv("IDMRG_FIG2_REAL_SMOKE" => "0") do
                @test_throws ArgumentError InfiniteCylinderDMRG._fig2_runner_optimization(
                    spec, override
                )
            end
        end
    end

    @testset "Fig. 2 momentum counting is sector-relative and uncertainty-gated" begin
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        counting = InfiniteCylinderDMRG._fig2_momentum_counting(
            spec, synthetic_counting_momentum()
        )
        @test counting.valid
        @test counting.counts == [1, 1, 2, 3, 5]
        @test counting.delta_ky == [0, -1, -2, -3, -4]
        @test counting.reference_raw_charge == 4
        @test counting.reference_ky_index == 2
        @test counting.minimum_cutoff_margin > 0.20

        shifted = InfiniteCylinderDMRG._fig2_momentum_counting(
            spec, synthetic_counting_momentum(; offset=3)
        )
        @test shifted.valid
        @test shifted.counts == counting.counts
        @test shifted.delta_ky == counting.delta_ky
        @test shifted.reference_ky_index == 5

        uncertain = InfiniteCylinderDMRG._fig2_momentum_counting(
            spec, synthetic_counting_momentum(; ambiguous=true)
        )
        @test !uncertain.valid
        @test isempty(uncertain.counts)
        @test occursin("digitization uncertainty", uncertain.reason)
    end

    @testset "Fig. 2 branch selection gates restart and sector correspondence" begin
        candidate(
            id,
            fidelity,
            sectors;
            restart_valid=true,
            momentum_valid=false,
            provenance_valid=true,
            requested_maxdim=4,
            achieved_maxlinkdim=requested_maxdim,
            checkpoint_maxlinkdim=achieved_maxlinkdim,
            converged=true,
            energy_per_site=-1.0,
        ) = (;
            candidate_id=id,
            sectors,
            provenance_valid,
            authoritative_energy_per_site=energy_per_site,
            evidence=(;
                requested_maxdim,
                achieved_maxlinkdim,
                checkpoint_maxlinkdim,
                converged,
                valid=true,
                energy_per_site,
                fidelity_valid=true,
                fidelity_to_previous=fidelity,
                restart_valid,
                momentum=(; valid=momentum_valid),
            ),
        )
        previous_sectors = Dict(0 => 0.7, 1 => 0.3)
        candidates = [
            candidate("provenance_bad", 1.0, Dict(1 => 0.7, 2 => 0.3);
                provenance_valid=false, momentum_valid=true),
            candidate("restart_bad", 0.999, Dict(1 => 0.7, 2 => 0.3);
                restart_valid=false, momentum_valid=true),
            candidate("checkpoint_bad", 0.9995, Dict(1 => 0.7, 2 => 0.3);
                checkpoint_maxlinkdim=2, momentum_valid=true),
            candidate("high_fidelity_large_sector_residual", 0.998,
                Dict(1 => 0.5, 2 => 0.5);
                momentum_valid=true),
            candidate("raw_momentum_invalid", 0.95, Dict(1 => 0.7, 2 => 0.3)),
        ]
        selected, basis = InfiniteCylinderDMRG._select_fig2_candidate(
            candidates,
            2;
            previous_sector_weights=previous_sectors,
            sector_tolerance=1e-6,
        )
        @test candidates[selected].candidate_id ==
            "high_fidelity_large_sector_residual"
        @test basis == "maximum_mixed_fidelity"
        @test candidates[selected].evidence.momentum.valid

        tied_candidates = [
            candidate("lower_energy_large_residual", 0.9,
                Dict(1 => 0.5, 2 => 0.5); energy_per_site=-2.0),
            candidate("higher_energy_small_residual", 0.9,
                Dict(1 => 0.7, 2 => 0.3); energy_per_site=-1.0),
        ]
        tied_selected, tied_basis =
            InfiniteCylinderDMRG._select_fig2_candidate(
                tied_candidates,
                2;
                previous_sector_weights=previous_sectors,
                sector_tolerance=1e-6,
            )
        @test tied_candidates[tied_selected].candidate_id ==
            "higher_energy_small_residual"
        @test tied_basis == "maximum_mixed_fidelity"

        invalid_candidates = [
            candidate("fidelity_bad", NaN, Dict(1 => 0.7, 2 => 0.3)),
            candidate("restart_bad", 0.999, Dict(1 => 0.7, 2 => 0.3);
                restart_valid=false),
            candidate("checkpoint_bad", 0.9995, Dict(1 => 0.7, 2 => 0.3);
                checkpoint_maxlinkdim=2),
        ]
        @test_throws WorkflowValidationError InfiniteCylinderDMRG._select_fig2_candidate(
            invalid_candidates,
            2;
            previous_sector_weights=previous_sectors,
            sector_tolerance=1e-6,
        )

        cap_candidates = [
            candidate(
                "underachieved_low_energy",
                NaN,
                Dict(0 => 1.0);
                achieved_maxlinkdim=2,
                energy_per_site=-2.0,
            ),
            candidate(
                "at_cap",
                NaN,
                Dict(0 => 1.0);
                energy_per_site=-1.0,
            ),
        ]
        cap_selected, cap_basis = InfiniteCylinderDMRG._select_fig2_candidate(
            cap_candidates, 1
        )
        @test cap_candidates[cap_selected].candidate_id == "at_cap"
        @test cap_basis == "minimum_zero_flux_energy"
    end

    @testset "Fig. 2 selection core ignores optional single-point failures" begin
        helper_name = :_fig2_selection_core_validity
        @test isdefined(InfiniteCylinderDMRG, helper_name)
        if isdefined(InfiniteCylinderDMRG, helper_name)
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            expected_pair = (1, spec.config.Ly)
            fixed_cut = (;
                cut_x=expected_pair[1],
                bond=expected_pair[2],
                entropy=0.5,
                raw_schmidt_polarization=0.25,
                levels=[(;
                    cut_x=expected_pair[1],
                    bond=expected_pair[2],
                    singular_value=1.0,
                    probability=1.0,
                    entanglement_energy=0.0,
                    raw_charge=0,
                    physical_charge=0.0,
                )],
                sectors=[(;
                    cut_x=expected_pair[1],
                    bond=expected_pair[2],
                    raw_charge=0,
                    physical_charge=0.0,
                    weight=1.0,
                )],
            )
            optional_cut = (;
                cut_x=2,
                bond=2 * spec.config.Ly,
                entropy=NaN,
                raw_schmidt_polarization=NaN,
                levels=NamedTuple[],
                sectors=NamedTuple[],
            )
            optional_reasons = [
                "neutral transfer failed: synthetic optional diagnostic",
                "entanglement cut 2 failed: synthetic optional diagnostic",
            ]
            single_point = (;
                optimization=(; converged=true),
                energy=(; per_site=-1.0),
                entanglements=[fixed_cut, optional_cut],
                valid=false,
                observable_reasons=optional_reasons,
            )
            core_validity = getproperty(InfiniteCylinderDMRG, helper_name)
            @test !single_point.valid
            @test core_validity(single_point, expected_pair)
            @test single_point.observable_reasons == optional_reasons
            @test !core_validity(
                merge(single_point, (; optimization=(; converged=false))),
                expected_pair,
            )
            @test !core_validity(
                merge(single_point, (; energy=(; per_site=NaN))),
                expected_pair,
            )
            @test !core_validity(
                merge(single_point, (; entanglements=[optional_cut])),
                expected_pair,
            )
        end
    end

    const FIG2_SYNTHETIC_CHARGE_SCALE = charge_scale(
        load_fig2_benchmark(FIG2_MANIFEST_PATH).config
    )
    function synthetic_raw_schmidt_polarization(sector_weights)
        normalized = InfiniteCylinderDMRG._normalized_sector_weights(
            sector_weights
        )
        return sum(
            charge * weight / FIG2_SYNTHETIC_CHARGE_SCALE for
            (charge, weight) in normalized
        )
    end
    function synthetic_sector_weights_for_polarization(polarization)
        value = Float64(polarization)
        @assert isfinite(value)
        lower = floor(Int, value)
        upper = ceil(Int, value)
        lower == upper && return Dict(
            FIG2_SYNTHETIC_CHARGE_SCALE * lower => 1.0
        )
        upper_weight = value - lower
        return Dict(
            FIG2_SYNTHETIC_CHARGE_SCALE * lower => 1 - upper_weight,
            FIG2_SYNTHETIC_CHARGE_SCALE * upper => upper_weight,
        )
    end

    function synthetic_fig2_convergence_rows(
        schedule;
        final_maxlinkdim=last(schedule),
        final_converged=true,
    )
        rows = NamedTuple[]
        for (stage, stage_maxdim) in enumerate(schedule)
            is_final = stage == length(schedule)
            maxlinkdim = is_final ? final_maxlinkdim : stage_maxdim
            push!(rows, (
                stage,
                iteration=1,
                maxlinkdim,
                delta_energy="missing",
                converged=false,
            ))
            if !is_final || final_converged
                push!(rows, (
                    stage,
                    iteration=2,
                    maxlinkdim,
                    delta_energy=1.0e-8,
                    converged=false,
                ))
                push!(rows, (
                    stage,
                    iteration=3,
                    maxlinkdim,
                    delta_energy=1.0e-8,
                    converged=true,
                ))
            end
        end
        return rows
    end

    function write_fake_candidate_files(
        directory,
        marker;
        requested_maxdim,
        achieved_maxlinkdim=requested_maxdim,
        checkpoint_maxlinkdim=achieved_maxlinkdim,
        summary_maxdim_schedule=
            InfiniteCylinderDMRG._fig2_maxdim_schedule(requested_maxdim),
        summary_converged=true,
        summary_valid=true,
        fixed_cut_valid=summary_converged,
        summary_energy_per_site=-1.0,
        convergence_energy_per_site=summary_energy_per_site,
        summary_optimization_overrides=Dict{String,Any}(),
        sector_weights=Dict(0 => 1.0),
        convergence_rows=synthetic_fig2_convergence_rows(
            summary_maxdim_schedule;
            final_maxlinkdim=achieved_maxlinkdim,
            final_converged=summary_converged,
        ),
    )
        mkpath(directory)
        open(joinpath(directory, "state.h5"), "w") do io
            write(
                io,
                "synthetic-checkpoint:maxlinkdim=$checkpoint_maxlinkdim:$marker\n",
            )
        end
        synthetic_state_path = joinpath(directory, "state.h5")
        open(joinpath(directory, "progress.toml"), "w") do io
            TOML.print(io, Dict(
                "format" => "fqahc_fig2_synthetic_progress_v1",
                "complete" => true,
                "event_count" => 3,
                "resume_count" => 0,
                "latest_maxlinkdim" => achieved_maxlinkdim,
                "final_state_sha256" =>
                    InfiniteCylinderDMRG._fig2_file_sha256(
                        synthetic_state_path
                    ),
            ); sorted=true)
        end
        summary_optimization = Dict{String,Any}(
            "maxdim_schedule" => summary_maxdim_schedule,
            "multisite_update_alg" => "sequential",
            "vumps_tol" => 1.0e-6,
            "energy_tol" => 1.0e-6,
            "energy_mismatch_tol" => 1.0e-6,
            "stable_iterations" => 2,
        )
        merge!(summary_optimization, summary_optimization_overrides)
        open(joinpath(directory, "summary.toml"), "w") do io
            TOML.print(io, Dict(
                "format" => "infinite_cylinder_vumps_text_v1",
                "algorithm" => "VUMPS",
                "valid" => summary_valid,
                "converged" => summary_converged,
                "energy" => Dict(
                    "per_site" => summary_energy_per_site,
                ),
                "optimization" => summary_optimization,
            ); sorted=true)
        end
        convergence = IOBuffer()
        println(
            convergence,
            "stage\titeration\tmaxlinkdim\tenergy_left\tenergy_right\tenergy_mismatch\tdelta_energy\teps_left\teps_right\tprecision_error\telapsed_seconds\tconverged",
        )
        default_means = zeros(Float64, length(convergence_rows))
        for row_index in 2:length(convergence_rows)
            row = convergence_rows[row_index]
            delta = hasproperty(row, :delta_energy) ? row.delta_energy : 0.0
            default_means[row_index] = default_means[row_index - 1] +
                (delta isa Real && isfinite(delta) ? abs(Float64(delta)) : 0.0)
        end
        if !isempty(default_means)
            default_means .+=
                convergence_energy_per_site * FIG2_SYNTHETIC_SITES_PER_CELL -
                last(default_means)
        end
        previous_written_mean = nothing
        for (row_index, row) in enumerate(convergence_rows)
            explicit_raw = hasproperty(row, :energy_left) &&
                hasproperty(row, :energy_right)
            requested_mismatch = hasproperty(row, :energy_mismatch) ?
                row.energy_mismatch : 0.0
            energy_left, energy_right = if explicit_raw
                (row.energy_left, row.energy_right)
            else
                raw_mismatch = requested_mismatch isa Real &&
                    isfinite(requested_mismatch) ?
                    abs(Float64(requested_mismatch)) : 0.0
                (
                    default_means[row_index] + raw_mismatch / 2,
                    default_means[row_index] - raw_mismatch / 2,
                )
            end
            written_mean = (energy_left + energy_right) / 2
            energy_mismatch = explicit_raw ? requested_mismatch :
                abs(energy_left - energy_right)
            requested_delta = hasproperty(row, :delta_energy) ?
                row.delta_energy : 0.0
            delta_energy = if row_index == 1 || explicit_raw ||
                    !(requested_delta isa Real)
                requested_delta
            else
                abs(written_mean - previous_written_mean)
            end
            println(
                convergence,
                join((
                    row.stage,
                    row.iteration,
                    row.maxlinkdim,
                    energy_left,
                    energy_right,
                    energy_mismatch,
                    delta_energy,
                    hasproperty(row, :eps_left) ? row.eps_left : 1.0e-8,
                    hasproperty(row, :eps_right) ? row.eps_right : 1.0e-8,
                    hasproperty(row, :precision_error) ?
                        row.precision_error : 1.0e-8,
                    hasproperty(row, :elapsed_seconds) ?
                        row.elapsed_seconds : 0.0,
                    row.converged,
                ), '\t'),
            )
            previous_written_mean = written_mean
        end
        normalized_sectors = InfiniteCylinderDMRG._normalized_sector_weights(
            sector_weights
        )
        entanglement_spectrum = IOBuffer()
        println(
            entanglement_spectrum,
            "cut_x\tbond\tlevel\tsingular_value\tprobability\tentanglement_energy\tqn\traw_charge\tphysical_charge\tvalid",
        )
        schmidt_sectors = IOBuffer()
        println(
            schmidt_sectors,
            "cut_x\tbond\tqn\traw_charge\tphysical_charge\tweight\tvalid",
        )
        for (level, charge) in enumerate(sort!(collect(keys(normalized_sectors))))
            probability = normalized_sectors[charge]
            singular_value = sqrt(probability)
            println(
                entanglement_spectrum,
                join((
                    1,
                    FIG2_SYNTHETIC_LY,
                    level,
                    singular_value,
                    singular_value^2,
                    -log(probability),
                    "Nf($charge)",
                    charge,
                    Float64(charge) / FIG2_SYNTHETIC_CHARGE_SCALE,
                    fixed_cut_valid,
                ), '\t'),
            )
            println(
                schmidt_sectors,
                join((
                    1,
                    FIG2_SYNTHETIC_LY,
                    "Nf($charge)",
                    charge,
                    Float64(charge) / FIG2_SYNTHETIC_CHARGE_SCALE,
                    probability,
                    fixed_cut_valid,
                ), '\t'),
            )
        end
        tables = Dict(
            "convergence.tsv" => String(take!(convergence)),
            "expansion.tsv" =>
                "stage\ttarget\tbefore\tafter\tprogressed\telapsed_seconds\n",
            "density.tsv" =>
                "site\tx\ty\tdensity\tvalid\n1\t0\t0\t0.5\ttrue\n",
            "entanglement_spectrum.tsv" =>
                String(take!(entanglement_spectrum)),
            "schmidt_sectors.tsv" =>
                String(take!(schmidt_sectors)),
        )
        for (filename, contents) in tables
            open(joinpath(directory, filename), "w") do io
                write(io, contents)
            end
        end
        return nothing
    end

    @testset "optional observables do not invalidate Fig. 2 selection core" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            optional_reasons = [
                "neutral transfer failed: synthetic optional diagnostic",
                "entanglement cut 2 failed: synthetic optional diagnostic",
            ]
            candidate_calls = String[]
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(args...) -> ["only"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    push!(candidate_calls, candidate_id)
                    write_fake_candidate_files(
                        candidate_directory,
                        "optional-observable-failure";
                        requested_maxdim=dimension,
                        summary_valid=false,
                        sector_weights=Dict(0 => 1.0),
                    )
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; candidate_id),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=-1.0,
                        raw_schmidt_polarization=0.0,
                        sector_weights=Dict(0 => 1.0),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                        reason=join(optional_reasons, "; "),
                    )
                end,
            )

            run = nothing
            run_error = try
                run = run_fig2_benchmark(
                    spec,
                    directory;
                    stage="optional_observable_core",
                    dimensions=[4],
                    fluxes=[0.0],
                    operations,
                )
                nothing
            catch error
                error
            end
            @test isnothing(run_error)
            if isnothing(run_error)
                @test candidate_calls == ["only"]
                @test only(run.selections).candidate_id == "only"
                ledger_path = joinpath(directory, "ledger.toml")
                ledger = TOML.parsefile(ledger_path)
                candidate_row = only(ledger["candidate"])
                candidate_path = joinpath(
                    directory, candidate_row["directory"], "candidate.toml"
                )
                summary_path = joinpath(
                    directory, candidate_row["directory"], "summary.toml"
                )
                metadata = TOML.parsefile(candidate_path)
                summary = TOML.parsefile(summary_path)
                @test summary["valid"] === false
                @test metadata["valid"] === true
                @test all(
                    reason -> occursin(reason, metadata["reason"]),
                    optional_reasons,
                )
                integrity = synthetic_fig2_acceptance_integrity_inputs(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                @test only(integrity.selections).candidate_id == "only"

                metadata["valid"] = false
                open(candidate_path, "w") do io
                    TOML.print(io, metadata; sorted=true)
                end
                candidate_row["checksums"]["candidate.toml"] =
                    InfiniteCylinderDMRG._fig2_file_sha256(candidate_path)
                open(ledger_path, "w") do io
                    TOML.print(io, ledger; sorted=true)
                end
                replay_error = try
                    synthetic_fig2_acceptance_integrity_inputs(
                        spec,
                        directory;
                        candidate_ids_provider=operations.candidate_ids,
                    )
                    nothing
                catch error
                    error
                end
                @test replay_error isa ArgumentError
                replay_message = isnothing(replay_error) ? "" :
                    lowercase(sprint(showerror, replay_error))
                @test occursin("selection-core validity", replay_message)
            end
        end
    end

    @testset "sector path uses physical endpoint translation" begin
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        sector_row(sectors; fidelity_valid=true, fidelity=0.99) = (;
            sector_weights=sectors,
            fidelity_valid,
            fidelity_to_previous=fidelity,
        )

        endpoint = InfiniteCylinderDMRG._fig2_sector_path_result(spec, [
            sector_row(Dict(0 => 1.0); fidelity_valid=false, fidelity=NaN),
            sector_row(Dict(18 => 1.0)),
        ])
        @test endpoint.valid
        @test endpoint.shift == 1
        @test get(endpoint, :raw_shift, nothing) == 18

        nondivisible = InfiniteCylinderDMRG._fig2_sector_path_result(spec, [
            sector_row(Dict(0 => 1.0); fidelity_valid=false, fidelity=NaN),
            sector_row(Dict(1 => 1.0)),
        ])
        @test !nondivisible.valid
        @test get(nondivisible, :raw_shift, nothing) == 1

        interpolation_rows = [
            (;
                raw=(phi_y=0.0, raw_schmidt_polarization=0.0),
                tracked=(valid=true, sector_distance_raw=0.2),
            ),
            (;
                raw=(phi_y=1.0, raw_schmidt_polarization=0.2),
                tracked=(valid=true, sector_distance_raw=0.2),
            ),
        ]
        bracket = try
            InfiniteCylinderDMRG._fig2_interpolation_bracket(
                interpolation_rows, 0.0, 1.0e-6
            )
        catch error
            error
        end
        @test bracket isa NamedTuple
    end

    @testset "persisted candidate replay derives selection evidence from raw tables" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(args...) -> ["only"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    first_point = point == 1
                    sectors = first_point ? Dict(0 => 1.0) : Dict(1 => 1.0)
                    write_fake_candidate_files(
                        candidate_directory,
                        "raw-replay-point=$point";
                        requested_maxdim=dimension,
                        summary_energy_per_site=-1.0,
                        sector_weights=sectors,
                    )
                    fidelity = first_point ? NaN : 0.99
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; point),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=-1.0,
                        raw_schmidt_polarization=
                            synthetic_raw_schmidt_polarization(sectors),
                        sector_weights=sectors,
                        fidelity_to_previous=fidelity,
                        fidelity_valid=!first_point,
                        mixed_fidelity=first_point ? nothing :
                            synthetic_mixed_fidelity(fidelity),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                    )
                end,
            )
            run_fig2_benchmark(
                spec,
                directory;
                stage="raw_replay_tdd",
                dimensions=[4],
                fluxes=[0.0, 6pi],
                operations,
            )
            ledger = TOML.parsefile(joinpath(directory, "ledger.toml"))
            candidates = Dict(
                Int(row["point"]) => row for row in ledger["candidate"]
            )
            audit_message(error) = isnothing(error) ? "" :
                lowercase(sprint(showerror, error))

            function artifact_audit_error(tamper::Function, point)
                row = candidates[point]
                return mktempdir() do scratch
                    candidate_directory = joinpath(scratch, "candidate")
                    cp(
                        joinpath(directory, row["directory"]),
                        candidate_directory;
                        force=true,
                    )
                    tamper(candidate_directory)
                    metadata = TOML.parsefile(joinpath(
                        candidate_directory, "candidate.toml"
                    ))
                    try
                        InfiniteCylinderDMRG._validate_fig2_candidate_artifacts(
                            spec,
                            candidate_directory,
                            4,
                            point,
                            Float64(row["phi_y"]),
                            "only",
                            metadata["generation_provenance"],
                            ;
                            progress_audit=synthetic_progress_audit,
                        )
                        nothing
                    catch error
                        error
                    end
                end
            end

            function rewrite_first_momentum_row!(tamper::Function, directory)
                path = joinpath(
                    directory, "momentum_entanglement_spectrum.tsv"
                )
                lines = readlines(path)
                fields = split(lines[2], '\t'; keepempty=true)
                tamper(fields)
                lines[2] = join(fields, '\t')
                open(path, "w") do io
                    foreach(line -> println(io, line), lines)
                end
            end

            missing_ky = artifact_audit_error(1) do candidate_directory
                rewrite_first_momentum_row!(candidate_directory) do fields
                    fields[5] = "missing"
                    fields[6] = "missing"
                    fields[7] = "missing"
                    fields[8] = "NaN"
                end
            end
            @test missing_ky isa ArgumentError
            @test occursin("momentum", audit_message(missing_ky))

            malformed_ky = artifact_audit_error(1) do candidate_directory
                rewrite_first_momentum_row!(candidate_directory) do fields
                    fields[5] = "not-an-integer"
                end
            end
            @test malformed_ky isa ArgumentError
            @test occursin("ky_index", audit_message(malformed_ky))

            excessive_phase_residual = artifact_audit_error(1) do candidate_directory
                rewrite_first_momentum_row!(candidate_directory) do fields
                    raw_charge = parse(Int, fields[4])
                    ky_index = parse(Int, fields[5])
                    offsets = InfiniteCylinderDMRG._fig2_momentum_offsets(
                        fields[18]
                    )
                    unshifted_index = mod(
                        ky_index - offsets[raw_charge], spec.config.Ny
                    )
                    expected_phase = 2pi * unshifted_index / spec.config.Ny
                    forged_phase = expected_phase + 1.0e-3
                    fields[7] = string(forged_phase)
                    fields[8] = string(abs(
                        cis(forged_phase) - cis(expected_phase)
                    ))
                end
            end
            @test excessive_phase_residual isa ArgumentError
            @test occursin(
                "phase residual", audit_message(excessive_phase_residual)
            )

            counting_mismatch = artifact_audit_error(1) do candidate_directory
                path = joinpath(
                    candidate_directory,
                    "momentum_entanglement_spectrum.tsv",
                )
                lines = readlines(path)
                fields = split(lines[4], '\t'; keepempty=true)
                fields[5] = "2"
                fields[6] = string(2pi * 2 / 6)
                fields[7] = string(2pi * 2 / 6)
                lines[4] = join(fields, '\t')
                open(path, "w") do io
                    foreach(line -> println(io, line), lines)
                end
            end
            @test counting_mismatch isa ArgumentError
            @test occursin("counting", audit_message(counting_mismatch))

            @test isnothing(artifact_audit_error(1) do candidate_directory
                _ = candidate_directory
            end)
            mixed_residual_forgery = artifact_audit_error(2) do candidate_directory
                metadata_path = joinpath(candidate_directory, "candidate.toml")
                metadata = TOML.parsefile(metadata_path)
                metadata["mixed_residual"] = 1.0e-4
                open(metadata_path, "w") do io
                    TOML.print(io, metadata; sorted=true)
                end
                table_path = joinpath(candidate_directory, "mixed_fidelity.tsv")
                lines = readlines(table_path)
                fields = split(lines[2], '\t'; keepempty=true)
                fields[8] = "0.0001"
                lines[2] = join(fields, '\t')
                open(table_path, "w") do io
                    foreach(line -> println(io, line), lines)
                end
            end
            @test mixed_residual_forgery isa ArgumentError
            @test occursin(
                "mixed", audit_message(mixed_residual_forgery)
            )

            energy_mismatch = artifact_audit_error(1) do candidate_directory
                summary_path = joinpath(candidate_directory, "summary.toml")
                summary = TOML.parsefile(summary_path)
                summary["energy"]["per_site"] = -0.5
                open(summary_path, "w") do io
                    TOML.print(io, summary; sorted=true)
                end
            end
            @test energy_mismatch isa ArgumentError
            @test occursin("energy", audit_message(energy_mismatch))
        end
    end

    function validate_production_convergence(path, schedule)
        return InfiniteCylinderDMRG._validate_fig2_convergence_tsv(
            path,
            schedule;
            optimization=(;
                vumps_tol=1.0e-6,
                energy_tol=1.0e-6,
                energy_mismatch_tol=1.0e-6,
                stable_iterations=2,
            ),
            energy_normalization_sites=36,
        )
    end

    @testset "convergence full-table audit accepts a valid stage schedule" begin
        mktempdir() do directory
            write_fake_candidate_files(
                directory,
                "valid-convergence-schedule";
                requested_maxdim=4,
                summary_maxdim_schedule=[4],
                sector_weights=Dict(0 => 1.0),
                convergence_rows=[
                    (
                        stage=1,
                        iteration=1,
                        maxlinkdim=2,
                        delta_energy="missing",
                        converged=false,
                    ),
                    (
                        stage=1,
                        iteration=2,
                        maxlinkdim=4,
                        delta_energy="missing",
                        converged=false,
                    ),
                    (
                        stage=1,
                        iteration=3,
                        maxlinkdim=4,
                        delta_energy=1.0e-8,
                        converged=false,
                    ),
                    (
                        stage=1,
                        iteration=4,
                        maxlinkdim=4,
                        delta_energy=1.0e-8,
                        converged=true,
                    ),
                ],
            )
            outcome = try
                validate_production_convergence(
                    joinpath(directory, "convergence.tsv"), [4]
                )
            catch error
                error
            end
            @test !(outcome isa Exception)
            if !(outcome isa Exception)
                @test outcome.stage == 1
                @test outcome.iteration == 4
                @test outcome.maxlinkdim == 4
                @test outcome.converged
            end
        end
    end

    @testset "convergence replay enforces the immutable optimization contract" begin
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        optimization_data = spec.data["optimization"]
        optimization = (
            vumps_tol=Float64(optimization_data["vumps_tol"]),
            energy_tol=Float64(optimization_data["energy_tol"]),
            energy_mismatch_tol=Float64(
                optimization_data["energy_mismatch_tol"]
            ),
            stable_iterations=Int(optimization_data["stable_iterations"]),
        )
        energy_normalization_sites = sites_per_cell(spec.config)

        function production_audit(path, schedule)
            return InfiniteCylinderDMRG._validate_fig2_convergence_tsv(
                path,
                schedule;
                optimization,
                energy_normalization_sites,
            )
        end

        function audit_outcome(rows; schedule=[4])
            return mktempdir() do directory
                write_fake_candidate_files(
                    directory,
                    "convergence-replay";
                    requested_maxdim=last(schedule),
                    summary_maxdim_schedule=schedule,
                    sector_weights=Dict(0 => 1.0),
                    convergence_rows=rows,
                )
                try
                    production_audit(
                        joinpath(directory, "convergence.tsv"), schedule
                    )
                catch error
                    error
                end
            end
        end

        stable_row(iteration; converged=false, kwargs...) = merge(
            (
                stage=1,
                iteration,
                maxlinkdim=4,
                delta_energy=1.0e-8,
                energy_mismatch=1.0e-8,
                precision_error=1.0e-8,
                converged,
            ),
            (; kwargs...),
        )
        first_row = (
            stage=1,
            iteration=1,
            maxlinkdim=2,
            delta_energy="missing",
            energy_mismatch=1.0e-8,
            precision_error=1.0e-8,
            converged=false,
        )
        activation_row = merge(
            stable_row(2), (; delta_energy="missing")
        )

        bad_delta = audit_outcome([
            first_row,
            stable_row(2),
            stable_row(3; delta_energy=0.25, converged=true),
        ])
        @test bad_delta isa ArgumentError

        bad_mismatch = audit_outcome([
            first_row,
            stable_row(2),
            stable_row(3; energy_mismatch=0.25, converged=true),
        ])
        @test bad_mismatch isa ArgumentError

        bad_precision = audit_outcome([
            first_row,
            stable_row(2),
            stable_row(3; precision_error=0.25, converged=true),
        ])
        @test bad_precision isa ArgumentError

        premature = audit_outcome([
            first_row,
            stable_row(2; converged=true),
        ])
        @test premature isa ArgumentError

        low_dimension_converged = audit_outcome([
            first_row,
            merge(stable_row(2), (; maxlinkdim=2)),
            merge(stable_row(3; converged=true), (; maxlinkdim=2)),
        ])
        @test low_dimension_converged isa ArgumentError

        late_reset = audit_outcome([
            first_row,
            activation_row,
            merge(stable_row(3), (; delta_energy="missing")),
            stable_row(4),
            stable_row(5; converged=true),
        ])
        @test late_reset isa ArgumentError

        unconverged_stage = audit_outcome(
            [
                merge(first_row, (; maxlinkdim=2)),
                merge(stable_row(2), (; maxlinkdim=2)),
                (
                    stage=2,
                    iteration=1,
                    maxlinkdim=4,
                    delta_energy=1.0e-8,
                    energy_mismatch=1.0e-8,
                    precision_error=1.0e-8,
                    converged=false,
                ),
                (
                    stage=2,
                    iteration=2,
                    maxlinkdim=4,
                    delta_energy=1.0e-8,
                    energy_mismatch=1.0e-8,
                    precision_error=1.0e-8,
                    converged=true,
                ),
            ];
            schedule=[2, 4],
        )
        @test unconverged_stage isa ArgumentError

        valid = audit_outcome([
            first_row,
            activation_row,
            stable_row(3),
            stable_row(4; converged=true),
        ])
        @test !(valid isa Exception)
        if !(valid isa Exception)
            @test valid.stage == 1
            @test valid.iteration == 4
            @test valid.converged
        end

        multi_stage = audit_outcome(
            [
                (
                    stage=1,
                    iteration=1,
                    maxlinkdim=2,
                    delta_energy="missing",
                    energy_mismatch=1.0e-8,
                    precision_error=1.0e-8,
                    converged=false,
                ),
                (
                    stage=1,
                    iteration=2,
                    maxlinkdim=2,
                    delta_energy=1.0e-8,
                    energy_mismatch=1.0e-8,
                    precision_error=1.0e-8,
                    converged=false,
                ),
                (
                    stage=1,
                    iteration=3,
                    maxlinkdim=2,
                    delta_energy=1.0e-8,
                    energy_mismatch=1.0e-8,
                    precision_error=1.0e-8,
                    converged=true,
                ),
                (
                    stage=2,
                    iteration=1,
                    maxlinkdim=4,
                    delta_energy="missing",
                    energy_mismatch=1.0e-8,
                    precision_error=1.0e-8,
                    converged=false,
                ),
                (
                    stage=2,
                    iteration=2,
                    maxlinkdim=4,
                    delta_energy=1.0e-8,
                    energy_mismatch=1.0e-8,
                    precision_error=1.0e-8,
                    converged=false,
                ),
                (
                    stage=2,
                    iteration=3,
                    maxlinkdim=4,
                    delta_energy=1.0e-8,
                    energy_mismatch=1.0e-8,
                    precision_error=1.0e-8,
                    converged=true,
                ),
            ];
            schedule=[2, 4],
        )
        @test !(multi_stage isa Exception)
        if !(multi_stage isa Exception)
            @test multi_stage.stage == 2
            @test multi_stage.iteration == 3
            @test multi_stage.maxlinkdim == 4
            @test multi_stage.converged
        end
    end

    @testset "convergence replay derives every diagnostic from raw columns" begin
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        optimization_data = spec.data["optimization"]
        optimization = (
            vumps_tol=Float64(optimization_data["vumps_tol"]),
            energy_tol=Float64(optimization_data["energy_tol"]),
            energy_mismatch_tol=Float64(
                optimization_data["energy_mismatch_tol"]
            ),
            stable_iterations=Int(optimization_data["stable_iterations"]),
        )
        energy_step = 2.0^-30
        mismatch_step = 2.0^-31
        eps_left = 2.0^-31
        eps_right = 2.0^-30
        target_energy = -1.0 * FIG2_SYNTHETIC_SITES_PER_CELL
        means = [target_energy - 2energy_step,
                 target_energy - energy_step,
                 target_energy]
        raw_energies = [(
            left=mean + mismatch_step / 2,
            right=mean - mismatch_step / 2,
        ) for mean in means]
        replayed_means = [
            (energy.left + energy.right) / 2 for energy in raw_energies
        ]
        baseline_rows = [
            (
                stage=1,
                iteration,
                maxlinkdim=4,
                energy_left=raw_energies[iteration].left,
                energy_right=raw_energies[iteration].right,
                energy_mismatch=abs(
                    raw_energies[iteration].left -
                        raw_energies[iteration].right
                ),
                delta_energy=iteration == 1 ? "missing" : abs(
                    replayed_means[iteration] -
                        replayed_means[iteration - 1]
                ),
                eps_left,
                eps_right,
                precision_error=max(eps_left, eps_right),
                converged=iteration == 3,
            ) for iteration in 1:3
        ]

        function audit_derived_rows(rows)
            return mktempdir() do directory
                write_fake_candidate_files(
                    directory,
                    "derived-convergence";
                    requested_maxdim=4,
                    summary_energy_per_site=-1.0,
                    sector_weights=Dict(0 => 1.0),
                    convergence_rows=rows,
                )
                try
                    InfiniteCylinderDMRG._validate_fig2_convergence_tsv(
                        joinpath(directory, "convergence.tsv"),
                        [4];
                        optimization,
                        energy_normalization_sites=
                            sites_per_cell(spec.config),
                    )
                catch error
                    error
                end
            end
        end

        baseline = audit_derived_rows(baseline_rows)
        @test !(baseline isa Exception)
        tamper_cases = (
            "small false energy mismatch" => rows ->
                (rows[3] = merge(rows[3], (; energy_mismatch=0.0))),
            "small false delta energy" => rows ->
                (rows[3] = merge(rows[3], (; delta_energy=0.0))),
            "small false precision error" => rows ->
                (rows[3] = merge(rows[3], (; precision_error=0.0))),
            "negative energy mismatch" => rows ->
                (rows[3] = merge(rows[3], (; energy_mismatch=-energy_step))),
            "negative delta energy" => rows ->
                (rows[3] = merge(rows[3], (; delta_energy=-energy_step))),
            "negative precision error" => rows ->
                (rows[3] = merge(rows[3], (; precision_error=-energy_step))),
            "negative eps" => rows ->
                (rows[3] = merge(rows[3], (; eps_left=-energy_step))),
        )
        for (label, tamper) in tamper_cases
            @testset "$label" begin
                rows = copy(baseline_rows)
                tamper(rows)
                @test audit_derived_rows(rows) isa ArgumentError
            end
        end
    end

    @testset "candidate ranking uses replayed convergence energy" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            raw_energies = Dict(
                "alpha" => -1.0,
                "beta" => -1.0000005,
            )
            reported_energies = merge(
                copy(raw_energies), Dict("alpha" => -1.00000075)
            )
            candidate_ids = (args...) -> ["alpha", "beta"]
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=candidate_ids,
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    energy = reported_energies[candidate_id]
                    sectors = Dict(0 => 1.0)
                    write_fake_candidate_files(
                        candidate_directory,
                        "near-degenerate-$candidate_id";
                        requested_maxdim=dimension,
                        summary_energy_per_site=energy,
                        convergence_energy_per_site=
                            raw_energies[candidate_id],
                        sector_weights=sectors,
                    )
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; candidate_id),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=energy,
                        raw_schmidt_polarization=0.0,
                        sector_weights=sectors,
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                    )
                end,
            )
            run = run_fig2_benchmark(
                spec,
                directory;
                stage="authoritative_energy",
                dimensions=[4],
                fluxes=[0.0],
                operations,
            )
            @test only(run.selections).candidate_id == "beta"

            ledger_path = joinpath(directory, "ledger.toml")
            ledger = TOML.parsefile(ledger_path)
            alpha_row = only(filter(
                row -> String(row["candidate_id"]) == "alpha",
                ledger["candidate"],
            ))
            alpha_directory = joinpath(directory, alpha_row["directory"])
            forged_energy = reported_energies["alpha"]
            summary_path = joinpath(alpha_directory, "summary.toml")
            summary = TOML.parsefile(summary_path)
            summary["energy"]["per_site"] = forged_energy
            open(summary_path, "w") do io
                TOML.print(io, summary; sorted=true)
            end
            metadata_path = joinpath(alpha_directory, "candidate.toml")
            metadata = TOML.parsefile(metadata_path)
            metadata["energy_per_site"] = forged_energy
            open(metadata_path, "w") do io
                TOML.print(io, metadata; sorted=true)
            end
            alpha_row["checksums"]["summary.toml"] =
                InfiniteCylinderDMRG._fig2_file_sha256(summary_path)
            alpha_row["checksums"]["candidate.toml"] =
                InfiniteCylinderDMRG._fig2_file_sha256(metadata_path)

            forged_alpha = InfiniteCylinderDMRG._persisted_fig2_candidate(
                spec,
                directory,
                alpha_row,
                4,
                1,
                0.0,
                "alpha",
            )
            ledger["selection"][1] = InfiniteCylinderDMRG._selection_dict(
                InfiniteCylinderDMRG._fig2_selection_from_candidate(
                    forged_alpha,
                    4,
                    1,
                    0.0,
                    "minimum_zero_flux_energy",
                )
            )
            open(ledger_path, "w") do io
                TOML.print(io, ledger; sorted=true)
            end

            replay_error = try
                InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                    spec,
                    directory;
                    checkpoint_audit=synthetic_checkpoint_audit,
                    progress_audit=synthetic_progress_audit,
                    candidate_ids_provider=candidate_ids,
                )
                nothing
            catch error
                error
            end
            @test replay_error isa ArgumentError
            replay_message = isnothing(replay_error) ? "" :
                lowercase(sprint(showerror, replay_error))
            @test occursin("selection replay winner mismatch", replay_message)
        end
    end

    @testset "convergence reset sentinel follows expansion activation" begin
        mktempdir() do directory
            write_fake_candidate_files(
                directory,
                "canonical-missing-first-delta";
                requested_maxdim=4,
                summary_maxdim_schedule=[4],
                sector_weights=Dict(0 => 1.0),
                convergence_rows=[
                    (
                        stage=1,
                        iteration=1,
                        maxlinkdim=2,
                        delta_energy="missing",
                        converged=false,
                    ),
                    (
                        stage=1,
                        iteration=2,
                        maxlinkdim=4,
                        delta_energy="missing",
                        converged=false,
                    ),
                    (
                        stage=1,
                        iteration=3,
                        maxlinkdim=4,
                        delta_energy=1.0e-8,
                        converged=false,
                    ),
                    (
                        stage=1,
                        iteration=4,
                        maxlinkdim=4,
                        delta_energy=1.0e-8,
                        converged=true,
                    ),
                ],
            )
            outcome = try
                validate_production_convergence(
                    joinpath(directory, "convergence.tsv"), [4]
                )
            catch error
                error
            end
            @test !(outcome isa Exception)
            if !(outcome isa Exception)
                @test outcome.iteration == 4
                @test outcome.converged
            end
        end

        function rejects_missing_at(row, column)
            mktempdir() do directory
                write_fake_candidate_files(
                    directory,
                    "misplaced-missing";
                    requested_maxdim=4,
                    summary_maxdim_schedule=[4],
                    sector_weights=Dict(0 => 1.0),
                    convergence_rows=[
                        (
                            stage=1,
                            iteration=1,
                            maxlinkdim=2,
                            delta_energy="missing",
                            converged=false,
                        ),
                        (
                            stage=1,
                            iteration=2,
                            maxlinkdim=4,
                            delta_energy="missing",
                            converged=false,
                        ),
                        (
                            stage=1,
                            iteration=3,
                            maxlinkdim=4,
                            delta_energy=1.0e-8,
                            converged=false,
                        ),
                        (
                            stage=1,
                            iteration=4,
                            maxlinkdim=4,
                            delta_energy=1.0e-8,
                            converged=true,
                        ),
                    ],
                )
                path = joinpath(directory, "convergence.tsv")
                lines = readlines(path)
                fields = split(lines[row + 1], '\t'; keepempty=true)
                fields[column] = "missing"
                lines[row + 1] = join(fields, '\t')
                open(path, "w") do io
                    foreach(line -> println(io, line), lines)
                end
                @test_throws ArgumentError begin
                    validate_production_convergence(path, [4])
                end
            end
        end
        rejects_missing_at(1, 4)
        rejects_missing_at(3, 7)
    end

    @testset "producer activation sentinels normalize atomically" begin
        function write_producer_convergence(
            directory;
            first_delta="Inf",
            second_delta="Inf",
        )
            write_fake_candidate_files(
                directory,
                "producer-convergence";
                requested_maxdim=4,
                summary_maxdim_schedule=[4],
                sector_weights=Dict(0 => 1.0),
                convergence_rows=[
                    (
                        stage=1,
                        iteration=1,
                        maxlinkdim=2,
                        delta_energy=first_delta,
                        converged=false,
                    ),
                    (
                        stage=1,
                        iteration=2,
                        maxlinkdim=4,
                        delta_energy=second_delta,
                        converged=false,
                    ),
                    (
                        stage=1,
                        iteration=3,
                        maxlinkdim=4,
                        delta_energy=0.25,
                        converged=true,
                    ),
                ],
            )
            return joinpath(directory, "convergence.tsv")
        end

        mktempdir() do directory
            path = write_producer_convergence(directory)
            normalization_error = try
                InfiniteCylinderDMRG._normalize_fig2_convergence_sentinel!(path)
                nothing
            catch error
                error
            end
            @test isnothing(normalization_error)
            if isnothing(normalization_error)
                lines = readlines(path)
                @test split(lines[2], '\t'; keepempty=true)[7] == "missing"
                @test split(lines[3], '\t'; keepempty=true)[7] == "missing"
                @test split(lines[4], '\t'; keepempty=true)[7] == "0.25"
                final = InfiniteCylinderDMRG._validate_fig2_convergence_tsv(path)
                @test final.iteration == 3
                @test final.converged
            end
        end

        function rejects_normalization(
            tamper; first_delta="Inf", second_delta="Inf"
        )
            mktempdir() do directory
                path = write_producer_convergence(
                    directory; first_delta, second_delta
                )
                tamper(path)
                original = read(path)
                @test_throws ArgumentError begin
                    InfiniteCylinderDMRG._normalize_fig2_convergence_sentinel!(path)
                end
                @test read(path) == original
            end
        end

        identity_tamper(path) = nothing
        rejects_normalization(identity_tamper; first_delta="0.0")
        rejects_normalization(identity_tamper; first_delta="missing")
        rejects_normalization(identity_tamper; second_delta="NaN")
        rejects_normalization() do path
            lines = readlines(path)
            lines[1] = "not\tthe\tproduction\theader"
            open(path, "w") do io
                foreach(line -> println(io, line), lines)
            end
        end
        rejects_normalization() do path
            lines = readlines(path)
            fields = split(lines[2], '\t'; keepempty=true)
            pop!(fields)
            lines[2] = join(fields, '\t')
            open(path, "w") do io
                foreach(line -> println(io, line), lines)
            end
        end
        rejects_normalization() do path
            lines = readlines(path)
            fields = split(lines[2], '\t'; keepempty=true)
            fields[4] = "NaN"
            lines[2] = join(fields, '\t')
            open(path, "w") do io
                foreach(line -> println(io, line), lines)
            end
        end
    end

    @testset "convergence full-table audit rejects invalid stage histories" begin
        function rejects_convergence(rows, schedule)
            mktempdir() do directory
                write_fake_candidate_files(
                    directory,
                    "invalid-convergence-history";
                    requested_maxdim=last(schedule),
                    summary_maxdim_schedule=schedule,
                    sector_weights=Dict(0 => 1.0),
                    convergence_rows=rows,
                )
                @test_throws ArgumentError begin
                    validate_production_convergence(
                        joinpath(directory, "convergence.tsv"), schedule
                    )
                end
            end
        end

        rejects_convergence([
            (stage=1, iteration=1, maxlinkdim=2, converged=false),
            (stage=1, iteration=1, maxlinkdim=4, converged=true),
        ], [4])
        rejects_convergence([
            (stage=1, iteration=2, maxlinkdim=2, converged=false),
            (stage=1, iteration=1, maxlinkdim=4, converged=true),
        ], [4])
        rejects_convergence([
            (stage=2, iteration=1, maxlinkdim=4, converged=true),
        ], [4])
        rejects_convergence([
            (stage=1, iteration=1, maxlinkdim=2, converged=false),
            (stage=3, iteration=1, maxlinkdim=8, converged=true),
        ], [2, 4, 8])
        rejects_convergence([
            (stage=1, iteration=1, maxlinkdim=5, converged=true),
        ], [4])

        function rejects_convergence_float(column, value)
            mktempdir() do directory
                write_fake_candidate_files(
                    directory,
                    "invalid-convergence-float";
                    requested_maxdim=4,
                    sector_weights=Dict(0 => 1.0),
                )
                path = joinpath(directory, "convergence.tsv")
                lines = readlines(path)
                fields = split(lines[2], '\t'; keepempty=true)
                fields[column] = value
                lines[2] = join(fields, '\t')
                open(path, "w") do io
                    foreach(line -> println(io, line), lines)
                end
                @test_throws ArgumentError begin
                    validate_production_convergence(path, [4])
                end
            end
        end

        rejects_convergence_float(4, "NaN")
        rejects_convergence_float(7, "Inf")
        rejects_convergence_float(10, "Inf")
        rejects_convergence_float(11, "-0.01")
    end

    @testset "Fig. 2 provenance commit is the actual repository HEAD" begin
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        actual_commit = synthetic_repository_commit()
        synthetic = synthetic_fig2_provenance(spec, pwd(), 0.0)
        @test synthetic["git_commit"] == actual_commit

        wrong_commit = actual_commit == repeat("a", 40) ?
            repeat("b", 40) : repeat("a", 40)
        forged = synthetic_fig2_provenance(
            spec,
            pwd(),
            0.0;
            git_commit=wrong_commit,
        )
        @test_throws ArgumentError InfiniteCylinderDMRG._fig2_generation_provenance(
            spec, forged
        )

        mktempdir() do directory
            InfiniteCylinderDMRG._write_fig2_toml(
                joinpath(directory, "provenance.toml"), forged
            )
            result = InfiniteCylinderDMRG._fig2_provenance_result(
                spec, directory
            )
            @test !result.valid
            @test occursin("commit", lowercase(result.reason))
        end

        for (key, invalid) in (
            ("blas_threads", 2),
            ("strided_threads", 2),
            ("blocksparse_threaded", Threads.nthreads() == 1),
        )
            tampered = deepcopy(synthetic)
            tampered[key] = invalid
            @test_throws ArgumentError InfiniteCylinderDMRG._fig2_generation_provenance(
                spec, tampered
            )
            mktempdir() do directory
                InfiniteCylinderDMRG._write_fig2_toml(
                    joinpath(directory, "provenance.toml"), tampered
                )
                result = InfiniteCylinderDMRG._fig2_provenance_result(
                    spec, directory
                )
                @test !result.valid
            end
        end
    end

    @testset "Fig. 2 repository HEAD resolver supports worktrees and packed refs" begin
        loose_commit = "0123456789abcdef0123456789abcdef01234567"
        packed_commit = "89abcdef0123456789abcdef0123456789abcdef"
        detached_commit = "fedcba9876543210fedcba9876543210fedcba98"

        mktempdir() do repository
            git_directory = joinpath(repository, ".git")
            reference_directory = joinpath(git_directory, "refs", "heads")
            mkpath(reference_directory)
            write(joinpath(git_directory, "HEAD"), "ref: refs/heads/DMRG\n")
            write(joinpath(reference_directory, "DMRG"), "$loose_commit\n")
            @test InfiniteCylinderDMRG._fig2_repository_commit(repository) ==
                loose_commit
        end

        mktempdir() do directory
            repository = joinpath(directory, "repository")
            git_directory = joinpath(directory, "worktrees", "benchmark")
            common_directory = joinpath(directory, "common")
            mkpath(repository)
            mkpath(git_directory)
            mkpath(joinpath(common_directory, "refs", "heads"))
            write(joinpath(repository, ".git"), "gitdir: ../worktrees/benchmark\n")
            write(joinpath(git_directory, "commondir"), "../../common\n")
            write(joinpath(git_directory, "HEAD"), "ref: refs/heads/DMRG\n")
            write(
                joinpath(common_directory, "refs", "heads", "DMRG"),
                "$loose_commit\n",
            )
            @test InfiniteCylinderDMRG._fig2_repository_commit(repository) ==
                loose_commit
        end

        mktempdir() do directory
            repository = joinpath(directory, "repository")
            git_directory = joinpath(directory, "worktrees", "benchmark")
            common_directory = joinpath(directory, "common")
            mkpath(repository)
            mkpath(git_directory)
            mkpath(common_directory)
            write(joinpath(repository, ".git"), "gitdir: ../worktrees/benchmark\n")
            write(joinpath(git_directory, "commondir"), "../../common\n")
            write(joinpath(git_directory, "HEAD"), "ref: refs/heads/DMRG\n")
            write(
                joinpath(common_directory, "packed-refs"),
                "# pack-refs with: peeled fully-peeled sorted\n" *
                "$packed_commit refs/heads/DMRG\n",
            )
            @test InfiniteCylinderDMRG._fig2_repository_commit(repository) ==
                packed_commit
        end

        mktempdir() do repository
            git_directory = joinpath(repository, ".git")
            mkpath(git_directory)
            write(joinpath(git_directory, "HEAD"), "$detached_commit\n")
            @test InfiniteCylinderDMRG._fig2_repository_commit(repository) ==
                detached_commit
        end
    end

    @testset "global two-process restart gate is recomputed fail closed" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            provenance = synthetic_fig2_provenance(spec, directory, 0.0)
            function gate_result(overrides=Dict{String,Any}())
                write_synthetic_restart_gate!(spec, directory; overrides)
                return InfiniteCylinderDMRG._fig2_restart_gate_result(
                    spec, directory, provenance
                )
            end

            valid = gate_result()
            @test valid.valid
            if valid.valid
                @test valid.measured["pre_maxlinkdim"] == 4
                @test valid.measured["post_maxlinkdim"] == 8
            end
            @test_throws InterruptException InfiniteCylinderDMRG._fig2_restart_gate_result(
                spec,
                directory,
                provenance;
                parse_gate=(args...) -> throw(InterruptException()),
            )
            ordinary_parse_failure =
                InfiniteCylinderDMRG._fig2_restart_gate_result(
                    spec,
                    directory,
                    provenance;
                    parse_gate=(args...) -> error("synthetic parse failure"),
                )
            @test !ordinary_parse_failure.valid
            @test occursin(
                "synthetic parse failure", ordinary_parse_failure.reason
            )

            evidence_controlled_energy = gate_result(Dict{String,Any}(
                "energy_absolute_error" => 1.0,
                "energy_tolerance" => 1.0,
            ))
            @test !evidence_controlled_energy.valid
            @test occursin(
                "energy", lowercase(evidence_controlled_energy.reason)
            )
            evidence_controlled_density = gate_result(Dict{String,Any}(
                "density_max_absolute_error" => 1.0,
                "density_tolerance" => 1.0,
            ))
            @test !evidence_controlled_density.valid
            @test occursin(
                "density", lowercase(evidence_controlled_density.reason)
            )

            @test spec.data["restart_gate"] == Dict{String,Any}(
                "format" => "fqahc_fig2_restart_gate_v1",
                "evidence_filename" => "restart_gate.toml",
                "backend_commit" =>
                    "765f2777703bc1138b009adbed1b97bde1973402",
                "energy_tolerance" => 1.0e-10,
                "density_tolerance" => 1.0e-10,
                "required_checks" => [
                    "load_completed",
                    "no_expansion_iteration_completed",
                    "controlled_expansion_completed",
                    "energy_reproduced",
                    "density_reproduced",
                    "schmidt_sectors_reproduced",
                    "entanglement_spectrum_reproduced",
                    "site_indices_reproduced",
                    "next_flux_iteration_completed",
                ],
            )
            function rejects_restart_manifest(tamper::Function)
                mktempdir() do scratch
                    data = TOML.parsefile(FIG2_MANIFEST_PATH)
                    tamper(data["restart_gate"])
                    path = joinpath(scratch, "fqahc_fig2.toml")
                    open(path, "w") do io
                        TOML.print(io, data; sorted=true)
                    end
                    @test_throws ArgumentError load_fig2_benchmark(path)
                end
            end
            for key in ("energy_tolerance", "density_tolerance")
                rejects_restart_manifest(
                    gate -> gate[key] = 1.0e-9
                )
                rejects_restart_manifest(
                    gate -> pop!(gate, key, nothing)
                )
            end

            for (label, removed_keys, reason_fragment) in (
                (
                    "missing energy pair",
                    ("energy_absolute_error", "energy_tolerance"),
                    "energy",
                ),
                (
                    "missing energy error",
                    ("energy_absolute_error",),
                    "energy",
                ),
                (
                    "missing energy tolerance",
                    ("energy_tolerance",),
                    "energy",
                ),
                (
                    "missing density pair",
                    ("density_max_absolute_error", "density_tolerance"),
                    "density",
                ),
                (
                    "missing density error",
                    ("density_max_absolute_error",),
                    "density",
                ),
                (
                    "missing density tolerance",
                    ("density_tolerance",),
                    "density",
                ),
            )
                @testset "$label" begin
                    gate = synthetic_restart_gate(spec, directory)
                    foreach(key -> delete!(gate, key), removed_keys)
                    InfiniteCylinderDMRG._write_fig2_toml(
                        joinpath(directory, "restart_gate.toml"), gate
                    )
                    result = InfiniteCylinderDMRG._fig2_restart_gate_result(
                        spec, directory, provenance
                    )
                    @test !result.valid
                    @test occursin(reason_fragment, lowercase(result.reason))
                end
            end

            tamper_cases = [
                ("format", Dict{String,Any}("format" => "wrong"), "format"),
                ("pass", Dict{String,Any}("pass" => false), "pass"),
                ("manifest", Dict{String,Any}(
                    "manifest_sha256" => repeat("0", 64),
                ), "manifest"),
                ("commit", Dict{String,Any}(
                    "git_commit" => repeat("0", 40),
                ), "commit"),
                ("Julia", Dict{String,Any}(
                    "julia_version" => "1.12.4",
                ), "julia"),
                ("project", Dict{String,Any}(
                    "project_manifest_sha256" => repeat("0", 64),
                ), "project"),
                ("backend", Dict{String,Any}(
                    "backend_commit" => repeat("0", 40),
                ), "backend"),
                ("save exit", Dict{String,Any}(
                    "save_process_exit_code" => 1,
                ), "save process"),
                ("resume exit", Dict{String,Any}(
                    "resume_process_exit_code" => 1,
                ), "resume process"),
                ("save identity", Dict{String,Any}(
                    "save_process_id" => "",
                ), "save process"),
                ("distinct identities", Dict{String,Any}(
                    "resume_process_id" => "save:101",
                ), "distinct"),
                ("pre maxlinkdim", Dict{String,Any}(
                    "pre_maxlinkdim" => 0,
                ), "pre_maxlinkdim"),
                ("bond growth", Dict{String,Any}(
                    "post_maxlinkdim" => 4,
                ), "post_maxlinkdim"),
                ("energy error", Dict{String,Any}(
                    "energy_absolute_error" => 2.0e-10,
                ), "energy"),
                ("density error", Dict{String,Any}(
                    "density_max_absolute_error" => 2.0e-10,
                ), "density"),
            ]
            for required_check in spec.data["restart_gate"]["required_checks"]
                push!(tamper_cases, (
                    required_check,
                    Dict{String,Any}(required_check => false),
                    required_check,
                ))
            end
            for (label, overrides, reason_fragment) in tamper_cases
                @testset "$label" begin
                    result = gate_result(overrides)
                    @test !result.valid
                    @test occursin(
                        lowercase(reason_fragment), lowercase(result.reason)
                    )
                end
            end
        end
    end

    @testset "candidate generation provenance survives cross-job resume" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            candidate_calls = String[]
            function operations_for(commit, job)
                return Fig2BenchmarkOperations(
                    checkpoint_audit=synthetic_checkpoint_audit,
                    progress_audit=synthetic_progress_audit,
                    provenance=(spec, output, runtime_seconds) ->
                        synthetic_fig2_provenance(
                            spec,
                            output,
                            runtime_seconds;
                            git_commit=commit,
                            pbs_job_id=job,
                        ),
                    candidate_ids=(args...) -> ["only"],
                    run_candidate=(spec, dimension, point, phi_y, candidate_id,
                                   previous_state, candidate_directory) -> begin
                        push!(candidate_calls, "$commit:$job")
                        write_fake_candidate_files(
                            candidate_directory,
                            job;
                            requested_maxdim=dimension,
                            sector_weights=Dict(0 => 1.0),
                        )
                        momentum = synthetic_counting_momentum()
                        counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                            spec, momentum
                        )
                        Fig2CandidateEvidence(;
                            state=(; commit, job),
                            requested_maxdim=dimension,
                            achieved_maxlinkdim=dimension,
                            converged=true,
                            valid=true,
                            energy_per_site=-1.0,
                            raw_schmidt_polarization=0.0,
                            sector_weights=Dict(0 => 1.0),
                            momentum,
                            momentum_counting=counting.counts,
                            momentum_counting_evidence=counting,
                            restart_valid=true,
                        )
                    end,
                )
            end

            first_commit = synthetic_repository_commit()
            wrong_commit = first_commit == repeat("b", 40) ?
                repeat("c", 40) : repeat("b", 40)
            run_fig2_benchmark(
                spec,
                directory;
                stage="provenance_resume",
                dimensions=[4],
                fluxes=[0.0],
                operations=operations_for(first_commit, "first.w003"),
            )
            @test candidate_calls == ["$first_commit:first.w003"]

            run_fig2_benchmark(
                spec,
                directory;
                stage="provenance_resume",
                dimensions=[4],
                fluxes=[0.0],
                operations=operations_for(first_commit, "second.w003"),
            )
            @test candidate_calls == ["$first_commit:first.w003"]

            ledger = TOML.parsefile(joinpath(directory, "ledger.toml"))
            candidate = only(ledger["candidate"])
            metadata = TOML.parsefile(joinpath(
                directory, candidate["directory"], "candidate.toml"
            ))
            @test metadata["format"] == "fqahc_fig2_candidate_v4"
            @test metadata["requested_maxdim"] == 4
            @test metadata["achieved_maxlinkdim"] == 4
            @test get(metadata, "checkpoint_maxlinkdim", nothing) == 4
            @test get(candidate, "checkpoint_maxlinkdim", nothing) == 4
            @test haskey(metadata, "generation_provenance")
            if haskey(metadata, "generation_provenance")
                generation = metadata["generation_provenance"]
                @test generation["git_commit"] == first_commit
                @test generation["pbs_job_id"] == "first.w003"
                @test generation["blas_threads"] == 1
                @test generation["strided_threads"] == 1
                @test generation["blocksparse_threaded"] ===
                    (Threads.nthreads() > 1)
                @test generation["benchmark_source_sha256"] ==
                    synthetic_fig2_provenance(
                        spec, directory, 0.0
                    )["benchmark_source_sha256"]
                @test generation["project_manifest_sha256"] ==
                    synthetic_fig2_provenance(
                        spec, directory, 0.0
                    )["project_manifest_sha256"]
            end
            @test haskey(candidate, "generation_provenance_sha256")
            root_provenance = TOML.parsefile(joinpath(directory, "provenance.toml"))
            @test get(root_provenance, "candidate_pbs_job_ids", String[]) ==
                ["first.w003"]

            @test_throws ArgumentError run_fig2_benchmark(
                spec,
                directory;
                stage="provenance_resume",
                dimensions=[4],
                fluxes=[0.0],
                operations=operations_for(wrong_commit, "third.w003"),
            )
        end
    end

    @testset "paper curve is recomputed from immutable digitization and raw pump" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            curve = spec.data["paper_curve_rule"]
            flux_units = collect(0.0:0.125:3.0)
            target_values = [
                0.000000, 0.039286, 0.080329, 0.122619, 0.165646,
                0.208617, 0.250794, 0.292290, 0.332370, 0.371882,
                0.410998, 0.450680, 0.490760, 0.532426, 0.575907,
                0.620918, 0.667290, 0.713152, 0.760601, 0.807029,
                0.851927, 0.889739, 0.925283, 0.961735, 1.000000,
            ]
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(args...) -> ["only"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    sectors = synthetic_sector_weights_for_polarization(
                        target_values[point]
                    )
                    write_fake_candidate_files(
                        candidate_directory,
                        "curve-$point";
                        requested_maxdim=dimension,
                        sector_weights=sectors,
                    )
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    fidelity = point == 1 ? NaN : 0.99
                    Fig2CandidateEvidence(;
                        state=(; point),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=-1.0,
                        raw_schmidt_polarization=target_values[point],
                        sector_weights=sectors,
                        fidelity_to_previous=fidelity,
                        fidelity_valid=point > 1,
                        mixed_fidelity=point == 1 ? nothing :
                            synthetic_mixed_fidelity(fidelity),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                    )
                end,
            )
            run_fig2_benchmark(
                spec,
                directory;
                stage="paper_curve",
                dimensions=[2000],
                fluxes=2pi .* flux_units,
                operations,
            )

            tracked_lines = readlines(joinpath(
                directory, "pump_sector_tracked.tsv"
            ))
            tracked_fields = split.(tracked_lines[2:end], '\t')
            sector_tolerance =
                spec.data["convergence"]["sector_residual_tol"]
            @test any(
                parse(Float64, fields[7]) > sector_tolerance &&
                    parse(Bool, fields[8]) for fields in tracked_fields
            )

            report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test report["paper_curve"]["pass"]
            comparison_path = joinpath(directory, "paper_curve_comparison.toml")
            @test isfile(comparison_path)
            if isfile(comparison_path)
                comparison = TOML.parsefile(comparison_path)
                @test comparison["source_vector_path_sha256"] ==
                    curve["source_vector_path_sha256"]
                points = comparison["point"]
                @test length(points) == 25
                if !isempty(points)
                    @test maximum(row["absolute_residual"] for row in points) <=
                        1e-12
                    @test all(row["pass"] for row in points)
                end
            end

            raw_path = joinpath(directory, "pump_raw.tsv")
            raw_before_tamper = read(raw_path)
            raw_lines = readlines(raw_path)
            fields = split(last(raw_lines), '\t')
            fields[5] = "1.3"
            fields[6] = "1.3"
            raw_lines[end] = join(fields, '\t')
            open(raw_path, "w") do io
                foreach(line -> println(io, line), raw_lines)
            end
            open(comparison_path, "w") do io
                TOML.print(io, Dict(
                    "format" => "fqahc_fig2_paper_curve_comparison_v1",
                    "source_pdf_sha256" => curve["source_pdf_sha256"],
                    "source_page" => curve["source_page"],
                    "source_panel" => curve["source_panel"],
                    "maximum_residual" => 0.0,
                    "tolerance" => spec.data["pump_tolerance"],
                    "pass" => true,
                ); sorted=true)
            end
            tampered = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test tampered["paper_curve"]["pass"]
            measured = tampered["paper_curve"]["measured"]
            @test measured isa Real
            if measured isa Real
                @test measured <= 1e-12
            end
            @test read(raw_path) == raw_before_tamper
        end
    end

    function fig2_reference_curve_value(rule, flux_units)
        reference_fluxes = Float64.(rule["flux_units_2pi"])
        reference_values = Float64.(rule["delta_q"])
        exact = findfirst(value -> isapprox(
            value, flux_units; atol=1e-12, rtol=0
        ), reference_fluxes)
        !isnothing(exact) && return reference_values[exact]
        right = findfirst(value -> value > flux_units, reference_fluxes)
        isnothing(right) && error("synthetic flux lies above the paper curve")
        right == 1 && error("synthetic flux lies below the paper curve")
        left = right - 1
        weight = (flux_units - reference_fluxes[left]) /
            (reference_fluxes[right] - reference_fluxes[left])
        return (1 - weight) * reference_values[left] +
            weight * reference_values[right]
    end

    function run_fig2_curve_grid(directory, flux_units; dimension=2000)
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        curve_rule = spec.data["paper_curve_rule"]
        operations = Fig2BenchmarkOperations(
            progress_audit=synthetic_progress_audit,
            provenance=synthetic_fig2_provenance,
            candidate_ids=(args...) -> ["only"],
            run_candidate=(spec, requested_dimension, point, phi_y, candidate_id,
                           previous_state, candidate_directory) -> begin
                polarization = fig2_reference_curve_value(
                    curve_rule, phi_y / (2pi)
                )
                sectors = synthetic_sector_weights_for_polarization(
                    polarization
                )
                write_fake_candidate_files(
                    candidate_directory,
                    "interpolation-D=$requested_dimension-point=$point";
                    requested_maxdim=requested_dimension,
                    sector_weights=sectors,
                )
                first_point = point == 1
                fidelity = first_point ? NaN : 0.99
                momentum = synthetic_counting_momentum()
                counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                    spec, momentum
                )
                Fig2CandidateEvidence(;
                    state=(; requested_dimension, point),
                    requested_maxdim=requested_dimension,
                    achieved_maxlinkdim=requested_dimension,
                    converged=true,
                    valid=true,
                    energy_per_site=-1.0,
                    raw_schmidt_polarization=polarization,
                    sector_weights=sectors,
                    fidelity_to_previous=fidelity,
                    fidelity_valid=!first_point,
                    mixed_fidelity=first_point ? nothing :
                        synthetic_mixed_fidelity(fidelity),
                    momentum,
                    momentum_counting=counting.counts,
                    momentum_counting_evidence=counting,
                    restart_valid=true,
                )
            end,
        )
        run = run_fig2_benchmark(
            spec,
            directory;
            stage="paper_curve_interpolation",
            dimensions=[dimension],
            fluxes=2pi .* Float64.(flux_units),
            operations,
        )
        candidate_ids_provider = operations.candidate_ids
        integrity = synthetic_fig2_acceptance_integrity_inputs(
            spec, directory; candidate_ids_provider
        )
        return (; spec, run, integrity, dimension, candidate_ids_provider)
    end

    function recompute_fig2_curve(
        fixture,
        comparison_filename,
    )
        curve_function = InfiniteCylinderDMRG._fig2_paper_curve_result
        if applicable(
            curve_function,
            fixture.spec,
            fixture.run.output,
            fixture.dimension,
            fixture.integrity.selections,
        )
            return curve_function(
                fixture.spec,
                fixture.run.output,
                fixture.dimension,
                fixture.integrity.selections;
                comparison_filename,
            )
        end
        return curve_function(
            fixture.spec,
            fixture.run.output,
            fixture.dimension;
            comparison_filename,
        )
    end

    @testset "deterministic paper-curve interpolation" begin
    @testset "paper curve linearly interpolates a sparse validated raw grid" begin
        mktempdir() do directory
            fixture = run_fig2_curve_grid(directory, collect(0.0:0.5:3.0))
            report = write_synthetic_fig2_acceptance_report!(
                fixture.spec,
                directory;
                candidate_ids_provider=fixture.candidate_ids_provider,
            )
            @test report["paper_curve"]["pass"]

            comparison_path = joinpath(directory, "paper_curve_comparison.toml")
            comparison = TOML.parsefile(comparison_path)
            @test comparison["format"] ==
                "fqahc_fig2_paper_curve_comparison_v3"
            @test comparison["sampling_policy"] ==
                "linear_interpolation_on_validated_raw_grid_without_extrapolation"
            @test comparison["pump_raw_sha256"] ==
                InfiniteCylinderDMRG._fig2_file_sha256(joinpath(
                    directory, "pump_raw.tsv"
                ))
            @test comparison["pump_sector_tracked_sha256"] ==
                InfiniteCylinderDMRG._fig2_file_sha256(joinpath(
                    directory, "pump_sector_tracked.tsv"
                ))
            @test length(comparison["point"]) == 25
            interpolated = comparison["point"][2]
            @test interpolated["sampling"] == "linear_interpolation"
            @test !interpolated["exact_match"]
            @test interpolated["left_point"] == 1
            @test interpolated["right_point"] == 2
            @test interpolated["left_candidate_id"] == "only"
            @test interpolated["right_candidate_id"] == "only"
            @test interpolated["left_phi_y"] == 0.0
            @test interpolated["right_phi_y"] == Float64(pi)
            @test interpolated["left_weight"] == 0.75
            @test interpolated["right_weight"] == 0.25
            @test interpolated["left_branch_valid"]
            @test interpolated["right_branch_valid"]
            @test interpolated["left_sector_residual"] == 0.0
            left_sectors = synthetic_sector_weights_for_polarization(
                fig2_reference_curve_value(
                    fixture.spec.data["paper_curve_rule"], 0.0
                )
            )
            right_sectors = synthetic_sector_weights_for_polarization(
                fig2_reference_curve_value(
                    fixture.spec.data["paper_curve_rule"], 0.5
                )
            )
            expected_right_residual = best_sector_shift(
                left_sectors, right_sectors
            ).distance
            @test interpolated["right_sector_residual"] ==
                expected_right_residual
            @test expected_right_residual >
                fixture.spec.data["convergence"]["sector_residual_tol"]
            @test interpolated["uncertainty"] ==
                fixture.spec.data["paper_curve_rule"]["uncertainty"][2]
        end
    end

    @testset "paper curve keeps exact brackets when refinement points are present" begin
        mktempdir() do directory
            reference = collect(0.0:0.125:3.0)
            refined = sort!(unique(vcat(reference, [0.0625, 0.1875, 1.3125])))
            fixture = run_fig2_curve_grid(directory, refined)
            report = write_synthetic_fig2_acceptance_report!(
                fixture.spec,
                directory;
                candidate_ids_provider=fixture.candidate_ids_provider,
            )
            @test report["paper_curve"]["pass"]
            comparison = TOML.parsefile(joinpath(
                directory, "paper_curve_comparison.toml"
            ))
            exact = comparison["point"][2]
            @test exact["sampling"] == "exact"
            @test exact["exact_match"]
            @test exact["left_point"] == exact["right_point"]
            @test exact["left_candidate_id"] == exact["right_candidate_id"]
            @test exact["left_phi_y"] == exact["right_phi_y"] == 0.25pi
            @test exact["left_weight"] == 1.0
            @test exact["right_weight"] == 0.0
        end
    end

    @testset "paper curve forbids extrapolation and requests refinement" begin
        mktempdir() do directory
            fixture = run_fig2_curve_grid(directory, collect(0.25:0.25:2.75))
            report = write_synthetic_fig2_acceptance_report!(
                fixture.spec,
                directory;
                candidate_ids_provider=fixture.candidate_ids_provider,
            )
            @test !report["paper_curve"]["pass"]
            reason = lowercase(report["paper_curve"]["reason"])
            @test occursin("extrapolat", reason)
            @test occursin("refinement", reason)
        end
    end

    @testset "paper curve strictly validates raw and tracked interpolation evidence" begin
        mktempdir() do directory
            fixture = run_fig2_curve_grid(directory, collect(0.0:0.125:3.0))
            raw_path = joinpath(directory, "pump_raw.tsv")
            tracked_path = joinpath(directory, "pump_sector_tracked.tsv")
            raw_lines = readlines(raw_path)
            tracked_lines = readlines(tracked_path)

            function rejected_curve(label, fragment)
                result = recompute_fig2_curve(
                    fixture, "paper_curve_rejected_$label.toml"
                )
                @test !result.pass
                @test occursin(fragment, lowercase(result.reason))
            end

            open(raw_path, "w") do io
                foreach(line -> println(io, line), vcat(raw_lines, raw_lines[2]))
            end
            rejected_curve("duplicate", "duplicate")

            open(raw_path, "w") do io
                swapped = copy(raw_lines)
                swapped[3], swapped[4] = swapped[4], swapped[3]
                foreach(line -> println(io, line), swapped)
            end
            rejected_curve("unsorted", "sorted")

            open(raw_path, "w") do io
                tampered = copy(raw_lines)
                fields = split(tampered[4], '\t')
                fields[6] = string(parse(Float64, fields[6]) + 0.25)
                tampered[4] = join(fields, '\t')
                foreach(line -> println(io, line), tampered)
            end
            rejected_curve("delta", "delta")

            open(raw_path, "w") do io
                tampered = copy(raw_lines)
                fields = split(tampered[4], '\t')
                fields[4] = "forged-candidate"
                tampered[4] = join(fields, '\t')
                foreach(line -> println(io, line), tampered)
            end
            rejected_curve("candidate", "candidate")

            open(raw_path, "w") do io
                foreach(line -> println(io, line), raw_lines)
            end
            open(tracked_path, "w") do io
                tampered = copy(tracked_lines)
                fields = split(tampered[4], '\t')
                fields[4] = "forged-tracked-candidate"
                tampered[4] = join(fields, '\t')
                foreach(line -> println(io, line), tampered)
            end
            rejected_curve("tracked_identity", "tracked")
            @test occursin(
                "identity",
                lowercase(recompute_fig2_curve(
                    fixture, "paper_curve_rejected_tracked_identity_detail.toml"
                ).reason),
            )
        end
    end

    @testset "paper curve rejects an invalid interpolation branch" begin
        mktempdir() do directory
            fixture = run_fig2_curve_grid(directory, collect(0.0:0.5:3.0))
            tracked_path = joinpath(directory, "pump_sector_tracked.tsv")
            tracked_lines = readlines(tracked_path)
            fields = split(tracked_lines[3], '\t')
            fields[8] = "false"
            fields[9] = "synthetic invalid branch"
            tracked_lines[3] = join(fields, '\t')
            open(tracked_path, "w") do io
                foreach(line -> println(io, line), tracked_lines)
            end
            result = recompute_fig2_curve(
                fixture, "paper_curve_invalid_branch.toml"
            )
            @test !result.pass
            reason = lowercase(result.reason)
            @test occursin("branch", reason)
            @test occursin("refinement", reason)
        end
    end

    @testset "acceptance rebuilds tampered derived interpolation tables" begin
        mktempdir() do directory
            fixture = run_fig2_curve_grid(directory, collect(0.0:0.5:3.0))
            first_report = write_synthetic_fig2_acceptance_report!(
                fixture.spec,
                directory;
                candidate_ids_provider=fixture.candidate_ids_provider,
            )
            @test first_report["paper_curve"]["pass"]
            raw_path = joinpath(directory, "pump_raw.tsv")
            tracked_path = joinpath(directory, "pump_sector_tracked.tsv")
            comparison_path = joinpath(directory, "paper_curve_comparison.toml")
            raw_before = read(raw_path)
            tracked_before = read(tracked_path)
            open(raw_path, "w") do io
                write(io, "tampered raw table\n")
            end
            open(tracked_path, "w") do io
                write(io, "tampered tracked table\n")
            end
            open(comparison_path, "w") do io
                TOML.print(io, Dict("pass" => true, "reason" => "forged"))
            end

            rebuilt = write_synthetic_fig2_acceptance_report!(
                fixture.spec,
                directory;
                candidate_ids_provider=fixture.candidate_ids_provider,
            )
            @test rebuilt["paper_curve"]["pass"]
            @test read(raw_path) == raw_before
            @test read(tracked_path) == tracked_before
            comparison = TOML.parsefile(comparison_path)
            @test comparison["format"] ==
                "fqahc_fig2_paper_curve_comparison_v3"
            @test comparison["pump_raw_sha256"] ==
                InfiniteCylinderDMRG._fig2_file_sha256(raw_path)
            @test comparison["pump_sector_tracked_sha256"] ==
                InfiniteCylinderDMRG._fig2_file_sha256(tracked_path)
        end
        end
    end

    @testset "selection replay requires complete contiguous flux-point keys" begin
        mktempdir() do directory
            fixture = run_fig2_curve_grid(directory, collect(0.0:0.5:3.0))
            ledger_path = joinpath(directory, "ledger.toml")
            ledger_bytes = read(ledger_path)

            function rejects_incomplete_replay(forged)
                open(ledger_path, "w") do io
                    TOML.print(io, forged; sorted=true)
                end
                report = write_synthetic_fig2_acceptance_report!(
                    fixture.spec,
                    directory;
                    candidate_ids_provider=fixture.candidate_ids_provider,
                )
                @test all(!row["pass"] for row in values(report))
                @test all(
                    occursin("selection replay", lowercase(row["reason"])) for
                    row in values(report)
                )
            end

            missing_selection = TOML.parse(String(copy(ledger_bytes)))
            filter!(
                row -> Int(row["point"]) != 4,
                missing_selection["selection"],
            )
            rejects_incomplete_replay(missing_selection)

            missing_group = TOML.parse(String(copy(ledger_bytes)))
            filter!(
                row -> Int(row["point"]) != 4,
                missing_group["selection"],
            )
            filter!(
                row -> Int(row["point"]) != 4,
                missing_group["candidate"],
            )
            rejects_incomplete_replay(missing_group)

            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end
        end
    end

    function fig2_argument_error_message(f::Function)
        exception = try
            f()
            nothing
        catch error
            error
        end
        @test exception isa ArgumentError
        return isnothing(exception) ? "" : sprint(showerror, exception)
    end

    function run_claimed_dimension_case(
        directory;
        dimension=3000,
        requested_maxdim=dimension,
        achieved_maxlinkdim=requested_maxdim,
        checkpoint_maxlinkdim=achieved_maxlinkdim,
        summary_maxdim_schedule=
            InfiniteCylinderDMRG._fig2_maxdim_schedule(dimension),
        summary_optimization_overrides=Dict{String,Any}(),
        convergence_rows=
            synthetic_fig2_convergence_rows(summary_maxdim_schedule),
    )
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        operations = Fig2BenchmarkOperations(
            progress_audit=synthetic_progress_audit,
            provenance=synthetic_fig2_provenance,
            candidate_ids=(args...) -> ["only"],
            run_candidate=(spec, dimension, point, phi_y, candidate_id,
                           previous_state, candidate_directory) -> begin
                write_fake_candidate_files(
                    candidate_directory,
                    "claimed-D";
                    requested_maxdim,
                    achieved_maxlinkdim,
                    summary_maxdim_schedule,
                    summary_optimization_overrides,
                    sector_weights=Dict(0 => 1.0),
                    convergence_rows,
                )
                momentum = synthetic_counting_momentum()
                counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                    spec, momentum
                )
                Fig2CandidateEvidence(;
                    state=(; dimension),
                    requested_maxdim,
                    achieved_maxlinkdim,
                    checkpoint_maxlinkdim,
                    converged=true,
                    valid=true,
                    energy_per_site=-1.0,
                    raw_schmidt_polarization=0.0,
                    sector_weights=Dict(0 => 1.0),
                    momentum,
                    momentum_counting=counting.counts,
                    momentum_counting_evidence=counting,
                    restart_valid=true,
                )
            end,
        )
        return run_fig2_benchmark(
            spec,
            directory;
            stage="claimed_dimension",
            dimensions=[dimension],
            fluxes=[0.0],
            operations,
        )
    end

    @testset "persisted candidate audit accepts only canonical maxdim suffixes" begin
        for schedule in ([32], [16, 32])
            mktempdir() do directory
                run = run_claimed_dimension_case(
                    directory;
                    dimension=32,
                    summary_maxdim_schedule=schedule,
                )
                @test run isa Fig2BenchmarkRun
                @test only(run.selections).dimension == 32
            end
        end

        for schedule in ([16], [8, 32], [12, 16, 32], [16, 32, 32])
            mktempdir() do directory
                message = fig2_argument_error_message() do
                    run_claimed_dimension_case(
                        directory;
                        dimension=32,
                        summary_maxdim_schedule=schedule,
                    )
                end
                @test occursin("maxdim schedule", message)
            end
        end
    end

    @testset "claimed bond dimension is proved by every candidate artifact" begin
        claimed_schedule = InfiniteCylinderDMRG._fig2_maxdim_schedule(3000)
        valid_claimed_rows = synthetic_fig2_convergence_rows(claimed_schedule)
        cases = [
            (
                label="requested D disagrees with ledger dimension",
                kwargs=(; requested_maxdim=2000, achieved_maxlinkdim=2000),
                fragments=("requested", "dimension"),
            ),
            (
                label="summary has the wrong maxdim schedule",
                kwargs=(; summary_maxdim_schedule=[2000]),
                fragments=("summary", "maxdim"),
            ),
            (
                label="summary update algorithm disagrees with the manifest",
                kwargs=(; summary_optimization_overrides=Dict(
                    "multisite_update_alg" => "parallel",
                )),
                fragments=("summary", "multisite_update_alg", "manifest"),
            ),
            (
                label="summary vumps tolerance disagrees with the manifest",
                kwargs=(; summary_optimization_overrides=Dict(
                    "vumps_tol" => 2.0e-6,
                )),
                fragments=("summary", "vumps_tol", "manifest"),
            ),
            (
                label="summary energy tolerance disagrees with the manifest",
                kwargs=(; summary_optimization_overrides=Dict(
                    "energy_tol" => 2.0e-6,
                )),
                fragments=("summary", "energy_tol", "manifest"),
            ),
            (
                label="summary energy mismatch tolerance disagrees with manifest",
                kwargs=(; summary_optimization_overrides=Dict(
                    "energy_mismatch_tol" => 2.0e-6,
                )),
                fragments=("summary", "energy_mismatch_tol", "manifest"),
            ),
            (
                label="summary stable iterations disagree with the manifest",
                kwargs=(; summary_optimization_overrides=Dict(
                    "stable_iterations" => 3,
                )),
                fragments=("summary", "stable_iterations", "manifest"),
            ),
            (
                label="restartable checkpoint disagrees with in-memory maxlinkdim",
                kwargs=(; checkpoint_maxlinkdim=2000),
                fragments=("checkpoint", "achieved"),
            ),
            (
                label="convergence table is empty",
                kwargs=(; convergence_rows=NamedTuple[]),
                fragments=("convergence.tsv", "data"),
            ),
            (
                label="final convergence flag disagrees with summary",
                kwargs=(; convergence_rows=valid_claimed_rows[1:(end - 1)]),
                fragments=("summary", "final"),
            ),
            (
                label="final convergence row has the wrong maxlinkdim",
                kwargs=(; convergence_rows=synthetic_fig2_convergence_rows(
                    claimed_schedule; final_maxlinkdim=2000
                )),
                fragments=("final", "maxlinkdim"),
            ),
            (
                label="convergence integer columns are parsed strictly",
                kwargs=(; convergence_rows=vcat(
                    valid_claimed_rows[1:(end - 3)],
                    [(
                        stage=length(claimed_schedule),
                        iteration=1,
                        maxlinkdim="3000.0",
                        converged=false,
                    )],
                )),
                fragments=("convergence.tsv", "maxlinkdim"),
            ),
        ]
        for case in cases
            mktempdir() do directory
                message = fig2_argument_error_message() do
                    run_claimed_dimension_case(directory; case.kwargs...)
                end
                @testset "$(case.label)" begin
                    for fragment in case.fragments
                        @test occursin(fragment, lowercase(message))
                    end
                end
            end
        end
    end

    function run_bond_convergence_case(
        directory;
        dimensions,
        fluxes,
        endpoint_pumps=Dict(dimension => 1.0 for dimension in dimensions),
        endpoint_energies=Dict(dimension => -1.0 for dimension in dimensions),
        transfer_residuals=Dict(dimension => 0.0 for dimension in dimensions),
        endpoint_sector_shifts=Dict(dimension => 1 for dimension in dimensions),
        early_plateau=false,
        early_plateau_overrides=Dict{String,Any}(),
    )
        spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
        curve_rule = spec.data["paper_curve_rule"]
        curve_fluxes = Float64.(curve_rule["flux_units_2pi"])
        curve_values = Float64.(curve_rule["delta_q"])
        curve_value(phi_y) = begin
            flux_units = phi_y / (2pi)
            index = findfirst(value -> isapprox(
                value, flux_units; atol=1e-10, rtol=0
            ), curve_fluxes)
            isnothing(index) && error("synthetic convergence flux is off the paper grid")
            return curve_values[index]
        end
        endpoint_sector_weights(polarization, physical_shift) = begin
            value = Float64(polarization)
            shift = Int(physical_shift)
            delta = value - shift
            iszero(delta) && return Dict(
                FIG2_SYNTHETIC_CHARGE_SCALE * shift => 1.0
            )
            remote_charge_distance = 1_000_000
            remote_weight = abs(delta) / remote_charge_distance
            @assert 0 < remote_weight < 1
            direction = delta > 0 ? 1 : -1
            remote_physical_charge =
                shift + direction * remote_charge_distance
            return Dict(
                FIG2_SYNTHETIC_CHARGE_SCALE * shift => 1 - remote_weight,
                FIG2_SYNTHETIC_CHARGE_SCALE * remote_physical_charge =>
                    remote_weight,
            )
        end
        operations = Fig2BenchmarkOperations(
            progress_audit=synthetic_progress_audit,
            provenance=synthetic_fig2_provenance,
            candidate_ids=(args...) -> ["only"],
            run_candidate=(spec, dimension, point, phi_y, candidate_id,
                           previous_state, candidate_directory) -> begin
                raw_polarization =
                    endpoint_pumps[dimension] * curve_value(phi_y)
                sectors = if point == length(fluxes)
                    endpoint_sector_weights(
                        raw_polarization,
                        endpoint_sector_shifts[dimension],
                    )
                else
                    synthetic_sector_weights_for_polarization(
                        raw_polarization
                    )
                end
                write_fake_candidate_files(
                    candidate_directory,
                    "D=$dimension,point=$point";
                    requested_maxdim=dimension,
                    summary_energy_per_site=endpoint_energies[dimension],
                    sector_weights=sectors,
                )
                first_point = point == 1
                fidelity = first_point ? NaN : 0.99
                momentum = synthetic_counting_momentum(;
                    transfer_residual=transfer_residuals[dimension]
                )
                counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                    spec, momentum
                )
                Fig2CandidateEvidence(;
                    state=(; dimension, point),
                    requested_maxdim=dimension,
                    achieved_maxlinkdim=dimension,
                    converged=true,
                    valid=true,
                    energy_per_site=endpoint_energies[dimension],
                    raw_schmidt_polarization=raw_polarization,
                    sector_weights=sectors,
                    fidelity_to_previous=fidelity,
                    fidelity_valid=!first_point,
                    mixed_fidelity=first_point ? nothing :
                        synthetic_mixed_fidelity(fidelity),
                    momentum,
                    momentum_counting=counting.counts,
                    momentum_counting_evidence=counting,
                    restart_valid=true,
                )
            end,
        )
        run_fig2_benchmark(
            spec,
            directory;
            stage="bond_convergence",
            dimensions,
            fluxes,
            operations,
        )
        if early_plateau
            evidence = Dict{String,Any}(
                "format" => "fqahc_fig2_early_plateau_v1",
                "manifest_sha256" => spec.sha256,
                "omit_dimension" => 3000,
                "lower_dimension" => 1000,
                "upper_dimension" => 2000,
                "reason" =>
                    "D=1000 and D=2000 independently satisfy every fixed acceptance gate",
            )
            merge!(evidence, early_plateau_overrides)
            open(joinpath(directory, "early_plateau.toml"), "w") do io
                TOML.print(io, evidence; sorted=true)
            end
        end
        return write_synthetic_fig2_acceptance_report!(
            spec,
            directory;
            candidate_ids_provider=operations.candidate_ids,
        )
    end

    full_curve_fluxes = 2pi .* collect(0.0:0.125:3.0)

    @testset "bond convergence requires declared progression and D3000" begin
        mktempdir() do directory
            report = run_bond_convergence_case(
                directory;
                dimensions=[1000, 2000],
                fluxes=[0.0, 6pi],
            )
            @test !report["bond_convergence"]["pass"]
            @test occursin(
                "256", lowercase(report["bond_convergence"]["reason"])
            )
        end
        mktempdir() do directory
            report = run_bond_convergence_case(
                directory;
                dimensions=[256, 512, 1000, 2000],
                fluxes=[0.0, 6pi],
            )
            @test !report["bond_convergence"]["pass"]
            @test occursin(
                "3000", lowercase(report["bond_convergence"]["reason"])
            )
            @test occursin(
                "early_plateau", lowercase(report["bond_convergence"]["reason"])
            )
        end
    end

    @testset "bond convergence rejects omitted D256 or D512" begin
        mktempdir() do directory
            report = run_bond_convergence_case(
                directory;
                dimensions=[512, 1000, 2000, 3000],
                fluxes=[0.0, 6pi],
            )
            @test !report["bond_convergence"]["pass"]
            @test occursin("256", report["bond_convergence"]["reason"])
        end
        mktempdir() do directory
            report = run_bond_convergence_case(
                directory;
                dimensions=[256, 1000, 2000, 3000],
                fluxes=[0.0, 6pi],
            )
            @test !report["bond_convergence"]["pass"]
            @test occursin("512", report["bond_convergence"]["reason"])
        end
    end

    @testset "D3000 convergence compares D2000 and D3000" begin
        mktempdir() do directory
            dimensions = [256, 512, 1000, 2000, 3000]
            energies = Dict(dimension => -1.0 for dimension in dimensions)
            energies[3000] = -1.0 + 5e-7
            report = run_bond_convergence_case(
                directory;
                dimensions,
                fluxes=full_curve_fluxes,
                endpoint_energies=energies,
            )
            @test report["bond_convergence"]["pass"]
            @test report["bond_convergence"]["measured"]["lower_dimension"] == 2000
            @test report["bond_convergence"]["measured"]["upper_dimension"] == 3000
            @test isfile(joinpath(
                directory, "paper_curve_comparison_D2000.toml"
            ))
            @test isfile(joinpath(
                directory, "paper_curve_comparison_D3000.toml"
            ))
        end
    end

    @testset "early plateau is explicit and recomputes both paper curves" begin
        mktempdir() do directory
            dimensions = [256, 512, 1000, 2000]
            energies = Dict(dimension => -1.0 for dimension in dimensions)
            energies[2000] = -1.0 + 5e-7
            report = run_bond_convergence_case(
                directory;
                dimensions,
                fluxes=full_curve_fluxes,
                endpoint_pumps=Dict(
                    256 => 1.0, 512 => 1.0, 1000 => 1.0, 2000 => 1.01
                ),
                endpoint_energies=energies,
                early_plateau=true,
            )
            @test report["bond_convergence"]["pass"]
            @test report["bond_convergence"]["evidence"] == "early_plateau.toml"
            @test report["bond_convergence"]["measured"]["lower_dimension"] == 1000
            @test report["bond_convergence"]["measured"]["upper_dimension"] == 2000
            @test isfile(joinpath(
                directory, "paper_curve_comparison_D1000.toml"
            ))
            @test isfile(joinpath(
                directory, "paper_curve_comparison_D2000.toml"
            ))
            @test haskey(
                report["bond_convergence"]["measured"], "early_plateau_sha256"
            )
            @test haskey(
                report["bond_convergence"]["measured"], "ledger_sha256"
            )
        end
        for pumps in (
            Dict(256 => 1.0, 512 => 1.0, 1000 => 1.03, 2000 => 1.015),
            Dict(256 => 1.0, 512 => 1.0, 1000 => 1.015, 2000 => 1.03),
        )
            mktempdir() do directory
                dimensions = [256, 512, 1000, 2000]
                report = run_bond_convergence_case(
                    directory;
                    dimensions,
                    fluxes=full_curve_fluxes,
                    endpoint_pumps=pumps,
                    early_plateau=true,
                )
                @test !report["bond_convergence"]["pass"]
                @test occursin(
                    "curve", lowercase(report["bond_convergence"]["reason"])
                )
            end
        end
        for overrides in (
            Dict{String,Any}("manifest_sha256" => repeat("0", 64)),
            Dict{String,Any}("reason" => ""),
        )
            mktempdir() do directory
                report = run_bond_convergence_case(
                    directory;
                    dimensions=[256, 512, 1000, 2000],
                    fluxes=[0.0, 6pi],
                    early_plateau=true,
                    early_plateau_overrides=overrides,
                )
                @test !report["bond_convergence"]["pass"]
                @test occursin(
                    haskey(overrides, "manifest_sha256") ? "manifest" : "reason",
                    lowercase(report["bond_convergence"]["reason"]),
                )
            end
        end
    end

    @testset "bond convergence requires complete high-D flux endpoints" begin
        mktempdir() do directory
            dimensions = [256, 512, 1000, 2000, 3000]
            report = run_bond_convergence_case(
                directory;
                dimensions,
                fluxes=[0.0],
            )
            @test !report["bond_convergence"]["pass"]
            @test occursin(
                "endpoint", lowercase(report["bond_convergence"]["reason"])
            )
        end
    end

    @testset "bond convergence gates pump, sector, and momentum residuals" begin
        mktempdir() do directory
            dimensions = [256, 512, 1000, 2000, 3000]
            pumps = Dict(dimension => 1.0 for dimension in dimensions)
            pumps[3000] = 1.03
            report = run_bond_convergence_case(
                directory;
                dimensions,
                fluxes=[0.0, 6pi],
                endpoint_pumps=pumps,
            )
            @test !report["bond_convergence"]["pass"]
            @test occursin(
                "pump", lowercase(report["bond_convergence"]["reason"])
            )
        end
        mktempdir() do directory
            dimensions = [256, 512, 1000, 2000, 3000]
            shifts = Dict(dimension => 1 for dimension in dimensions)
            shifts[3000] = 0
            report = run_bond_convergence_case(
                directory;
                dimensions,
                fluxes=[0.0, 6pi],
                endpoint_sector_shifts=shifts,
            )
            @test !report["bond_convergence"]["pass"]
            @test occursin(
                "sector", lowercase(report["bond_convergence"]["reason"])
            )
        end
        mktempdir() do directory
            dimensions = [256, 512, 1000, 2000, 3000]
            residuals = Dict(dimension => 0.0 for dimension in dimensions)
            residuals[3000] = 1.0e-3
            message = lowercase(fig2_argument_error_message() do
                run_bond_convergence_case(
                    directory;
                    dimensions,
                    fluxes=[0.0, 6pi],
                    transfer_residuals=residuals,
                )
            end)
            @test occursin("momentum", message)
            @test occursin("transfer", message)
            @test occursin("tolerance", message)
        end
    end

    @testset "Fig. 2 candidate completion rejects unsafe IDs and audits bad restart" begin
        @test_throws ArgumentError InfiniteCylinderDMRG._candidate_relative_directory(
            4, 1, "a/b"
        )
        @test_throws ArgumentError InfiniteCylinderDMRG._candidate_relative_directory(
            4, 1, "a?b"
        )

        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(args...) -> ["not_restartable"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    write_fake_candidate_files(
                        candidate_directory,
                        "restart=false";
                        requested_maxdim=dimension,
                        checkpoint_maxlinkdim=0,
                        sector_weights=Dict(0 => 1.0),
                    )
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; point),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        checkpoint_maxlinkdim=0,
                        converged=true,
                        valid=true,
                        energy_per_site=-1.0,
                        raw_schmidt_polarization=0.0,
                        sector_weights=Dict(0 => 1.0),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=false,
                    )
                end,
            )
            restart_error = try
                run_fig2_benchmark(
                    spec,
                    directory;
                    stage="restart_gate",
                    dimensions=[4],
                    fluxes=[0.0],
                    operations,
                )
                nothing
            catch error
                error
            end
            @test restart_error isa WorkflowValidationError
            ledger_path = joinpath(directory, "ledger.toml")
            @test isfile(ledger_path)
            metadata_path = joinpath(
                directory,
                "D_0004",
                "phi_0000",
                "candidate_not_restartable",
                "candidate.toml",
            )
            @test isfile(metadata_path)
            if isfile(ledger_path) && isfile(metadata_path)
                ledger = TOML.parsefile(ledger_path)
                @test length(ledger["candidate"]) == 1
                @test isempty(ledger["selection"])
                @test only(ledger["candidate"])["complete"]
                @test only(ledger["candidate"])["checkpoint_maxlinkdim"] == 0
                metadata = TOML.parsefile(metadata_path)
                @test metadata["restart_valid"] === false
                @test metadata["checkpoint_maxlinkdim"] == 0
                low_dimension_error = try
                    InfiniteCylinderDMRG._validate_persisted_fig2_candidate_files(
                        spec,
                        directory,
                        only(ledger["candidate"]),
                        synthetic_fig2_provenance(spec, directory, 0.0);
                        checkpoint_audit=(args...) -> (
                            restart_valid=true,
                            checkpoint_maxlinkdim=1,
                        ),
                        progress_audit=synthetic_progress_audit,
                    )
                    nothing
                catch error
                    error
                end
                @test low_dimension_error isa ArgumentError
                @test occursin(
                    "checkpoint",
                    lowercase(sprint(showerror, low_dimension_error)),
                )
            end
            resumed_error = try
                run_fig2_benchmark(
                    spec,
                    directory;
                    stage="restart_gate",
                    dimensions=[4],
                    fluxes=[0.0],
                    operations,
                )
                nothing
            catch error
                error
            end
            @test resumed_error isa WorkflowValidationError
            resumed_ledger = TOML.parsefile(ledger_path)
            @test length(resumed_ledger["candidate"]) == 1
            InfiniteCylinderDMRG._write_fig2_provenance!(
                spec,
                directory,
                time_ns(),
                synthetic_fig2_provenance(spec, directory, 0.0),
                resumed_ledger,
            )
            incomplete_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test all(!row["pass"] for row in values(incomplete_report))
            @test all(
                occursin("selection replay", lowercase(row["reason"])) for
                row in values(incomplete_report)
            )
        end
    end

    @testset "complete failed candidates are retained but never selected" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            calls = String[]
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(args...) -> ["failed_low_energy", "good"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    push!(calls, candidate_id)
                    failed = candidate_id == "failed_low_energy"
                    achieved = failed ? 2 : dimension
                    converged = !failed
                    sectors = Dict(0 => 1.0)
                    write_fake_candidate_files(
                        candidate_directory,
                        candidate_id;
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=achieved,
                        summary_converged=converged,
                        summary_energy_per_site=failed ? -2.0 : -1.0,
                        sector_weights=sectors,
                        convergence_rows=converged ? [
                            (
                                stage=1,
                                iteration=1,
                                maxlinkdim=achieved,
                                delta_energy="missing",
                                converged=false,
                            ),
                            (
                                stage=1,
                                iteration=2,
                                maxlinkdim=achieved,
                                delta_energy=1.0e-8,
                                converged=false,
                            ),
                            (
                                stage=1,
                                iteration=3,
                                maxlinkdim=achieved,
                                delta_energy=1.0e-8,
                                converged=true,
                            ),
                        ] : [(
                            stage=1,
                            iteration=1,
                            maxlinkdim=achieved,
                            delta_energy="missing",
                            converged=false,
                        )],
                    )
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; candidate_id),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=achieved,
                        converged,
                        valid=converged,
                        energy_per_site=failed ? -2.0 : -1.0,
                        raw_schmidt_polarization=0.0,
                        sector_weights=sectors,
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                        reason=failed ?
                            "synthetic iteration limit" : "synthetic converged",
                    )
                end,
            )

            run = nothing
            run_error = try
                run = run_fig2_benchmark(
                    spec,
                    directory;
                    stage="retain_failed",
                    dimensions=[4],
                    fluxes=[0.0],
                    operations,
                )
                nothing
            catch error
                error
            end
            @test isnothing(run_error)
            ledger_path = joinpath(directory, "ledger.toml")
            if isnothing(run_error)
                @test calls == ["failed_low_energy", "good"]
                ledger = TOML.parsefile(ledger_path)
                @test length(ledger["candidate"]) == 2
                @test all(row["complete"] for row in ledger["candidate"])
                @test only(ledger["selection"])["candidate_id"] == "good"
                failed_row = only(filter(
                    row -> row["candidate_id"] == "failed_low_energy",
                    ledger["candidate"],
                ))
                @test failed_row["requested_maxdim"] == 4
                @test failed_row["achieved_maxlinkdim"] == 2
                failed_metadata = TOML.parsefile(joinpath(
                    directory, failed_row["directory"], "candidate.toml"
                ))
                @test failed_metadata["requested_maxdim"] == 4
                @test failed_metadata["achieved_maxlinkdim"] == 2
                @test failed_metadata["converged"] === false
                @test failed_metadata["valid"] === false
                @test occursin(
                    "synthetic iteration limit", failed_metadata["reason"]
                )
                failed_candidate = InfiniteCylinderDMRG._persisted_fig2_candidate(
                    spec,
                    directory,
                    failed_row,
                    4,
                    1,
                    0.0,
                    "failed_low_energy",
                )
                @test !InfiniteCylinderDMRG._fig2_candidate_selection_eligible(
                    failed_candidate
                )

                resumed = run_fig2_benchmark(
                    spec,
                    directory;
                    stage="retain_failed",
                    dimensions=[4],
                    fluxes=[0.0],
                    operations,
                )
                @test calls == ["failed_low_energy", "good"]
                @test only(resumed.selections).candidate_id == "good"
                integrity = synthetic_fig2_acceptance_integrity_inputs(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                @test only(integrity.selections).candidate_id == "good"
            else
                @test run_error isa ArgumentError
                @test occursin(
                    "achieved maxlinkdim",
                    lowercase(sprint(showerror, run_error)),
                )
                @test !isfile(ledger_path)
            end
        end
    end

    @testset "segmented resume replays the persisted winner before loading" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            candidate_calls = NamedTuple[]
            load_calls = NamedTuple[]
            candidate_ids = (spec, dimension, point, previous) ->
                point == 1 ? ["alpha", "beta"] : ["warm", "cold"]
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=candidate_ids,
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    push!(candidate_calls, (; point, candidate_id, previous_state))
                    first_point = point == 1
                    energy = first_point ?
                        (candidate_id == "alpha" ? -1.0 : -1.1) : -1.0
                    sectors = Dict(0 => 1.0)
                    write_fake_candidate_files(
                        candidate_directory,
                        "segmented-point=$point,candidate=$candidate_id";
                        requested_maxdim=dimension,
                        summary_energy_per_site=energy,
                        sector_weights=sectors,
                    )
                    fidelity = first_point ? NaN :
                        (candidate_id == "warm" ? 0.9 : 0.8)
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; point, candidate_id),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=energy,
                        raw_schmidt_polarization=0.0,
                        sector_weights=sectors,
                        fidelity_to_previous=fidelity,
                        fidelity_valid=!first_point,
                        mixed_fidelity=first_point ? nothing :
                            synthetic_mixed_fidelity(fidelity),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                    )
                end,
                load_state=(spec, dimension, point, selection, checkpoint) -> begin
                    push!(load_calls, (; point, candidate_id=selection.candidate_id))
                    return (; point, candidate_id=selection.candidate_id)
                end,
            )
            first_run = run_fig2_benchmark(
                spec,
                directory;
                stage="segmented_resume",
                dimensions=[4],
                fluxes=[0.0],
                operations,
            )
            @test only(first_run.selections).candidate_id == "beta"

            ledger_path = joinpath(directory, "ledger.toml")
            ledger = TOML.parsefile(ledger_path)
            alpha_row = only(filter(
                row -> String(row["candidate_id"]) == "alpha",
                ledger["candidate"],
            ))
            alpha = InfiniteCylinderDMRG._persisted_fig2_candidate(
                spec,
                directory,
                alpha_row,
                4,
                1,
                0.0,
                "alpha",
            )
            ledger["selection"][1] = InfiniteCylinderDMRG._selection_dict(
                InfiniteCylinderDMRG._fig2_selection_from_candidate(
                    alpha,
                    4,
                    1,
                    0.0,
                    "minimum_zero_flux_energy",
                )
            )
            open(ledger_path, "w") do io
                TOML.print(io, ledger; sorted=true)
            end
            empty!(candidate_calls)
            empty!(load_calls)

            resume_error = try
                run_fig2_benchmark(
                    spec,
                    directory;
                    stage="segmented_resume",
                    dimensions=[4],
                    fluxes=[0.0, 6pi],
                    operations,
                )
                nothing
            catch error
                error
            end
            @test resume_error isa ArgumentError
            resume_message = isnothing(resume_error) ? "" :
                lowercase(sprint(showerror, resume_error))
            @test occursin("selection replay winner mismatch", resume_message)
            @test isempty(load_calls)
            @test isempty(candidate_calls)
        end
    end

    @testset "atomic candidate ledger preserves raw pump evidence" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            function assert_integrity_failure(report)
                @test sort!(collect(keys(report))) ==
                    sort!(copy(spec.acceptance_rows))
                @test all(!row["pass"] for row in values(report))
                @test all(
                    row["evidence"] == "acceptance_integrity_failure.toml" for
                    row in values(report)
                )
                @test all(
                    occursin("integrity", lowercase(row["reason"])) for
                    row in values(report)
                )
                evidence = TOML.parsefile(joinpath(
                    directory, "acceptance_integrity_failure.toml"
                ))
                @test evidence["pass"] === false
                @test evidence["manifest_sha256"] == spec.sha256
                return nothing
            end
            calls = NamedTuple[]
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(spec, dimension, point, previous) ->
                    point == 1 ? ["alpha", "beta"] : ["warm", "cold"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    push!(calls, (; dimension, point, phi_y, candidate_id,
                                  previous_state))
                    first_point = point == 1
                    energy = first_point ?
                        (candidate_id == "alpha" ? -1.0 : -1.1) :
                        (candidate_id == "warm" ? -1.08 : -1.2)
                    fidelity = first_point ? NaN :
                        (candidate_id == "warm" ? 0.97 : 0.65)
                    raw_polarization = first_point ?
                        (candidate_id == "alpha" ? 0.20 : 0.25) :
                        (candidate_id == "warm" ? 1.25 : 9.0)
                    sectors = synthetic_sector_weights_for_polarization(
                        raw_polarization
                    )
                    write_fake_candidate_files(
                        candidate_directory,
                        "D=$dimension,point=$point,candidate=$candidate_id";
                        requested_maxdim=dimension,
                        summary_energy_per_site=energy,
                        sector_weights=sectors,
                    )
                    momentum_valid = !(point == 2 && candidate_id == "warm")
                    momentum = momentum_valid ?
                        synthetic_counting_momentum() : synthetic_momentum(; valid=false)
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; point, candidate_id),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=energy,
                        raw_schmidt_polarization=raw_polarization,
                        sector_weights=sectors,
                        fidelity_to_previous=fidelity,
                        fidelity_valid=!first_point,
                        mixed_fidelity=first_point ? nothing :
                            synthetic_mixed_fidelity(fidelity),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                        reason="synthetic candidate",
                    )
                end,
            )
            run = run_fig2_benchmark(
                spec,
                directory;
                stage="test",
                dimensions=[4],
                fluxes=[0.0, 6pi],
                operations,
            )

            @test run isa Fig2BenchmarkRun
            provenance_path = joinpath(directory, "provenance.toml")
            @test isfile(provenance_path)
            if isfile(provenance_path)
                provenance = TOML.parsefile(provenance_path)
                @test provenance["format"] == "fqahc_fig2_provenance_v2"
                @test provenance["manifest_sha256"] == spec.sha256
                @test occursin(r"^[0-9a-f]{40}$", provenance["git_commit"])
                @test provenance["git_tree_clean"]
                @test provenance["julia_version"] == "1.12.5"
                @test !isempty(provenance["pbs_job_id"])
                @test provenance["threads"] == Threads.nthreads()
                @test provenance["blas_threads"] == 1
                @test provenance["strided_threads"] == 1
                @test provenance["blocksparse_threaded"] ===
                    (Threads.nthreads() > 1)
                @test provenance["runtime_seconds"] >= 0
                @test occursin(
                    r"^[0-9a-f]{64}$",
                    provenance["project_manifest_sha256"],
                )
                @test occursin(
                    r"^[0-9a-f]{64}$",
                    provenance["benchmark_source_sha256"],
                )
            end
            @test [row.candidate_id for row in run.selections] == ["beta", "warm"]
            @test [row.selection_basis for row in run.selections] ==
                ["minimum_zero_flux_energy", "maximum_mixed_fidelity"]
            @test length(calls) == 4
            @test all(isnothing(call.previous_state) for call in calls[1:2])
            @test all(call.previous_state == (; point=1, candidate_id="beta") for call in calls[3:4])

            ledger = TOML.parsefile(joinpath(directory, "ledger.toml"))
            @test ledger["format"] == "fqahc_fig2_ledger_v3"
            @test length(ledger["candidate"]) == 4
            @test length(ledger["selection"]) == 2
            @test all(row["complete"] for row in ledger["candidate"])
            @test [row["candidate_id"] for row in ledger["selection"]] ==
                ["beta", "warm"]
            @test all(row["state_sha256"] != "" for row in ledger["candidate"])
            @test all(
                row["progress_complete"] === true for row in ledger["candidate"]
            )
            @test all(
                row["progress_event_count"] == 3 for row in ledger["candidate"]
            )
            @test all(
                row["progress_resume_count"] == 0 for row in ledger["candidate"]
            )
            @test all(
                row["progress_latest_maxlinkdim"] == row["dimension"] for
                row in ledger["candidate"]
            )
            @test all(
                row["progress_final_state_sha256"] == row["state_sha256"] for
                row in ledger["candidate"]
            )
            @test all(
                row["progress_sha256"] == row["checksums"]["progress.toml"] for
                row in ledger["candidate"]
            )

            warm_row = only(filter(
                row -> row["point"] == 2 && row["candidate_id"] == "warm",
                ledger["candidate"],
            ))
            warm_directory = joinpath(directory, warm_row["directory"])
            @test first(readlines(joinpath(warm_directory, "mixed_fidelity.tsv"))) ==
                "fidelity_to_previous\tmixed_value_real\tmixed_value_imag\tself_value1_real\tself_value1_imag\tself_value2_real\tself_value2_imag\tmixed_residual\tself_residual1\tself_residual2\tconverged\tvalid\treason\treference_present\treference_dimension\treference_point\treference_candidate_id\treference_directory\treference_state_sha256"
            @test first(readlines(joinpath(
                warm_directory, "momentum_entanglement_spectrum.tsv"
            ))) ==
                "level\tprobability\tentanglement_energy\traw_charge\tky_index\tky\teigenphase\tphase_residual\tunitary_residual\tpower_residual\tcommutator_residual\ttransfer_residual\tvalidation_valid\tvalidation_reason\tmomentum_valid\tmomentum_reason\tmomentum_convention\tsector_ky_offsets"
            warm_metadata = TOML.parsefile(joinpath(warm_directory, "candidate.toml"))
            for key in (
                "mixed_value_real",
                "mixed_value_imag",
                "self_value1_real",
                "self_value1_imag",
                "self_value2_real",
                "self_value2_imag",
                "mixed_residual",
                "self_residual1",
                "self_residual2",
                "mixed_converged",
                "mixed_valid",
                "mixed_reason",
                "momentum_unitary_residual",
                "momentum_power_residual",
                "momentum_commutator_residual",
                "momentum_transfer_residual",
                "momentum_validation_valid",
                "momentum_validation_reason",
                "momentum_convention",
                "momentum_sector_ky_offsets",
                "momentum_counting_valid",
                "momentum_counting_reason",
                "momentum_counting_delta_ky",
                "momentum_counting_cutoffs",
                "momentum_counting_has_reference",
                "momentum_counting_reference_raw_charge",
                "momentum_counting_reference_level",
                "momentum_counting_reference_ky_index",
                "momentum_counting_minimum_cutoff_margin",
            )
                @test haskey(warm_metadata, key)
            end
            @test warm_metadata["mixed_residual"] == 1.0e-10
            @test warm_metadata["self_residual1"] == 2.0e-10
            @test warm_metadata["self_residual2"] == 3.0e-10

            first_selection = first(ledger["selection"])
            for key in (
                "momentum_unitary_residual",
                "momentum_power_residual",
                "momentum_commutator_residual",
                "momentum_transfer_residual",
                "momentum_validation_valid",
                "momentum_validation_reason",
                "momentum_valid",
                "momentum_reason",
                "momentum_convention",
                "momentum_sector_ky_offsets",
                "momentum_counting_valid",
                "momentum_counting_reason",
                "momentum_counting_delta_ky",
                "momentum_counting_cutoffs",
                "momentum_counting_has_reference",
                "momentum_counting_reference_raw_charge",
                "momentum_counting_reference_level",
                "momentum_counting_reference_ky_index",
                "momentum_counting_minimum_cutoff_margin",
            )
                @test haskey(first_selection, key)
            end
            @test first_selection["momentum_unitary_residual"] == 0.0
            @test first_selection["momentum_sector_ky_offsets"] ==
                Dict("3" => 0, "4" => 0)

            for candidate in ledger["candidate"]
                candidate_directory = joinpath(directory, candidate["directory"])
                for filename in (
                    "state.h5",
                    "summary.toml",
                    "convergence.tsv",
                    "expansion.tsv",
                    "density.tsv",
                    "entanglement_spectrum.tsv",
                    "schmidt_sectors.tsv",
                    "progress.toml",
                    "momentum_entanglement_spectrum.tsv",
                    "mixed_fidelity.tsv",
                    "candidate.toml",
                )
                    @test isfile(joinpath(candidate_directory, filename))
                end
                metadata = TOML.parsefile(joinpath(
                    candidate_directory, "candidate.toml"
                ))
                for key in (
                    "progress_complete",
                    "progress_event_count",
                    "progress_resume_count",
                    "progress_latest_maxlinkdim",
                    "progress_final_state_sha256",
                    "progress_sha256",
                )
                    @test metadata[key] == candidate[key]
                end
                @test metadata["progress_sha256"] ==
                    candidate["checksums"]["progress.toml"]
                @test metadata["progress_final_state_sha256"] ==
                    candidate["checksums"]["state.h5"]
            end

            raw_path = joinpath(directory, "pump_raw.tsv")
            tracked_path = joinpath(directory, "pump_sector_tracked.tsv")
            raw_before = read(raw_path)
            report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            raw_after = read(raw_path)
            @test raw_after == raw_before
            raw_text = String(raw_after)
            @test occursin("0.25", raw_text)
            @test occursin("1.25", raw_text)
            @test !occursin("unwrap", lowercase(raw_text))
            @test report["endpoint_pump"]["pass"]
            @test report["sector_shift"]["pass"]
            @test !report["momentum_counting"]["pass"]
            @test occursin(
                "invalid momentum",
                lowercase(report["momentum_counting"]["reason"]),
            )
            @test report["paper_curve"]["pass"]

            FIG2_INTERRUPT_PROBES[:acceptance] = try
                write_fig2_acceptance_report!(
                    spec,
                    directory;
                    checkpoint_audit=(args...) -> throw(InterruptException()),
                    progress_audit=synthetic_progress_audit,
                    candidate_ids_provider=operations.candidate_ids,
                )
            catch error
                error
            end
            @testset "Fig. 2 workflow propagates InterruptException" begin
                @test FIG2_INTERRUPT_PROBES[:checkpoint_restart] isa
                    InterruptException
                @test FIG2_INTERRUPT_PROBES[:checkpoint_audit_load] isa
                    InterruptException
                @test FIG2_INTERRUPT_PROBES[:checkpoint_audit_measure] isa
                    InterruptException
                @test FIG2_INTERRUPT_PROBES[:integrity_digest] isa
                    InterruptException
                @test FIG2_INTERRUPT_PROBES[:integrity_digest_ordinary] ==
                    "unreadable"
                @test FIG2_INTERRUPT_PROBES[:acceptance] isa InterruptException
            end

            @testset "Fig. 2 validation wrappers preserve interrupts" begin
                contexts = (
                    "candidate sector metadata is incomplete or invalid",
                    "candidate identity fields are missing",
                    "candidate audit metadata is incomplete or invalid",
                    "persisted selection audit fields are incomplete",
                    "selection replay failed at dimension 4 point 1",
                )
                for context in contexts
                    interrupt = try
                        InfiniteCylinderDMRG._fig2_argument_error(context) do
                            throw(InterruptException())
                        end
                    catch error
                        error
                    end
                    @test interrupt isa InterruptException

                    ordinary = try
                        InfiniteCylinderDMRG._fig2_argument_error(context) do
                            error("synthetic validation failure")
                        end
                    catch error
                        error
                    end
                    @test ordinary isa ArgumentError
                end
            end

            curve_path = joinpath(directory, "paper_curve_comparison.toml")
            malicious_curve = Dict(
                "format" => "fqahc_fig2_paper_curve_comparison_v1",
                "source_pdf_sha256" =>
                    "23f0691570b906ae7a9394cef2e3dd6eebd4e36d79e882a4bbd2f17824cc89a0",
                "source_page" => 3,
                "source_panel" => "d",
                "maximum_residual" => 0.20,
                "tolerance" => 0.05,
                "pass" => true,
                "reason" => "untrusted external pass flag",
            )
            open(curve_path, "w") do io
                TOML.print(io, malicious_curve; sorted=true)
            end
            malicious_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test malicious_report["paper_curve"]["pass"]
            @test malicious_report["paper_curve"]["measured"] < 0.02
            @test malicious_report["paper_curve"] == report["paper_curve"]
            @test occursin(
                "raw",
                lowercase(malicious_report["paper_curve"]["reason"]),
            )
            regenerated_curve = TOML.parsefile(curve_path)
            @test regenerated_curve["pass"]
            @test regenerated_curve["source_vector_path_sha256"] ==
                spec.data["paper_curve_rule"]["source_vector_path_sha256"]
            @test regenerated_curve["pump_raw_sha256"] ==
                InfiniteCylinderDMRG._fig2_file_sha256(raw_path)
            @test regenerated_curve["pump_sector_tracked_sha256"] ==
                InfiniteCylinderDMRG._fig2_file_sha256(tracked_path)
            @test regenerated_curve["sampling_policy"] ==
                "linear_interpolation_on_validated_raw_grid_without_extrapolation"

            valid_curve = copy(malicious_curve)
            valid_curve["maximum_residual"] = 0.01
            valid_curve["pass"] = false
            open(curve_path, "w") do io
                TOML.print(io, valid_curve; sorted=true)
            end
            valid_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test valid_report["paper_curve"]["pass"]
            @test valid_report["paper_curve"] == report["paper_curve"]
            @test !report["bond_convergence"]["pass"]
            @test !report["restart"]["pass"]
            @test report["restart"]["evidence"] == "restart_gate.toml"
            @test occursin(
                "missing", lowercase(report["restart"]["reason"])
            )
            write_synthetic_restart_gate!(spec, directory)
            restart_valid_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test restart_valid_report["restart"]["pass"]
            @test restart_valid_report["restart"]["evidence"] ==
                "restart_gate.toml"
            @test occursin(
                "two-process",
                lowercase(restart_valid_report["restart"]["reason"]),
            )
            write_synthetic_restart_gate!(
                spec,
                directory;
                overrides=Dict{String,Any}(
                    "next_flux_iteration_completed" => false,
                ),
            )
            restart_tampered_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test !restart_tampered_report["restart"]["pass"]
            @test occursin(
                "next_flux_iteration_completed",
                restart_tampered_report["restart"]["reason"],
            )
            write_synthetic_restart_gate!(spec, directory)
            @test report["provenance"]["pass"]
            provenance_bytes = read(provenance_path)
            tampered_provenance = TOML.parsefile(provenance_path)
            tampered_provenance["git_tree_clean"] = false
            open(provenance_path, "w") do io
                TOML.print(io, tampered_provenance; sorted=true)
            end
            tampered_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            assert_integrity_failure(tampered_report)
            @test occursin(
                "provenance",
                lowercase(tampered_report["provenance"]["reason"]),
            )
            open(provenance_path, "w") do io
                write(io, provenance_bytes)
            end
            tracked_path = joinpath(directory, "pump_sector_tracked.tsv")
            tracked_before = read(tracked_path)
            selected_keys = Set(
                (Int(row["dimension"]), Int(row["point"]), String(row["candidate_id"]))
                for row in ledger["selection"]
            )
            unselected = first(filter(
                row -> (
                    Int(row["dimension"]),
                    Int(row["point"]),
                    String(row["candidate_id"]),
                ) ∉ selected_keys,
                ledger["candidate"],
            ))
            state_path = joinpath(
                directory, unselected["directory"], "state.h5"
            )
            state_bytes = read(state_path)
            open(state_path, "a") do io
                write(io, "tampered state")
            end
            state_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            assert_integrity_failure(state_report)
            @test occursin(
                "checksum", lowercase(state_report["restart"]["reason"])
            )
            open(state_path, "w") do io
                write(io, state_bytes)
            end

            ledger_path = joinpath(directory, "ledger.toml")
            ledger_bytes = read(ledger_path)
            progress_candidate = unselected
            progress_directory = joinpath(
                directory, progress_candidate["directory"]
            )
            progress_path = joinpath(progress_directory, "progress.toml")
            progress_metadata_path = joinpath(
                progress_directory, "candidate.toml"
            )
            progress_bytes = read(progress_path)
            progress_metadata_bytes = read(progress_metadata_path)
            forged_final_state_sha256 = repeat("0", 64)
            @test forged_final_state_sha256 !=
                progress_candidate["state_sha256"]

            forged_progress = TOML.parse(String(copy(progress_bytes)))
            forged_progress["final_state_sha256"] =
                forged_final_state_sha256
            open(progress_path, "w") do io
                TOML.print(io, forged_progress; sorted=true)
            end
            forged_progress_sha256 =
                InfiniteCylinderDMRG._fig2_file_sha256(progress_path)

            forged_metadata = TOML.parse(
                String(copy(progress_metadata_bytes))
            )
            forged_metadata["progress_final_state_sha256"] =
                forged_final_state_sha256
            forged_metadata["progress_sha256"] = forged_progress_sha256
            open(progress_metadata_path, "w") do io
                TOML.print(io, forged_metadata; sorted=true)
            end
            forged_metadata_sha256 =
                InfiniteCylinderDMRG._fig2_file_sha256(progress_metadata_path)

            forged_ledger = TOML.parse(String(copy(ledger_bytes)))
            forged_row = only(filter(
                row -> Int(row["dimension"]) ==
                        Int(progress_candidate["dimension"]) &&
                    Int(row["point"]) == Int(progress_candidate["point"]) &&
                    String(row["candidate_id"]) ==
                        String(progress_candidate["candidate_id"]),
                forged_ledger["candidate"],
            ))
            forged_row["progress_final_state_sha256"] =
                forged_final_state_sha256
            forged_row["progress_sha256"] = forged_progress_sha256
            forged_row["checksums"]["progress.toml"] =
                forged_progress_sha256
            forged_row["checksums"]["candidate.toml"] =
                forged_metadata_sha256
            open(ledger_path, "w") do io
                TOML.print(io, forged_ledger; sorted=true)
            end

            forged_progress_report =
                write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
            assert_integrity_failure(forged_progress_report)
            @test occursin(
                "final-state",
                lowercase(forged_progress_report["restart"]["reason"]),
            )
            open(progress_path, "w") do io
                write(io, progress_bytes)
            end
            open(progress_metadata_path, "w") do io
                write(io, progress_metadata_bytes)
            end
            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end

            fixed_candidate_ids = (spec, dimension, point, previous_state) ->
                point == 1 ? ["alpha", "beta"] : ["warm", "cold"]

            deleted_candidate_ledger = TOML.parse(String(copy(ledger_bytes)))
            filter!(
                row -> !(Int(row["dimension"]) == 4 &&
                    Int(row["point"]) == 1 &&
                    String(row["candidate_id"]) == "alpha"),
                deleted_candidate_ledger["candidate"],
            )
            open(ledger_path, "w") do io
                TOML.print(io, deleted_candidate_ledger; sorted=true)
            end
            deleted_candidate_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            assert_integrity_failure(deleted_candidate_report)
            deleted_candidate_message = fig2_argument_error_message() do
                InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                    spec,
                    directory;
                    checkpoint_audit=synthetic_checkpoint_audit,
                    progress_audit=synthetic_progress_audit,
                    candidate_ids_provider=fixed_candidate_ids,
                )
            end
            @test occursin("candidate", lowercase(deleted_candidate_message))
            @test occursin("set", lowercase(deleted_candidate_message))
            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end

            extra_candidate_ledger = TOML.parse(String(copy(ledger_bytes)))
            source_candidate = only(filter(
                row -> Int(row["dimension"]) == 4 &&
                    Int(row["point"]) == 1 &&
                    String(row["candidate_id"]) == "alpha",
                extra_candidate_ledger["candidate"],
            ))
            extra_candidate = deepcopy(source_candidate)
            extra_candidate_id = "extra_without_winner_change"
            extra_relative_directory =
                InfiniteCylinderDMRG._candidate_relative_directory(
                    4, 1, extra_candidate_id
                )
            source_directory = joinpath(directory, source_candidate["directory"])
            extra_directory = joinpath(directory, extra_relative_directory)
            cp(source_directory, extra_directory; force=true)
            extra_metadata_path = joinpath(extra_directory, "candidate.toml")
            extra_metadata = TOML.parsefile(extra_metadata_path)
            extra_metadata["candidate_id"] = extra_candidate_id
            open(extra_metadata_path, "w") do io
                TOML.print(io, extra_metadata; sorted=true)
            end
            extra_checksums = Dict(
                filename => InfiniteCylinderDMRG._fig2_file_sha256(
                    joinpath(extra_directory, filename)
                ) for filename in InfiniteCylinderDMRG.FIG2_REQUIRED_CANDIDATE_FILES
            )
            extra_candidate["candidate_id"] = extra_candidate_id
            extra_candidate["directory"] = extra_relative_directory
            extra_candidate["checksums"] = extra_checksums
            extra_candidate["state_sha256"] = extra_checksums["state.h5"]
            push!(extra_candidate_ledger["candidate"], extra_candidate)
            open(ledger_path, "w") do io
                TOML.print(io, extra_candidate_ledger; sorted=true)
            end
            extra_candidate_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            assert_integrity_failure(extra_candidate_report)
            extra_candidate_message = fig2_argument_error_message() do
                InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                    spec,
                    directory;
                    checkpoint_audit=synthetic_checkpoint_audit,
                    progress_audit=synthetic_progress_audit,
                    candidate_ids_provider=fixed_candidate_ids,
                )
            end
            @test occursin("candidate", lowercase(extra_candidate_message))
            @test occursin("set", lowercase(extra_candidate_message))
            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end

            legacy_ledger = TOML.parse(String(copy(ledger_bytes)))
            legacy_ledger["format"] = "fqahc_fig2_ledger_v2"
            open(ledger_path, "w") do io
                TOML.print(io, legacy_ledger; sorted=true)
            end
            legacy_error = try
                run_fig2_benchmark(
                    spec,
                    directory;
                    stage="test",
                    dimensions=[4],
                    fluxes=[0.0, 6pi],
                    operations,
                )
                nothing
            catch error
                error
            end
            @test legacy_error isa ArgumentError
            legacy_message = isnothing(legacy_error) ? "" :
                lowercase(sprint(showerror, legacy_error))
            @test occursin(
                "ledger format",
                legacy_message,
            )
            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end
            selected_row = first(ledger["selection"])
            selected_candidate = only(filter(
                row ->
                    Int(row["dimension"]) == Int(selected_row["dimension"]) &&
                    Int(row["point"]) == Int(selected_row["point"]) &&
                    String(row["candidate_id"]) ==
                        String(selected_row["candidate_id"]),
                ledger["candidate"],
            ))

            @testset "Schmidt replay is pinned to the physical cut" begin
                candidate_directory = joinpath(
                    directory, selected_candidate["directory"]
                )
                sector_path = joinpath(
                    candidate_directory, "schmidt_sectors.tsv"
                )
                spectrum_path = joinpath(
                    candidate_directory, "entanglement_spectrum.tsv"
                )
                metadata_path = joinpath(candidate_directory, "candidate.toml")
                sector_bytes = read(sector_path)
                spectrum_bytes = read(spectrum_path)
                metadata_bytes = read(metadata_path)
                wrong_weights = synthetic_sector_weights_for_polarization(0.5)
                wrong_pair = (2, 2 * spec.config.Ly)

                sector_lines = readlines(sector_path)
                forged_sector_lines = String[first(sector_lines)]
                for line in sector_lines[2:end]
                    fields = split(line, '\t'; keepempty=true)
                    charge = parse(Int, fields[4])
                    fields[1] = string(wrong_pair[1])
                    fields[2] = string(wrong_pair[2])
                    fields[6] = string(wrong_weights[charge])
                    push!(forged_sector_lines, join(fields, '\t'))
                end
                append!(forged_sector_lines, sector_lines[2:end])
                open(sector_path, "w") do io
                    foreach(line -> println(io, line), forged_sector_lines)
                end

                spectrum_lines = readlines(spectrum_path)
                forged_spectrum_lines = String[first(spectrum_lines)]
                for line in spectrum_lines[2:end]
                    fields = split(line, '\t'; keepempty=true)
                    charge = parse(Int, fields[8])
                    probability = wrong_weights[charge]
                    fields[1] = string(wrong_pair[1])
                    fields[2] = string(wrong_pair[2])
                    fields[4] = string(sqrt(probability))
                    fields[5] = string(probability)
                    fields[6] = string(-log(probability))
                    push!(forged_spectrum_lines, join(fields, '\t'))
                end
                append!(forged_spectrum_lines, spectrum_lines[2:end])
                open(spectrum_path, "w") do io
                    foreach(line -> println(io, line), forged_spectrum_lines)
                end

                metadata = TOML.parsefile(metadata_path)
                metadata["raw_schmidt_polarization"] = 0.5
                metadata["sector_weights"] = Dict(
                    string(charge) => weight for
                    (charge, weight) in wrong_weights
                )
                open(metadata_path, "w") do io
                    TOML.print(io, metadata; sorted=true)
                end

                forged = TOML.parse(String(copy(ledger_bytes)))
                forged_candidate = only(filter(
                    row -> String(row["directory"]) ==
                        String(selected_candidate["directory"]),
                    forged["candidate"],
                ))
                for (filename, path) in (
                    "schmidt_sectors.tsv" => sector_path,
                    "entanglement_spectrum.tsv" => spectrum_path,
                    "candidate.toml" => metadata_path,
                )
                    forged_candidate["checksums"][filename] =
                        InfiniteCylinderDMRG._fig2_file_sha256(path)
                end
                forged_selection = only(filter(
                    row -> Int(row["dimension"]) ==
                            Int(selected_candidate["dimension"]) &&
                        Int(row["point"]) == Int(selected_candidate["point"]),
                    forged["selection"],
                ))
                forged_selection["raw_schmidt_polarization"] = 0.5
                forged_selection["sector_weights"] = Dict(
                    string(charge) => weight for
                    (charge, weight) in wrong_weights
                )
                open(ledger_path, "w") do io
                    TOML.print(io, forged; sorted=true)
                end

                replay_error = try
                    InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                        spec,
                        directory;
                        checkpoint_audit=synthetic_checkpoint_audit,
                        progress_audit=synthetic_progress_audit,
                        candidate_ids_provider=operations.candidate_ids,
                    )
                    nothing
                catch error
                    error
                end
                @test replay_error isa ArgumentError
                replay_message = isnothing(replay_error) ? "" :
                    lowercase(sprint(showerror, replay_error))
                @test occursin("cut", replay_message)
                @test occursin("bond", replay_message)
                forged_report = write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                @test all(
                    row["evidence"] == "acceptance_integrity_failure.toml" for
                    row in values(forged_report)
                )

                open(sector_path, "w") do io
                    write(io, sector_bytes)
                end
                open(spectrum_path, "w") do io
                    write(io, spectrum_bytes)
                end
                open(metadata_path, "w") do io
                    write(io, metadata_bytes)
                end
                open(ledger_path, "w") do io
                    write(io, ledger_bytes)
                end
                open(raw_path, "w") do io
                    write(io, raw_before)
                end
                open(tracked_path, "w") do io
                    write(io, tracked_before)
                end
            end

            @testset "candidate energy is replayed from convergence evidence" begin
                candidate_directory = joinpath(
                    directory, selected_candidate["directory"]
                )
                summary_path = joinpath(candidate_directory, "summary.toml")
                metadata_path = joinpath(candidate_directory, "candidate.toml")
                summary_bytes = read(summary_path)
                metadata_bytes = read(metadata_path)
                forged_energy = -2.0

                summary = TOML.parsefile(summary_path)
                summary["energy"]["per_site"] = forged_energy
                open(summary_path, "w") do io
                    TOML.print(io, summary; sorted=true)
                end
                metadata = TOML.parsefile(metadata_path)
                metadata["energy_per_site"] = forged_energy
                open(metadata_path, "w") do io
                    TOML.print(io, metadata; sorted=true)
                end

                forged = TOML.parse(String(copy(ledger_bytes)))
                forged_candidate = only(filter(
                    row -> String(row["directory"]) ==
                        String(selected_candidate["directory"]),
                    forged["candidate"],
                ))
                forged_candidate["checksums"]["summary.toml"] =
                    InfiniteCylinderDMRG._fig2_file_sha256(summary_path)
                forged_candidate["checksums"]["candidate.toml"] =
                    InfiniteCylinderDMRG._fig2_file_sha256(metadata_path)
                forged_selection = only(filter(
                    row -> Int(row["dimension"]) ==
                            Int(selected_candidate["dimension"]) &&
                        Int(row["point"]) == Int(selected_candidate["point"]),
                    forged["selection"],
                ))
                forged_selection["energy_per_site"] = forged_energy
                open(ledger_path, "w") do io
                    TOML.print(io, forged; sorted=true)
                end

                replay_error = try
                    InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                        spec,
                        directory;
                        checkpoint_audit=synthetic_checkpoint_audit,
                        progress_audit=synthetic_progress_audit,
                        candidate_ids_provider=operations.candidate_ids,
                    )
                    nothing
                catch error
                    error
                end
                @test replay_error isa ArgumentError
                replay_message = isnothing(replay_error) ? "" :
                    lowercase(sprint(showerror, replay_error))
                @test occursin("energy", replay_message)
                @test occursin("convergence", replay_message)
                forged_report = write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                @test all(
                    row["evidence"] == "acceptance_integrity_failure.toml" for
                    row in values(forged_report)
                )

                open(summary_path, "w") do io
                    write(io, summary_bytes)
                end
                open(metadata_path, "w") do io
                    write(io, metadata_bytes)
                end
                open(ledger_path, "w") do io
                    write(io, ledger_bytes)
                end
                open(raw_path, "w") do io
                    write(io, raw_before)
                end
                open(tracked_path, "w") do io
                    write(io, tracked_before)
                end
            end

            @testset "coordinated final raw energy forgery is rejected" begin
                candidate_directory = joinpath(
                    directory, selected_candidate["directory"]
                )
                convergence_path = joinpath(
                    candidate_directory, "convergence.tsv"
                )
                summary_path = joinpath(candidate_directory, "summary.toml")
                metadata_path = joinpath(candidate_directory, "candidate.toml")
                convergence_bytes = read(convergence_path)
                summary_bytes = read(summary_path)
                metadata_bytes = read(metadata_path)
                forged_energy = -2.0

                convergence_lines = readlines(convergence_path)
                final_fields = split(
                    last(convergence_lines), '\t'; keepempty=true
                )
                forged_raw_energy =
                    forged_energy * sites_per_cell(spec.config)
                final_fields[4] = string(forged_raw_energy)
                final_fields[5] = string(forged_raw_energy)
                convergence_lines[end] = join(final_fields, '\t')
                open(convergence_path, "w") do io
                    foreach(line -> println(io, line), convergence_lines)
                end
                summary = TOML.parsefile(summary_path)
                summary["energy"]["per_site"] = forged_energy
                open(summary_path, "w") do io
                    TOML.print(io, summary; sorted=true)
                end
                metadata = TOML.parsefile(metadata_path)
                metadata["energy_per_site"] = forged_energy
                open(metadata_path, "w") do io
                    TOML.print(io, metadata; sorted=true)
                end

                forged = TOML.parse(String(copy(ledger_bytes)))
                forged_candidate = only(filter(
                    row -> String(row["directory"]) ==
                        String(selected_candidate["directory"]),
                    forged["candidate"],
                ))
                for (filename, path) in (
                    "convergence.tsv" => convergence_path,
                    "summary.toml" => summary_path,
                    "candidate.toml" => metadata_path,
                )
                    forged_candidate["checksums"][filename] =
                        InfiniteCylinderDMRG._fig2_file_sha256(path)
                end
                forged_selection = only(filter(
                    row -> Int(row["dimension"]) ==
                            Int(selected_candidate["dimension"]) &&
                        Int(row["point"]) == Int(selected_candidate["point"]),
                    forged["selection"],
                ))
                forged_selection["energy_per_site"] = forged_energy
                open(ledger_path, "w") do io
                    TOML.print(io, forged; sorted=true)
                end

                replay_error = try
                    InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                        spec,
                        directory;
                        checkpoint_audit=synthetic_checkpoint_audit,
                        progress_audit=synthetic_progress_audit,
                        candidate_ids_provider=operations.candidate_ids,
                    )
                    nothing
                catch error
                    error
                end
                @test replay_error isa ArgumentError
                replay_message = isnothing(replay_error) ? "" :
                    lowercase(sprint(showerror, replay_error))
                @test occursin("convergence", replay_message)
                @test occursin("delta_energy", replay_message)
                forged_report = write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                @test all(
                    row["evidence"] == "acceptance_integrity_failure.toml" for
                    row in values(forged_report)
                )

                open(convergence_path, "w") do io
                    write(io, convergence_bytes)
                end
                open(summary_path, "w") do io
                    write(io, summary_bytes)
                end
                open(metadata_path, "w") do io
                    write(io, metadata_bytes)
                end
                open(ledger_path, "w") do io
                    write(io, ledger_bytes)
                end
                open(raw_path, "w") do io
                    write(io, raw_before)
                end
                open(tracked_path, "w") do io
                    write(io, tracked_before)
                end
            end

            @testset "charge entanglement energy is replayed from probability" begin
                candidate_directory = joinpath(
                    directory, selected_candidate["directory"]
                )
                spectrum_path = joinpath(
                    candidate_directory, "entanglement_spectrum.tsv"
                )
                spectrum_bytes = read(spectrum_path)
                spectrum_lines = readlines(spectrum_path)
                fields = split(spectrum_lines[2], '\t'; keepempty=true)
                @test parse(Float64, fields[5]) > 0
                fields[6] = string(parse(Float64, fields[6]) + 0.5)
                spectrum_lines[2] = join(fields, '\t')
                open(spectrum_path, "w") do io
                    foreach(line -> println(io, line), spectrum_lines)
                end

                forged = TOML.parse(String(copy(ledger_bytes)))
                forged_candidate = only(filter(
                    row -> String(row["directory"]) ==
                        String(selected_candidate["directory"]),
                    forged["candidate"],
                ))
                forged_candidate["checksums"]["entanglement_spectrum.tsv"] =
                    InfiniteCylinderDMRG._fig2_file_sha256(spectrum_path)
                open(ledger_path, "w") do io
                    TOML.print(io, forged; sorted=true)
                end

                replay_error = try
                    InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                        spec,
                        directory;
                        checkpoint_audit=synthetic_checkpoint_audit,
                        progress_audit=synthetic_progress_audit,
                        candidate_ids_provider=operations.candidate_ids,
                    )
                    nothing
                catch error
                    error
                end
                @test replay_error isa ArgumentError
                replay_message = isnothing(replay_error) ? "" :
                    lowercase(sprint(showerror, replay_error))
                @test occursin("entanglement_energy", replay_message)
                @test occursin("probability", replay_message)
                forged_report = write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                @test all(
                    row["evidence"] == "acceptance_integrity_failure.toml" for
                    row in values(forged_report)
                )

                open(spectrum_path, "w") do io
                    write(io, spectrum_bytes)
                end
                open(ledger_path, "w") do io
                    write(io, ledger_bytes)
                end
                open(raw_path, "w") do io
                    write(io, raw_before)
                end
                open(tracked_path, "w") do io
                    write(io, tracked_before)
                end
            end

            @testset "mixed fidelity is bound to the prior selected checkpoint" begin
                reference_columns = [
                    "reference_present",
                    "reference_dimension",
                    "reference_point",
                    "reference_candidate_id",
                    "reference_directory",
                    "reference_state_sha256",
                ]
                metadata_keys = ["mixed_$column" for column in reference_columns]
                warm_candidate = only(filter(
                    row -> Int(row["point"]) == 2 &&
                        String(row["candidate_id"]) == "warm",
                    ledger["candidate"],
                ))
                alpha_candidate = only(filter(
                    row -> Int(row["point"]) == 1 &&
                        String(row["candidate_id"]) == "alpha",
                    ledger["candidate"],
                ))
                warm_directory = joinpath(
                    directory, warm_candidate["directory"]
                )
                mixed_path = joinpath(warm_directory, "mixed_fidelity.tsv")
                metadata_path = joinpath(warm_directory, "candidate.toml")
                mixed_bytes = read(mixed_path)
                metadata_bytes = read(metadata_path)
                mixed_lines = readlines(mixed_path)
                header_fields = split(first(mixed_lines), '\t'; keepempty=true)
                @test length(header_fields) == 19
                @test length(header_fields) >= length(reference_columns) &&
                    header_fields[(end - length(reference_columns) + 1):end] ==
                        reference_columns

                warm_metadata = TOML.parsefile(metadata_path)
                @test all(haskey(warm_metadata, key) for key in metadata_keys)
                point1_directory = joinpath(
                    directory,
                    selected_candidate["directory"],
                )
                point1_metadata = TOML.parsefile(joinpath(
                    point1_directory, "candidate.toml"
                ))
                @test all(haskey(point1_metadata, key) for key in metadata_keys)
                @test (
                    point1_metadata["mixed_reference_present"],
                    point1_metadata["mixed_reference_dimension"],
                    point1_metadata["mixed_reference_point"],
                    point1_metadata["mixed_reference_candidate_id"],
                    point1_metadata["mixed_reference_directory"],
                    point1_metadata["mixed_reference_state_sha256"],
                ) == (false, 0, 0, "absent", "absent", "absent")
                point1_mixed_fields = split(
                    readlines(joinpath(point1_directory, "mixed_fidelity.tsv"))[2],
                    '\t';
                    keepempty=true,
                )
                @test point1_mixed_fields[14:19] ==
                    ["false", "0", "0", "absent", "absent", "absent"]
                @test String(selected_candidate["candidate_id"]) == "beta"
                @test (
                    warm_metadata["mixed_reference_present"],
                    warm_metadata["mixed_reference_dimension"],
                    warm_metadata["mixed_reference_point"],
                    warm_metadata["mixed_reference_candidate_id"],
                    warm_metadata["mixed_reference_directory"],
                    warm_metadata["mixed_reference_state_sha256"],
                ) == (
                    true,
                    4,
                    1,
                    "beta",
                    String(selected_candidate["directory"]),
                    String(selected_candidate["state_sha256"]),
                )

                wrong_reference = Dict(
                    "reference_present" => true,
                    "reference_dimension" => Int(alpha_candidate["dimension"]),
                    "reference_point" => Int(alpha_candidate["point"]),
                    "reference_candidate_id" =>
                        String(alpha_candidate["candidate_id"]),
                    "reference_directory" => String(alpha_candidate["directory"]),
                    "reference_state_sha256" =>
                        String(alpha_candidate["state_sha256"]),
                )
                for column in reference_columns
                    warm_metadata["mixed_$column"] = wrong_reference[column]
                end
                open(metadata_path, "w") do io
                    TOML.print(io, warm_metadata; sorted=true)
                end

                if all(column in header_fields for column in reference_columns)
                    indices = Dict(
                        column => only(findall(==(column), header_fields)) for
                        column in reference_columns
                    )
                    fields = split(mixed_lines[2], '\t'; keepempty=true)
                    for column in reference_columns
                        fields[indices[column]] = string(wrong_reference[column])
                    end
                    fields[indices["reference_present"]] = "true"
                    mixed_lines[2] = join(fields, '\t')
                    open(mixed_path, "w") do io
                        foreach(line -> println(io, line), mixed_lines)
                    end
                end

                forged = TOML.parse(String(copy(ledger_bytes)))
                forged_candidate = only(filter(
                    row -> String(row["directory"]) ==
                        String(warm_candidate["directory"]),
                    forged["candidate"],
                ))
                forged_candidate["checksums"]["mixed_fidelity.tsv"] =
                    InfiniteCylinderDMRG._fig2_file_sha256(mixed_path)
                forged_candidate["checksums"]["candidate.toml"] =
                    InfiniteCylinderDMRG._fig2_file_sha256(metadata_path)
                open(ledger_path, "w") do io
                    TOML.print(io, forged; sorted=true)
                end

                replay_error = try
                    InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                        spec,
                        directory;
                        checkpoint_audit=synthetic_checkpoint_audit,
                        progress_audit=synthetic_progress_audit,
                        candidate_ids_provider=operations.candidate_ids,
                    )
                    nothing
                catch error
                    error
                end
                @test replay_error isa ArgumentError
                replay_message = isnothing(replay_error) ? "" :
                    lowercase(sprint(showerror, replay_error))
                @test occursin("reference", replay_message)
                @test occursin("selected", replay_message)
                forged_report = write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                @test all(
                    row["evidence"] == "acceptance_integrity_failure.toml" for
                    row in values(forged_report)
                )

                open(mixed_path, "w") do io
                    write(io, mixed_bytes)
                end
                open(metadata_path, "w") do io
                    write(io, metadata_bytes)
                end
                open(ledger_path, "w") do io
                    write(io, ledger_bytes)
                end
                open(raw_path, "w") do io
                    write(io, raw_before)
                end
                open(tracked_path, "w") do io
                    write(io, tracked_before)
                end
            end

            @testset "transfer residual replay uses the producer tolerance" begin
                candidate_directory = joinpath(
                    directory, selected_candidate["directory"]
                )
                momentum_path = joinpath(
                    candidate_directory,
                    "momentum_entanglement_spectrum.tsv",
                )
                metadata_path = joinpath(candidate_directory, "candidate.toml")
                momentum_bytes = read(momentum_path)
                metadata_bytes = read(metadata_path)

                momentum_lines = readlines(momentum_path)
                for line_index in 2:length(momentum_lines)
                    fields = split(
                        momentum_lines[line_index], '\t'; keepempty=true
                    )
                    fields[12] = "5.0e-7"
                    momentum_lines[line_index] = join(fields, '\t')
                end
                open(momentum_path, "w") do io
                    foreach(line -> println(io, line), momentum_lines)
                end

                metadata = TOML.parsefile(metadata_path)
                metadata["momentum_transfer_residual"] = 5.0e-7
                open(metadata_path, "w") do io
                    TOML.print(io, metadata; sorted=true)
                end

                forged = TOML.parse(String(copy(ledger_bytes)))
                forged_candidate = only(filter(
                    row -> String(row["directory"]) ==
                        String(selected_candidate["directory"]),
                    forged["candidate"],
                ))
                forged_candidate["checksums"][
                    "momentum_entanglement_spectrum.tsv"
                ] = InfiniteCylinderDMRG._fig2_file_sha256(momentum_path)
                forged_candidate["checksums"]["candidate.toml"] =
                    InfiniteCylinderDMRG._fig2_file_sha256(metadata_path)
                forged_selection = only(filter(
                    row -> Int(row["dimension"]) ==
                            Int(selected_candidate["dimension"]) &&
                        Int(row["point"]) == Int(selected_candidate["point"]),
                    forged["selection"],
                ))
                forged_selection["momentum_transfer_residual"] = 5.0e-7
                open(ledger_path, "w") do io
                    TOML.print(io, forged; sorted=true)
                end

                forged_endpoint =
                    InfiniteCylinderDMRG._selection_from_dict(forged_selection)
                @test !InfiniteCylinderDMRG._fig2_momentum_endpoint_valid(
                    spec, forged_endpoint
                )

                replay_error = try
                    InfiniteCylinderDMRG._fig2_acceptance_integrity_inputs(
                        spec,
                        directory;
                        checkpoint_audit=synthetic_checkpoint_audit,
                        progress_audit=synthetic_progress_audit,
                        candidate_ids_provider=operations.candidate_ids,
                    )
                    nothing
                catch error
                    error
                end
                @test replay_error isa ArgumentError
                replay_message = isnothing(replay_error) ? "" :
                    lowercase(sprint(showerror, replay_error))
                @test occursin("transfer", replay_message)
                @test occursin("tolerance", replay_message)

                forged_report = write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                @test all(!row["pass"] for row in values(forged_report))
                @test all(
                    row["evidence"] == "acceptance_integrity_failure.toml" for
                    row in values(forged_report)
                )

                open(momentum_path, "w") do io
                    write(io, momentum_bytes)
                end
                open(metadata_path, "w") do io
                    write(io, metadata_bytes)
                end
                open(ledger_path, "w") do io
                    write(io, ledger_bytes)
                end
            end

            physical_charge_error = mktempdir() do scratch
                candidate_directory = joinpath(scratch, "candidate")
                cp(
                    joinpath(directory, selected_candidate["directory"]),
                    candidate_directory;
                    force=true,
                )
                spectrum_path = joinpath(
                    candidate_directory, "entanglement_spectrum.tsv"
                )
                lines = readlines(spectrum_path)
                fields = split(lines[2], '\t'; keepempty=true)
                fields[9] = string(parse(Float64, fields[9]) + 0.25)
                lines[2] = join(fields, '\t')
                open(spectrum_path, "w") do io
                    foreach(line -> println(io, line), lines)
                end
                metadata = TOML.parsefile(joinpath(
                    candidate_directory, "candidate.toml"
                ))
                try
                    InfiniteCylinderDMRG._validate_fig2_candidate_artifacts(
                        spec,
                        candidate_directory,
                        Int(selected_candidate["dimension"]),
                        Int(selected_candidate["point"]),
                        Float64(selected_candidate["phi_y"]),
                        String(selected_candidate["candidate_id"]),
                        metadata["generation_provenance"],
                        ;
                        progress_audit=synthetic_progress_audit,
                    )
                    nothing
                catch error
                    error
                end
            end
            @test physical_charge_error isa ArgumentError
            physical_charge_message = isnothing(physical_charge_error) ? "" :
                lowercase(sprint(showerror, physical_charge_error))
            @test occursin("physical_charge", physical_charge_message)
            @test occursin("entanglement_spectrum", physical_charge_message)

            function forge_candidate_table(tamper::Function, filename)
                table_path = joinpath(
                    directory, selected_candidate["directory"], filename
                )
                table_bytes = read(table_path)
                lines = readlines(table_path)
                fields = split(lines[2], '\t'; keepempty=true)
                tamper(fields)
                lines[2] = join(fields, '\t')
                open(table_path, "w") do io
                    foreach(line -> println(io, line), lines)
                end

                forged = TOML.parse(String(copy(ledger_bytes)))
                ledger_candidate = only(filter(
                    row -> String(row["directory"]) ==
                        String(selected_candidate["directory"]),
                    forged["candidate"],
                ))
                ledger_candidate["checksums"][filename] =
                    InfiniteCylinderDMRG._fig2_file_sha256(table_path)
                open(ledger_path, "w") do io
                    TOML.print(io, forged; sorted=true)
                end
                forged_report = write_synthetic_fig2_acceptance_report!(
                    spec,
                    directory;
                    candidate_ids_provider=operations.candidate_ids,
                )
                assert_integrity_failure(forged_report)
                @test all(
                    occursin(filename, lowercase(row["reason"])) for
                    row in values(forged_report)
                )

                open(table_path, "w") do io
                    write(io, table_bytes)
                end
                open(ledger_path, "w") do io
                    write(io, ledger_bytes)
                end
                return nothing
            end

            forge_candidate_table("schmidt_sectors.tsv") do fields
                fields[4] = "99"
                fields[6] = "0.6"
            end
            forge_candidate_table("entanglement_spectrum.tsv") do fields
                fields[5] = "0.125"
                fields[8] = "99"
            end

            candidate_path = joinpath(
                directory, selected_candidate["directory"], "candidate.toml"
            )
            candidate_bytes = read(candidate_path)
            current_generation_provenance = synthetic_fig2_provenance(
                spec, directory, 0.0
            )

            checkpoint_path = joinpath(
                directory, selected_candidate["directory"], "state.h5"
            )
            checkpoint_bytes = read(checkpoint_path)
            progress_path = joinpath(
                directory, selected_candidate["directory"], "progress.toml"
            )
            progress_bytes = read(progress_path)
            open(checkpoint_path, "w") do io
                write(io, "corrupt-but-checksummed-checkpoint\n")
            end
            checkpoint_forged_row = deepcopy(selected_candidate)
            checkpoint_sha256 =
                InfiniteCylinderDMRG._fig2_file_sha256(checkpoint_path)
            checkpoint_progress = TOML.parsefile(progress_path)
            checkpoint_progress["final_state_sha256"] = checkpoint_sha256
            open(progress_path, "w") do io
                TOML.print(io, checkpoint_progress; sorted=true)
            end
            checkpoint_progress_sha256 =
                InfiniteCylinderDMRG._fig2_file_sha256(progress_path)
            checkpoint_metadata = TOML.parsefile(candidate_path)
            checkpoint_metadata["progress_final_state_sha256"] =
                checkpoint_sha256
            checkpoint_metadata["progress_sha256"] =
                checkpoint_progress_sha256
            open(candidate_path, "w") do io
                TOML.print(io, checkpoint_metadata; sorted=true)
            end
            checkpoint_forged_row["checksums"]["state.h5"] =
                checkpoint_sha256
            checkpoint_forged_row["state_sha256"] = checkpoint_sha256
            checkpoint_forged_row["checksums"]["progress.toml"] =
                checkpoint_progress_sha256
            checkpoint_forged_row["progress_sha256"] =
                checkpoint_progress_sha256
            checkpoint_forged_row["progress_final_state_sha256"] =
                checkpoint_sha256
            checkpoint_forged_row["checksums"]["candidate.toml"] =
                InfiniteCylinderDMRG._fig2_file_sha256(candidate_path)
            checkpoint_forgery_error = try
                InfiniteCylinderDMRG._validate_persisted_fig2_candidate_files(
                    spec,
                    directory,
                    checkpoint_forged_row,
                    current_generation_provenance,
                    ;
                    progress_audit=synthetic_progress_audit,
                )
                nothing
            catch error
                error
            end
            @test checkpoint_forgery_error isa ArgumentError
            @test occursin(
                "checkpoint",
                lowercase(sprint(showerror, checkpoint_forgery_error)),
            )
            open(checkpoint_path, "w") do io
                write(io, checkpoint_bytes)
            end
            open(progress_path, "w") do io
                write(io, progress_bytes)
            end
            open(candidate_path, "w") do io
                write(io, candidate_bytes)
            end

            low_dimension_error = try
                InfiniteCylinderDMRG._validate_persisted_fig2_candidate_files(
                    spec,
                    directory,
                    selected_candidate,
                    current_generation_provenance;
                    checkpoint_audit=(args...) -> (
                        restart_valid=true,
                        checkpoint_maxlinkdim=1,
                    ),
                    progress_audit=synthetic_progress_audit,
                )
                nothing
            catch error
                error
            end
            @test low_dimension_error isa ArgumentError
            @test occursin(
                "checkpoint",
                lowercase(sprint(showerror, low_dimension_error)),
            )

            forged_checkpoint_metadata = TOML.parsefile(candidate_path)
            forged_checkpoint_metadata["checkpoint_maxlinkdim"] -= 1
            open(candidate_path, "w") do io
                TOML.print(io, forged_checkpoint_metadata; sorted=true)
            end
            metadata_forged_row = deepcopy(selected_candidate)
            metadata_forged_row["checksums"]["candidate.toml"] =
                InfiniteCylinderDMRG._fig2_file_sha256(candidate_path)
            metadata_checkpoint_error = fig2_argument_error_message() do
                InfiniteCylinderDMRG._validate_persisted_fig2_candidate_files(
                    spec,
                    directory,
                    metadata_forged_row,
                    current_generation_provenance,
                    ;
                    progress_audit=synthetic_progress_audit,
                )
            end
            @test occursin("checkpoint", lowercase(metadata_checkpoint_error))
            open(candidate_path, "w") do io
                write(io, candidate_bytes)
            end

            ledger_forged_row = deepcopy(selected_candidate)
            ledger_forged_row["checkpoint_maxlinkdim"] -= 1
            ledger_checkpoint_error = fig2_argument_error_message() do
                InfiniteCylinderDMRG._validate_persisted_fig2_candidate_files(
                    spec,
                    directory,
                    ledger_forged_row,
                    current_generation_provenance,
                    ;
                    progress_audit=synthetic_progress_audit,
                )
            end
            @test occursin("checkpoint", lowercase(ledger_checkpoint_error))

            candidate_metadata = TOML.parsefile(candidate_path)
            candidate_metadata["raw_schmidt_polarization"] += 0.5
            open(candidate_path, "w") do io
                TOML.print(io, candidate_metadata; sorted=true)
            end
            forged_ledger = TOML.parsefile(ledger_path)
            forged_candidate = only(filter(
                row -> String(row["directory"]) ==
                    String(selected_candidate["directory"]),
                forged_ledger["candidate"],
            ))
            forged_candidate["checksums"]["candidate.toml"] =
                InfiniteCylinderDMRG._fig2_file_sha256(candidate_path)
            forged_selection = only(filter(
                row -> Int(row["dimension"]) ==
                        Int(selected_candidate["dimension"]) &&
                    Int(row["point"]) == Int(selected_candidate["point"]),
                forged_ledger["selection"],
            ))
            forged_selection["raw_schmidt_polarization"] += 0.5
            open(ledger_path, "w") do io
                TOML.print(io, forged_ledger; sorted=true)
            end
            metadata_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            assert_integrity_failure(metadata_report)
            @test all(
                occursin("polarization", lowercase(row["reason"])) for
                row in values(metadata_report)
            )
            @test all(
                occursin("entanglement_spectrum", lowercase(row["reason"])) for
                row in values(metadata_report)
            )
            open(candidate_path, "w") do io
                write(io, candidate_bytes)
            end
            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end
            open(raw_path, "w") do io
                write(io, raw_before)
            end
            open(tracked_path, "w") do io
                write(io, tracked_before)
            end

            function forge_selection_winner(candidate_id, point)
                forged = TOML.parse(String(copy(ledger_bytes)))
                candidate_row = only(filter(
                    row -> Int(row["dimension"]) == 4 &&
                        Int(row["point"]) == point &&
                        String(row["candidate_id"]) == candidate_id,
                    forged["candidate"],
                ))
                phi_y = Float64(candidate_row["phi_y"])
                candidate = InfiniteCylinderDMRG._persisted_fig2_candidate(
                    spec,
                    directory,
                    candidate_row,
                    4,
                    point,
                    phi_y,
                    candidate_id,
                )
                basis = point == 1 ?
                    "minimum_zero_flux_energy" : "maximum_mixed_fidelity"
                replacement = InfiniteCylinderDMRG._selection_dict(
                    InfiniteCylinderDMRG._fig2_selection_from_candidate(
                        candidate,
                        4,
                        point,
                        phi_y,
                        basis,
                    )
                )
                selection_index = only(findall(
                    row -> Int(row["dimension"]) == 4 &&
                        Int(row["point"]) == point,
                    forged["selection"],
                ))
                forged["selection"][selection_index] = replacement
                return forged
            end

            forged_zero_flux = forge_selection_winner("alpha", 1)
            open(ledger_path, "w") do io
                TOML.print(io, forged_zero_flux; sorted=true)
            end
            zero_flux_replay_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            assert_integrity_failure(zero_flux_replay_report)
            @test all(
                occursin(
                    "candidate mixed reference does not match the exact prior selected checkpoint",
                    lowercase(row["reason"]),
                ) for
                row in values(zero_flux_replay_report)
            )

            forged_continuation = forge_selection_winner("cold", 2)
            open(ledger_path, "w") do io
                TOML.print(io, forged_continuation; sorted=true)
            end
            continuation_replay_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            assert_integrity_failure(continuation_replay_report)
            @test all(
                occursin("selection replay", lowercase(row["reason"])) for
                row in values(continuation_replay_report)
            )
            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end

            forged_ledger = TOML.parse(String(copy(ledger_bytes)))
            forged_ledger["selection"][1]["candidate_id"] = "alpha"
            open(ledger_path, "w") do io
                TOML.print(io, forged_ledger; sorted=true)
            end
            association_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            assert_integrity_failure(association_report)
            @test occursin(
                "selection", lowercase(association_report["sector_shift"]["reason"])
            )
            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end

            open(raw_path, "w") do io
                write(io, "tampered raw pump table\n")
            end
            open(tracked_path, "w") do io
                write(io, "tampered tracked pump table\n")
            end
            rebuilt_report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test rebuilt_report["endpoint_pump"] == report["endpoint_pump"]
            @test rebuilt_report["sector_shift"] == report["sector_shift"]
            @test read(raw_path) == raw_before
            @test read(tracked_path) == tracked_before

            candidate_call_count = length(calls)
            fully_resumed = run_fig2_benchmark(
                spec,
                directory;
                stage="test",
                dimensions=[4],
                fluxes=[0.0, 6pi],
                operations,
            )
            @test length(fully_resumed.selections) == 2
            @test length(calls) == candidate_call_count
            @test read(raw_path) == raw_before
            @test read(tracked_path) == tracked_before
            @test sort!(collect(keys(report))) == sort!(TOML.parsefile(
                FIG2_MANIFEST_PATH
            )["acceptance"]["rows"])
            @test isfile(joinpath(directory, "acceptance.toml"))
            @test isfile(joinpath(directory, "pump_sector_tracked.tsv"))
            @test isempty(filter(
                name -> occursin(".tmp", name),
                readdir(directory),
            ))
        end
    end

    @testset "nonempty ledger resumes from the selected checkpoint" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            candidate_calls = NamedTuple[]
            load_calls = NamedTuple[]
            operations = Fig2BenchmarkOperations(
                checkpoint_audit=synthetic_checkpoint_audit,
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(spec, dimension, point, previous) ->
                    point == 1 ? ["alpha", "beta"] : ["warm", "cold"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    push!(candidate_calls, (; point, candidate_id, previous_state))
                    first_point = point == 1
                    warm = candidate_id == "warm"
                    energy = first_point ?
                        (candidate_id == "alpha" ? -1.0 : -1.1) :
                        (warm ? -1.08 : -1.2)
                    raw_polarization = first_point ?
                        (candidate_id == "alpha" ? 0.20 : 0.25) :
                        (warm ? 1.25 : 9.0)
                    sectors = synthetic_sector_weights_for_polarization(
                        raw_polarization
                    )
                    write_fake_candidate_files(
                        candidate_directory,
                        "point=$point,candidate=$candidate_id";
                        requested_maxdim=dimension,
                        summary_energy_per_site=energy,
                        sector_weights=sectors,
                    )
                    fidelity = first_point ? NaN : (warm ? 0.97 : 0.65)
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; source=:fresh, point, candidate_id),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=energy,
                        raw_schmidt_polarization=raw_polarization,
                        sector_weights=sectors,
                        fidelity_to_previous=fidelity,
                        fidelity_valid=!first_point,
                        mixed_fidelity=first_point ? nothing :
                            synthetic_mixed_fidelity(fidelity),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                    )
                end,
                load_state=(spec, dimension, point, selection, checkpoint) -> begin
                    push!(load_calls, (; dimension, point,
                                      candidate_id=selection.candidate_id,
                                      checkpoint))
                    return (; source=:checkpoint, point,
                            candidate_id=selection.candidate_id)
                end,
            )

            first_run = run_fig2_benchmark(
                spec,
                directory;
                stage="resume",
                dimensions=[4],
                fluxes=[0.0],
                operations,
            )
            @test [row.candidate_id for row in first_run.selections] == ["beta"]
            @test length(candidate_calls) == 2

            empty!(candidate_calls)
            resumed = run_fig2_benchmark(
                spec,
                directory;
                stage="resume",
                dimensions=[4],
                fluxes=[0.0, 6pi],
                operations,
            )
            @test [row.candidate_id for row in resumed.selections] == ["beta", "warm"]
            @test [(call.point, call.candidate_id) for call in candidate_calls] ==
                [(2, "warm"), (2, "cold")]
            @test all(
                call.previous_state ==
                    (; source=:checkpoint, point=1, candidate_id="beta")
                for call in candidate_calls
            )
            @test length(load_calls) == 1
            @test load_calls[1].candidate_id == "beta"
            @test endswith(load_calls[1].checkpoint, joinpath(
                "D_0004", "phi_0000", "candidate_beta", "state.h5"
            ))

            ledger = TOML.parsefile(joinpath(directory, "ledger.toml"))
            @test length(ledger["candidate"]) == 4
            @test length(ledger["selection"]) == 2

            ledger_path = joinpath(directory, "ledger.toml")
            ledger_bytes = read(ledger_path)
            ledger["selection"][1]["momentum_unitary_residual"] = 0.125
            open(ledger_path, "w") do io
                TOML.print(io, ledger; sorted=true)
            end
            selection_error = fig2_argument_error_message() do
                run_fig2_benchmark(
                    spec,
                    directory;
                    stage="resume",
                    dimensions=[4],
                    fluxes=[0.0, 6pi],
                    operations,
                )
            end
            @test occursin("selection momentum", lowercase(selection_error))
            open(ledger_path, "w") do io
                write(io, ledger_bytes)
            end
            ledger = TOML.parsefile(ledger_path)

            first_candidate = ledger["candidate"][1]
            first_directory = joinpath(directory, first_candidate["directory"])
            density_path = joinpath(first_directory, "density.tsv")
            density_bytes = read(density_path)
            density_lines = readlines(density_path)
            open(density_path, "w") do io
                println(io, "bad_density_header")
                foreach(line -> println(io, line), density_lines[2:end])
            end
            density_error = fig2_argument_error_message() do
                run_fig2_benchmark(
                    spec,
                    directory;
                    stage="resume",
                    dimensions=[4],
                    fluxes=[0.0, 6pi],
                    operations,
                )
            end
            @test occursin("density.tsv header", density_error)
            open(density_path, "w") do io
                write(io, density_bytes)
            end

            candidate_path = joinpath(first_directory, "candidate.toml")
            candidate_bytes = read(candidate_path)
            candidate_metadata = TOML.parsefile(candidate_path)
            candidate_metadata["candidate_id"] = "wrong_identity"
            open(candidate_path, "w") do io
                TOML.print(io, candidate_metadata; sorted=true)
            end
            identity_error = fig2_argument_error_message() do
                run_fig2_benchmark(
                    spec,
                    directory;
                    stage="resume",
                    dimensions=[4],
                    fluxes=[0.0, 6pi],
                    operations,
                )
            end
            @test occursin("candidate identity", identity_error)
            open(candidate_path, "w") do io
                write(io, candidate_bytes)
            end

            @test_throws ArgumentError run_fig2_benchmark(
                spec,
                directory;
                stage="resume",
                dimensions=[8],
                fluxes=[0.0, 6pi],
                operations,
            )

            tampered_state = joinpath(
                directory, ledger["candidate"][1]["directory"], "state.h5"
            )
            open(tampered_state, "a") do io
                write(io, "tampered\n")
            end
            @test_throws ArgumentError run_fig2_benchmark(
                spec,
                directory;
                stage="resume",
                dimensions=[4],
                fluxes=[0.0, 6pi],
                operations,
            )
        end
    end

    @testset "manifest copy is immutable across resume" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            operations = Fig2BenchmarkOperations(
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(args...) -> String[],
                run_candidate=(args...) -> error("must not run"),
            )
            run_fig2_benchmark(
                spec,
                directory;
                stage="empty",
                dimensions=Int[],
                fluxes=Float64[],
                operations,
            )
            copied = joinpath(directory, "manifest.toml")
            @test read(copied) == read(FIG2_MANIFEST_PATH)
            open(copied, "a") do io
                write(io, "\n# changed\n")
            end
            @test_throws ArgumentError run_fig2_benchmark(
                spec,
                directory;
                stage="empty",
                dimensions=Int[],
                fluxes=Float64[],
                operations,
            )
        end
    end

    @testset "Fig. 2 benchmark CLI parsing" begin
        settings = parse_fig2_benchmark_args([
            "--manifest=$FIG2_MANIFEST_PATH",
            "--stage=pilot_D32",
            "--output=out/fqahc_fig2",
            "--dimensions=32,64,128",
            "--flux_units_2pi=0,0.5,1,1.5,2,2.5,3",
            "--threads=24",
        ])
        @test settings isa Fig2BenchmarkSettings
        @test settings.manifest == FIG2_MANIFEST_PATH
        @test settings.stage == "pilot_D32"
        @test settings.output == "out/fqahc_fig2"
        @test settings.dimensions == [32, 64, 128]
        @test settings.fluxes == 2pi .* [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        @test settings.threads == 24

        valid = [
            "--manifest=$FIG2_MANIFEST_PATH",
            "--stage=pilot",
            "--output=out",
            "--dimensions=32,64",
            "--flux_units_2pi=0,3",
        ]
        @test_throws ArgumentError parse_fig2_benchmark_args(valid[2:end])
        @test_throws ArgumentError parse_fig2_benchmark_args([
            valid; "--unknown=value"
        ])
        @test_throws ArgumentError parse_fig2_benchmark_args(replace.(
            valid, "--dimensions=32,64" => "--dimensions=64,32"
        ))
        @test_throws ArgumentError parse_fig2_benchmark_args(replace.(
            valid, "--flux_units_2pi=0,3" => "--flux_units_2pi=0.5,3"
        ))
        @test_throws ArgumentError parse_fig2_benchmark_args(replace.(
            valid, "--flux_units_2pi=0,3" => "--flux_units_2pi=0,2.5"
        ))
        @test_throws ArgumentError parse_fig2_benchmark_args([
            valid; "--threads=0"
        ])
    end

    @testset "guarded Fig. 2 benchmark entry point" begin
        benchmark_bin = joinpath(
            @__DIR__, "..", "bin", "run_fig2_benchmark.jl"
        )
        @test include(benchmark_bin) === nothing
        @test isdefined(Main, :run_fig2_benchmark_main)
    end

    @testset "Fig. 2 benchmark entry point routes the audited workflow" begin
        mktempdir() do directory
            configured_threads = Int[]
            operations = Fig2BenchmarkOperations(
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(args...) -> ["only"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    raw_polarization = point - 1.0
                    sectors = synthetic_sector_weights_for_polarization(
                        raw_polarization
                    )
                    write_fake_candidate_files(
                        candidate_directory,
                        "point=$point";
                        requested_maxdim=dimension,
                        sector_weights=sectors,
                    )
                    fidelity = point == 1 ? NaN : 0.99
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; point),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=-1.0,
                        raw_schmidt_polarization=raw_polarization,
                        sector_weights=sectors,
                        fidelity_to_previous=fidelity,
                        fidelity_valid=point > 1,
                        mixed_fidelity=point == 1 ? nothing :
                            synthetic_mixed_fidelity(fidelity),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                    )
                end,
            )
            result = run_fig2_benchmark_main([
                "--manifest=$FIG2_MANIFEST_PATH",
                "--stage=cli_test",
                "--output=$directory",
                "--dimensions=4",
                "--flux_units_2pi=0,3",
                "--threads=3",
            ];
                operations,
                configure_threads=threads -> push!(configured_threads, threads),
            )
            @test configured_threads == [3]
            @test result.run isa Fig2BenchmarkRun
            @test length(result.run.selections) == 2
            @test all(!row["pass"] for row in values(result.report))
            @test all(
                occursin("checkpoint", lowercase(row["reason"])) for
                row in values(result.report)
            )
            @test isfile(joinpath(directory, "ledger.toml"))
            @test isfile(joinpath(directory, "acceptance.toml"))
        end
    end

    @testset "sector tracking retains adjacent residuals as diagnostics" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            operations = Fig2BenchmarkOperations(
                progress_audit=synthetic_progress_audit,
                provenance=synthetic_fig2_provenance,
                candidate_ids=(args...) -> ["only"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    raw_polarization = point - 1.0
                    sectors = point == 1 ? Dict(0 => 1.0) : Dict(
                        0 => 0.5,
                        2 * FIG2_SYNTHETIC_CHARGE_SCALE => 0.5,
                    )
                    write_fake_candidate_files(
                        candidate_directory,
                        "point=$point";
                        requested_maxdim=dimension,
                        sector_weights=sectors,
                    )
                    fidelity = point == 1 ? NaN : 0.99
                    momentum = synthetic_counting_momentum()
                    counting = InfiniteCylinderDMRG._fig2_momentum_counting(
                        spec, momentum
                    )
                    Fig2CandidateEvidence(;
                        state=(; point),
                        requested_maxdim=dimension,
                        achieved_maxlinkdim=dimension,
                        converged=true,
                        valid=true,
                        energy_per_site=-1.0,
                        raw_schmidt_polarization=raw_polarization,
                        sector_weights=sectors,
                        fidelity_to_previous=fidelity,
                        fidelity_valid=point > 1,
                        mixed_fidelity=point == 1 ? nothing :
                            synthetic_mixed_fidelity(fidelity),
                        momentum,
                        momentum_counting=counting.counts,
                        momentum_counting_evidence=counting,
                        restart_valid=true,
                    )
                end,
            )
            run = run_fig2_benchmark(
                spec,
                directory;
                stage="sector_residual",
                dimensions=[4],
                fluxes=[0.0, 6pi],
                operations,
            )
            @test length(run.selections) == 2
            ledger = TOML.parsefile(joinpath(directory, "ledger.toml"))
            @test length(ledger["candidate"]) == 2
            @test length(ledger["selection"]) == 2
            @test ledger["candidate"][end]["point"] == 2
            @test ledger["candidate"][end]["complete"]
            tracked_fields = split(last(readlines(joinpath(
                directory, "pump_sector_tracked.tsv"
            ))), '\t')
            @test parse(Bool, tracked_fields[8])
            @test parse(Float64, tracked_fields[7]) >
                spec.data["convergence"]["sector_residual_tol"]
            report = write_synthetic_fig2_acceptance_report!(
                spec,
                directory;
                candidate_ids_provider=operations.candidate_ids,
            )
            @test !report["sector_shift"]["pass"]
            @test occursin(
                "endpoint", lowercase(report["sector_shift"]["reason"])
            )
        end
    end
end
