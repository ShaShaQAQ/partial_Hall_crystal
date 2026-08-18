using Test
using TOML
using InfiniteCylinderDMRG

const FIG2_MANIFEST_PATH = normpath(joinpath(
    @__DIR__, "..", "benchmarks", "fqahc_fig2.toml"
))

@testset "immutable Fig. 2 benchmark manifest" begin
    manifest = TOML.parsefile(FIG2_MANIFEST_PATH)
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
    @test manifest["pump_tolerance"] == 0.05
    @test manifest["counting"] == [1, 1, 2, 3, 5]
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

    function write_fake_candidate_files(directory, marker)
        mkpath(directory)
        for filename in (
            "state.h5",
            "summary.toml",
            "convergence.tsv",
            "expansion.tsv",
            "density.tsv",
            "entanglement_spectrum.tsv",
            "schmidt_sectors.tsv",
        )
            open(joinpath(directory, filename), "w") do io
                write(io, "$filename:$marker\n")
            end
        end
        return nothing
    end

    @testset "atomic candidate ledger preserves raw pump evidence" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            calls = NamedTuple[]
            operations = Fig2BenchmarkOperations(
                candidate_ids=(spec, dimension, point, previous) ->
                    point == 1 ? ["alpha", "beta"] : ["warm", "cold"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    push!(calls, (; dimension, point, phi_y, candidate_id,
                                  previous_state))
                    write_fake_candidate_files(
                        candidate_directory,
                        "D=$dimension,point=$point,candidate=$candidate_id",
                    )
                    first_point = point == 1
                    beta_or_warm = candidate_id in ("beta", "warm")
                    energy = first_point ?
                        (candidate_id == "alpha" ? -1.0 : -1.1) :
                        (candidate_id == "warm" ? -1.08 : -1.2)
                    fidelity = first_point ? NaN :
                        (candidate_id == "warm" ? 0.97 : 0.65)
                    raw_polarization = first_point ?
                        (candidate_id == "alpha" ? 0.20 : 0.25) :
                        (candidate_id == "warm" ? 1.25 : 9.0)
                    sectors = first_point ?
                        Dict(0 => 0.7, 1 => 0.3) :
                        Dict(1 => 0.7, 2 => 0.3)
                    momentum_valid = !(point == 2 && candidate_id == "warm")
                    Fig2CandidateEvidence(;
                        state=(; point, candidate_id),
                        converged=true,
                        valid=true,
                        energy_per_site=energy,
                        raw_schmidt_polarization=raw_polarization,
                        sector_weights=sectors,
                        fidelity_to_previous=fidelity,
                        fidelity_valid=!first_point,
                        momentum=synthetic_momentum(; valid=momentum_valid),
                        momentum_counting=momentum_valid ? [1, 1, 2, 3, 5] : Int[],
                        restart_valid=beta_or_warm,
                        provenance_valid=true,
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
            @test [row.candidate_id for row in run.selections] == ["beta", "warm"]
            @test [row.selection_basis for row in run.selections] ==
                ["minimum_zero_flux_energy", "maximum_mixed_fidelity"]
            @test length(calls) == 4
            @test all(isnothing(call.previous_state) for call in calls[1:2])
            @test all(call.previous_state == (; point=1, candidate_id="beta") for call in calls[3:4])

            ledger = TOML.parsefile(joinpath(directory, "ledger.toml"))
            @test ledger["format"] == "fqahc_fig2_ledger_v1"
            @test length(ledger["candidate"]) == 4
            @test length(ledger["selection"]) == 2
            @test all(row["complete"] for row in ledger["candidate"])
            @test [row["candidate_id"] for row in ledger["selection"]] ==
                ["beta", "warm"]
            @test all(row["state_sha256"] != "" for row in ledger["candidate"])

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
                    "momentum_entanglement_spectrum.tsv",
                    "mixed_fidelity.tsv",
                    "candidate.toml",
                )
                    @test isfile(joinpath(candidate_directory, filename))
                end
            end

            raw_path = joinpath(directory, "pump_raw.tsv")
            raw_before = read(raw_path)
            report = write_fig2_acceptance_report!(spec, directory)
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
            @test !report["paper_curve"]["pass"]
            @test !report["bond_convergence"]["pass"]
            @test report["restart"]["pass"]
            @test report["provenance"]["pass"]
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
                candidate_ids=(spec, dimension, point, previous) ->
                    point == 1 ? ["alpha", "beta"] : ["warm", "cold"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    push!(candidate_calls, (; point, candidate_id, previous_state))
                    write_fake_candidate_files(
                        candidate_directory,
                        "point=$point,candidate=$candidate_id",
                    )
                    first_point = point == 1
                    warm = candidate_id == "warm"
                    Fig2CandidateEvidence(;
                        state=(; source=:fresh, point, candidate_id),
                        converged=true,
                        valid=true,
                        energy_per_site=first_point ?
                            (candidate_id == "alpha" ? -1.0 : -1.1) :
                            (warm ? -1.08 : -1.2),
                        raw_schmidt_polarization=first_point ?
                            (candidate_id == "alpha" ? 0.20 : 0.25) :
                            (warm ? 1.25 : 9.0),
                        sector_weights=first_point ?
                            Dict(0 => 0.7, 1 => 0.3) :
                            Dict(1 => 0.7, 2 => 0.3),
                        fidelity_to_previous=first_point ? NaN :
                            (warm ? 0.97 : 0.65),
                        fidelity_valid=!first_point,
                        momentum=synthetic_momentum(),
                        momentum_counting=[1, 1, 2, 3, 5],
                        restart_valid=true,
                        provenance_valid=true,
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
                candidate_ids=(args...) -> ["only"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    write_fake_candidate_files(candidate_directory, "point=$point")
                    Fig2CandidateEvidence(;
                        state=(; point),
                        converged=true,
                        valid=true,
                        energy_per_site=-1.0,
                        raw_schmidt_polarization=point - 1.0,
                        sector_weights=point == 1 ? Dict(0 => 1.0) : Dict(1 => 1.0),
                        fidelity_to_previous=point == 1 ? NaN : 0.99,
                        fidelity_valid=point > 1,
                        momentum=synthetic_momentum(),
                        momentum_counting=[1, 1, 2, 3, 5],
                        restart_valid=true,
                        provenance_valid=true,
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
            @test result.report["endpoint_pump"]["pass"]
            @test isfile(joinpath(directory, "ledger.toml"))
            @test isfile(joinpath(directory, "acceptance.toml"))
        end
    end

    @testset "sector tracking requires a valid correspondence residual" begin
        mktempdir() do directory
            spec = load_fig2_benchmark(FIG2_MANIFEST_PATH)
            operations = Fig2BenchmarkOperations(
                candidate_ids=(args...) -> ["only"],
                run_candidate=(spec, dimension, point, phi_y, candidate_id,
                               previous_state, candidate_directory) -> begin
                    write_fake_candidate_files(candidate_directory, "point=$point")
                    Fig2CandidateEvidence(;
                        state=(; point),
                        converged=true,
                        valid=true,
                        energy_per_site=-1.0,
                        raw_schmidt_polarization=point - 1.0,
                        sector_weights=point == 1 ?
                            Dict(0 => 1.0) : Dict(0 => 0.5, 1 => 0.5),
                        fidelity_to_previous=point == 1 ? NaN : 0.99,
                        fidelity_valid=point > 1,
                        momentum=synthetic_momentum(),
                        momentum_counting=[1, 1, 2, 3, 5],
                        restart_valid=true,
                        provenance_valid=true,
                    )
                end,
            )
            run_fig2_benchmark(
                spec,
                directory;
                stage="sector_residual",
                dimensions=[4],
                fluxes=[0.0, 6pi],
                operations,
            )
            endpoint = split(readlines(joinpath(
                directory, "pump_sector_tracked.tsv"
            ))[end], '\t')
            @test endpoint[6] == "0"
            @test endpoint[8] == "false"
            report = write_fig2_acceptance_report!(spec, directory)
            @test !report["sector_shift"]["pass"]
        end
    end
end
