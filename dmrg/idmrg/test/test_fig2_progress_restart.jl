using Test
using TOML
using ITensorInfiniteMPS
using InfiniteCylinderDMRG

struct Fig2SyntheticProgressFailure <: Exception end

function fig2_progress_exception(f::Function)
    try
        f()
        return nothing
    catch error
        return error
    end
end

@testset "Fig. 2 candidate progress is atomic, chained, and resumable" begin
    manifest = joinpath(
        @__DIR__, "..", "benchmarks", "fqahc_fig2.toml"
    )
    spec = load_fig2_benchmark(manifest)
    _, _, psi = initial_infinite_mps(spec.config)
    dimension = 4
    point = 1
    phi_y = 0.0
    candidate_id = "cold_progress_probe"

    iteration_record = VUMPSRecord(
        1, 1, 1, -1.0, -1.0, 0.0, Inf,
        1.0e-4, 1.0e-4, 1.0e-4, 0.25, false,
    )
    first_event = VUMPSProgressEvent(
        1,
        :iteration,
        1,
        1,
        4,
        psi,
        nothing,
        iteration_record,
    )

    mktempdir() do directory
        first = InfiniteCylinderDMRG._fig2_persist_progress_event!(
            spec,
            directory,
            dimension,
            point,
            phi_y,
            candidate_id,
            [4],
            first_event;
            resume_count=0,
            canonicalize_state=identity,
        )
        pointer_path = joinpath(directory, ".progress", "latest.toml")
        @test isfile(pointer_path)
        pointer = TOML.parsefile(pointer_path)
        @test pointer["convergence_state_policy"] ==
            "reset_on_resume"
        @test pointer["progress_generations_to_keep"] == 2
        @test pointer["git_commit"] ==
            InfiniteCylinderDMRG._fig2_repository_commit()
        @test pointer["benchmark_source_sha256"] ==
            InfiniteCylinderDMRG._fig2_file_sha256(joinpath(
                @__DIR__, "..", "src", "Fig2Benchmark.jl"
            ))
        @test first.event_sequence == 1
        @test first.resume_count == 0
        @test first.verified_event_count == 1
        @test isfile(joinpath(directory, first.state_path))
        @test isfile(joinpath(directory, first.event_path))

        loaded = InfiniteCylinderDMRG._fig2_load_progress(
            spec,
            directory,
            dimension,
            point,
            phi_y,
            candidate_id,
        )
        @test loaded.state isa InfiniteCanonicalMPS
        @test loaded.event_sequence == 1
        @test loaded.resume_count == 0
        @test loaded.next_resume_count == 1
        @test loaded.verified_event_count == 1
        @test maximum(link_dimensions(loaded.state)) == first.maxlinkdim

        expansion_record = SubspaceExpansionRecord(
            1, 4, [1], [2], true, 0.5
        )
        resumed_event = VUMPSProgressEvent(
            1,
            :expansion,
            1,
            1,
            4,
            psi,
            expansion_record,
            nothing,
        )
        second = InfiniteCylinderDMRG._fig2_persist_progress_event!(
            spec,
            directory,
            dimension,
            point,
            phi_y,
            candidate_id,
            [4],
            resumed_event;
            resume_count=loaded.next_resume_count,
            canonicalize_state=identity,
        )
        @test second.event_sequence == 2
        @test second.runner_event_sequence == 1
        @test second.resume_count == 1
        @test second.verified_event_count == 2
        @test second.state_path != first.state_path
        @test second.event_path != first.event_path
        @test isfile(joinpath(directory, first.state_path))
        @test isfile(joinpath(directory, second.state_path))

        pointer_before_failure = read(pointer_path)
        second_iteration_record = VUMPSRecord(
            1, 2, 1, -1.1, -1.1, 0.0, 0.1,
            1.0e-5, 1.0e-5, 1.0e-5, 0.25, false,
        )
        failed_event = VUMPSProgressEvent(
            2,
            :iteration,
            1,
            2,
            4,
            psi,
            nothing,
            second_iteration_record,
        )
        failure = fig2_progress_exception() do
            InfiniteCylinderDMRG._fig2_persist_progress_event!(
                spec,
                directory,
                dimension,
                point,
                phi_y,
                candidate_id,
                [4],
                failed_event;
                resume_count=1,
                canonicalize_state=identity,
                save_state=(args...) -> throw(Fig2SyntheticProgressFailure()),
            )
        end
        @test failure isa Fig2SyntheticProgressFailure
        @test read(pointer_path) == pointer_before_failure
        @test InfiniteCylinderDMRG._fig2_load_progress(
            spec,
            directory,
            dimension,
            point,
            phi_y,
            candidate_id,
        ).event_sequence == 2

        third = InfiniteCylinderDMRG._fig2_persist_progress_event!(
            spec,
            directory,
            dimension,
            point,
            phi_y,
            candidate_id,
            [4],
            failed_event;
            resume_count=1,
            canonicalize_state=identity,
        )
        @test third.event_sequence == 3
        @test third.verified_event_count == 3
        @test !isfile(joinpath(directory, first.state_path))
        @test isfile(joinpath(directory, second.state_path))
        @test isfile(joinpath(directory, third.state_path))

        mismatch = fig2_progress_exception() do
            InfiniteCylinderDMRG._fig2_load_progress(
                spec,
                directory,
                dimension,
                point,
                phi_y,
                "different_candidate",
            )
        end
        @test mismatch isa ArgumentError
        @test occursin("candidate", lowercase(sprint(showerror, mismatch)))

        first_event_path = joinpath(directory, first.event_path)
        first_event_bytes = read(first_event_path)
        open(first_event_path, "a") do io
            write(io, "# tampered\n")
        end
        chain_error = fig2_progress_exception() do
            InfiniteCylinderDMRG._fig2_load_progress(
                spec,
                directory,
                dimension,
                point,
                phi_y,
                candidate_id,
            )
        end
        @test chain_error isa ArgumentError
        @test occursin("checksum", lowercase(sprint(showerror, chain_error)))
        open(first_event_path, "w") do io
            write(io, first_event_bytes)
        end

        latest_state_path = joinpath(directory, third.state_path)
        latest_state_bytes = read(latest_state_path)
        open(latest_state_path, "a") do io
            write(io, "tampered")
        end
        state_error = fig2_progress_exception() do
            InfiniteCylinderDMRG._fig2_load_progress(
                spec,
                directory,
                dimension,
                point,
                phi_y,
                candidate_id,
            )
        end
        @test state_error isa ArgumentError
        @test occursin("state", lowercase(sprint(showerror, state_error)))
        @test occursin("checksum", lowercase(sprint(showerror, state_error)))
        open(latest_state_path, "w") do io
            write(io, latest_state_bytes)
        end

        final_checkpoint = joinpath(directory, "state.h5")
        save_checkpoint(final_checkpoint, psi, spec.config)
        finalized = InfiniteCylinderDMRG._fig2_finalize_progress!(
            spec,
            directory,
            dimension,
            point,
            phi_y,
            candidate_id,
            final_checkpoint,
        )
        @test finalized.complete
        @test finalized.event_count == 3
        @test finalized.resume_count == 1
        @test finalized.final_state_sha256 ==
            InfiniteCylinderDMRG._fig2_file_sha256(final_checkpoint)
        progress_path = joinpath(directory, "progress.toml")
        @test isfile(progress_path)
        @test TOML.parsefile(progress_path)["convergence_state_policy"] ==
            "reset_on_resume"

        audit = InfiniteCylinderDMRG._fig2_validate_progress_artifact(
            spec,
            directory,
            dimension,
            point,
            phi_y,
            candidate_id,
            final_checkpoint,
        )
        @test audit.complete
        @test audit.event_count == 3
        @test audit.resume_count == 1
        @test audit.latest_maxlinkdim == maximum(link_dimensions(psi))
    end
end
