using Test
using TOML
using InfiniteCylinderDMRG

import MPSKit

const MPSKIT_REFINEMENT_MANIFEST_PATH = normpath(joinpath(
    @__DIR__, "..", "benchmarks", "fqahc_fig2.toml"
))

@testset "immutable Fig. 2 same-dimension refinement contract" begin
    manifest = TOML.parsefile(MPSKIT_REFINEMENT_MANIFEST_PATH)
    @test manifest["format"] == "fqahc_fig2_benchmark_v6"
    @test get(manifest["optimization"], "max_refinement_chunks", 0) == 8

    spec = load_fig2_benchmark(MPSKIT_REFINEMENT_MANIFEST_PATH)
    @test spec.data["optimization"]["max_refinement_chunks"] == 8
end

struct MPSKitRefinementInterrupted <: Exception end

function mpskit_refinement_fixture_record(
    state,
    hamiltonian;
    stage,
    iteration,
    kind,
    requested_maxdim,
    residual,
)
    environments = MPSKit.environments(state, hamiltonian, state)
    energy = ComplexF64(
        MPSKit.expectation_value(state, hamiltonian, environments) /
        length(state),
    )
    dimensions = InfiniteCylinderDMRG._mpskit_solver_link_dimensions(state)
    return MPSKitSolverStageRecord(
        stage,
        iteration,
        kind,
        requested_maxdim,
        maximum(dimensions),
        1.0e-9,
        kind == :idmrg2_vumps ? 1 : 0,
        1,
        real(energy),
        imag(energy),
        kind == :idmrg2_vumps ? Float64(residual) : 0.0,
        Float64(residual),
        Float64(residual),
        0.01,
    )
end

function mpskit_refinement_fixture_result(
    state,
    hamiltonian,
    records;
    converged,
)
    environments = MPSKit.environments(state, hamiltonian, state)
    final = last(records)
    dimensions = InfiniteCylinderDMRG._mpskit_solver_link_dimensions(state)
    return MPSKitSolverResult(
        state,
        environments,
        records,
        final.energy_per_site,
        final.energy_imaginary,
        final.vumps_residual,
        final.recomputed_galerkin_residual,
        dimensions,
        converged,
    )
end

@testset "MPSKit refinement progress preserves iteration and kind" begin
    record = MPSKitSolverStageRecord(
        1,
        3,
        :vumps_refinement,
        1,
        1,
        1.0e-9,
        0,
        1,
        -0.25,
        0.0,
        0.0,
        4.0e-7,
        4.0e-7,
        0.01,
    )
    event = InfiniteCylinderDMRG.MPSKitFig2ProgressEvent(1, nothing, record)
    @test event.kind == :vumps_refinement
    @test event.stage == 1
    @test event.iteration == 3
    payload = InfiniteCylinderDMRG._fig2_progress_event_payload(event)
    stored = payload["mpskit_stage_record"]
    @test stored["iteration"] == 3
    @test stored["kind"] == "vumps_refinement"
    restored = InfiniteCylinderDMRG._mpskit_fig2_record_from_progress(stored)
    @test restored.stage == 1
    @test restored.iteration == 3
    @test restored.kind == :vumps_refinement
end

@testset "MPSKit candidate resumes same-D refinement from the next chunk" begin
    spec = load_fig2_benchmark(MPSKIT_REFINEMENT_MANIFEST_PATH)
    candidate_id = first(fig2_initial_candidates(spec.config)).id
    growth_calls = Ref(0)
    refinement_calls = Ref(0)

    function nonconverged_growth(
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
        _ = (
            cutoff,
            idmrg_maxiter,
            vumps_maxiter,
            galerkin_tol,
            energy_imag_tol,
            verbosity,
        )
        growth_calls[] += 1
        @test maxdim_schedule == [1]
        record = mpskit_refinement_fixture_record(
            initial_state,
            hamiltonian;
            stage=1,
            iteration=1,
            kind=:idmrg2_vumps,
            requested_maxdim=1,
            residual=2.0e-4,
        )
        result = mpskit_refinement_fixture_result(
            initial_state,
            hamiltonian,
            [record];
            converged=false,
        )
        progress_callback(initial_state, result.environments, [record])
        return result
    end

    function interrupt_then_refine(
        hamiltonian,
        initial_state;
        stage,
        start_iteration,
        requested_maxdim,
        cutoff,
        vumps_maxiter,
        galerkin_tol,
        energy_imag_tol=1e-10,
        max_chunks,
        stable_iterations,
        verbosity=0,
        progress_callback=(args...) -> nothing,
    )
        _ = (
            cutoff,
            vumps_maxiter,
            galerkin_tol,
            energy_imag_tol,
            verbosity,
        )
        refinement_calls[] += 1
        @test stage == 1
        @test requested_maxdim == 1
        if refinement_calls[] == 1
            @test start_iteration == 2
            @test max_chunks == 8
            @test stable_iterations == 2
            record = mpskit_refinement_fixture_record(
                initial_state,
                hamiltonian;
                stage,
                iteration=2,
                kind=:vumps_refinement,
                requested_maxdim,
                residual=5.0e-7,
            )
            result = mpskit_refinement_fixture_result(
                initial_state,
                hamiltonian,
                [record];
                converged=false,
            )
            progress_callback(initial_state, result.environments, [record])
            throw(MPSKitRefinementInterrupted())
        end
        @test start_iteration == 3
        @test max_chunks == 7
        @test stable_iterations == 1
        record = mpskit_refinement_fixture_record(
            initial_state,
            hamiltonian;
            stage,
            iteration=3,
            kind=:vumps_refinement,
            requested_maxdim,
            residual=4.0e-7,
        )
        result = mpskit_refinement_fixture_result(
            initial_state,
            hamiltonian,
            [record];
            converged=true,
        )
        progress_callback(initial_state, result.environments, [record])
        return result
    end

    operations = mpskit_fig2_operations(
        spec;
        solver=nonconverged_growth,
        refiner=interrupt_then_refine,
    )
    mktempdir() do directory
        arguments = (
            spec,
            1,
            1,
            0.0,
            candidate_id,
            nothing,
            directory,
        )
        interrupted = try
            operations.run_candidate(arguments...)
            nothing
        catch error
            error
        end
        @test interrupted isa MPSKitRefinementInterrupted
        first_progress = InfiniteCylinderDMRG._mpskit_fig2_load_progress(
            spec,
            directory,
            1,
            1,
            0.0,
            candidate_id,
        )
        @test first_progress.event_sequence == 2
        @test getproperty.(first_progress.records, :iteration) == [1, 2]
        @test getproperty.(first_progress.records, :kind) ==
            [:idmrg2_vumps, :vumps_refinement]

        evidence = operations.run_candidate(arguments...)
        @test evidence.converged
        @test evidence.valid
        @test growth_calls[] == 1
        @test refinement_calls[] == 2
        progress = InfiniteCylinderDMRG._mpskit_fig2_load_progress(
            spec,
            directory,
            1,
            1,
            0.0,
            candidate_id,
        )
        @test progress.event_sequence == 3
        @test progress.resume_count == 1
        @test getproperty.(progress.records, :stage) == [1, 1, 1]
        @test getproperty.(progress.records, :iteration) == [1, 2, 3]
        @test getproperty.(progress.records, :kind) == [
            :idmrg2_vumps,
            :vumps_refinement,
            :vumps_refinement,
        ]
        @test isfile(joinpath(directory, "state.h5"))
        @test TOML.parsefile(joinpath(directory, "summary.toml"))[
            "converged"
        ]
        @test length(readlines(joinpath(directory, "convergence.tsv"))) == 4
    end
end

@testset "MPSKit candidate saves an invalid state after refinement budget" begin
    spec = load_fig2_benchmark(MPSKIT_REFINEMENT_MANIFEST_PATH)
    candidate_id = first(fig2_initial_candidates(spec.config)).id

    function nonconverged_growth(
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
        _ = (
            maxdim_schedule,
            cutoff,
            idmrg_maxiter,
            vumps_maxiter,
            galerkin_tol,
            energy_imag_tol,
            verbosity,
        )
        record = mpskit_refinement_fixture_record(
            initial_state,
            hamiltonian;
            stage=1,
            iteration=1,
            kind=:idmrg2_vumps,
            requested_maxdim=1,
            residual=2.0e-4,
        )
        result = mpskit_refinement_fixture_result(
            initial_state,
            hamiltonian,
            [record];
            converged=false,
        )
        progress_callback(initial_state, result.environments, [record])
        return result
    end

    function exhaust_refinement(
        hamiltonian,
        initial_state;
        stage,
        start_iteration,
        requested_maxdim,
        cutoff,
        vumps_maxiter,
        galerkin_tol,
        energy_imag_tol=1e-10,
        max_chunks,
        stable_iterations,
        verbosity=0,
        progress_callback=(args...) -> nothing,
    )
        _ = (
            cutoff,
            vumps_maxiter,
            galerkin_tol,
            energy_imag_tol,
            stable_iterations,
            verbosity,
        )
        @test max_chunks == 8
        records = MPSKitSolverStageRecord[]
        for chunk in 1:max_chunks
            push!(records, mpskit_refinement_fixture_record(
                initial_state,
                hamiltonian;
                stage,
                iteration=start_iteration + chunk - 1,
                kind=:vumps_refinement,
                requested_maxdim,
                residual=1.0e-4,
            ))
            result = mpskit_refinement_fixture_result(
                initial_state,
                hamiltonian,
                copy(records);
                converged=false,
            )
            progress_callback(initial_state, result.environments, copy(records))
        end
        return mpskit_refinement_fixture_result(
            initial_state,
            hamiltonian,
            records;
            converged=false,
        )
    end

    operations = mpskit_fig2_operations(
        spec;
        solver=nonconverged_growth,
        refiner=exhaust_refinement,
    )
    mktempdir() do directory
        evidence = operations.run_candidate(
            spec,
            1,
            1,
            0.0,
            candidate_id,
            nothing,
            directory,
        )
        @test !evidence.converged
        @test !evidence.valid
        @test isfile(joinpath(directory, "state.h5"))
        summary = TOML.parsefile(joinpath(directory, "summary.toml"))
        @test !summary["converged"]
        progress = InfiniteCylinderDMRG._mpskit_fig2_load_progress(
            spec,
            directory,
            1,
            1,
            0.0,
            candidate_id,
        )
        @test progress.event_sequence == 9
        @test last(progress.records).iteration == 9
        @test last(progress.records).kind == :vumps_refinement
        @test length(readlines(joinpath(directory, "convergence.tsv"))) == 10
    end
end
