function _mpskit_fig2_checkpoint_save(path, state, config, metadata)
    return save_mpskit_checkpoint(path, state, config, metadata)
end

function _mpskit_fig2_checkpoint_load(path, config)
    return load_mpskit_checkpoint(path, config).state
end

function _mpskit_fig2_state_maxlinkdim(state)
    dimensions = _mpskit_solver_link_dimensions(state)
    isempty(dimensions) && error("MPSKit Fig. 2 state has no virtual bonds")
    return maximum(dimensions)
end

function _mpskit_fig2_state_is_valid(state, config)
    try
        mpskit_state_space_fingerprints(state, config)
        return true
    catch error
        error isa InterruptException && rethrow()
        return false
    end
end

struct MPSKitFig2ProgressEvent{S}
    sequence::Int
    kind::Symbol
    stage::Int
    iteration::Int
    target::Int
    psi::S
    record::MPSKitSolverStageRecord

    function MPSKitFig2ProgressEvent(
        sequence::Integer,
        state,
        record::MPSKitSolverStageRecord,
    )
        sequence > 0 || throw(
            ArgumentError("MPSKit Fig. 2 progress sequence must be positive"),
        )
        return new{typeof(state)}(
            Int(sequence),
            :mpskit_stage,
            record.stage,
            1,
            record.requested_maxdim,
            state,
            record,
        )
    end
end

function _fig2_progress_event_payload(event::MPSKitFig2ProgressEvent)
    record = event.record
    return Dict{String,Any}(
        "mpskit_stage_record" => Dict(
            "stage" => record.stage,
            "requested_maxdim" => record.requested_maxdim,
            "actual_maxdim" => record.actual_maxdim,
            "cutoff" => record.cutoff,
            "idmrg_maxiter" => record.idmrg_maxiter,
            "vumps_maxiter" => record.vumps_maxiter,
            "energy_per_site" => record.energy_per_site,
            "energy_imaginary" => record.energy_imaginary,
            "idmrg_diagnostic" => record.idmrg_diagnostic,
            "vumps_residual" => record.vumps_residual,
            "recomputed_galerkin_residual" =>
                record.recomputed_galerkin_residual,
            "elapsed_seconds" => record.elapsed_seconds,
        ),
    )
end

function _mpskit_fig2_progress_metadata(event::MPSKitFig2ProgressEvent)
    record = event.record
    return (
        artifact="fig2_progress_generation",
        sequence=event.sequence,
        stage=record.stage,
        requested_maxdim=record.requested_maxdim,
        actual_maxdim=record.actual_maxdim,
        recomputed_galerkin_residual=record.recomputed_galerkin_residual,
    )
end

function _mpskit_fig2_persist_progress_event!(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    maxdim_schedule,
    event::MPSKitFig2ProgressEvent;
    resume_count::Integer,
)
    config = with_flux(spec.config, Float64(phi_y))
    return _fig2_persist_progress_event!(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        maxdim_schedule,
        event;
        resume_count,
        canonicalize_state=identity,
        save_state=(path, state, checkpoint_config) ->
            _mpskit_fig2_checkpoint_save(
                path,
                state,
                checkpoint_config,
                _mpskit_fig2_progress_metadata(event),
            ),
        state_is_valid=state -> _mpskit_fig2_state_is_valid(state, config),
        state_maxlinkdim=_mpskit_fig2_state_maxlinkdim,
        continuation_builder=identity,
    )
end

function _mpskit_fig2_finalize_progress!(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    final_checkpoint,
)
    config = with_flux(spec.config, Float64(phi_y))
    return _fig2_finalize_progress!(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        final_checkpoint;
        checkpoint_filename="state.h5",
        load_state=_mpskit_fig2_checkpoint_load,
        state_is_valid=state -> _mpskit_fig2_state_is_valid(state, config),
        state_maxlinkdim=_mpskit_fig2_state_maxlinkdim,
    )
end

function _mpskit_fig2_progress_audit(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    final_checkpoint,
)
    config = with_flux(spec.config, Float64(phi_y))
    return _fig2_validate_progress_artifact(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        final_checkpoint;
        load_state=_mpskit_fig2_checkpoint_load,
        state_is_valid=state -> _mpskit_fig2_state_is_valid(state, config),
        state_maxlinkdim=_mpskit_fig2_state_maxlinkdim,
    )
end

function _mpskit_fig2_load_state(
    spec,
    dimension,
    point,
    selection,
    checkpoint,
)
    _ = (dimension, point, selection)
    _fig2_validated_snapshot(spec)
    config = with_flux(spec.config, selection.phi_y)
    return _mpskit_fig2_checkpoint_load(checkpoint, config)
end

function _mpskit_fig2_checkpoint_audit(spec, checkpoint, phi_y)
    _fig2_validated_snapshot(spec)
    config = with_flux(spec.config, phi_y)
    return _fig2_checkpoint_audit(
        checkpoint,
        config;
        loader=_mpskit_fig2_checkpoint_load,
        maxlinkdim_of=_mpskit_fig2_state_maxlinkdim,
    )
end

function _mpskit_fig2_run_candidate(args...)
    throw(ArgumentError(
        "MPSKit Fig. 2 candidate adapter is not configured yet",
    ))
end

function mpskit_fig2_operations(spec::Fig2BenchmarkSpec)
    snapshot = _fig2_validated_snapshot(spec)
    backend = _fig2_required(snapshot, "backend")
    _fig2_required(backend, "id") == FIG2_PRODUCTION_BACKEND_ID || throw(
        ArgumentError("Fig. 2 manifest does not select the MPSKit backend"),
    )
    return Fig2BenchmarkOperations(
        backend_id=FIG2_PRODUCTION_BACKEND_ID,
        checkpoint_filename="state.h5",
        checkpoint_save=_mpskit_fig2_checkpoint_save,
        checkpoint_load=_mpskit_fig2_checkpoint_load,
        state_maxlinkdim=_mpskit_fig2_state_maxlinkdim,
        candidate_ids=_default_fig2_candidate_ids,
        run_candidate=_mpskit_fig2_run_candidate,
        load_state=_mpskit_fig2_load_state,
        provenance=_default_fig2_provenance,
        checkpoint_audit=_mpskit_fig2_checkpoint_audit,
        progress_audit=_mpskit_fig2_progress_audit,
    )
end

function fig2_operations_for_backend(spec::Fig2BenchmarkSpec)
    snapshot = _fig2_validated_snapshot(spec)
    backend = _fig2_required(snapshot, "backend")
    backend_id = _fig2_required(backend, "id")
    backend_id == FIG2_PRODUCTION_BACKEND_ID || throw(
        ArgumentError("unsupported Fig. 2 production backend id"),
    )
    return mpskit_fig2_operations(spec)
end
