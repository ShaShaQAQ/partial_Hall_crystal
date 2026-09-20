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
        progress_audit=_default_fig2_progress_audit,
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
