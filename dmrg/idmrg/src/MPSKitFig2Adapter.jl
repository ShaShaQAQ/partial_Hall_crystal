const MPSKIT_FIG2_ALGORITHM = "MPSKit.IDMRG2+VUMPS"

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

function _mpskit_fig2_record_from_progress(data)
    data isa AbstractDict || throw(
        ArgumentError("MPSKit progress stage record must be a table"),
    )
    return MPSKitSolverStageRecord(
        Int(data["stage"]),
        Int(data["requested_maxdim"]),
        Int(data["actual_maxdim"]),
        Float64(data["cutoff"]),
        Int(data["idmrg_maxiter"]),
        Int(data["vumps_maxiter"]),
        Float64(data["energy_per_site"]),
        Float64(data["energy_imaginary"]),
        Float64(data["idmrg_diagnostic"]),
        Float64(data["vumps_residual"]),
        Float64(data["recomputed_galerkin_residual"]),
        Float64(data["elapsed_seconds"]),
    )
end

function _mpskit_fig2_load_progress(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
)
    config = with_flux(spec.config, Float64(phi_y))
    progress = _fig2_load_progress(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id;
        load_state=_mpskit_fig2_checkpoint_load,
        state_is_valid=state -> _mpskit_fig2_state_is_valid(state, config),
        state_maxlinkdim=_mpskit_fig2_state_maxlinkdim,
    )
    isnothing(progress) && return nothing
    records = MPSKitSolverStageRecord[]
    for event in progress.events
        haskey(event, "mpskit_stage_record") || throw(ArgumentError(
            "MPSKit progress event is missing its stage record",
        ))
        push!(
            records,
            _mpskit_fig2_record_from_progress(
                event["mpskit_stage_record"],
            ),
        )
    end
    getproperty.(records, :stage) == collect(eachindex(records)) || throw(
        ArgumentError("MPSKit progress stages are not contiguous"),
    )
    return merge(progress, (; records))
end

function _mpskit_fig2_global_record(record, stage_offset::Integer)
    return MPSKitSolverStageRecord(
        record.stage + stage_offset,
        record.requested_maxdim,
        record.actual_maxdim,
        record.cutoff,
        record.idmrg_maxiter,
        record.vumps_maxiter,
        record.energy_per_site,
        record.energy_imaginary,
        record.idmrg_diagnostic,
        record.vumps_residual,
        record.recomputed_galerkin_residual,
        record.elapsed_seconds,
    )
end

function _mpskit_fig2_convergence_rows(result, config, tolerance)
    previous_energy = nothing
    return [
        let
            energy = record.energy_per_site * sites_per_cell(config)
            delta = isnothing(previous_energy) ? "missing" :
                abs(energy - previous_energy)
            previous_energy = energy
            stage_converged =
                isfinite(record.energy_per_site) &&
                isfinite(record.energy_imaginary) &&
                abs(record.energy_imaginary) <= 1e-10 &&
                isfinite(record.recomputed_galerkin_residual) &&
                record.recomputed_galerkin_residual <= tolerance
            (
                record.stage,
                1,
                record.actual_maxdim,
                energy,
                energy,
                0.0,
                delta,
                record.idmrg_diagnostic,
                record.vumps_residual,
                record.recomputed_galerkin_residual,
                record.elapsed_seconds,
                stage_converged,
            )
        end for record in result.records
    ]
end

function _mpskit_fig2_summary(
    spec,
    config,
    result,
    energy,
    density,
    entanglements,
    transfer,
    maxdim_schedule,
)
    snapshot = _fig2_validated_snapshot(spec)
    optimization = snapshot["optimization"]
    backend = mpskit_backend_provenance()
    fixed = only(filter(data -> data.cut_x == 1, entanglements))
    valid = result.converged &&
        isfinite(energy.per_site) &&
        all(isfinite, density.densities) &&
        !isempty(fixed.levels) &&
        !isempty(fixed.sectors)
    return Dict{String,Any}(
        "format" => TEXT_OUTPUT_FORMAT,
        "algorithm" => MPSKIT_FIG2_ALGORITHM,
        "valid" => valid,
        "converged" => result.converged,
        "reason" => result.converged ?
            "recomputed MPSKit Galerkin gate passed" :
            "recomputed MPSKit Galerkin gate failed",
        "backend" => Dict(
            "id" => backend.backend,
            "mpskit_commit" => backend.mpskit_commit,
            "tensorkittensors_commit" => backend.tensorkittensors_commit,
            "blocktensorkit_commit" => backend.blocktensorkit_commit,
        ),
        "configuration" => Dict(
            "Ly" => config.Ly,
            "x_period" => config.x_period,
            "filling_num" => config.filling_num,
            "filling_den" => config.filling_den,
            "phi_y" => config.phi_y,
            "sites_per_cell" => sites_per_cell(config),
            "particles_per_cell" => particles_per_cell(config),
            "charge_scale" => charge_scale(config),
            "signature" => configuration_signature(config),
        ),
        "energy" => Dict(
            "valid" => isfinite(energy.per_site),
            "per_cell" => energy.per_cell,
            "per_x" => energy.per_x,
            "per_unit_cell" => energy.per_unit_cell,
            "per_site" => energy.per_site,
            "imaginary_per_cell" => energy.imaginary_per_cell,
        ),
        "observables" => Dict(
            "density_valid" => all(isfinite, density.densities),
            "entanglement_valid" =>
                all(data -> !isempty(data.levels), entanglements),
            "transfer_valid" => transfer.valid,
            "transfer_reason" => transfer.reason,
            "neutral_transfer_ratio" => transfer.ratio,
            "correlation_length_neutral_cell" => transfer.xi_cell,
            "correlation_length_neutral_x" => transfer.xi_x,
        ),
        "entanglement" => [
            Dict(
                "cut_x" => data.cut_x,
                "bond" => data.bond,
                "valid" => result.converged && !isempty(data.levels),
                "entropy" => data.entropy,
                "raw_schmidt_polarization" => data.mean_physical_charge,
            ) for data in entanglements
        ],
        "optimization" => Dict(
            "maxdim_schedule" => Int.(maxdim_schedule),
            "cutoff" => Float64(optimization["cutoff"]),
            "vumps_tol" => Float64(optimization["vumps_tol"]),
            "energy_tol" => Float64(optimization["energy_tol"]),
            "energy_mismatch_tol" =>
                Float64(optimization["energy_mismatch_tol"]),
            "stable_iterations" => Int(optimization["stable_iterations"]),
            "max_iterations" => Int(optimization["max_iterations"]),
            "multisite_update_alg" =>
                String(optimization["multisite_update_alg"]),
            "solver_tolerance_policy" =>
                String(optimization["solver_tolerance_policy"]),
            "recomputed_galerkin_gate" =>
                Float64(optimization["vumps_tol"]),
        ),
    )
end

function _mpskit_fig2_write_solver_outputs!(
    spec,
    config,
    directory,
    result,
    energy,
    density,
    entanglements,
    transfer,
    maxdim_schedule,
)
    mkpath(directory)
    summary = _mpskit_fig2_summary(
        spec,
        config,
        result,
        energy,
        density,
        entanglements,
        transfer,
        maxdim_schedule,
    )
    tolerance = Float64(spec.data["optimization"]["vumps_tol"])
    converged_rows = _mpskit_fig2_convergence_rows(
        result,
        config,
        tolerance,
    )
    contents = Dict{String,String}(
        "summary.toml" => _render_summary(summary),
        "convergence.tsv" => _render_tsv(
            FIG2_ARTIFACT_HEADERS["convergence.tsv"],
            converged_rows,
        ),
        "expansion.tsv" => _render_tsv(
            FIG2_ARTIFACT_HEADERS["expansion.tsv"],
            Tuple[],
        ),
        "density.tsv" => _render_tsv(
            FIG2_ARTIFACT_HEADERS["density.tsv"],
            (
                (
                    row.site,
                    row.x,
                    row.y,
                    row.density,
                    result.converged && isfinite(row.density),
                ) for row in density.rows
            ),
        ),
        "entanglement_spectrum.tsv" => _render_tsv(
            FIG2_ARTIFACT_HEADERS["entanglement_spectrum.tsv"],
            (
                (
                    level.cut_x,
                    level.bond,
                    level.level,
                    level.singular_value,
                    level.probability,
                    level.entanglement_energy,
                    level.sector,
                    level.raw_charge,
                    level.physical_charge,
                    result.converged,
                ) for data in entanglements for level in data.levels
            ),
        ),
        "schmidt_sectors.tsv" => _render_tsv(
            FIG2_ARTIFACT_HEADERS["schmidt_sectors.tsv"],
            (
                (
                    sector.cut_x,
                    sector.bond,
                    sector.sector,
                    sector.raw_charge,
                    sector.physical_charge,
                    sector.weight,
                    result.converged,
                ) for data in entanglements for sector in data.sectors
            ),
        ),
    )
    for (filename, bytes) in contents
        _write_fig2_bytes(joinpath(directory, filename), bytes)
    end
    return summary
end

function _validate_mpskit_fig2_convergence_tsv(
    path,
    summary_schedule;
    optimization,
    energy_normalization_sites,
)
    contract = _fig2_convergence_optimization(optimization)
    schedule = Int.(summary_schedule)
    parsed = _parse_fig2_convergence_tsv(path)
    length(parsed) == length(schedule) || throw(ArgumentError(
        "MPSKit convergence.tsv must contain exactly one audited row per stage",
    ))
    previous_energy = nothing
    for (stage, row) in enumerate(parsed)
        row.stage == stage && row.iteration == 1 || throw(ArgumentError(
            "MPSKit convergence.tsv stages must be contiguous single rows",
        ))
        0 < row.maxlinkdim <= schedule[stage] || throw(ArgumentError(
            "MPSKit convergence.tsv actual maxlinkdim exceeds its stage cap",
        ))
        _fig2_same_float(row.energy_left, row.energy_right) &&
            iszero(row.energy_mismatch) || throw(ArgumentError(
            "MPSKit convergence.tsv real energy copies are inconsistent",
        ))
        energy = (row.energy_left + row.energy_right) / 2
        if stage == 1
            ismissing(row.delta_energy) || throw(ArgumentError(
                "first MPSKit convergence row must use the missing delta sentinel",
            ))
        else
            !ismissing(row.delta_energy) && _fig2_same_float(
                row.delta_energy,
                abs(energy - something(previous_energy)),
            ) || throw(ArgumentError(
                "MPSKit convergence delta energy disagrees with adjacent stages",
            ))
        end
        all(
            value -> isfinite(value) && value >= 0,
            (
                row.eps_left,
                row.eps_right,
                row.precision_error,
                row.elapsed_seconds,
            ),
        ) || throw(ArgumentError(
            "MPSKit convergence diagnostics must be finite and nonnegative",
        ))
        expected_converged =
            row.precision_error <= contract.vumps_tol &&
            row.energy_mismatch / energy_normalization_sites <
                contract.energy_mismatch_tol
        row.converged == expected_converged || throw(ArgumentError(
            "MPSKit convergence flag disagrees with the recomputed Galerkin gate",
        ))
        previous_energy = energy
    end
    return last(parsed)
end

function _mpskit_fig2_candidate_state(
    spec,
    config,
    candidate_id,
    previous_state,
)
    if candidate_id == "warm"
        isnothing(previous_state) && throw(
            ArgumentError("warm MPSKit candidate needs the previous state"),
        )
        return previous_state
    end
    occupied_sites = _fig2_candidate_pattern(spec, candidate_id)
    return mpskit_product_state(config, occupied_sites)
end

function _mpskit_fig2_selection_core_validity(result, energy, entanglement)
    return result.converged &&
        isfinite(energy.per_site) &&
        isfinite(entanglement.mean_physical_charge) &&
        !isempty(entanglement.levels) &&
        !isempty(entanglement.sectors) &&
        all(level -> isfinite(level.singular_value) &&
            isfinite(level.probability), entanglement.levels) &&
        all(sector -> isfinite(sector.weight), entanglement.sectors)
end

function _mpskit_fig2_run_candidate(
    spec,
    dimension,
    point,
    phi_y,
    candidate_id,
    previous_state,
    candidate_directory;
    solver=run_mpskit_idmrg,
    build_hamiltonian=mpskit_infinite_hamiltonian,
)
    snapshot = _fig2_validated_snapshot(spec)
    config = with_flux(spec.config, phi_y)
    full_schedule = _fig2_maxdim_schedule(dimension)
    progress = _mpskit_fig2_load_progress(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
    )
    state = isnothing(progress) ? _mpskit_fig2_candidate_state(
        spec,
        config,
        candidate_id,
        previous_state,
    ) : progress.state
    prior_records = isnothing(progress) ? MPSKitSolverStageRecord[] :
        progress.records
    completed_stages = length(prior_records)
    completed_stages <= length(full_schedule) || throw(ArgumentError(
        "MPSKit progress has more stages than the requested schedule",
    ))
    getproperty.(prior_records, :requested_maxdim) ==
        full_schedule[1:completed_stages] || throw(ArgumentError(
        "MPSKit progress schedule disagrees with the requested schedule",
    ))
    remaining_schedule = full_schedule[(completed_stages + 1):end]
    hamiltonian = build_hamiltonian(config, spec.model)
    optimization = snapshot["optimization"]
    resume_count = isnothing(progress) ? 0 : progress.next_resume_count

    result = if isempty(remaining_schedule)
        environments = MPSKit.environments(state, hamiltonian, state)
        record = last(prior_records)
        dimensions = _mpskit_solver_link_dimensions(state)
        MPSKitSolverResult(
            state,
            environments,
            prior_records,
            record.energy_per_site,
            record.energy_imaginary,
            record.vumps_residual,
            record.recomputed_galerkin_residual,
            dimensions,
            mpskit_solver_converged(;
                idmrg_diagnostic=record.idmrg_diagnostic,
                recomputed_galerkin_residual=
                    record.recomputed_galerkin_residual,
                energy_per_site=record.energy_per_site,
                energy_imaginary=record.energy_imaginary,
                final_stage_reached=true,
                galerkin_tol=Float64(optimization["vumps_tol"]),
                energy_imag_tol=1e-10,
            ),
        )
    else
        callback_sequence = Ref(0)
        local_result = solver(
            hamiltonian,
            state;
            maxdim_schedule=remaining_schedule,
            cutoff=Float64(optimization["cutoff"]),
            idmrg_maxiter=Int(optimization["max_iterations"]),
            vumps_maxiter=Int(optimization["max_iterations"]),
            galerkin_tol=Float64(optimization["vumps_tol"]),
            energy_imag_tol=1e-10,
            progress_callback=(stage_state, environments, records) -> begin
                _ = environments
                callback_sequence[] += 1
                global_record = _mpskit_fig2_global_record(
                    last(records),
                    completed_stages,
                )
                event = MPSKitFig2ProgressEvent(
                    callback_sequence[],
                    stage_state,
                    global_record,
                )
                _mpskit_fig2_persist_progress_event!(
                    spec,
                    candidate_directory,
                    dimension,
                    point,
                    phi_y,
                    candidate_id,
                    full_schedule,
                    event;
                    resume_count,
                )
            end,
        )
        global_records = [
            prior_records;
            [_mpskit_fig2_global_record(record, completed_stages) for
             record in local_result.records]
        ]
        MPSKitSolverResult(
            local_result.state,
            local_result.environments,
            global_records,
            local_result.energy_per_site,
            local_result.energy_imaginary,
            local_result.galerkin_residual,
            local_result.recomputed_galerkin_residual,
            local_result.link_dimensions,
            local_result.converged,
        )
    end

    isempty(result.records) && error("MPSKit Fig. 2 solver produced no stages")
    checkpoint = joinpath(candidate_directory, "state.h5")
    checkpoint_metadata = (
        artifact="fig2_candidate",
        candidate_id=String(candidate_id),
        dimension=Int(dimension),
        point=Int(point),
        completed_stages=length(result.records),
        requested_maxdim=Int(dimension),
        achieved_maxlinkdim=maximum(result.link_dimensions),
        recomputed_galerkin_residual=result.recomputed_galerkin_residual,
    )
    _mpskit_fig2_checkpoint_save(
        checkpoint,
        result.state,
        config,
        checkpoint_metadata,
    )
    _mpskit_fig2_finalize_progress!(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        checkpoint,
    )

    energy = mpskit_energy_data(
        result.state,
        hamiltonian,
        config;
        environments=result.environments,
    )
    density = mpskit_density_data(result.state, config)
    entanglements = mpskit_entanglement_data(result.state, config)
    transfer = mpskit_neutral_transfer_data(
        result.state,
        config;
        tol=Float64(optimization["transfer_tol"]),
        seed=_fig2_candidate_seed(dimension, point, candidate_id),
    )
    _mpskit_fig2_write_solver_outputs!(
        spec,
        config,
        candidate_directory,
        result,
        energy,
        density,
        entanglements,
        transfer,
        full_schedule,
    )

    fixed = only(filter(data -> data.cut_x == 1, entanglements))
    sector_weights = Dict(
        sector.raw_charge => sector.weight for sector in fixed.sectors
    )
    momentum_result = mpskit_momentum_entanglement_data(
        result.state,
        config;
        cut_x=1,
        tol=Float64(optimization["transfer_tol"]),
        residual_tol=Float64(optimization["momentum_residual_tol"]),
        seed=_fig2_candidate_seed(dimension, point, candidate_id),
    )
    momentum = momentum_result.spectrum
    momentum_counting = _fig2_momentum_counting(spec, momentum)
    fidelity = isnothing(previous_state) ? nothing :
        mpskit_mixed_transfer_fidelity(
            previous_state,
            result.state,
            config;
            tol=Float64(optimization["transfer_tol"]),
            seed=_fig2_candidate_seed(dimension, point, candidate_id),
        )
    checkpoint_audit = _mpskit_fig2_checkpoint_audit(
        spec,
        checkpoint,
        phi_y,
    )
    achieved_maxlinkdim = maximum(result.link_dimensions)
    restart_valid = checkpoint_audit.restart_valid &&
        checkpoint_audit.checkpoint_maxlinkdim == achieved_maxlinkdim
    reasons = String[]
    result.converged || push!(reasons, "recomputed Galerkin gate failed")
    transfer.valid || push!(reasons, transfer.reason)
    momentum.valid || push!(reasons, momentum.reason)
    !isnothing(fidelity) && !fidelity.valid && push!(reasons, fidelity.reason)
    restart_valid || push!(reasons, "MPSKit checkpoint restart audit failed")
    return Fig2CandidateEvidence(;
        state=result.state,
        requested_maxdim=Int(dimension),
        achieved_maxlinkdim,
        checkpoint_maxlinkdim=checkpoint_audit.checkpoint_maxlinkdim,
        converged=result.converged,
        valid=_mpskit_fig2_selection_core_validity(
            result,
            energy,
            fixed,
        ),
        energy_per_site=energy.per_site,
        raw_schmidt_polarization=fixed.mean_physical_charge,
        sector_weights,
        fidelity_to_previous=isnothing(fidelity) ? NaN :
            fidelity.fidelity_cell,
        fidelity_valid=!isnothing(fidelity) && fidelity.valid,
        mixed_fidelity=fidelity,
        momentum,
        momentum_counting=momentum_counting.valid ?
            momentum_counting.counts : Int[],
        momentum_counting_evidence=momentum_counting,
        restart_valid,
        reason=join(reasons, "; "),
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

function _mpskit_fig2_provenance(
    spec,
    output,
    runtime_seconds;
    base_provenance=_default_fig2_provenance,
)
    provenance = Dict{String,Any}(
        base_provenance(spec, output, runtime_seconds)
    )
    backend = mpskit_backend_provenance()
    adapter_source = abspath(@__FILE__)
    merge!(provenance, Dict{String,Any}(
        "backend_id" => backend.backend,
        "mpskit_commit" => backend.mpskit_commit,
        "tensorkittensors_commit" => backend.tensorkittensors_commit,
        "blocktensorkit_commit" => backend.blocktensorkit_commit,
        "backend_adapter_source" => adapter_source,
        "backend_adapter_source_sha256" =>
            _fig2_file_sha256(adapter_source),
    ))
    return provenance
end

function mpskit_fig2_operations(
    spec::Fig2BenchmarkSpec;
    solver=run_mpskit_idmrg,
    provenance=_mpskit_fig2_provenance,
    candidate_ids=_default_fig2_candidate_ids,
    build_hamiltonian=mpskit_infinite_hamiltonian,
)
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
        candidate_ids=candidate_ids,
        run_candidate=(args...) -> _mpskit_fig2_run_candidate(
            args...;
            solver,
            build_hamiltonian,
        ),
        load_state=_mpskit_fig2_load_state,
        provenance=provenance,
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
