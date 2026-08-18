const FIG2_MANIFEST_FORMAT = "fqahc_fig2_benchmark_v1"
const FIG2_LEDGER_FORMAT = "fqahc_fig2_ledger_v1"
const FIG2_REQUIRED_CANDIDATE_FILES = (
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

struct Fig2BenchmarkSpec
    source_path::String
    source_bytes::Vector{UInt8}
    sha256::String
    data::Dict{String,Any}
    config::InfiniteCylinderConfig
    model::CylinderModelParams
    acceptance_rows::Vector{String}
end

Base.@kwdef struct Fig2CandidateEvidence
    state::Any
    converged::Bool
    valid::Bool
    energy_per_site::Float64
    raw_schmidt_polarization::Float64
    sector_weights::Dict{Int,Float64}
    fidelity_to_previous::Float64=NaN
    fidelity_valid::Bool=false
    momentum::MomentumEntanglementData
    momentum_counting::Vector{Int}=Int[]
    restart_valid::Bool=false
    provenance_valid::Bool=false
    reason::String=""
end

function _default_fig2_candidate_ids end
function _default_fig2_run_candidate end
function _default_fig2_load_state end

Base.@kwdef struct Fig2BenchmarkOperations
    candidate_ids::Function=_default_fig2_candidate_ids
    run_candidate::Function=_default_fig2_run_candidate
    load_state::Function=_default_fig2_load_state
end

struct Fig2Selection
    dimension::Int
    point::Int
    phi_y::Float64
    candidate_id::String
    selection_basis::String
    directory::String
    energy_per_site::Float64
    raw_schmidt_polarization::Float64
    fidelity_to_previous::Float64
    fidelity_valid::Bool
    sector_weights::Dict{Int,Float64}
    momentum_valid::Bool
    momentum_reason::String
    momentum_counting::Vector{Int}
    restart_valid::Bool
    provenance_valid::Bool
end

struct Fig2BenchmarkRun
    spec::Fig2BenchmarkSpec
    stage::String
    output::String
    selections::Vector{Fig2Selection}
end

function _fig2_required(data, key::String)
    haskey(data, key) || throw(ArgumentError("Fig. 2 manifest is missing $key"))
    return data[key]
end

function _fig2_sha256(bytes::AbstractVector{UInt8})
    return bytes2hex(SHA.sha256(bytes))
end

_fig2_file_sha256(path::AbstractString) = _fig2_sha256(read(path))

function _fig2_number(data, key::String)
    value = _fig2_required(data, key)
    value isa Real && !(value isa Bool) || throw(
        ArgumentError("Fig. 2 manifest $key must be a real number")
    )
    converted = Float64(value)
    isfinite(converted) || throw(
        ArgumentError("Fig. 2 manifest $key must be finite")
    )
    return converted
end

function _fig2_integer(data, key::String)
    value = _fig2_required(data, key)
    value isa Integer && !(value isa Bool) || throw(
        ArgumentError("Fig. 2 manifest $key must be an integer")
    )
    return Int(value)
end

function _validate_fig2_manifest(data)
    _fig2_required(data, "format") == FIG2_MANIFEST_FORMAT || throw(
        ArgumentError("unsupported Fig. 2 benchmark manifest format")
    )
    _fig2_required(data, "geometry") == "paper_straight" || throw(
        ArgumentError("Fig. 2 benchmark requires paper_straight geometry")
    )
    _fig2_integer(data, "Ny") == 6 || throw(
        ArgumentError("Fig. 2 benchmark requires Ny=6")
    )
    _fig2_integer(data, "x_period") == 3 || throw(
        ArgumentError("Fig. 2 benchmark requires x_period=3")
    )
    (_fig2_integer(data, "filling_num"), _fig2_integer(data, "filling_den")) ==
        (7, 9) || throw(ArgumentError("Fig. 2 benchmark requires filling 7/9"))
    (_fig2_number(data, "t1"), _fig2_number(data, "t3"),
     _fig2_number(data, "V1"), _fig2_number(data, "V2"),
     _fig2_number(data, "V3")) == (1.0, 0.2, 10.0, 2.0, 2.0) || throw(
        ArgumentError("Fig. 2 benchmark model parameters do not match the paper")
    )
    _fig2_number(data, "flux_start") == 0.0 || throw(
        ArgumentError("Fig. 2 benchmark flux must start at zero")
    )
    _fig2_number(data, "flux_stop") == 6pi || throw(
        ArgumentError("Fig. 2 benchmark flux must stop at 6pi")
    )
    _fig2_required(data, "raw_pump_policy") ==
        "no_offset_no_unwrap_no_sign_flip_no_branch_translation" || throw(
        ArgumentError("Fig. 2 manifest must preserve raw pump data")
    )
    counting = _fig2_required(data, "counting")
    counting == [1, 1, 2, 3, 5] || throw(
        ArgumentError("Fig. 2 benchmark counting target is not fixed")
    )
    acceptance = _fig2_required(data, "acceptance")
    rows = String.(_fig2_required(acceptance, "rows"))
    rows == [
        "endpoint_pump",
        "sector_shift",
        "momentum_counting",
        "paper_curve",
        "bond_convergence",
        "restart",
        "provenance",
    ] || throw(ArgumentError("Fig. 2 acceptance rows are not the fixed contract"))
    return rows
end

function load_fig2_benchmark(path::AbstractString)
    source_path = abspath(path)
    isfile(source_path) || throw(ArgumentError("Fig. 2 manifest does not exist"))
    source_bytes = read(source_path)
    parsed = TOML.parse(String(copy(source_bytes)))
    data = Dict{String,Any}(parsed)
    rows = _validate_fig2_manifest(data)
    config = InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=_fig2_integer(data, "Ny"),
        x_period=_fig2_integer(data, "x_period"),
        filling_num=_fig2_integer(data, "filling_num"),
        filling_den=_fig2_integer(data, "filling_den"),
        phi_y=_fig2_number(data, "flux_start"),
    )
    model = CylinderModelParams(;
        t1=_fig2_number(data, "t1"),
        t3=_fig2_number(data, "t3"),
        V1=_fig2_number(data, "V1"),
        V2=_fig2_number(data, "V2"),
        V3=_fig2_number(data, "V3"),
    )
    return Fig2BenchmarkSpec(
        source_path,
        copy(source_bytes),
        _fig2_sha256(source_bytes),
        deepcopy(data),
        config,
        model,
        rows,
    )
end

function _write_fig2_bytes(path::AbstractString, bytes)
    return _atomic_replace(path) do temporary
        open(temporary, "w") do io
            write(io, bytes)
        end
    end
end

function _write_fig2_toml(path::AbstractString, data)
    return _write_fig2_bytes(path, _render_summary(data))
end

function _ensure_fig2_manifest_copy(spec::Fig2BenchmarkSpec, output)
    path = joinpath(output, "manifest.toml")
    if isfile(path)
        read(path) == spec.source_bytes || throw(
            ArgumentError("existing benchmark manifest copy differs byte-for-byte")
        )
    else
        _write_fig2_bytes(path, spec.source_bytes)
    end
    return path
end

function _new_fig2_ledger(spec, stage)
    return Dict{String,Any}(
        "format" => FIG2_LEDGER_FORMAT,
        "manifest_sha256" => spec.sha256,
        "manifest_path" => "manifest.toml",
        "stage" => String(stage),
        "candidate" => Any[],
        "selection" => Any[],
    )
end

function _load_fig2_ledger(spec, output, stage)
    path = joinpath(output, "ledger.toml")
    if !isfile(path)
        return _new_fig2_ledger(spec, stage)
    end
    ledger = Dict{String,Any}(TOML.parsefile(path))
    get(ledger, "format", "") == FIG2_LEDGER_FORMAT || throw(
        ArgumentError("unsupported Fig. 2 ledger format")
    )
    get(ledger, "manifest_sha256", "") == spec.sha256 || throw(
        ArgumentError("ledger manifest checksum does not match")
    )
    get(ledger, "stage", "") == String(stage) || throw(
        ArgumentError("ledger stage does not match requested stage")
    )
    haskey(ledger, "candidate") || (ledger["candidate"] = Any[])
    haskey(ledger, "selection") || (ledger["selection"] = Any[])
    return ledger
end

function _candidate_relative_directory(dimension, point, candidate_id)
    safe_id = replace(String(candidate_id), r"[^A-Za-z0-9_.-]" => "_")
    isempty(safe_id) && throw(ArgumentError("candidate ID has no safe filename characters"))
    return joinpath(
        "D_$(lpad(dimension, 4, '0'))",
        "phi_$(lpad(point - 1, 4, '0'))",
        "candidate_$safe_id",
    )
end

function _validate_fig2_evidence(evidence::Fig2CandidateEvidence, point)
    isfinite(evidence.energy_per_site) || throw(
        ArgumentError("candidate energy per site must be finite")
    )
    isfinite(evidence.raw_schmidt_polarization) || throw(
        ArgumentError("candidate raw Schmidt polarization must be finite")
    )
    normalized = _normalized_sector_weights(evidence.sector_weights)
    if point == 1
        evidence.fidelity_valid && throw(
            ArgumentError("zero-flux candidate must not claim previous-state fidelity")
        )
    elseif evidence.fidelity_valid
        isfinite(evidence.fidelity_to_previous) &&
            0 <= evidence.fidelity_to_previous <= 1 || throw(
            ArgumentError("valid mixed fidelity must lie in [0, 1]")
        )
    end
    if evidence.momentum.valid
        evidence.momentum.validation.valid || throw(
            ArgumentError("valid momentum data has invalid validation metadata")
        )
        all(level -> !ismissing(level.ky_index), evidence.momentum.levels) || throw(
            ArgumentError("valid momentum data contains missing labels")
        )
    else
        isempty(evidence.momentum_counting) || throw(
            ArgumentError("invalid momentum data must not report counting")
        )
        all(level -> ismissing(level.ky_index), evidence.momentum.levels) || throw(
            ArgumentError("invalid momentum data contains fabricated labels")
        )
    end
    return normalized
end

function _render_fig2_momentum(momentum)
    header = "level\tprobability\tentanglement_energy\traw_charge\tky_index\tky\teigenphase\tphase_residual\tvalid"
    rows = (
        (
            level.level,
            level.probability,
            level.entanglement_energy,
            level.raw_charge,
            level.ky_index,
            level.ky,
            level.eigenphase,
            level.phase_residual,
            momentum.valid,
        ) for level in momentum.levels
    )
    return _render_tsv(header, rows)
end

function _complete_fig2_candidate!(
    spec,
    output,
    ledger,
    dimension,
    point,
    phi_y,
    candidate_id,
    evidence,
    relative_directory,
)
    directory = joinpath(output, relative_directory)
    for filename in FIG2_REQUIRED_CANDIDATE_FILES[1:7]
        path = joinpath(directory, filename)
        isfile(path) && filesize(path) > 0 || throw(
            ArgumentError("candidate output is missing required file $filename")
        )
    end
    normalized_sectors = _validate_fig2_evidence(evidence, point)
    _write_fig2_bytes(
        joinpath(directory, "momentum_entanglement_spectrum.tsv"),
        _render_fig2_momentum(evidence.momentum),
    )
    _write_fig2_bytes(
        joinpath(directory, "mixed_fidelity.tsv"),
        _render_tsv(
            "fidelity_to_previous\tvalid",
            [(evidence.fidelity_to_previous, evidence.fidelity_valid)],
        ),
    )
    candidate_metadata = Dict{String,Any}(
        "format" => "fqahc_fig2_candidate_v1",
        "dimension" => dimension,
        "point" => point,
        "phi_y" => phi_y,
        "candidate_id" => String(candidate_id),
        "converged" => evidence.converged,
        "valid" => evidence.valid,
        "energy_per_site" => evidence.energy_per_site,
        "raw_schmidt_polarization" => evidence.raw_schmidt_polarization,
        "sector_weights" => Dict(string(charge) => weight for (charge, weight) in normalized_sectors),
        "fidelity_to_previous" => evidence.fidelity_to_previous,
        "fidelity_valid" => evidence.fidelity_valid,
        "momentum_valid" => evidence.momentum.valid,
        "momentum_reason" => evidence.momentum.reason,
        "momentum_counting" => copy(evidence.momentum_counting),
        "momentum_convention" => string(evidence.momentum.momentum_convention),
        "restart_valid" => evidence.restart_valid,
        "provenance_valid" => evidence.provenance_valid,
        "reason" => evidence.reason,
        "manifest_sha256" => spec.sha256,
    )
    _write_fig2_toml(joinpath(directory, "candidate.toml"), candidate_metadata)

    checksums = Dict{String,String}()
    for filename in FIG2_REQUIRED_CANDIDATE_FILES
        path = joinpath(directory, filename)
        isfile(path) && filesize(path) > 0 || throw(
            ArgumentError("candidate output is missing required file $filename")
        )
        checksums[filename] = _fig2_file_sha256(path)
    end
    row = Dict{String,Any}(
        "complete" => true,
        "dimension" => dimension,
        "point" => point,
        "phi_y" => phi_y,
        "candidate_id" => String(candidate_id),
        "directory" => relative_directory,
        "state_sha256" => checksums["state.h5"],
        "checksums" => checksums,
    )
    push!(ledger["candidate"], row)
    _write_fig2_toml(joinpath(output, "ledger.toml"), ledger)
    return normalized_sectors
end

function _select_fig2_candidate(candidates, point)
    eligible = findall(candidate ->
        candidate.evidence.converged && candidate.evidence.valid &&
        isfinite(candidate.evidence.energy_per_site), candidates)
    isempty(eligible) && throw(WorkflowValidationError(
        "no converged valid Fig. 2 candidate is selectable"
    ))
    if point == 1
        index = first(sort(eligible; by=i -> (
            candidates[i].evidence.energy_per_site,
            candidates[i].candidate_id,
        )))
        return index, "minimum_zero_flux_energy"
    end
    fidelity_eligible = filter(i ->
        candidates[i].evidence.fidelity_valid &&
        isfinite(candidates[i].evidence.fidelity_to_previous) &&
        0 <= candidates[i].evidence.fidelity_to_previous <= 1,
        eligible,
    )
    isempty(fidelity_eligible) && throw(WorkflowValidationError(
        "no candidate has valid mixed fidelity to the previous selected state"
    ))
    index = first(sort(fidelity_eligible; by=i -> (
        -candidates[i].evidence.fidelity_to_previous,
        candidates[i].evidence.energy_per_site,
        candidates[i].candidate_id,
    )))
    return index, "maximum_mixed_fidelity"
end

function _selection_dict(selection::Fig2Selection)
    return Dict{String,Any}(
        "dimension" => selection.dimension,
        "point" => selection.point,
        "phi_y" => selection.phi_y,
        "candidate_id" => selection.candidate_id,
        "selection_basis" => selection.selection_basis,
        "directory" => selection.directory,
        "energy_per_site" => selection.energy_per_site,
        "raw_schmidt_polarization" => selection.raw_schmidt_polarization,
        "fidelity_to_previous" => selection.fidelity_to_previous,
        "fidelity_valid" => selection.fidelity_valid,
        "sector_weights" => Dict(string(k) => v for (k, v) in selection.sector_weights),
        "momentum_valid" => selection.momentum_valid,
        "momentum_reason" => selection.momentum_reason,
        "momentum_counting" => selection.momentum_counting,
        "restart_valid" => selection.restart_valid,
        "provenance_valid" => selection.provenance_valid,
    )
end

function _selection_from_dict(row)
    return Fig2Selection(
        Int(row["dimension"]),
        Int(row["point"]),
        Float64(row["phi_y"]),
        String(row["candidate_id"]),
        String(row["selection_basis"]),
        String(row["directory"]),
        Float64(row["energy_per_site"]),
        Float64(row["raw_schmidt_polarization"]),
        Float64(row["fidelity_to_previous"]),
        Bool(row["fidelity_valid"]),
        Dict(parse(Int, k) => Float64(v) for (k, v) in row["sector_weights"]),
        Bool(row["momentum_valid"]),
        String(row["momentum_reason"]),
        Int.(row["momentum_counting"]),
        Bool(row["restart_valid"]),
        Bool(row["provenance_valid"]),
    )
end

function _fig2_candidate_key(row)
    return (Int(row["dimension"]), Int(row["point"]), String(row["candidate_id"]))
end

function _fig2_selection_key(row)
    return (Int(row["dimension"]), Int(row["point"]))
end

function _fig2_unique_rows(rows, key_function, description)
    indexed = Dict{Any,Any}()
    for row in rows
        key = key_function(row)
        haskey(indexed, key) && throw(
            ArgumentError("duplicate $description row for $(repr(key))")
        )
        indexed[key] = row
    end
    return indexed
end

function _validate_persisted_fig2_candidate_files(root, row)
    get(row, "complete", false) === true || throw(
        ArgumentError("persisted Fig. 2 candidate is not marked complete")
    )
    dimension, point, candidate_id = _fig2_candidate_key(row)
    relative_directory = String(row["directory"])
    relative_directory == _candidate_relative_directory(
        dimension, point, candidate_id
    ) || throw(
        ArgumentError("persisted candidate directory does not match its identity")
    )
    checksums = get(row, "checksums", nothing)
    checksums isa AbstractDict || throw(
        ArgumentError("persisted candidate checksums are missing")
    )
    directory = joinpath(root, relative_directory)
    for filename in FIG2_REQUIRED_CANDIDATE_FILES
        haskey(checksums, filename) || throw(
            ArgumentError("persisted candidate checksum is missing for $filename")
        )
        path = joinpath(directory, filename)
        isfile(path) && filesize(path) > 0 || throw(
            ArgumentError("persisted candidate output is missing $filename")
        )
        _fig2_file_sha256(path) == String(checksums[filename]) || throw(
            ArgumentError("persisted candidate checksum mismatch for $filename")
        )
    end
    get(row, "state_sha256", "") == String(checksums["state.h5"]) || throw(
        ArgumentError("persisted candidate state checksum is inconsistent")
    )
    return nothing
end

function _validate_fig2_ledger_schedule(rows, dims, phis, description)
    allowed_dimensions = Set(dims)
    for row in rows
        dimension = Int(row["dimension"])
        point = Int(row["point"])
        dimension in allowed_dimensions || throw(
            ArgumentError("persisted $description dimension is not in the requested schedule")
        )
        1 <= point <= length(phis) || throw(
            ArgumentError("persisted $description point is not in the requested flux schedule")
        )
        Float64(row["phi_y"]) == phis[point] || throw(
            ArgumentError("persisted $description flux does not match the requested schedule")
        )
    end
    return nothing
end

function _persisted_fig2_candidate(
    spec,
    root,
    row,
    dimension,
    point,
    phi_y,
    candidate_id,
)
    get(row, "complete", false) === true || throw(
        ArgumentError("persisted Fig. 2 candidate is not marked complete")
    )
    relative_directory = String(row["directory"])
    directory = joinpath(root, relative_directory)
    metadata_path = joinpath(directory, "candidate.toml")
    isfile(metadata_path) || throw(
        ArgumentError("persisted Fig. 2 candidate metadata is missing")
    )
    metadata = TOML.parsefile(metadata_path)
    get(metadata, "format", "") == "fqahc_fig2_candidate_v1" || throw(
        ArgumentError("unsupported persisted Fig. 2 candidate format")
    )
    get(metadata, "manifest_sha256", "") == spec.sha256 || throw(
        ArgumentError("persisted candidate manifest checksum does not match")
    )
    (Int(metadata["dimension"]), Int(metadata["point"]),
     Float64(metadata["phi_y"]), String(metadata["candidate_id"])) ==
        (dimension, point, phi_y, candidate_id) || throw(
        ArgumentError("persisted candidate identity does not match the requested schedule")
    )
    sectors = _normalized_sector_weights(Dict(
        parse(Int, charge) => Float64(weight) for
        (charge, weight) in metadata["sector_weights"]
    ))
    evidence = (
        state=nothing,
        converged=Bool(metadata["converged"]),
        valid=Bool(metadata["valid"]),
        energy_per_site=Float64(metadata["energy_per_site"]),
        raw_schmidt_polarization=Float64(metadata["raw_schmidt_polarization"]),
        sector_weights=sectors,
        fidelity_to_previous=Float64(metadata["fidelity_to_previous"]),
        fidelity_valid=Bool(metadata["fidelity_valid"]),
        momentum=(
            valid=Bool(metadata["momentum_valid"]),
            reason=String(metadata["momentum_reason"]),
        ),
        momentum_counting=Int.(metadata["momentum_counting"]),
        restart_valid=Bool(metadata["restart_valid"]),
        provenance_valid=Bool(metadata["provenance_valid"]),
    )
    return (; candidate_id, relative_directory, evidence, sectors)
end

function _fig2_selection_from_candidate(
    candidate,
    dimension,
    point,
    phi_y,
    selection_basis,
)
    evidence = candidate.evidence
    return Fig2Selection(
        dimension,
        point,
        phi_y,
        candidate.candidate_id,
        selection_basis,
        candidate.relative_directory,
        evidence.energy_per_site,
        evidence.raw_schmidt_polarization,
        evidence.fidelity_to_previous,
        evidence.fidelity_valid,
        candidate.sectors,
        evidence.momentum.valid,
        evidence.momentum.reason,
        copy(evidence.momentum_counting),
        evidence.restart_valid,
        evidence.provenance_valid,
    )
end

function _write_fig2_pump_tables(spec, output, selections)
    ordered = sort(selections; by=row -> (row.dimension, row.point))
    sector_tolerance = _fig2_number(
        spec.data["convergence"], "sector_residual_tol"
    )
    reference = Dict{Int,Float64}()
    raw_rows = NamedTuple[]
    tracked_rows = NamedTuple[]
    previous = Dict{Int,Fig2Selection}()
    cumulative_shift = Dict{Int,Int}()
    path_valid = Dict{Int,Bool}()
    for row in ordered
        get!(reference, row.dimension, row.raw_schmidt_polarization)
        push!(raw_rows, (
            dimension=row.dimension,
            point=row.point,
            phi_y=row.phi_y,
            candidate_id=row.candidate_id,
            raw_schmidt_polarization=row.raw_schmidt_polarization,
            delta_raw_schmidt_polarization=
                row.raw_schmidt_polarization - reference[row.dimension],
        ))
        if !haskey(previous, row.dimension)
            cumulative_shift[row.dimension] = 0
            path_valid[row.dimension] = true
            push!(tracked_rows, (
                dimension=row.dimension,
                point=row.point,
                phi_y=row.phi_y,
                candidate_id=row.candidate_id,
                alignment_shift_to_previous_raw=0,
                cumulative_sector_shift=0,
                sector_distance_raw=0.0,
                valid=true,
                reason="reference charge-sector labeling",
            ))
        else
            prior = previous[row.dimension]
            alignment = best_sector_shift(prior.sector_weights, row.sector_weights)
            step_valid = row.fidelity_valid &&
                isfinite(row.fidelity_to_previous) &&
                alignment.distance <= sector_tolerance
            valid = path_valid[row.dimension] && step_valid
            valid && (cumulative_shift[row.dimension] -= alignment.shift)
            path_valid[row.dimension] = valid
            push!(tracked_rows, (
                dimension=row.dimension,
                point=row.point,
                phi_y=row.phi_y,
                candidate_id=row.candidate_id,
                alignment_shift_to_previous_raw=alignment.shift,
                cumulative_sector_shift=cumulative_shift[row.dimension],
                sector_distance_raw=alignment.distance,
                valid,
                reason=if valid
                    "derived from mixed-fidelity-supported sector correspondence"
                elseif !(row.fidelity_valid && isfinite(row.fidelity_to_previous))
                    "invalid because mixed fidelity is unavailable"
                elseif alignment.distance > sector_tolerance
                    "invalid because sector correspondence residual exceeds tolerance"
                else
                    "invalid because an earlier sector correspondence was invalid"
                end,
            ))
        end
        previous[row.dimension] = row
    end
    _write_fig2_bytes(
        joinpath(output, "pump_raw.tsv"),
        _render_tsv(
            "dimension\tpoint\tphi_y\tcandidate_id\traw_schmidt_polarization\tdelta_raw_schmidt_polarization",
            (Tuple(values(row)) for row in raw_rows),
        ),
    )
    _write_fig2_bytes(
        joinpath(output, "pump_sector_tracked.tsv"),
        _render_tsv(
            "dimension\tpoint\tphi_y\tcandidate_id\talignment_shift_to_previous_raw\tcumulative_sector_shift\tsector_distance_raw\tvalid\treason",
            (Tuple(values(row)) for row in tracked_rows),
        ),
    )
    return nothing
end

function _validate_fig2_schedule(dimensions, fluxes)
    dims = Int[]
    for dimension in dimensions
        dimension isa Integer && !(dimension isa Bool) && dimension > 0 || throw(
            ArgumentError("Fig. 2 dimensions must be positive integers")
        )
        push!(dims, Int(dimension))
    end
    length(unique(dims)) == length(dims) || throw(
        ArgumentError("Fig. 2 dimensions must be unique")
    )
    phis = Float64[]
    for flux in fluxes
        flux isa Real && !(flux isa Bool) && isfinite(flux) || throw(
            ArgumentError("Fig. 2 flux values must be finite real numbers")
        )
        push!(phis, Float64(flux))
    end
    issorted(phis) || throw(ArgumentError("Fig. 2 flux values must be sorted"))
    length(unique(phis)) == length(phis) || throw(
        ArgumentError("Fig. 2 flux values must be unique")
    )
    isempty(dims) == isempty(phis) || throw(
        ArgumentError("dimensions and fluxes must both be empty or both nonempty")
    )
    return dims, phis
end

function run_fig2_benchmark(
    spec::Fig2BenchmarkSpec,
    output::AbstractString;
    stage::AbstractString,
    dimensions,
    fluxes,
    operations::Fig2BenchmarkOperations=Fig2BenchmarkOperations(),
)
    isempty(strip(stage)) && throw(ArgumentError("benchmark stage must not be empty"))
    dims, phis = _validate_fig2_schedule(dimensions, fluxes)
    root = abspath(output)
    mkpath(root)
    _ensure_fig2_manifest_copy(spec, root)
    ledger = _load_fig2_ledger(spec, root, stage)
    persisted_candidates = _fig2_unique_rows(
        ledger["candidate"], _fig2_candidate_key, "candidate"
    )
    foreach(
        row -> _validate_persisted_fig2_candidate_files(root, row),
        values(persisted_candidates),
    )
    persisted_selections = _fig2_unique_rows(
        ledger["selection"], _fig2_selection_key, "selection"
    )
    _validate_fig2_ledger_schedule(
        values(persisted_candidates), dims, phis, "candidate"
    )
    _validate_fig2_ledger_schedule(
        values(persisted_selections), dims, phis, "selection"
    )
    selections = Fig2Selection[]
    for dimension in dims
        previous_state = nothing
        previous_selection = nothing
        for (point, phi_y) in enumerate(phis)
            selection_key = (dimension, point)
            if haskey(persisted_selections, selection_key)
                selection = _selection_from_dict(persisted_selections[selection_key])
                selection.phi_y == phi_y || throw(
                    ArgumentError("persisted selection flux does not match the requested schedule")
                )
                haskey(
                    persisted_candidates,
                    (dimension, point, selection.candidate_id),
                ) || throw(
                    ArgumentError("persisted selection has no complete candidate row")
                )
                push!(selections, selection)
                previous_selection = selection
                previous_state = nothing
                continue
            end
            if isnothing(previous_state) && !isnothing(previous_selection)
                checkpoint = joinpath(
                    root, previous_selection.directory, "state.h5"
                )
                previous_state = operations.load_state(
                    spec,
                    dimension,
                    previous_selection.point,
                    previous_selection,
                    checkpoint,
                )
            end
            candidate_ids = String.(operations.candidate_ids(
                spec, dimension, point, previous_state
            ))
            isempty(candidate_ids) && throw(
                WorkflowValidationError("candidate list is empty at point $point")
            )
            length(unique(candidate_ids)) == length(candidate_ids) || throw(
                ArgumentError("candidate IDs must be unique at each point")
            )
            candidates = NamedTuple[]
            for candidate_id in candidate_ids
                candidate_key = (dimension, point, candidate_id)
                if haskey(persisted_candidates, candidate_key)
                    push!(
                        candidates,
                        _persisted_fig2_candidate(
                            spec,
                            root,
                            persisted_candidates[candidate_key],
                            dimension,
                            point,
                            phi_y,
                            candidate_id,
                        ),
                    )
                    continue
                end
                relative_directory = _candidate_relative_directory(
                    dimension, point, candidate_id
                )
                candidate_directory = joinpath(root, relative_directory)
                evidence = operations.run_candidate(
                    spec,
                    dimension,
                    point,
                    phi_y,
                    candidate_id,
                    previous_state,
                    candidate_directory,
                )
                evidence isa Fig2CandidateEvidence || throw(
                    ArgumentError("candidate operation must return Fig2CandidateEvidence")
                )
                sectors = _complete_fig2_candidate!(
                    spec,
                    root,
                    ledger,
                    dimension,
                    point,
                    phi_y,
                    candidate_id,
                    evidence,
                    relative_directory,
                )
                push!(candidates, (;
                    candidate_id,
                    relative_directory,
                    evidence,
                    sectors,
                ))
                persisted_candidates[candidate_key] = last(ledger["candidate"])
            end
            selected_index, selection_basis = _select_fig2_candidate(candidates, point)
            selected = candidates[selected_index]
            selection = _fig2_selection_from_candidate(
                selected,
                dimension,
                point,
                phi_y,
                selection_basis,
            )
            push!(selections, selection)
            push!(ledger["selection"], _selection_dict(selection))
            persisted_selections[selection_key] = last(ledger["selection"])
            _write_fig2_toml(joinpath(root, "ledger.toml"), ledger)
            _write_fig2_pump_tables(spec, root, selections)
            previous_selection = selection
            previous_state = selected.evidence.state
            if isnothing(previous_state)
                checkpoint = joinpath(root, selection.directory, "state.h5")
                previous_state = operations.load_state(
                    spec, dimension, point, selection, checkpoint
                )
            end
        end
    end
    isempty(dims) && begin
        _write_fig2_toml(joinpath(root, "ledger.toml"), ledger)
        _write_fig2_pump_tables(spec, root, selections)
    end
    return Fig2BenchmarkRun(spec, String(stage), root, selections)
end

function _acceptance_row(pass, measured, tolerance, evidence, reason)
    return Dict{String,Any}(
        "pass" => Bool(pass),
        "measured" => measured,
        "tolerance" => tolerance,
        "evidence" => String(evidence),
        "reason" => String(reason),
    )
end

function _latest_dimension_selections(ledger)
    rows = [_selection_from_dict(row) for row in ledger["selection"]]
    isempty(rows) && return Fig2Selection[]
    dimension = maximum(row.dimension for row in rows)
    return sort!(filter(row -> row.dimension == dimension, rows); by=row -> row.point)
end

function _read_tracked_endpoint(path, dimension)
    lines = readlines(path)
    length(lines) >= 2 || return nothing
    fields = split.(lines[2:end], '\t')
    selected = filter(row -> parse(Int, row[1]) == dimension, fields)
    isempty(selected) && return nothing
    row = last(selected)
    return (
        shift=parse(Int, row[6]),
        distance=parse(Float64, row[7]),
        valid=parse(Bool, row[8]),
    )
end

function write_fig2_acceptance_report!(
    spec::Fig2BenchmarkSpec,
    output::AbstractString,
)
    root = abspath(output)
    _ensure_fig2_manifest_copy(spec, root)
    ledger = _load_fig2_ledger(
        spec,
        root,
        TOML.parsefile(joinpath(root, "ledger.toml"))["stage"],
    )
    selected = _latest_dimension_selections(ledger)
    report = Dict{String,Any}()
    if length(selected) >= 2
        measured = last(selected).raw_schmidt_polarization -
            first(selected).raw_schmidt_polarization
        expected = _fig2_number(spec.data, "expected_endpoint_pump")
        tolerance = _fig2_number(spec.data, "pump_tolerance")
        pass = abs(measured - expected) <= tolerance
        report["endpoint_pump"] = _acceptance_row(
            pass,
            measured,
            tolerance,
            "pump_raw.tsv",
            pass ? "raw endpoint pump matches the paper target" :
                "raw endpoint pump misses the paper target",
        )
    else
        report["endpoint_pump"] = _acceptance_row(
            false, "missing", _fig2_number(spec.data, "pump_tolerance"),
            "pump_raw.tsv", "insufficient raw pump points",
        )
    end

    tracked = isempty(selected) ? nothing : _read_tracked_endpoint(
        joinpath(root, "pump_sector_tracked.tsv"), last(selected).dimension
    )
    expected_shift = _fig2_integer(spec.data, "expected_sector_shift")
    sector_tolerance = _fig2_number(spec.data["convergence"], "sector_residual_tol")
    sector_pass = !isnothing(tracked) && tracked.valid &&
        tracked.shift == expected_shift && tracked.distance <= sector_tolerance
    report["sector_shift"] = _acceptance_row(
        sector_pass,
        isnothing(tracked) ? "missing" : tracked.shift,
        sector_tolerance,
        "pump_sector_tracked.tsv",
        sector_pass ? "overlap-supported sector shift matches one charge" :
            "sector shift or correspondence residual is invalid",
    )

    target_counting = Int.(spec.data["counting"])
    endpoint = isempty(selected) ? nothing : last(selected)
    momentum_pass = !isnothing(endpoint) && endpoint.momentum_valid &&
        endpoint.momentum_counting == target_counting
    momentum_reason = if isnothing(endpoint)
        "missing momentum evidence"
    elseif !endpoint.momentum_valid
        "invalid momentum evidence: $(endpoint.momentum_reason)"
    elseif endpoint.momentum_counting != target_counting
        "validated momentum counting does not match the fixed target"
    else
        "validated momentum counting matches the fixed target"
    end
    report["momentum_counting"] = _acceptance_row(
        momentum_pass,
        isnothing(endpoint) ? Int[] : endpoint.momentum_counting,
        0,
        isnothing(endpoint) ? "ledger.toml" :
            joinpath(endpoint.directory, "momentum_entanglement_spectrum.tsv"),
        momentum_reason,
    )

    curve_path = joinpath(root, "paper_curve_comparison.toml")
    curve_pass = false
    curve_measured = "missing"
    curve_reason = "paper curve digitization/comparison evidence is missing"
    if isfile(curve_path)
        curve = TOML.parsefile(curve_path)
        curve_pass = get(curve, "pass", false) === true
        curve_measured = get(curve, "maximum_residual", "missing")
        curve_reason = get(curve, "reason", "paper curve comparison supplied")
    end
    report["paper_curve"] = _acceptance_row(
        curve_pass, curve_measured, _fig2_number(spec.data, "pump_tolerance"),
        "paper_curve_comparison.toml", curve_reason,
    )

    selection_rows = [_selection_from_dict(row) for row in ledger["selection"]]
    dimensions = sort!(unique(row.dimension for row in selection_rows))
    convergence_pass = length(dimensions) >= 2
    convergence_measured = "missing"
    if convergence_pass
        previous = filter(row -> row.dimension == dimensions[end - 1], selection_rows)
        current = filter(row -> row.dimension == dimensions[end], selection_rows)
        if isempty(previous) || isempty(current)
            convergence_pass = false
        else
            convergence_measured = abs(last(current).energy_per_site - last(previous).energy_per_site)
            convergence_pass &= convergence_measured <=
                _fig2_number(spec.data["convergence"], "energy_per_site_change_tol")
        end
    end
    report["bond_convergence"] = _acceptance_row(
        convergence_pass,
        convergence_measured,
        _fig2_number(spec.data["convergence"], "energy_per_site_change_tol"),
        "ledger.toml",
        convergence_pass ? "successive bond dimensions satisfy energy tolerance" :
            "at least two converged bond dimensions are required",
    )

    restart_pass = !isempty(selected) && all(row -> row.restart_valid, selected)
    report["restart"] = _acceptance_row(
        restart_pass,
        count(row -> row.restart_valid, selected),
        length(selected),
        "ledger.toml",
        restart_pass ? "every selected checkpoint passed restart validation" :
            "at least one selected checkpoint lacks restart validation",
    )
    provenance_pass = !isempty(selected) &&
        all(row -> row.provenance_valid, selected) &&
        _fig2_file_sha256(joinpath(root, "manifest.toml")) == spec.sha256
    report["provenance"] = _acceptance_row(
        provenance_pass,
        provenance_pass ? spec.sha256 : "invalid",
        spec.sha256,
        "manifest.toml",
        provenance_pass ? "manifest checksum and candidate provenance are valid" :
            "manifest checksum or candidate provenance is invalid",
    )

    sort!(collect(keys(report))) == sort!(copy(spec.acceptance_rows)) || error(
        "acceptance report rows drifted from the immutable manifest"
    )
    _write_fig2_toml(joinpath(root, "acceptance.toml"), report)
    return report
end

function _default_fig2_candidate_ids(spec, dimension, point, previous_state)
    candidates = fig2_initial_candidates(spec.config)
    point == 1 && return [candidate.id for candidate in candidates]
    cold_count = Int(get(spec.data["pilot"], "cold_candidates", 2))
    bounded = first(candidates, min(cold_count, length(candidates)))
    return ["warm"; ["cold_$(candidate.id)" for candidate in bounded]]
end

function _default_fig2_load_state(spec, dimension, point, selection, checkpoint)
    config = with_flux(spec.config, selection.phi_y)
    return load_checkpoint(checkpoint, config)
end

function _fig2_candidate_pattern(spec, candidate_id)
    id = startswith(candidate_id, "cold_") ? candidate_id[6:end] : candidate_id
    matches = filter(candidate -> candidate.id == id, fig2_initial_candidates(spec.config))
    length(matches) == 1 || throw(ArgumentError("unknown Fig. 2 candidate ID $candidate_id"))
    return only(matches).occupied_sites
end

function _default_fig2_run_candidate(
    spec,
    dimension,
    point,
    phi_y,
    candidate_id,
    previous_state,
    candidate_directory,
)
    config = with_flux(spec.config, phi_y)
    if candidate_id == "warm"
        isnothing(previous_state) && throw(ArgumentError("warm candidate needs previous state"))
        psi = previous_state
        sites = siteinds(only, psi.AL)
        occupied_sites = default_occupied_sites(config)
    else
        occupied_sites = _fig2_candidate_pattern(spec, candidate_id)
        sites, _, psi = initial_infinite_mps(config; occupied_sites)
    end
    optimization = spec.data["optimization"]
    settings = SinglePointSettings(
        config,
        spec.model,
        [dimension],
        Float64(optimization["cutoff"]),
        Float64(optimization["vumps_tol"]),
        Float64(optimization["energy_tol"]),
        Float64(optimization["energy_mismatch_tol"]),
        1e-12,
        Float64(optimization["transfer_tol"]),
        4,
        Int(optimization["max_iterations"]),
        Int(optimization["stable_iterations"]),
        Base.Threads.nthreads(),
        candidate_directory,
        joinpath(candidate_directory, "state.h5"),
        nothing,
        collect(occupied_sites),
        0,
        true,
    )
    H = build_infinite_mpo(config, spec.model, sites)
    result = _run_prepared_point(
        settings,
        H,
        psi,
        SinglePointOperations();
        transfer_rng=Random.Xoshiro(_derived_seed(0, :fig2, dimension, point, candidate_id)),
    )
    entanglement = first(result.entanglements)
    sectors = Dict(sector.raw_charge => sector.weight for sector in entanglement.sectors)
    momentum = momentum_entanglement_data(
        result.optimization.psi,
        config;
        cut_x=1,
        residual_tol=Float64(optimization["momentum_residual_tol"]),
        transfer_tol=Float64(optimization["transfer_tol"]),
        rng_seed=_derived_seed(0, :momentum, dimension, point, candidate_id),
    )
    fidelity = if isnothing(previous_state)
        nothing
    else
        mixed_transfer_fidelity(
            previous_state,
            result.optimization.psi,
            config;
            tol=Float64(optimization["transfer_tol"]),
            rng=Random.Xoshiro(_derived_seed(0, :fidelity, dimension, point, candidate_id)),
        )
    end
    restart_valid = try
        load_checkpoint(joinpath(candidate_directory, "state.h5"), config)
        true
    catch
        false
    end
    reasons = [result.optimization.reason; result.observable_reasons]
    !isnothing(fidelity) && !fidelity.valid && push!(reasons, fidelity.reason)
    return Fig2CandidateEvidence(;
        state=result.optimization.psi,
        converged=result.optimization.converged,
        valid=result.valid,
        energy_per_site=result.energy.per_site,
        raw_schmidt_polarization=entanglement.raw_schmidt_polarization,
        sector_weights=sectors,
        fidelity_to_previous=isnothing(fidelity) ? NaN : fidelity.fidelity_cell,
        fidelity_valid=!isnothing(fidelity) && fidelity.valid,
        momentum,
        momentum_counting=Int[],
        restart_valid,
        provenance_valid=true,
        reason=join(reasons, "; "),
    )
end
