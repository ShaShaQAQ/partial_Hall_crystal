const MPSKIT_CHECKPOINT_FORMAT = "mpskit_idmrg_checkpoint_v1"

struct MPSKitCheckpointData{S,M,P}
    state::S
    metadata::M
    provenance::P
end

function _mpskit_sector_descriptor(sector)
    charge_descriptor = if applicable(length, sector) && length(sector) == 2
        u1_charge = try
            TensorKit.charge(sector[2])
        catch
            nothing
        end
        isnothing(u1_charge) ? "" : ";u1=$(repr(u1_charge))"
    else
        ""
    end
    return "sector=$(repr(string(sector)))$charge_descriptor"
end

function _mpskit_space_descriptor(space)
    sector_rows = String[]
    for sector in TensorKit.sectors(space)
        push!(
            sector_rows,
            "$(_mpskit_sector_descriptor(sector));dim=$(TensorKit.dim(space, sector))",
        )
    end
    sort!(sector_rows)
    return join(
        [
            "dual=$(TensorKit.isdual(space))",
            "dim=$(TensorKit.dim(space))",
            sector_rows...,
        ],
        "|",
    )
end

function _mpskit_spaces_fingerprint(kind::AbstractString, spaces)
    io = IOBuffer()
    println(io, "mpskit_space_fingerprint_v1")
    println(io, "kind=", kind)
    println(io, "count=", length(spaces))
    for (site, space) in enumerate(spaces)
        println(io, "site=", site, ";", _mpskit_space_descriptor(space))
    end
    return bytes2hex(SHA.sha256(take!(io)))
end

function mpskit_state_space_fingerprints(
    state,
    c::InfiniteCylinderConfig,
)
    length(state) == sites_per_cell(c) || throw(
        ArgumentError("state and configuration unit-cell lengths do not match"),
    )
    physical_spaces = [MPSKit.physicalspace(state, site) for site in 1:length(state)]
    virtual_spaces = [
        MPSKit.right_virtualspace(state, site) for site in 1:length(state)
    ]
    expected_physical = mpskit_physical_spaces(c)
    physical_spaces == expected_physical || throw(
        ArgumentError("state does not use the configuration physical spaces"),
    )
    return (
        physical=_mpskit_spaces_fingerprint("physical", physical_spaces),
        virtual=_mpskit_spaces_fingerprint("right_virtual", virtual_spaces),
    )
end

function _mpskit_checkpoint_configuration(c::InfiniteCylinderConfig)
    return (
        geometry=string(c.geometry),
        Ly=c.Ly,
        Ny=c.Ny,
        x_period=c.x_period,
        filling_num=c.filling_num,
        filling_den=c.filling_den,
        phi_y=c.phi_y,
        sites_per_cell=sites_per_cell(c),
        particles_per_cell=particles_per_cell(c),
        charge_scale=charge_scale(c),
        configuration_signature=configuration_signature(c),
    )
end

function _mpskit_checkpoint_provenance(
    c::InfiniteCylinderConfig,
    fingerprints,
)
    backend = mpskit_backend_provenance()
    return (
        format=MPSKIT_CHECKPOINT_FORMAT,
        complete=true,
        backend=backend.backend,
        blocktensorkit_commit=backend.blocktensorkit_commit,
        mpskit_commit=backend.mpskit_commit,
        tensorkittensors_commit=backend.tensorkittensors_commit,
        configuration_signature=configuration_signature(c),
        physical_space_fingerprint=fingerprints.physical,
        virtual_space_fingerprint=fingerprints.virtual,
        julia_version=string(VERSION),
        blocktensorkit_version=backend.blocktensorkit_version,
        mpskit_version=backend.mpskit_version,
        tensorkit_version=backend.tensorkit_version,
        tensorkittensors_version=backend.tensorkittensors_version,
        jld2_version=backend.jld2_version,
    )
end

function _mpskit_checkpoint_write(
    path::AbstractString,
    state,
    c::InfiniteCylinderConfig,
    metadata,
    fingerprints,
)
    configuration = _mpskit_checkpoint_configuration(c)
    provenance = _mpskit_checkpoint_provenance(c, fingerprints)
    JLD2.jldopen(path, "w") do file
        for (name, value) in pairs(provenance)
            file[string(name)] = value
        end
        for (name, value) in pairs(configuration)
            file["configuration/$(string(name))"] = value
        end
        file["metadata"] = metadata
        file["state"] = state
    end
    return nothing
end

function _mpskit_checkpoint_missing(name::AbstractString)
    return CheckpointFormatError("MPSKit checkpoint is missing required key $name")
end

function _mpskit_checkpoint_read_value(file, name::AbstractString, expected_type::Type)
    haskey(file, name) || throw(_mpskit_checkpoint_missing(name))
    value = file[name]
    value isa expected_type &&
        (expected_type === Bool || !(value isa Bool && expected_type <: Integer)) || throw(
        CheckpointFormatError(
            "MPSKit checkpoint key $name has type $(typeof(value)); expected $expected_type",
        ),
    )
    return value
end

function _validate_mpskit_checkpoint_configuration(file, c::InfiniteCylinderConfig)
    requested = _mpskit_checkpoint_configuration(c)
    for (name, requested_value) in pairs(requested)
        key = "configuration/$(string(name))"
        haskey(file, key) || throw(_mpskit_checkpoint_missing(key))
        stored_value = file[key]
        isequal(stored_value, requested_value) || throw(
            CheckpointCompatibilityError(
                "$(string(name)) is $(repr(stored_value)) in the MPSKit checkpoint but $(repr(requested_value)) was requested",
            ),
        )
    end
    return nothing
end

function _validate_mpskit_checkpoint_backend(file)
    expected = mpskit_backend_provenance()
    checks = (
        backend=(expected.backend, "backend"),
        blocktensorkit_commit=(
            expected.blocktensorkit_commit,
            "BlockTensorKit commit",
        ),
        mpskit_commit=(expected.mpskit_commit, "MPSKit commit"),
        tensorkittensors_commit=(
            expected.tensorkittensors_commit,
            "TensorKitTensors commit",
        ),
    )
    for (name, (expected_value, label)) in pairs(checks)
        stored_value = _mpskit_checkpoint_read_value(
            file,
            string(name),
            AbstractString,
        )
        stored_value == expected_value || throw(
            CheckpointCompatibilityError(
                "$label is $(repr(stored_value)) in the checkpoint but $(repr(expected_value)) is required",
            ),
        )
    end
    return nothing
end

function _load_mpskit_checkpoint_file(file, c::InfiniteCylinderConfig)
    format = _mpskit_checkpoint_read_value(file, "format", AbstractString)
    format == MPSKIT_CHECKPOINT_FORMAT || throw(
        CheckpointFormatError(
            "unsupported MPSKit checkpoint format $(repr(format))",
        ),
    )
    complete = _mpskit_checkpoint_read_value(file, "complete", Bool)
    complete || throw(CheckpointFormatError("MPSKit checkpoint is not complete"))
    _validate_mpskit_checkpoint_backend(file)
    _validate_mpskit_checkpoint_configuration(file, c)

    haskey(file, "metadata") || throw(_mpskit_checkpoint_missing("metadata"))
    haskey(file, "state") || throw(_mpskit_checkpoint_missing("state"))
    metadata = file["metadata"]
    state = file["state"]
    fingerprints = try
        mpskit_state_space_fingerprints(state, c)
    catch error
        error isa InterruptException && rethrow()
        throw(
            CheckpointFormatError(
                "invalid MPSKit state spaces: $(sprint(showerror, error))",
            ),
        )
    end
    stored_physical = _mpskit_checkpoint_read_value(
        file,
        "physical_space_fingerprint",
        AbstractString,
    )
    stored_virtual = _mpskit_checkpoint_read_value(
        file,
        "virtual_space_fingerprint",
        AbstractString,
    )
    fingerprints.physical == stored_physical || throw(
        CheckpointFormatError("MPSKit physical-space fingerprint is inconsistent"),
    )
    fingerprints.virtual == stored_virtual || throw(
        CheckpointFormatError("MPSKit virtual-space fingerprint is inconsistent"),
    )
    provenance = _mpskit_checkpoint_provenance(c, fingerprints)
    return MPSKitCheckpointData(state, metadata, provenance)
end

function _mpskit_checkpoint_format_error(context::AbstractString, error)
    error isa InterruptException && rethrow()
    error isa Union{CheckpointFormatError,CheckpointCompatibilityError} && rethrow()
    return CheckpointFormatError(
        "$context: $(sprint(showerror, error)); expected a $MPSKIT_CHECKPOINT_FORMAT MPSKit checkpoint",
    )
end

function load_mpskit_checkpoint(
    path::AbstractString,
    c::InfiniteCylinderConfig,
)
    isfile(path) || throw(
        CheckpointFormatError("MPSKit checkpoint file does not exist: $path"),
    )
    try
        return JLD2.jldopen(path, "r") do file
            _load_mpskit_checkpoint_file(file, c)
        end
    catch error
        throw(_mpskit_checkpoint_format_error("could not load MPSKit checkpoint", error))
    end
end

function save_mpskit_checkpoint(
    path::AbstractString,
    state,
    c::InfiniteCylinderConfig,
    metadata,
)
    fingerprints = mpskit_state_space_fingerprints(state, c)
    return _atomic_replace(path) do temporary
        try
            _mpskit_checkpoint_write(
                temporary,
                state,
                c,
                metadata,
                fingerprints,
            )
            reopened = JLD2.jldopen(temporary, "r") do file
                _load_mpskit_checkpoint_file(file, c)
            end
            reopened.metadata == metadata || throw(
                CheckpointFormatError(
                    "MPSKit checkpoint metadata changed during atomic verification",
                ),
            )
        catch error
            throw(
                _mpskit_checkpoint_format_error(
                    "could not save and verify MPSKit checkpoint",
                    error,
                ),
            )
        end
    end
end

function _mpskit_restart_gate_value(data, key::String, description::String)
    haskey(data, key) || throw(ArgumentError(
        "$description evidence is missing $key",
    ))
    return data[key]
end

function _mpskit_restart_gate_bool(data, key::String, description::String)
    value = _mpskit_restart_gate_value(data, key, description)
    value isa Bool || throw(ArgumentError(
        "$description evidence $key must be Boolean",
    ))
    return value
end

function _mpskit_restart_gate_integer(data, key::String, description::String)
    value = _mpskit_restart_gate_value(data, key, description)
    value isa Integer && !(value isa Bool) || throw(ArgumentError(
        "$description evidence $key must be an integer",
    ))
    return Int(value)
end

function _mpskit_restart_gate_number(data, key::String, description::String)
    value = _mpskit_restart_gate_value(data, key, description)
    value isa Real && !(value isa Bool) || throw(ArgumentError(
        "$description evidence $key must be a real number",
    ))
    converted = Float64(value)
    isfinite(converted) || throw(ArgumentError(
        "$description evidence $key must be finite",
    ))
    return converted
end

function _mpskit_restart_gate_string(data, key::String, description::String)
    value = _mpskit_restart_gate_value(data, key, description)
    value isa AbstractString && !isempty(value) || throw(ArgumentError(
        "$description evidence $key must be a nonempty string",
    ))
    return String(value)
end

function _mpskit_restart_gate_vector(data, key::String, description::String)
    value = _mpskit_restart_gate_value(data, key, description)
    value isa AbstractVector || throw(ArgumentError(
        "$description evidence $key must be an array",
    ))
    converted = Float64[]
    for entry in value
        entry isa Real && !(entry isa Bool) || throw(ArgumentError(
            "$description evidence $key entries must be real numbers",
        ))
        number = Float64(entry)
        isfinite(number) || throw(ArgumentError(
            "$description evidence $key entries must be finite",
        ))
        push!(converted, number)
    end
    return converted
end

function _mpskit_restart_gate_sector_weights(
    data,
    key::String,
    description::String,
)
    rows = _mpskit_restart_gate_value(data, key, description)
    rows isa AbstractVector || throw(ArgumentError(
        "$description evidence $key must be an array of tables",
    ))
    weights = Dict{String,Float64}()
    for row in rows
        row isa AbstractDict || throw(ArgumentError(
            "$description evidence $key rows must be tables",
        ))
        sector = _mpskit_restart_gate_string(row, "sector", description)
        haskey(weights, sector) && throw(ArgumentError(
            "$description evidence $key has duplicate sector $sector",
        ))
        weights[sector] = _mpskit_restart_gate_number(
            row,
            "weight",
            description,
        )
    end
    isempty(weights) && throw(ArgumentError(
        "$description evidence $key must not be empty",
    ))
    return weights
end

function _mpskit_restart_gate_spectrum(
    data,
    key::String,
    description::String,
)
    rows = _mpskit_restart_gate_value(data, key, description)
    rows isa AbstractVector || throw(ArgumentError(
        "$description evidence $key must be an array of tables",
    ))
    spectrum = Dict{String,Vector{Float64}}()
    for row in rows
        row isa AbstractDict || throw(ArgumentError(
            "$description evidence $key rows must be tables",
        ))
        sector = _mpskit_restart_gate_string(row, "sector", description)
        haskey(spectrum, sector) && throw(ArgumentError(
            "$description evidence $key has duplicate sector $sector",
        ))
        spectrum[sector] = _mpskit_restart_gate_vector(
            row,
            "singular_values",
            description,
        )
    end
    isempty(spectrum) && throw(ArgumentError(
        "$description evidence $key must not be empty",
    ))
    return spectrum
end

function _mpskit_restart_gate_max_error(first_values, second_values)
    length(first_values) == length(second_values) || return Inf
    isempty(first_values) && return 0.0
    return maximum(abs.(first_values .- second_values))
end

function _mpskit_restart_gate_dict_error(first_values, second_values)
    keys(first_values) == keys(second_values) || return Inf
    errors = Float64[]
    for key in keys(first_values)
        first_value = first_values[key]
        second_value = second_values[key]
        if first_value isa AbstractVector
            push!(
                errors,
                _mpskit_restart_gate_max_error(first_value, second_value),
            )
        else
            push!(errors, abs(first_value - second_value))
        end
    end
    return isempty(errors) ? 0.0 : maximum(errors)
end

function write_fig2_restart_gate!(
    spec::Fig2BenchmarkSpec,
    output::AbstractString;
    before_path::AbstractString,
    after_path::AbstractString,
    save_process_exit_code::Integer,
    resume_process_exit_code::Integer,
)
    before_file = abspath(before_path)
    after_file = abspath(after_path)
    before_file != after_file || throw(ArgumentError(
        "restart-gate save and resume evidence paths must be distinct",
    ))
    isfile(before_file) || throw(ArgumentError(
        "restart-gate save evidence does not exist: $before_file",
    ))
    isfile(after_file) || throw(ArgumentError(
        "restart-gate resume evidence does not exist: $after_file",
    ))
    save_process_exit_code isa Bool && throw(ArgumentError(
        "save process exit code must be an integer",
    ))
    resume_process_exit_code isa Bool && throw(ArgumentError(
        "resume process exit code must be an integer",
    ))

    snapshot = _fig2_validated_snapshot(spec)
    contract = snapshot["restart_gate"]
    before = TOML.parsefile(before_file)
    after = TOML.parsefile(after_file)
    save_process_id = _mpskit_restart_gate_string(
        before,
        "process_id",
        "save-process",
    )
    resume_process_id = _mpskit_restart_gate_string(
        after,
        "process_id",
        "resume-process",
    )
    save_process_id != resume_process_id || throw(ArgumentError(
        "restart-gate save and resume process IDs must be distinct",
    ))

    pre_maxlinkdim = _mpskit_restart_gate_integer(
        after,
        "pre_maxlinkdim",
        "resume-process",
    )
    no_expansion_maxlinkdim = _mpskit_restart_gate_integer(
        after,
        "no_expansion_post_maxlinkdim",
        "resume-process",
    )
    post_maxlinkdim = _mpskit_restart_gate_integer(
        after,
        "post_maxlinkdim",
        "resume-process",
    )
    saved_maxlinkdim = _mpskit_restart_gate_integer(
        before,
        "maxlinkdim",
        "save-process",
    )
    energy_error = abs(
        _mpskit_restart_gate_number(
            after,
            "pre_energy_per_site",
            "resume-process",
        ) -
        _mpskit_restart_gate_number(
            before,
            "energy_per_site",
            "save-process",
        ),
    )
    density_error = _mpskit_restart_gate_max_error(
        _mpskit_restart_gate_vector(
            after,
            "pre_densities",
            "resume-process",
        ),
        _mpskit_restart_gate_vector(
            before,
            "densities",
            "save-process",
        ),
    )
    sector_error = _mpskit_restart_gate_dict_error(
        _mpskit_restart_gate_sector_weights(
            after,
            "pre_sector_weights",
            "resume-process",
        ),
        _mpskit_restart_gate_sector_weights(
            before,
            "sector_weights",
            "save-process",
        ),
    )
    spectrum_error = _mpskit_restart_gate_dict_error(
        _mpskit_restart_gate_spectrum(
            after,
            "pre_entanglement_spectrum",
            "resume-process",
        ),
        _mpskit_restart_gate_spectrum(
            before,
            "entanglement_spectrum",
            "save-process",
        ),
    )
    energy_tolerance = Float64(contract["energy_tolerance"])
    density_tolerance = Float64(contract["density_tolerance"])

    evidence_directory = dirname(before_file)
    dirname(after_file) == evidence_directory || throw(ArgumentError(
        "restart-gate save and resume evidence must share a directory",
    ))
    checkpoint = joinpath(evidence_directory, "state.h5")
    resumed_checkpoint = joinpath(evidence_directory, "state-resumed.h5")
    next_flux_checkpoint = joinpath(evidence_directory, "state-next-flux.h5")
    all(isfile, (checkpoint, resumed_checkpoint, next_flux_checkpoint)) ||
        throw(ArgumentError(
            "restart-gate checkpoint evidence is incomplete",
        ))
    checkpoint_sha256 = _fig2_file_sha256(checkpoint)
    resumed_checkpoint_sha256 = _fig2_file_sha256(resumed_checkpoint)
    next_flux_checkpoint_sha256 = _fig2_file_sha256(next_flux_checkpoint)
    checkpoint_integrity =
        checkpoint_sha256 == _mpskit_restart_gate_string(
            before,
            "checkpoint_sha256",
            "save-process",
        ) == _mpskit_restart_gate_string(
            after,
            "original_checkpoint_sha256",
            "resume-process",
        ) &&
        resumed_checkpoint_sha256 == _mpskit_restart_gate_string(
            after,
            "resumed_checkpoint_sha256",
            "resume-process",
        ) &&
        next_flux_checkpoint_sha256 == _mpskit_restart_gate_string(
            after,
            "next_flux_checkpoint_sha256",
            "resume-process",
        )

    site_spaces_reproduced =
        _mpskit_restart_gate_string(
            before,
            "configuration_signature",
            "save-process",
        ) == _mpskit_restart_gate_string(
            after,
            "configuration_signature",
            "resume-process",
        ) &&
        _mpskit_restart_gate_string(
            before,
            "physical_space_fingerprint",
            "save-process",
        ) == _mpskit_restart_gate_string(
            after,
            "physical_space_fingerprint",
            "resume-process",
        ) &&
        _mpskit_restart_gate_string(
            before,
            "virtual_space_fingerprint",
            "save-process",
        ) == _mpskit_restart_gate_string(
            after,
            "virtual_space_fingerprint",
            "resume-process",
        )
    next_flux_phi_y = _mpskit_restart_gate_number(
        after,
        "next_flux_phi_y",
        "resume-process",
    )
    checks = Dict{String,Bool}(
        "load_completed" =>
            _mpskit_restart_gate_bool(
                after,
                "load_completed",
                "resume-process",
            ) && checkpoint_integrity,
        "no_expansion_iteration_completed" =>
            _mpskit_restart_gate_bool(
                after,
                "no_expansion_iteration_completed",
                "resume-process",
            ) &&
            _mpskit_restart_gate_integer(
                after,
                "refinement_iterations",
                "resume-process",
            ) == 1 &&
            no_expansion_maxlinkdim == pre_maxlinkdim,
        "controlled_expansion_completed" =>
            _mpskit_restart_gate_bool(
                after,
                "controlled_expansion_completed",
                "resume-process",
            ) &&
            saved_maxlinkdim == pre_maxlinkdim == 4 &&
            post_maxlinkdim == 8,
        "energy_reproduced" => energy_error <= energy_tolerance,
        "density_reproduced" => density_error <= density_tolerance,
        "schmidt_sectors_reproduced" => sector_error <= density_tolerance,
        "entanglement_spectrum_reproduced" =>
            spectrum_error <= density_tolerance,
        "site_indices_reproduced" => site_spaces_reproduced,
        "next_flux_iteration_completed" =>
            _mpskit_restart_gate_bool(
                after,
                "next_flux_iteration_completed",
                "resume-process",
            ) &&
            _mpskit_restart_gate_integer(
                after,
                "next_flux_iterations",
                "resume-process",
            ) == 1 &&
            next_flux_phi_y == 0.1,
    )
    required_checks = String.(contract["required_checks"])
    Set(keys(checks)) == Set(required_checks) || error(
        "restart-gate writer and manifest required checks disagree",
    )
    pass = Int(save_process_exit_code) == 0 &&
        Int(resume_process_exit_code) == 0 &&
        all(checks[check] for check in required_checks)

    active_project, project_manifest = _fig2_project_manifest()
    backend = mpskit_backend_provenance()
    gate = Dict{String,Any}(
        "format" => String(contract["format"]),
        "pass" => pass,
        "manifest_sha256" => spec.sha256,
        "git_commit" => _fig2_repository_commit(),
        "julia_version" => string(VERSION),
        "project_manifest_sha256" => _fig2_file_sha256(project_manifest),
        "backend_commit" => backend.mpskit_commit,
        "save_process_exit_code" => Int(save_process_exit_code),
        "resume_process_exit_code" => Int(resume_process_exit_code),
        "save_process_id" => "save:$save_process_id",
        "resume_process_id" => "resume:$resume_process_id",
        "pre_maxlinkdim" => pre_maxlinkdim,
        "post_maxlinkdim" => post_maxlinkdim,
        "next_flux_phi_y" => next_flux_phi_y,
        "energy_absolute_error" => energy_error,
        "energy_tolerance" => energy_tolerance,
        "density_max_absolute_error" => density_error,
        "density_tolerance" => density_tolerance,
        "schmidt_sector_max_absolute_error" => sector_error,
        "entanglement_spectrum_max_absolute_error" => spectrum_error,
        "save_evidence_sha256" => _fig2_file_sha256(before_file),
        "resume_evidence_sha256" => _fig2_file_sha256(after_file),
        "checkpoint_sha256" => checkpoint_sha256,
        "resumed_checkpoint_sha256" => resumed_checkpoint_sha256,
        "next_flux_checkpoint_sha256" => next_flux_checkpoint_sha256,
        "active_project" => active_project,
        "project_manifest" => project_manifest,
    )
    merge!(gate, checks)
    root = abspath(output)
    mkpath(root)
    gate_path = joinpath(root, String(contract["evidence_filename"]))
    _write_fig2_toml(gate_path, gate)
    return TOML.parsefile(gate_path)
end
