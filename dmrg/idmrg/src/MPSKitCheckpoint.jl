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
