const CHECKPOINT_FORMAT = "infinite_cylinder_vumps_v1"
const TEXT_OUTPUT_FORMAT = "infinite_cylinder_vumps_text_v1"
const ITENSOR_INFINITE_MPS_COMMIT = "765f2777703bc1138b009adbed1b97bde1973402"

struct CheckpointFormatError <: Exception
    message::String
end

function Base.showerror(io::IO, error::CheckpointFormatError)
    return print(io, "invalid checkpoint format: ", error.message)
end

struct CheckpointCompatibilityError <: Exception
    message::String
end

function Base.showerror(io::IO, error::CheckpointCompatibilityError)
    return print(io, "checkpoint configuration mismatch: ", error.message)
end

const EXPECTED_CHECKPOINT_IO_ERRORS = Union{
    HDF5.API.H5Error,
    Base.IOError,
    SystemError,
    EOFError,
    ArgumentError,
    KeyError,
    BoundsError,
}

function _checkpoint_format_error(context::AbstractString, error)
    error isa InterruptException && rethrow()
    error isa EXPECTED_CHECKPOINT_IO_ERRORS || rethrow()
    return CheckpointFormatError("$context: $(sprint(showerror, error))")
end

function _atomic_replace(writer::Function, path::AbstractString)
    destination = abspath(path)
    mkpath(dirname(destination))
    temporary, handle = mktemp(dirname(destination); cleanup=false)
    close(handle)
    committed = false
    try
        writer(temporary)
        mv(temporary, destination; force=true)
        committed = true
    finally
        committed || rm(temporary; force=true)
    end
    return String(path)
end

function _isometry_residual(tensor, retained_index)
    gram = tensor * dag(prime(tensor, retained_index))
    identity = ITensors.denseblocks(delta(inds(gram)...))
    return norm(gram - identity)
end

function _checkpoint_metadata(c::InfiniteCylinderConfig)
    return (
        format=CHECKPOINT_FORMAT,
        backend_commit=ITENSOR_INFINITE_MPS_COMMIT,
        Ly=c.Ly,
        x_period=c.x_period,
        filling_num=c.filling_num,
        filling_den=c.filling_den,
        phi_y=c.phi_y,
        sites_per_cell=sites_per_cell(c),
        particles_per_cell=particles_per_cell(c),
        charge_scale=charge_scale(c),
        configuration_signature=_configuration_signature(c),
    )
end

function _validate_checkpoint_state(
    psi::InfiniteCanonicalMPS, c::InfiniteCylinderConfig
)
    expected_sites = sites_per_cell(c)
    nsites(psi) == expected_sites || throw(
        CheckpointFormatError("state cell size does not match checkpoint metadata")
    )
    all(component -> nsites(component) == expected_sites, (psi.AL, psi.C, psi.AR)) ||
        throw(CheckpointFormatError("canonical state components have unequal cell sizes"))
    all(
        component -> translator(component) === translator(psi.AL),
        (psi.C, psi.AR),
    ) || throw(
        CheckpointFormatError("canonical state components use different translators")
    )

    signature = _configuration_signature(c)
    left_sites = try
        siteinds(only, psi.AL)
    catch error
        throw(_checkpoint_format_error("left-canonical tensors have invalid site indices", error))
    end
    right_sites = try
        siteinds(only, psi.AR)
    catch error
        throw(_checkpoint_format_error("right-canonical tensors have invalid site indices", error))
    end
    all(hastags(site, signature) for site in left_sites) || throw(
        CheckpointFormatError("state site indices do not carry the configuration signature")
    )
    left_sites == right_sites || throw(
        CheckpointFormatError("left- and right-canonical site indices do not match")
    )

    try
        validate_centered_site_charges(c, left_sites)
        for site in 1:expected_sites
            order(psi.AL[site]) == 3 ||
                throw(ArgumentError("AL[$site] is not a rank-three site tensor"))
            order(psi.AR[site]) == 3 ||
                throw(ArgumentError("AR[$site] is not a rank-three site tensor"))
            order(psi.C[site]) == 2 ||
                throw(ArgumentError("C[$site] is not a rank-two bond tensor"))
            length(commoninds(psi.C[site], psi.AL[site])) == 1 || throw(
                ArgumentError("C[$site] is not connected to AL[$site]")
            )
            length(commoninds(psi.C[site], psi.AR[site])) == 1 || throw(
                ArgumentError("C[$site] is not connected to AR[$site]")
            )
            left_link = only(commoninds(psi.AL[site - 1], psi.AL[site]))
            right_link = only(commoninds(psi.AL[site], psi.AL[site + 1]))
            left_link in inds(psi.C[site - 1]) || throw(
                ArgumentError("C[$(site - 1)] is disconnected from AL[$site]")
            )
            right_link in inds(psi.C[site]) || throw(
                ArgumentError("C[$site] is disconnected from the AL link chain")
            )

            right_left_link = only(commoninds(psi.AR[site - 1], psi.AR[site]))
            right_right_link = only(commoninds(psi.AR[site], psi.AR[site + 1]))
            right_left_link in inds(psi.C[site - 1]) || throw(
                ArgumentError("C[$(site - 1)] is disconnected from the AR link chain")
            )
            right_right_link in inds(psi.C[site]) || throw(
                ArgumentError("C[$site] is disconnected from AR[$site]")
            )

            _isometry_residual(psi.AL[site], right_link) <= 1e-10 || throw(
                ArgumentError("AL[$site] is not left-isometric to tolerance 1e-10")
            )
            _isometry_residual(psi.AR[site], right_left_link) <= 1e-10 || throw(
                ArgumentError("AR[$site] is not right-isometric to tolerance 1e-10")
            )

            left_center = psi.AL[site] * psi.C[site]
            right_center = psi.C[site - 1] * psi.AR[site]
            center_norms = (norm(left_center), norm(right_center), norm(psi.C[site]))
            all(value -> isfinite(value) && value > 0, center_norms) || throw(
                ArgumentError("canonical center tensors must have finite positive norm")
            )
            center_scale = max(norm(left_center), norm(right_center), 1.0)
            norm(left_center - right_center) <= 1e-10 * center_scale || throw(
                ArgumentError(
                    "AL[$site] * C[$site] does not match C[$(site - 1)] * AR[$site] to relative tolerance 1e-10"
                ),
            )
        end
    catch error
        error isa CheckpointFormatError && rethrow()
        throw(_checkpoint_format_error("state does not form a valid canonical cell", error))
    end
    return nothing
end

function save_checkpoint(
    path::AbstractString,
    psi::InfiniteCanonicalMPS,
    c::InfiniteCylinderConfig,
)
    _validate_checkpoint_state(psi, c)
    metadata = _checkpoint_metadata(c)
    return _atomic_replace(path) do temporary
        h5open(temporary, "w") do file
            attrs = attributes(file)
            for (name, value) in pairs(metadata)
                attrs[string(name)] = value
            end
            write(file, "state", psi)
        end
    end
end

function _read_checkpoint_attribute(attrs, name::String, expected_type::Type)
    haskey(attrs, name) ||
        throw(CheckpointFormatError("missing required attribute $name"))
    value = try
        read(attrs[name])
    catch error
        throw(_checkpoint_format_error("could not read attribute $name", error))
    end
    value isa expected_type && !(value isa Bool && expected_type <: Integer) || throw(
        CheckpointFormatError(
            "attribute $name has type $(typeof(value)); expected $expected_type"
        ),
    )
    return value
end

function _read_checkpoint_metadata(file)
    attrs = attributes(file)
    return (
        format=_read_checkpoint_attribute(attrs, "format", AbstractString),
        backend_commit=_read_checkpoint_attribute(
            attrs, "backend_commit", AbstractString
        ),
        Ly=_read_checkpoint_attribute(attrs, "Ly", Integer),
        x_period=_read_checkpoint_attribute(attrs, "x_period", Integer),
        filling_num=_read_checkpoint_attribute(attrs, "filling_num", Integer),
        filling_den=_read_checkpoint_attribute(attrs, "filling_den", Integer),
        phi_y=_read_checkpoint_attribute(attrs, "phi_y", Float64),
        sites_per_cell=_read_checkpoint_attribute(
            attrs, "sites_per_cell", Integer
        ),
        particles_per_cell=_read_checkpoint_attribute(
            attrs, "particles_per_cell", Integer
        ),
        charge_scale=_read_checkpoint_attribute(attrs, "charge_scale", Integer),
        configuration_signature=_read_checkpoint_attribute(
            attrs, "configuration_signature", AbstractString
        ),
    )
end

function _stored_checkpoint_configuration(metadata)
    metadata.format == CHECKPOINT_FORMAT || throw(
        CheckpointFormatError("unsupported format $(repr(metadata.format))")
    )
    metadata.backend_commit == ITENSOR_INFINITE_MPS_COMMIT || throw(
        CheckpointFormatError(
            "unsupported backend commit $(repr(metadata.backend_commit))"
        ),
    )
    stored = try
        InfiniteCylinderConfig(
            Int(metadata.Ly),
            Int(metadata.x_period),
            Int(metadata.filling_num),
            Int(metadata.filling_den),
            metadata.phi_y,
        )
    catch error
        throw(_checkpoint_format_error("stored configuration is invalid", error))
    end
    derived = _checkpoint_metadata(stored)
    for name in (
        :sites_per_cell,
        :particles_per_cell,
        :charge_scale,
        :configuration_signature,
    )
        getproperty(metadata, name) == getproperty(derived, name) || throw(
            CheckpointFormatError("attribute $name is inconsistent with the configuration")
        )
    end
    return stored
end

function _validate_checkpoint_compatibility(
    stored::InfiniteCylinderConfig, requested::InfiniteCylinderConfig
)
    for name in (:Ly, :x_period, :filling_num, :filling_den, :phi_y)
        stored_value = getproperty(stored, name)
        requested_value = getproperty(requested, name)
        isequal(stored_value, requested_value) || throw(
            CheckpointCompatibilityError(
                "$name is $(repr(stored_value)) in the checkpoint but $(repr(requested_value)) was requested"
            ),
        )
    end
    return nothing
end

function load_checkpoint(path::AbstractString, c::InfiniteCylinderConfig)
    isfile(path) ||
        throw(CheckpointFormatError("checkpoint file does not exist: $path"))
    try
        return h5open(path, "r") do file
            metadata = _read_checkpoint_metadata(file)
            stored = _stored_checkpoint_configuration(metadata)
            _validate_checkpoint_compatibility(stored, c)
            haskey(file, "state") ||
                throw(CheckpointFormatError("missing canonical state group"))
            psi = try
                read(file, "state", InfiniteCanonicalMPS)
            catch error
                throw(_checkpoint_format_error("could not read canonical state", error))
            end
            _validate_checkpoint_state(psi, c)
            return psi
        end
    catch error
        error isa Union{CheckpointFormatError,CheckpointCompatibilityError} &&
            rethrow()
        throw(_checkpoint_format_error("could not open HDF5 checkpoint", error))
    end
end

_output_value(value::AbstractFloat) = string(value)
_output_value(value) = string(value)

function _render_tsv(header::AbstractString, rows)
    io = IOBuffer()
    println(io, header)
    for row in rows
        println(io, join(_output_value.(row), '\t'))
    end
    return String(take!(io))
end

function _toml_value(value::NamedTuple)
    return Dict(string(name) => _toml_value(item) for (name, item) in pairs(value))
end

function _toml_value(value::AbstractDict)
    return Dict(string(name) => _toml_value(item) for (name, item) in pairs(value))
end

_toml_value(value::Tuple) = [_toml_value(item) for item in value]
_toml_value(value::AbstractVector) = [_toml_value(item) for item in value]
_toml_value(value) = value

function _render_summary(summary)
    io = IOBuffer()
    TOML.print(io, summary; sorted=true)
    return String(take!(io))
end

function _validated_output_shapes(
    c::InfiniteCylinderConfig,
    result::VUMPSResult,
    densities,
    entanglements,
    transfer::NeutralTransferData,
)
    _validate_checkpoint_state(result.psi, c)
    length(densities) == sites_per_cell(c) || throw(
        ArgumentError("density output must contain one row per reference-cell site")
    )
    for (site, row) in enumerate(densities)
        row.site == site && (row.x, row.y) == site_coordinates(c, site) || throw(
            ArgumentError("density output row $site has inconsistent coordinates")
        )
    end
    all(data -> 1 <= data.cut_x <= c.x_period, entanglements) || throw(
        ArgumentError("entanglement output contains a cut outside the reference cell")
    )
    length(transfer.eigenvalues) == length(transfer.phases) ==
        length(transfer.residual_norms) || throw(
        ArgumentError("transfer output arrays must have equal lengths")
    )
    0 <= transfer.converged <= length(transfer.eigenvalues) || throw(
        ArgumentError("transfer converged count is outside the spectrum")
    )
    return nothing
end

function _summary_data(
    c::InfiniteCylinderConfig,
    result::VUMPSResult,
    energy::EnergyData,
    densities,
    entanglements,
    transfer::NeutralTransferData,
    optimization,
    dependencies,
)
    energy_valid = result.converged && all(
        isfinite,
        (energy.per_cell, energy.per_x, energy.per_unit_cell, energy.per_site),
    )
    density_valid = result.converged && all(row -> isfinite(row.density), densities)
    entanglement_valid = result.converged && !isempty(entanglements) && all(
        data ->
            isfinite(data.entropy) &&
            isfinite(data.raw_schmidt_polarization) &&
            !isempty(data.levels) &&
            !isempty(data.sectors) &&
            all(
                level -> all(
                    isfinite,
                    (
                        level.singular_value,
                        level.probability,
                        level.entanglement_energy,
                        level.physical_charge,
                    ),
                ),
                data.levels,
            ) &&
            all(
                sector -> all(isfinite, (sector.physical_charge, sector.weight)),
                data.sectors,
            ),
        entanglements,
    )
    transfer_valid =
        result.converged &&
        transfer.valid &&
        transfer.converged >= 2 &&
        all(isfinite, transfer.eigenvalues) &&
        all(isfinite, transfer.phases) &&
        all(isfinite, transfer.residual_norms) &&
        all(isfinite, (transfer.ratio, transfer.xi_cell, transfer.xi_x))
    valid = energy_valid && density_valid && entanglement_valid && transfer_valid
    maybe(validity, value) = validity ? value : NaN

    dependency_table = _toml_value(dependencies)
    dependency_table["itensor_infinite_mps_commit"] = ITENSOR_INFINITE_MPS_COMMIT
    entanglement_summary = [
        Dict(
            "cut_x" => data.cut_x,
            "bond" => data.bond,
            "valid" => entanglement_valid,
            "entropy" => maybe(entanglement_valid, data.entropy),
            "raw_schmidt_polarization" => maybe(
                entanglement_valid, data.raw_schmidt_polarization
            ),
        ) for data in entanglements
    ]
    return Dict{String,Any}(
        "format" => TEXT_OUTPUT_FORMAT,
        "algorithm" => "VUMPS",
        "valid" => valid,
        "converged" => result.converged,
        "reason" => result.reason,
        "configuration" => Dict(
            "Ly" => c.Ly,
            "x_period" => c.x_period,
            "filling_num" => c.filling_num,
            "filling_den" => c.filling_den,
            "phi_y" => c.phi_y,
            "sites_per_cell" => sites_per_cell(c),
            "particles_per_cell" => particles_per_cell(c),
            "charge_scale" => charge_scale(c),
            "signature" => _configuration_signature(c),
        ),
        "energy" => Dict(
            "valid" => energy_valid,
            "per_cell" => maybe(energy_valid, energy.per_cell),
            "per_x" => maybe(energy_valid, energy.per_x),
            "per_unit_cell" => maybe(energy_valid, energy.per_unit_cell),
            "per_site" => maybe(energy_valid, energy.per_site),
        ),
        "observables" => Dict(
            "density_valid" => density_valid,
            "entanglement_valid" => entanglement_valid,
            "transfer_valid" => transfer_valid,
            "transfer_reason" => transfer.reason,
            "neutral_transfer_ratio" => maybe(transfer_valid, transfer.ratio),
            "correlation_length_neutral_cell" => maybe(transfer_valid, transfer.xi_cell),
            "correlation_length_neutral_x" => maybe(transfer_valid, transfer.xi_x),
        ),
        "entanglement" => entanglement_summary,
        "optimization" => _toml_value(optimization),
        "dependencies" => dependency_table,
    )
end

function write_output_files(
    directory::AbstractString,
    c::InfiniteCylinderConfig,
    result::VUMPSResult,
    energy::EnergyData,
    densities::AbstractVector{<:DensityRow},
    entanglements::AbstractVector{<:EntanglementData},
    transfer::NeutralTransferData;
    optimization=(;),
    dependencies=(;),
)
    _validated_output_shapes(c, result, densities, entanglements, transfer)
    summary = _summary_data(
        c,
        result,
        energy,
        densities,
        entanglements,
        transfer,
        optimization,
        dependencies,
    )
    transfer_valid = summary["observables"]["transfer_valid"]
    density_valid = summary["observables"]["density_valid"]
    entanglement_valid = summary["observables"]["entanglement_valid"]
    measured(validity, value) = validity ? value : NaN

    contents = [
        "summary.toml" => _render_summary(summary),
        "convergence.tsv" => _render_tsv(
            "stage\titeration\tmaxlinkdim\tenergy_left\tenergy_right\tenergy_mismatch\tdelta_energy\teps_left\teps_right\tprecision_error\telapsed_seconds\tconverged",
            (
                (
                    record.stage,
                    record.iteration,
                    record.maxlinkdim,
                    record.energy_left,
                    record.energy_right,
                    record.energy_mismatch,
                    record.delta_energy,
                    record.eps_left,
                    record.eps_right,
                    record.precision_error,
                    record.elapsed_seconds,
                    record.converged,
                ) for record in result.records
            ),
        ),
        "density.tsv" => _render_tsv(
            "site\tx\ty\tdensity",
            (
                (row.site, row.x, row.y, measured(density_valid, row.density)) for
                row in densities
            ),
        ),
        "entanglement_spectrum.tsv" => _render_tsv(
            "cut_x\tbond\tlevel\tsingular_value\tprobability\tentanglement_energy\tqn\traw_charge\tphysical_charge",
            (
                (
                    level.cut_x,
                    level.bond,
                    level.level,
                    measured(entanglement_valid, level.singular_value),
                    measured(entanglement_valid, level.probability),
                    measured(entanglement_valid, level.entanglement_energy),
                    level.qn,
                    level.raw_charge,
                    measured(entanglement_valid, level.physical_charge),
                ) for data in entanglements for level in data.levels
            ),
        ),
        "schmidt_sectors.tsv" => _render_tsv(
            "cut_x\tbond\tqn\traw_charge\tphysical_charge\tweight",
            (
                (
                    sector.cut_x,
                    sector.bond,
                    sector.qn,
                    sector.raw_charge,
                    measured(entanglement_valid, sector.physical_charge),
                    measured(entanglement_valid, sector.weight),
                ) for data in entanglements for sector in data.sectors
            ),
        ),
        "transfer_spectrum.tsv" => _render_tsv(
            "level\teigenvalue_real\teigenvalue_imag\tmagnitude\tphase\tresidual_norm\tconverged\tvalid",
            (
                (
                    level,
                    real(value),
                    imag(value),
                    abs(value),
                    transfer.phases[level],
                    transfer.residual_norms[level],
                    transfer_valid && level <= transfer.converged,
                    transfer_valid,
                ) for (level, value) in enumerate(transfer.eigenvalues)
            ),
        ),
    ]

    mkpath(directory)
    for (filename, content) in contents
        _atomic_replace(joinpath(directory, filename)) do temporary
            open(temporary, "w") do io
                write(io, content)
            end
        end
    end
    return String(directory)
end
