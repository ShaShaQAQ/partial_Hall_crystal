const SINGLE_POINT_KEYS = Set([
    "geometry",
    "Ly",
    "Ny",
    "x_period",
    "filling_num",
    "filling_den",
    "phi_y",
    "maxdim",
    "t1",
    "t3",
    "V1",
    "V2",
    "V3",
    "cutoff",
    "vumps_tol",
    "energy_tol",
    "mismatch_tol",
    "imaginary_tol",
    "transfer_tol",
    "transfer_neigs",
    "maxiter",
    "stability",
    "threads",
    "output",
    "checkpoint",
    "load",
    "occupied_sites",
    "seed",
    "allow_nonconverged",
])

const SINGLE_POINT_REQUIRED_KEYS = (
    "x_period",
    "filling_num",
    "filling_den",
    "phi_y",
    "maxdim",
    "output",
)

function _configuration_from_options(options)
    geometry = Symbol(get(options, "geometry", "legacy_sheared"))
    geometry in SUPPORTED_GEOMETRIES || throw(
        ArgumentError("--geometry must be legacy_sheared or paper_straight")
    )
    common = (
        x_period=_parse_int_option(options["x_period"], "x_period"),
        filling_num=_parse_int_option(options["filling_num"], "filling_num"),
        filling_den=_parse_int_option(options["filling_den"], "filling_den"),
        phi_y=_parse_float_option(options["phi_y"], "phi_y"),
    )
    if geometry == :legacy_sheared
        haskey(options, "Ny") && throw(
            ArgumentError("--Ny requires --geometry=paper_straight")
        )
        Ly = _parse_int_option(_required(options, "Ly"), "Ly")
        return InfiniteCylinderConfig(; geometry, Ly, common...)
    end
    haskey(options, "Ly") && throw(
        ArgumentError("--Ly cannot be combined with --geometry=paper_straight")
    )
    Ny = _parse_int_option(_required(options, "Ny"), "Ny")
    return InfiniteCylinderConfig(; geometry, Ny, common...)
end

struct SinglePointSettings
    config::InfiniteCylinderConfig
    model::CylinderModelParams
    maxdim_schedule::Vector{Int}
    cutoff::Float64
    vumps_tol::Float64
    energy_tol::Float64
    energy_mismatch_tol::Float64
    imaginary_tol::Float64
    transfer_tol::Float64
    transfer_neigs::Int
    max_iterations::Int
    stable_iterations::Int
    threads::Int
    output::String
    checkpoint::String
    load::Union{Nothing,String}
    occupied_sites::Vector{Int}
    seed::Int
    allow_nonconverged::Bool
end

function _parse_key_value_args(args, allowed_keys)
    parsed = Dict{String,String}()
    for argument in args
        startswith(argument, "--") && occursin('=', argument) || throw(
            ArgumentError("arguments must use --key=value syntax: $(repr(argument))")
        )
        key, value = split(argument[3:end], '='; limit=2)
        isempty(key) && throw(ArgumentError("option key must not be empty"))
        key in allowed_keys || throw(ArgumentError("unknown option --$key"))
        haskey(parsed, key) && throw(ArgumentError("duplicate option --$key"))
        parsed[key] = value
    end
    return parsed
end

function _required(options, key::String)
    haskey(options, key) || throw(ArgumentError("missing required option --$key"))
    return options[key]
end

function _parse_int_option(value::AbstractString, key::String)
    return try
        parse(Int, value)
    catch
        throw(ArgumentError("--$key must be an integer"))
    end
end

function _parse_float_option(value::AbstractString, key::String)
    parsed = try
        parse(Float64, value)
    catch
        throw(ArgumentError("--$key must be a number"))
    end
    isfinite(parsed) || throw(ArgumentError("--$key must be finite"))
    return parsed
end

function _parse_bool_option(value::AbstractString, key::String)
    value == "true" && return true
    value == "false" && return false
    throw(ArgumentError("--$key must be true or false"))
end

function _parse_maxdim_schedule(value::AbstractString)
    isempty(value) && throw(ArgumentError("--maxdim must not be empty"))
    fields = split(value, ','; keepempty=true)
    any(isempty, fields) && throw(
        ArgumentError("--maxdim must be a comma-separated list of positive integers")
    )
    targets = [_parse_int_option(field, "maxdim") for field in fields]
    all(>(0), targets) || throw(
        ArgumentError("--maxdim targets must be positive integers")
    )
    all(targets[index] < targets[index + 1] for index in 1:(length(targets) - 1)) ||
        throw(ArgumentError("--maxdim must be strictly increasing"))
    return targets
end

function _positive_float(options, key::String, default::Float64)
    value = haskey(options, key) ? _parse_float_option(options[key], key) : default
    value > 0 || throw(ArgumentError("--$key must be positive"))
    return value
end

function _nonnegative_float(options, key::String, default::Float64)
    value = haskey(options, key) ? _parse_float_option(options[key], key) : default
    value >= 0 || throw(ArgumentError("--$key must be nonnegative"))
    return value
end

function _positive_int(options, key::String, default::Int)
    value = haskey(options, key) ? _parse_int_option(options[key], key) : default
    value > 0 || throw(ArgumentError("--$key must be positive"))
    return value
end

function _nonempty_path(value::AbstractString, key::String)
    isempty(strip(value)) && throw(ArgumentError("--$key must not be empty"))
    return String(value)
end

const FIG2_BENCHMARK_KEYS = Set([
    "manifest",
    "stage",
    "output",
    "dimensions",
    "flux_units_2pi",
    "threads",
])

const FIG2_BENCHMARK_REQUIRED_KEYS = (
    "manifest",
    "stage",
    "output",
    "dimensions",
    "flux_units_2pi",
)

struct Fig2BenchmarkSettings
    manifest::String
    stage::String
    output::String
    dimensions::Vector{Int}
    fluxes::Vector{Float64}
    threads::Int
end

function _parse_fig2_dimensions(value::AbstractString)
    isempty(value) && throw(ArgumentError("--dimensions must not be empty"))
    fields = split(value, ','; keepempty=true)
    any(isempty, fields) && throw(
        ArgumentError("--dimensions must be a comma-separated list of positive integers")
    )
    dimensions = [_parse_int_option(field, "dimensions") for field in fields]
    all(>(0), dimensions) || throw(
        ArgumentError("--dimensions entries must be positive")
    )
    all(
        dimensions[index] < dimensions[index + 1] for
        index in 1:(length(dimensions) - 1)
    ) || throw(ArgumentError("--dimensions must be strictly increasing"))
    return dimensions
end

function _parse_fig2_fluxes(value::AbstractString)
    isempty(value) && throw(ArgumentError("--flux_units_2pi must not be empty"))
    fields = split(value, ','; keepempty=true)
    any(isempty, fields) && throw(
        ArgumentError("--flux_units_2pi must be a comma-separated list of numbers")
    )
    units = [_parse_float_option(field, "flux_units_2pi") for field in fields]
    length(units) >= 2 || throw(
        ArgumentError("--flux_units_2pi must contain both 0 and 3")
    )
    all(
        units[index] < units[index + 1] for index in 1:(length(units) - 1)
    ) || throw(ArgumentError("--flux_units_2pi must be strictly increasing"))
    first(units) == 0.0 || throw(
        ArgumentError("--flux_units_2pi must start at 0")
    )
    last(units) == 3.0 || throw(
        ArgumentError("--flux_units_2pi must stop at 3")
    )
    return 2pi .* units
end

function parse_fig2_benchmark_args(args=ARGS)
    options = _parse_key_value_args(args, FIG2_BENCHMARK_KEYS)
    foreach(key -> _required(options, key), FIG2_BENCHMARK_REQUIRED_KEYS)
    return Fig2BenchmarkSettings(
        _nonempty_path(options["manifest"], "manifest"),
        _nonempty_path(options["stage"], "stage"),
        _nonempty_path(options["output"], "output"),
        _parse_fig2_dimensions(options["dimensions"]),
        _parse_fig2_fluxes(options["flux_units_2pi"]),
        _positive_int(options, "threads", 1),
    )
end

function _parse_occupied_pattern(
    value::AbstractString,
    c::InfiniteCylinderConfig,
    key::String,
)
    isempty(value) && throw(ArgumentError("--$key must not be empty"))
    fields = split(value, ','; keepempty=true)
    any(isempty, fields) && throw(
        ArgumentError("--$key must be a comma-separated list of site integers")
    )
    pattern = [_parse_int_option(field, key) for field in fields]
    length(pattern) == particles_per_cell(c) || throw(
        ArgumentError("--$key particle count does not match the cylinder configuration")
    )
    all(site -> 1 <= site <= sites_per_cell(c), pattern) || throw(
        ArgumentError("--$key site lies outside the reference cell")
    )
    length(unique(pattern)) == length(pattern) || throw(
        ArgumentError("--$key contains a duplicate occupied site")
    )
    return pattern
end

function parse_single_point_args(args=ARGS)
    options = _parse_key_value_args(args, SINGLE_POINT_KEYS)
    foreach(key -> _required(options, key), SINGLE_POINT_REQUIRED_KEYS)

    config = _configuration_from_options(options)
    model_values = [
        haskey(options, key) ? _parse_float_option(options[key], key) : default for
        (key, default) in zip(
            ("t1", "t3", "V1", "V2", "V3"),
            (1.0, 0.2, 1.0, 0.0, 0.0),
        )
    ]
    model = CylinderModelParams(model_values...)
    maxdim_schedule = _parse_maxdim_schedule(options["maxdim"])
    cutoff = _positive_float(options, "cutoff", 1e-8)
    vumps_tol = _positive_float(options, "vumps_tol", 1e-5)
    derived_energy_tol = 10vumps_tol
    isfinite(derived_energy_tol) || throw(
        ArgumentError("--energy_tol default derived from --vumps_tol must be finite")
    )
    energy_tol = _positive_float(options, "energy_tol", derived_energy_tol)
    mismatch_tol = _positive_float(options, "mismatch_tol", derived_energy_tol)
    imaginary_tol = _nonnegative_float(options, "imaginary_tol", 1e-12)
    transfer_tol = _positive_float(options, "transfer_tol", 1e-8)
    transfer_neigs = _positive_int(options, "transfer_neigs", 4)
    transfer_neigs >= 2 || throw(ArgumentError("--transfer_neigs must be at least two"))
    max_iterations = _positive_int(options, "maxiter", 50)
    stable_iterations = _positive_int(options, "stability", 2)
    stable_iterations <= max_iterations || throw(
        ArgumentError("--stability must not exceed --maxiter")
    )
    threads = _positive_int(options, "threads", 1)
    output = _nonempty_path(options["output"], "output")
    checkpoint = haskey(options, "checkpoint") ?
        _nonempty_path(options["checkpoint"], "checkpoint") :
        joinpath(output, "state.h5")
    load = haskey(options, "load") ? _nonempty_path(options["load"], "load") : nothing
    occupied_sites = haskey(options, "occupied_sites") ?
        _parse_occupied_pattern(options["occupied_sites"], config, "occupied_sites") :
        default_occupied_sites(config)
    seed = haskey(options, "seed") ? _parse_int_option(options["seed"], "seed") : 0
    seed >= 0 || throw(ArgumentError("--seed must be nonnegative"))
    allow_nonconverged = haskey(options, "allow_nonconverged") ?
        _parse_bool_option(options["allow_nonconverged"], "allow_nonconverged") : false

    return SinglePointSettings(
        config,
        model,
        maxdim_schedule,
        cutoff,
        vumps_tol,
        energy_tol,
        mismatch_tol,
        imaginary_tol,
        transfer_tol,
        transfer_neigs,
        max_iterations,
        stable_iterations,
        threads,
        output,
        checkpoint,
        load,
        occupied_sites,
        seed,
        allow_nonconverged,
    )
end

const FLUX_SCAN_ONLY_KEYS = Set([
    "phi_start",
    "phi_stop",
    "phi_steps",
    "branch_mode",
    "cold_patterns",
    "fidelity_drop_tol",
    "charge_jump_tol",
    "sector_tol",
    "spectrum_tol",
    "spectrum_levels",
])
const FLUX_SCAN_KEYS = union(setdiff(SINGLE_POINT_KEYS, Set(["phi_y", "checkpoint"])), FLUX_SCAN_ONLY_KEYS)

struct FluxScanSettings
    point::SinglePointSettings
    phi_stop::Float64
    phi_steps::Int
    branch_mode::Symbol
    cold_patterns::Vector{Vector{Int}}
    fidelity_drop_tol::Float64
    charge_jump_tol::Float64
    sector_tol::Float64
    spectrum_tol::Float64
    spectrum_levels::Int
end

const BRANCH_CUT_X = 1

function _scan_metadata(settings::FluxScanSettings)
    c = settings.point.config
    dependencies = _toml_value(_runtime_dependencies())
    dependencies["itensor_infinite_mps_commit"] = ITENSOR_INFINITE_MPS_COMMIT
    return Dict{String,Any}(
        "format" => "infinite_cylinder_vumps_scan_v1",
        "algorithm" => "VUMPS",
        "configuration" => Dict(
            "Ly" => c.Ly,
            "x_period" => c.x_period,
            "filling_num" => c.filling_num,
            "filling_den" => c.filling_den,
            "sites_per_cell" => sites_per_cell(c),
            "unit_cells_per_cell" => unit_cells_per_cell(c),
            "particles_per_cell" => particles_per_cell(c),
            "physical_site_density" => Float64(physical_site_density(c)),
            "charge_scale" => charge_scale(c),
            "signature" => _configuration_signature(c),
        ),
        "scan" => Dict(
            "phi_start" => c.phi_y,
            "phi_stop" => settings.phi_stop,
            "phi_steps" => settings.phi_steps,
        ),
        "branch" => Dict(
            "mode" => string(settings.branch_mode),
            "cut_x" => BRANCH_CUT_X,
            "fidelity_drop_tol" => settings.fidelity_drop_tol,
            "charge_jump_tol" => settings.charge_jump_tol,
            "sector_tol" => settings.sector_tol,
            "spectrum_tol" => settings.spectrum_tol,
            "spectrum_levels" => settings.spectrum_levels,
        ),
        "initialization" => Dict(
            "occupied_sites" => copy(settings.point.occupied_sites),
            "cold_patterns" => copy.(settings.cold_patterns),
            "seed" => settings.point.seed,
        ),
        "optimization" => _toml_value(_optimization_metadata(settings.point)),
        "dependencies" => dependencies,
    )
end

function _parse_cold_patterns(value::AbstractString, c::InfiniteCylinderConfig)
    isempty(value) && return Vector{Int}[]
    patterns = Vector{Int}[]
    for raw_pattern in split(value, ';'; keepempty=true)
        isempty(raw_pattern) && throw(ArgumentError("--cold_patterns contains an empty pattern"))
        raw_sites = split(raw_pattern, ','; keepempty=true)
        any(isempty, raw_sites) && throw(
            ArgumentError("--cold_patterns must contain comma-separated site integers")
        )
        pattern = [_parse_int_option(site, "cold_patterns") for site in raw_sites]
        length(pattern) == particles_per_cell(c) || throw(
            ArgumentError("cold pattern particle count does not match the cylinder configuration")
        )
        all(site -> 1 <= site <= sites_per_cell(c), pattern) || throw(
            ArgumentError("cold pattern site lies outside the reference cell")
        )
        length(unique(pattern)) == length(pattern) || throw(
            ArgumentError("cold pattern contains a duplicate occupied site")
        )
        push!(patterns, pattern)
    end
    return patterns
end

function parse_flux_scan_args(args=ARGS)
    options = _parse_key_value_args(args, FLUX_SCAN_KEYS)
    for key in ("phi_start", "phi_stop", "phi_steps")
        _required(options, key)
    end

    point_arguments = String[]
    for (key, value) in options
        key in FLUX_SCAN_ONLY_KEYS && continue
        push!(point_arguments, "--$key=$value")
    end
    push!(point_arguments, "--phi_y=$(options["phi_start"])")
    point = parse_single_point_args(point_arguments)
    phi_stop = _parse_float_option(options["phi_stop"], "phi_stop")
    phi_steps = _positive_int(options, "phi_steps", 0)
    point.config.phi_y != phi_stop || phi_steps == 1 || throw(
        ArgumentError("--phi_steps must be one when start and stop are equal")
    )
    point.config.phi_y == phi_stop || phi_steps >= 2 || throw(
        ArgumentError("--phi_steps must be at least two for unequal start and stop")
    )

    branch_mode = Symbol(get(options, "branch_mode", "ground"))
    branch_mode in (:ground, :adiabatic) || throw(
        ArgumentError("--branch_mode must be ground or adiabatic")
    )
    cold_patterns = _parse_cold_patterns(get(options, "cold_patterns", ""), point.config)
    fidelity_drop_tol = _nonnegative_float(options, "fidelity_drop_tol", 1e-3)
    charge_jump_tol = _nonnegative_float(options, "charge_jump_tol", 0.1)
    sector_tol = _nonnegative_float(options, "sector_tol", 0.1)
    spectrum_tol = _nonnegative_float(options, "spectrum_tol", 0.5)
    spectrum_levels = _positive_int(options, "spectrum_levels", 20)
    return FluxScanSettings(
        point,
        phi_stop,
        phi_steps,
        branch_mode,
        cold_patterns,
        fidelity_drop_tol,
        charge_jump_tol,
        sector_tol,
        spectrum_tol,
        spectrum_levels,
    )
end

function flux_grid(phi_start::Real, phi_stop::Real, points::Integer)
    start = Float64(phi_start)
    stop = Float64(phi_stop)
    all(isfinite, (start, stop)) || throw(ArgumentError("flux endpoints must be finite"))
    !(points isa Bool) && points > 0 || throw(ArgumentError("flux point count must be positive"))
    count = Int(points)
    start == stop || count >= 2 || throw(
        ArgumentError("unequal flux endpoints require at least two points")
    )
    start != stop || count == 1 || throw(
        ArgumentError("equal flux endpoints require exactly one point")
    )
    count == 1 && return [start]
    return collect(range(start, stop; length=count))
end

flux_grid(settings::FluxScanSettings) =
    flux_grid(settings.point.config.phi_y, settings.phi_stop, settings.phi_steps)

struct FluxCandidateData
    label::String
    cut_x::Int
    energy_raw::Float64
    raw_schmidt_polarization::Float64
    sector_weights::Dict{Int,Float64}
    spectrum_levels::Vector{NamedTuple}
    converged::Bool
    valid::Bool
    energy_valid::Bool
    polarization_valid::Bool
    sector_valid::Bool
    spectrum_valid::Bool
    fidelity_to_previous::Float64
    fidelity_valid::Bool
    reason::String
end

function FluxCandidateData(
    label::AbstractString,
    energy_raw::Real,
    raw_schmidt_polarization::Real,
    sector_weights::AbstractDict,
    spectrum_levels;
    cut_x::Integer=1,
    converged::Bool=true,
    valid::Bool=true,
    energy_valid::Bool=isfinite(energy_raw),
    polarization_valid::Bool=isfinite(raw_schmidt_polarization),
    sector_valid::Bool=!isempty(sector_weights),
    spectrum_valid::Bool=!isempty(spectrum_levels),
    fidelity_to_previous::Real=NaN,
    fidelity_valid::Bool=false,
    reason::AbstractString=valid ? "valid candidate" : "invalid candidate",
)
    isempty(strip(label)) && throw(ArgumentError("flux candidate label must not be empty"))
    !(cut_x isa Bool) && cut_x > 0 || throw(
        ArgumentError("flux candidate cut_x must be positive")
    )
    weights = Dict{Int,Float64}()
    for (charge, weight) in sector_weights
        charge isa Integer && !(charge isa Bool) || throw(
            ArgumentError("sector charges must be integers")
        )
        weights[Int(charge)] = Float64(weight)
    end
    levels = NamedTuple[level for level in spectrum_levels]
    return FluxCandidateData(
        String(label),
        Int(cut_x),
        Float64(energy_raw),
        Float64(raw_schmidt_polarization),
        weights,
        levels,
        converged,
        valid,
        energy_valid,
        polarization_valid,
        sector_valid,
        spectrum_valid,
        Float64(fidelity_to_previous),
        fidelity_valid,
        String(reason),
    )
end

function select_flux_candidate_data(candidates; mode)
    isempty(candidates) && throw(ArgumentError("candidate list must not be empty"))
    mode in (:ground, :adiabatic) || throw(ArgumentError("unsupported branch mode: $mode"))
    selected = 0
    best = mode === :ground ? Inf : -Inf
    for (index, candidate) in pairs(candidates)
        value = if mode === :ground
            candidate.converged && candidate.energy_valid &&
                isfinite(candidate.energy_raw) || continue
            candidate.energy_raw
        else
            candidate.converged || continue
            candidate.fidelity_valid && isfinite(candidate.fidelity_to_previous) || continue
            0 <= candidate.fidelity_to_previous <= 1 || continue
            candidate.fidelity_to_previous
        end
        if (mode === :ground && value < best) || (mode === :adiabatic && value > best)
            selected = index
            best = value
        end
    end
    iszero(selected) && throw(ArgumentError("no selectable candidate is available"))
    return selected
end

function _fallback_flux_candidate(candidates; mode)
    isempty(candidates) && throw(ArgumentError("candidate list must not be empty"))
    mode in (:ground, :adiabatic) || throw(ArgumentError("unsupported branch mode: $mode"))
    valid_energy_indices = findall(
        candidate -> candidate.energy_valid && isfinite(candidate.energy_raw),
        candidates,
    )
    if !isempty(valid_energy_indices)
        selected = first(valid_energy_indices)
        selected_energy = candidates[selected].energy_raw
        for index in Iterators.drop(valid_energy_indices, 1)
            energy = candidates[index].energy_raw
            if energy < selected_energy
                selected = index
                selected_energy = energy
            end
        end
        return selected
    end
    return firstindex(candidates)
end

function flux_scan_row(
    point::Integer,
    phi_y::Real,
    candidate::FluxCandidateData,
    selected_candidate::Integer;
    reference_polarization::Real,
    x_period::Integer,
    selection_mode::Symbol,
    ground_energy::Real,
    ground_energy_valid::Bool,
    selection_valid::Bool=true,
    selection_reason::AbstractString="selected by requested branch mode",
)
    !(x_period isa Bool) && x_period > 0 || throw(
        ArgumentError("x_period must be positive")
    )
    selection_mode in (:ground, :adiabatic) || throw(
        ArgumentError("selection_mode must be ground or adiabatic")
    )
    polarization_valid = candidate.converged && candidate.polarization_valid &&
        isfinite(candidate.raw_schmidt_polarization) && isfinite(reference_polarization)
    delta = polarization_valid ?
        candidate.raw_schmidt_polarization - Float64(reference_polarization) : NaN
    energy_delta_valid = candidate.converged && candidate.energy_valid &&
        isfinite(candidate.energy_raw) && ground_energy_valid && isfinite(ground_energy)
    energy_delta = energy_delta_valid ?
        candidate.energy_raw - Float64(ground_energy) : NaN
    fidelity_valid = candidate.converged && candidate.fidelity_valid &&
        isfinite(candidate.fidelity_to_previous) &&
        0 <= candidate.fidelity_to_previous <= 1
    fidelity_density_valid = fidelity_valid
    fidelity_density_divergent = fidelity_valid && iszero(candidate.fidelity_to_previous)
    fidelity_density = fidelity_valid ?
        (fidelity_density_divergent ? Inf :
         -log(candidate.fidelity_to_previous) / Int(x_period)) : NaN
    return (
        point=Int(point),
        phi_y=Float64(phi_y),
        selected_candidate=Int(selected_candidate),
        selected_label=candidate.label,
        selection_mode=string(selection_mode),
        cut_x=candidate.cut_x,
        selection_valid,
        selection_reason=String(selection_reason),
        converged=candidate.converged,
        valid=candidate.valid,
        energy_raw=candidate.energy_raw,
        energy_valid=candidate.energy_valid && isfinite(candidate.energy_raw),
        energy_delta_to_ground_raw=energy_delta,
        energy_delta_to_ground_valid=energy_delta_valid,
        raw_schmidt_polarization=candidate.raw_schmidt_polarization,
        polarization_valid,
        delta_raw_schmidt_polarization=delta,
        fidelity_raw=candidate.fidelity_to_previous,
        fidelity_valid,
        fidelity_density_x_raw=fidelity_density,
        fidelity_density_x_valid=fidelity_density_valid,
        fidelity_density_x_divergent=fidelity_density_divergent,
        reason=candidate.reason,
    )
end

function _same_flux_sector_gauge_rows(point::Integer, phi_y::Real, candidates)
    cold = filter(candidate -> startswith(candidate.label, "cold_"), candidates)
    rows = NamedTuple[]
    for first_index in 1:(length(cold) - 1), second_index in (first_index + 1):length(cold)
        first_candidate = cold[first_index]
        second_candidate = cold[second_index]
        alignment = try
            first_candidate.cut_x == second_candidate.cut_x || throw(
                ArgumentError("candidate branch cuts do not match")
            )
            first_candidate.sector_valid && second_candidate.sector_valid ||
                throw(ArgumentError("candidate sector table is invalid"))
            best_sector_shift(
                first_candidate.sector_weights,
                second_candidate.sector_weights,
            )
        catch error
            error isa InterruptException && rethrow()
            push!(
                rows,
                (
                    point=Int(point),
                    phi_y=Float64(phi_y),
                    cut_x=first_candidate.cut_x,
                    first_label=first_candidate.label,
                    second_label=second_candidate.label,
                    sector_shift_raw=0,
                    sector_distance_raw=NaN,
                    valid=false,
                    reason="same-flux sector alignment failed: $(sprint(showerror, error))",
                ),
            )
            continue
        end
        push!(
            rows,
            (
                point=Int(point),
                phi_y=Float64(phi_y),
                cut_x=first_candidate.cut_x,
                first_label=first_candidate.label,
                second_label=second_candidate.label,
                sector_shift_raw=alignment.shift,
                sector_distance_raw=alignment.distance,
                valid=true,
                reason="diagnostic shift is added to second candidate raw-QN keys only for comparison",
            ),
        )
    end
    return rows
end

function _energy_crossing(previous, current)
    previous_by_label = Dict(candidate.label => candidate for candidate in previous)
    current_by_label = Dict(candidate.label => candidate for candidate in current)
    labels = sort!(collect(intersect(keys(previous_by_label), keys(current_by_label))))
    length(labels) >= 2 || return (raw=false, valid=false)
    compared = false
    for left in 1:(length(labels) - 1), right in (left + 1):length(labels)
        previous_left = previous_by_label[labels[left]]
        previous_right = previous_by_label[labels[right]]
        current_left = current_by_label[labels[left]]
        current_right = current_by_label[labels[right]]
        all(candidate -> candidate.energy_valid && isfinite(candidate.energy_raw),
            (previous_left, previous_right, current_left, current_right)) || continue
        compared = true
        old_order = sign(previous_left.energy_raw - previous_right.energy_raw)
        new_order = sign(current_left.energy_raw - current_right.energy_raw)
        old_order != 0 && new_order != 0 && old_order != new_order &&
            return (raw=true, valid=true)
    end
    return (raw=false, valid=compared)
end

function flux_branch_event(
    previous_row,
    current_row,
    previous_candidates,
    current_candidates;
    fidelity_drop_tol,
    charge_jump_tol,
    sector_tol,
    spectrum_tol,
    spectrum_levels,
)
    previous = previous_candidates[previous_row.selected_candidate]
    current = current_candidates[current_row.selected_candidate]
    fidelity_valid = current.fidelity_valid && isfinite(current.fidelity_to_previous) &&
        0 <= current.fidelity_to_previous <= 1
    charge_valid = previous_row.polarization_valid && current_row.polarization_valid
    charge_step = charge_valid ?
        current.raw_schmidt_polarization - previous.raw_schmidt_polarization : NaN
    sector_valid = previous.sector_valid && current.sector_valid
    sector_distance = sector_valid ?
        bhattacharyya_distance(previous.sector_weights, current.sector_weights) : NaN
    sector_shift = sector_valid ?
        best_sector_shift(previous.sector_weights, current.sector_weights).shift : 0
    spectrum = previous.spectrum_valid && current.spectrum_valid ?
        charge_matched_spectrum_distance(
            previous.spectrum_levels,
            current.spectrum_levels;
            nlevels=spectrum_levels,
        ) : (distance=NaN, valid=false)
    crossing = _energy_crossing(previous_candidates, current_candidates)
    flags_valid = fidelity_valid && charge_valid && sector_valid && spectrum.valid && crossing.valid
    return (
        from_point=previous_row.point,
        to_point=current_row.point,
        from_phi_y=previous_row.phi_y,
        to_phi_y=current_row.phi_y,
        fidelity_raw=current.fidelity_to_previous,
        fidelity_valid,
        charge_step_raw=charge_step,
        charge_step_valid=charge_valid,
        sector_distance_raw=sector_distance,
        sector_distance_valid=sector_valid,
        sector_shift_raw=sector_shift,
        sector_shift_valid=sector_valid,
        spectrum_distance_raw=spectrum.distance,
        spectrum_distance_valid=spectrum.valid,
        energy_crossing_raw=crossing.raw,
        energy_crossing_valid=crossing.valid,
        fidelity_drop_tol=Float64(fidelity_drop_tol),
        charge_jump_tol=Float64(charge_jump_tol),
        sector_tol=Float64(sector_tol),
        spectrum_tol=Float64(spectrum_tol),
        spectrum_levels=Int(spectrum_levels),
        flag_fidelity=fidelity_valid && 1 - current.fidelity_to_previous > fidelity_drop_tol,
        flag_charge=charge_valid && abs(charge_step) > charge_jump_tol,
        flag_sector=sector_valid && sector_distance > sector_tol,
        flag_spectrum=spectrum.valid && spectrum.distance > spectrum_tol,
        flag_energy_crossing=crossing.valid && crossing.raw,
        flags_valid,
    )
end

struct WorkflowValidationError <: Exception
    message::String
end

Base.showerror(io::IO, error::WorkflowValidationError) = print(io, error.message)

Base.@kwdef struct SinglePointOperations
    configure_threads::Function = configure_cli_threads
    initialize::Function = (config, occupied_sites) -> begin
        sites, _, psi = initial_infinite_mps(config; occupied_sites)
        (; sites, psi)
    end
    load_state::Function = load_checkpoint
    state_sites::Function = psi -> siteinds(only, psi.AL)
    build_hamiltonian::Function = build_infinite_mpo
    optimize::Function = _default_optimize
    energy::Function = (psi, H, config, settings) ->
        energy_data(psi, H, config; imaginary_tol=settings.imaginary_tol)
    density::Function = density_data
    entanglement::Function = (psi, config, cut_x) ->
        entanglement_data(psi, config; cut_x)
    transfer::Function = (psi, config, settings, rng) -> neutral_transfer_data(
        psi,
        config;
        neigs=settings.transfer_neigs,
        tol=settings.transfer_tol,
        rng,
    )
    fidelity::Function = (previous, current, config, settings, rng) ->
        mixed_transfer_fidelity(
            previous,
            current,
            config;
            tol=settings.transfer_tol,
            rng,
        )
    write_outputs::Function = _default_write_outputs
    save_state::Function = save_checkpoint
    cold_state::Function = _default_cold_state
end

struct SinglePointRunResult
    settings::SinglePointSettings
    optimization::VUMPSResult
    energy::EnergyData
    densities::Vector{DensityRow}
    entanglements::Vector{EntanglementData}
    transfer::NeutralTransferData
    valid::Bool
    observable_reasons::Vector{String}
end

function configure_cli_threads(threads::Integer)
    !(threads isa Bool) && threads > 0 || throw(ArgumentError("threads must be positive"))
    requested = Int(threads)
    if requested == 1
        LinearAlgebra.BLAS.set_num_threads(1)
        ITensors.Strided.disable_threads()
        ITensors.disable_threaded_blocksparse()
    else
        requested <= Base.Threads.nthreads() || throw(
            ArgumentError(
                "requested $requested ITensor threads but Julia started with only $(Base.Threads.nthreads()) threads",
            )
        )
        LinearAlgebra.BLAS.set_num_threads(1)
        ITensors.Strided.disable_threads()
        ITensors.enable_threaded_blocksparse()
    end
    return (
        requested,
        julia=Base.Threads.nthreads(),
        blas=LinearAlgebra.BLAS.get_num_threads(),
        strided=ITensors.Strided.get_num_threads(),
        blocksparse=ITensors.using_threaded_blocksparse(),
    )
end

function _default_optimize(H, psi, settings::SinglePointSettings)
    return run_vumps(
        H,
        psi;
        maxdim_schedule=settings.maxdim_schedule,
        cutoff=settings.cutoff,
        max_iterations=settings.max_iterations,
        vumps_tol=settings.vumps_tol,
        energy_tol=settings.energy_tol,
        energy_mismatch_tol=settings.energy_mismatch_tol,
        stable_iterations=settings.stable_iterations,
        imaginary_tol=settings.imaginary_tol,
        canonical_seed=_single_point_canonical_seed(settings.seed),
    )
end

function _optimization_metadata(settings::SinglePointSettings)
    return (
        model=(
            t1=settings.model.t1,
            t3=settings.model.t3,
            V1=settings.model.V1,
            V2=settings.model.V2,
            V3=settings.model.V3,
        ),
        maxdim_schedule=settings.maxdim_schedule,
        cutoff=settings.cutoff,
        vumps_tol=settings.vumps_tol,
        energy_tol=settings.energy_tol,
        energy_mismatch_tol=settings.energy_mismatch_tol,
        imaginary_tol=settings.imaginary_tol,
        transfer_tol=settings.transfer_tol,
        transfer_neigs=settings.transfer_neigs,
        max_iterations=settings.max_iterations,
        stable_iterations=settings.stable_iterations,
        threads=settings.threads,
        occupied_sites=settings.occupied_sites,
        seed=settings.seed,
        allow_nonconverged=settings.allow_nonconverged,
    )
end

function _runtime_dependencies()
    return (
        julia_version=string(VERSION),
        itensor_infinite_mps_version=string(Base.pkgversion(ITensorInfiniteMPS)),
        itensor_mps_version=string(Base.pkgversion(ITensorMPS)),
        itensors_version=string(Base.pkgversion(ITensors)),
        krylovkit_version=string(Base.pkgversion(KrylovKit)),
        hdf5_version=string(Base.pkgversion(HDF5)),
    )
end

function _seed_word(value::Integer)
    value >= 0 || throw(ArgumentError("seed stream integers must be nonnegative"))
    return UInt64(value)
end

function _seed_word(value::Union{Symbol,AbstractString})
    word = UInt64(0xcbf29ce484222325)
    for byte in codeunits(string(value))
        word = xor(word, UInt64(byte)) * UInt64(0x100000001b3)
    end
    return word
end

function _derived_seed(seed::Integer, tags...)
    word = xor(_seed_word(seed), UInt64(0x9e3779b97f4a7c15))
    for tag in tags
        word = xor(
            word,
            _seed_word(tag) + UInt64(0x9e3779b97f4a7c15) + (word << 6) + (word >> 2),
        )
    end
    return word
end

_single_point_canonical_seed(seed::Integer) =
    _derived_seed(seed, :vumps_canonicalization)

_flux_candidate_seed(seed::Integer, point::Integer, label::AbstractString)::Int =
    Int(mod(
        _derived_seed(seed, :flux_candidate, point, label),
        UInt64(typemax(Int)),
    ))

function _default_write_outputs(
    directory,
    config,
    result,
    energy,
    densities,
    entanglements,
    transfer,
    settings,
)
    return write_output_files(
        directory,
        config,
        result,
        energy,
        densities,
        entanglements,
        transfer;
        optimization=_optimization_metadata(settings),
        dependencies=_runtime_dependencies(),
    )
end

function _default_cold_state(sites, config::InfiniteCylinderConfig, pattern)
    occupied = Set(pattern)
    initstate(n::Integer) = mod1(n, sites_per_cell(config)) in occupied ? "Occ" : "Emp"
    return InfMPS(sites, initstate)
end

function _with_phi(
    settings::SinglePointSettings,
    phi_y::Real;
    output=settings.output,
    checkpoint=settings.checkpoint,
    seed=settings.seed,
)
    config = with_flux(settings.config, phi_y)
    return SinglePointSettings(
        config,
        settings.model,
        copy(settings.maxdim_schedule),
        settings.cutoff,
        settings.vumps_tol,
        settings.energy_tol,
        settings.energy_mismatch_tol,
        settings.imaginary_tol,
        settings.transfer_tol,
        settings.transfer_neigs,
        settings.max_iterations,
        settings.stable_iterations,
        settings.threads,
        String(output),
        String(checkpoint),
        nothing,
        copy(settings.occupied_sites),
        seed,
        settings.allow_nonconverged,
    )
end

function _observable_or_fallback(compute::Function, fallback::Function, name::String, reasons)
    try
        return compute()
    catch error
        error isa InterruptException && rethrow()
        push!(reasons, "$name failed: $(sprint(showerror, error))")
        return fallback()
    end
end

function _single_point_observables(
    result::VUMPSResult,
    H,
    settings::SinglePointSettings,
    operations::SinglePointOperations,
    transfer_rng::AbstractRNG,
)
    config = settings.config
    reasons = String[]
    energy = _observable_or_fallback(
        () -> operations.energy(result.psi, H, config, settings),
        () -> EnergyData(NaN, NaN, NaN, NaN),
        "energy",
        reasons,
    )
    densities = _observable_or_fallback(
        () -> Vector{DensityRow}(operations.density(result.psi, config)),
        () -> [
            DensityRow(site, site_coordinates(config, site)..., NaN) for
            site in 1:sites_per_cell(config)
        ],
        "density",
        reasons,
    )
    entanglements = EntanglementData[]
    for cut_x in 1:config.x_period
        data = _observable_or_fallback(
            () -> operations.entanglement(result.psi, config, cut_x),
            () -> EntanglementData(cut_x, cut_x * config.Ly, NaN, NaN, EntanglementLevel[], EntanglementSector[]),
            "entanglement cut $cut_x",
            reasons,
        )
        push!(entanglements, data)
    end
    transfer = _observable_or_fallback(
        () -> operations.transfer(
            result.psi,
            config,
            settings,
            transfer_rng,
        ),
        () -> _invalid_neutral_transfer(
            ComplexF64[],
            Float64[],
            0,
            "neutral transfer observable failed",
        ),
        "neutral transfer",
        reasons,
    )
    return (; energy, densities, entanglements, transfer, reasons)
end

function _required_observables_valid(result::VUMPSResult, observables)
    energy_valid = all(isfinite, (
        observables.energy.per_cell,
        observables.energy.per_x,
        observables.energy.per_unit_cell,
        observables.energy.per_site,
    ))
    density_valid = length(observables.densities) == nsites(result.psi) &&
        all(row -> isfinite(row.density), observables.densities)
    entanglement_valid = !isempty(observables.entanglements) &&
        all(data -> isfinite(data.entropy) &&
                    isfinite(data.raw_schmidt_polarization) &&
                    !isempty(data.levels) && !isempty(data.sectors), observables.entanglements)
    transfer = observables.transfer
    transfer_valid = transfer.valid && transfer.converged >= 2 &&
        all(isfinite, (transfer.ratio, transfer.xi_cell, transfer.xi_x))
    return result.converged && energy_valid && density_valid && entanglement_valid && transfer_valid
end

function _run_prepared_point(
    settings::SinglePointSettings,
    H,
    psi,
    operations::SinglePointOperations,
    ;
    transfer_rng::AbstractRNG=Random.Xoshiro(
        _derived_seed(settings.seed, :transfer),
    ),
)
    optimization = operations.optimize(H, psi, settings)
    observables = _single_point_observables(
        optimization,
        H,
        settings,
        operations,
        transfer_rng,
    )
    valid = _required_observables_valid(optimization, observables)
    operations.write_outputs(
        settings.output,
        settings.config,
        optimization,
        observables.energy,
        observables.densities,
        observables.entanglements,
        observables.transfer,
        settings,
    )
    operations.save_state(settings.checkpoint, optimization.psi, settings.config)
    return SinglePointRunResult(
        settings,
        optimization,
        observables.energy,
        observables.densities,
        observables.entanglements,
        observables.transfer,
        valid,
        observables.reasons,
    )
end

function _enforce_workflow_validity(result::SinglePointRunResult)
    result.valid || result.settings.allow_nonconverged || throw(
        WorkflowValidationError(
            result.optimization.converged ?
            "required single-point observable is invalid; raw output and checkpoint were preserved" :
            "VUMPS is nonconverged; raw output and checkpoint were preserved: $(result.optimization.reason)",
        )
    )
    return result
end

function run_single_point(
    settings::SinglePointSettings;
    operations::SinglePointOperations=SinglePointOperations(),
)
    operations.configure_threads(settings.threads)
    if isnothing(settings.load)
        state = operations.initialize(settings.config, settings.occupied_sites)
        sites = state.sites
        psi = state.psi
    else
        psi = operations.load_state(settings.load, settings.config)
        sites = operations.state_sites(psi)
    end
    H = operations.build_hamiltonian(settings.config, settings.model, sites)
    return _enforce_workflow_validity(
        _run_prepared_point(settings, H, psi, operations)
    )
end

const SCAN_SUMMARY_HEADER =
    "point\tphi_y\tselected_candidate\tselected_label\tselection_mode\tcut_x\tselection_valid\tselection_reason\tconverged\tvalid\tenergy_raw\tenergy_valid\tenergy_delta_to_ground_raw\tenergy_delta_to_ground_valid\traw_schmidt_polarization\tpolarization_valid\tdelta_raw_schmidt_polarization\tfidelity_raw\tfidelity_valid\tfidelity_density_x_raw\tfidelity_density_x_valid\tfidelity_density_x_divergent\treason"
const BRANCH_EVENTS_HEADER =
    "from_point\tto_point\tfrom_phi_y\tto_phi_y\tfidelity_raw\tfidelity_valid\tcharge_step_raw\tcharge_step_valid\tsector_distance_raw\tsector_distance_valid\tsector_shift_raw\tsector_shift_valid\tspectrum_distance_raw\tspectrum_distance_valid\tenergy_crossing_raw\tenergy_crossing_valid\tfidelity_drop_tol\tcharge_jump_tol\tsector_tol\tspectrum_tol\tspectrum_levels\tflag_fidelity\tflag_charge\tflag_sector\tflag_spectrum\tflag_energy_crossing\tflags_valid"
const SECTOR_GAUGE_HEADER =
    "point\tphi_y\tcut_x\tfirst_label\tsecond_label\tsector_shift_raw\tsector_distance_raw\tvalid\treason"

struct FluxScanRunResult
    settings::FluxScanSettings
    rows::Vector{NamedTuple}
    events::Vector{NamedTuple}
    sector_gauge::Vector{NamedTuple}
    candidates::Vector{Vector{FluxCandidateData}}
    selected_results::Vector{SinglePointRunResult}
end

function _candidate_data(result::SinglePointRunResult, label; fidelity=nothing)
    entanglement = isempty(result.entanglements) ? nothing : first(result.entanglements)
    sector_weights = isnothing(entanglement) ? Dict{Int,Float64}() :
        Dict(sector.raw_charge => sector.weight for sector in entanglement.sectors)
    spectrum = isnothing(entanglement) ? NamedTuple[] : [
        (
            raw_charge=level.raw_charge,
            entanglement_energy=level.entanglement_energy,
            probability=level.probability,
        ) for level in entanglement.levels
    ]
    fidelity_raw = isnothing(fidelity) ? NaN : fidelity.fidelity_cell
    fidelity_valid = !isnothing(fidelity) && fidelity.valid && isfinite(fidelity_raw)
    reasons = String[result.optimization.reason]
    append!(reasons, result.observable_reasons)
    if !isnothing(fidelity) && !fidelity.valid
        push!(reasons, "mixed fidelity invalid: $(fidelity.reason)")
    end
    return FluxCandidateData(
        label,
        result.energy.per_site,
        isnothing(entanglement) ? NaN : entanglement.raw_schmidt_polarization,
        sector_weights,
        spectrum;
        cut_x=isnothing(entanglement) ? BRANCH_CUT_X : entanglement.cut_x,
        converged=result.optimization.converged,
        valid=result.valid,
        energy_valid=all(isfinite, (
            result.energy.per_cell,
            result.energy.per_x,
            result.energy.per_unit_cell,
            result.energy.per_site,
        )),
        polarization_valid=!isnothing(entanglement) &&
            isfinite(entanglement.raw_schmidt_polarization),
        sector_valid=!isempty(sector_weights),
        spectrum_valid=!isempty(spectrum) && all(level -> isfinite(level.entanglement_energy), spectrum),
        fidelity_to_previous=fidelity_raw,
        fidelity_valid,
        reason=join(reasons, "; "),
    )
end

function _safe_fidelity(
    operations::SinglePointOperations,
    previous,
    current,
    config,
    settings,
    rng,
)
    try
        return operations.fidelity(previous, current, config, settings, rng)
    catch error
        error isa InterruptException && rethrow()
        return _invalid_mixed_transfer(
            ComplexF64(NaN),
            ComplexF64(NaN),
            ComplexF64(NaN),
            NaN,
            NaN,
            NaN,
            "mixed fidelity failed: $(sprint(showerror, error))",
        )
    end
end

function _write_scan_table(path, header, rows)
    _atomic_replace(path) do temporary
        open(temporary, "w") do io
            println(io, header)
            for row in rows
                println(io, join(_output_value.(values(row)), '\t'))
            end
        end
    end
    return path
end

function _write_scan_metadata(path, settings::FluxScanSettings)
    _atomic_replace(path) do temporary
        open(temporary, "w") do io
            write(io, _render_summary(_scan_metadata(settings)))
        end
    end
    return path
end

function _ground_energy(candidates)
    valid = filter(
        candidate -> candidate.converged && candidate.energy_valid &&
            isfinite(candidate.energy_raw),
        candidates,
    )
    isempty(valid) && return (value=NaN, valid=false)
    return (value=minimum(candidate.energy_raw for candidate in valid), valid=true)
end

function run_flux_scan(
    settings::FluxScanSettings;
    operations::SinglePointOperations=SinglePointOperations(),
)
    operations.configure_threads(settings.point.threads)
    if isnothing(settings.point.load)
        initial = operations.initialize(
            settings.point.config,
            settings.point.occupied_sites,
        )
        shared_sites = initial.sites
        warm_state = initial.psi
    else
        warm_state = operations.load_state(
            settings.point.load,
            settings.point.config,
        )
        shared_sites = operations.state_sites(warm_state)
    end
    previous_selected_state = nothing
    reference_polarization = NaN
    rows = NamedTuple[]
    events = NamedTuple[]
    sector_gauge = NamedTuple[]
    all_candidates = Vector{FluxCandidateData}[]
    selected_results = SinglePointRunResult[]

    for (point, phi_y) in enumerate(flux_grid(settings))
        config = with_flux(settings.point.config, phi_y)
        H = operations.build_hamiltonian(config, settings.point.model, shared_sites)
        point_directory = joinpath(settings.point.output, "phi_$(lpad(point - 1, 3, '0'))")
        labels = ["warm"; ["cold_$(lpad(index, 2, '0'))" for index in eachindex(settings.cold_patterns)]]
        states = Any[warm_state]
        append!(states, [
            operations.cold_state(shared_sites, config, pattern) for
            pattern in settings.cold_patterns
        ])
        point_results = SinglePointRunResult[]
        point_candidates = FluxCandidateData[]
        for (label, state) in zip(labels, states)
            candidate_directory = joinpath(point_directory, "candidate_$label")
            point_settings = _with_phi(
                settings.point,
                phi_y;
                output=candidate_directory,
                checkpoint=joinpath(candidate_directory, "state.h5"),
                seed=_flux_candidate_seed(settings.point.seed, point, label),
            )
            candidate_result = _run_prepared_point(
                point_settings,
                H,
                state,
                operations;
                transfer_rng=Random.Xoshiro(
                    _derived_seed(settings.point.seed, :transfer, point, label),
                ),
            )
            push!(point_results, candidate_result)
            fidelity = isnothing(previous_selected_state) ? nothing :
                _safe_fidelity(
                    operations,
                    previous_selected_state,
                    candidate_result.optimization.psi,
                    config,
                    point_settings,
                    Random.Xoshiro(
                        _derived_seed(settings.point.seed, :fidelity, point, label),
                    ),
                )
            push!(point_candidates, _candidate_data(candidate_result, label; fidelity))
        end

        mode = isnothing(previous_selected_state) ? :ground : settings.branch_mode
        selected_index = 0
        selection_valid = true
        selection_reason = "selected by requested branch mode"
        try
            selected_index = select_flux_candidate_data(point_candidates; mode)
        catch error
            error isa ArgumentError || rethrow()
            selected_index = _fallback_flux_candidate(point_candidates; mode)
            selection_valid = false
            selection_reason = "raw fallback after selection failure: $(sprint(showerror, error))"
        end
        selected_result = point_results[selected_index]
        selected_candidate = point_candidates[selected_index]
        ground_energy = _ground_energy(point_candidates)
        if point == 1
            reference_polarization = selected_candidate.raw_schmidt_polarization
        end
        row = flux_scan_row(
            point,
            phi_y,
            selected_candidate,
            selected_index;
            reference_polarization,
            x_period=config.x_period,
            selection_mode=mode,
            ground_energy=ground_energy.value,
            ground_energy_valid=ground_energy.valid,
            selection_valid,
            selection_reason,
        )
        if !isempty(rows)
            push!(
                events,
                flux_branch_event(
                    rows[end],
                    row,
                    all_candidates[end],
                    point_candidates;
                    fidelity_drop_tol=settings.fidelity_drop_tol,
                    charge_jump_tol=settings.charge_jump_tol,
                    sector_tol=settings.sector_tol,
                    spectrum_tol=settings.spectrum_tol,
                    spectrum_levels=settings.spectrum_levels,
                ),
            )
        end
        push!(rows, row)
        append!(
            sector_gauge,
            _same_flux_sector_gauge_rows(point, phi_y, point_candidates),
        )
        push!(all_candidates, point_candidates)
        push!(selected_results, selected_result)
        warm_state = selected_result.optimization.psi
        previous_selected_state = warm_state
    end

    mkpath(settings.point.output)
    _write_scan_table(joinpath(settings.point.output, "scan_summary.tsv"), SCAN_SUMMARY_HEADER, rows)
    _write_scan_table(joinpath(settings.point.output, "branch_events.tsv"), BRANCH_EVENTS_HEADER, events)
    _write_scan_table(
        joinpath(settings.point.output, "sector_gauge.tsv"),
        SECTOR_GAUGE_HEADER,
        sector_gauge,
    )
    _write_scan_metadata(
        joinpath(settings.point.output, "scan_metadata.toml"),
        settings,
    )
    scan = FluxScanRunResult(
        settings,
        rows,
        events,
        sector_gauge,
        all_candidates,
        selected_results,
    )
    all(result -> result.valid, selected_results) &&
        all(row -> row.selection_valid, rows) ||
        settings.point.allow_nonconverged || throw(
        WorkflowValidationError(
            "flux scan contains a nonconverged or invalid selected point; raw candidate outputs and scan tables were preserved",
        )
    )
    return scan
end
