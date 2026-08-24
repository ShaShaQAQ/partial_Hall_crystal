const FIG2_MANIFEST_FORMAT = "fqahc_fig2_benchmark_v4"
const FIG2_LEDGER_FORMAT = "fqahc_fig2_ledger_v3"
const FIG2_CANDIDATE_FORMAT = "fqahc_fig2_candidate_v4"
const FIG2_PROGRESS_EVENT_FORMAT = "fqahc_fig2_progress_event_v1"
const FIG2_PROGRESS_POINTER_FORMAT = "fqahc_fig2_progress_pointer_v1"
const FIG2_PROGRESS_SUMMARY_FORMAT = "fqahc_fig2_progress_summary_v1"
const FIG2_PROGRESS_CONVERGENCE_STATE_POLICY = "reset_on_resume"
const FIG2_MIXED_REFERENCE_ABSENT = "absent"
const FIG2_PAPER_CURVE_FLUX_UNITS_2PI = collect(0.0:0.125:3.0)
const FIG2_PAPER_CURVE_DELTA_Q = Float64[
    0.000000, 0.039286, 0.080329, 0.122619, 0.165646,
    0.208617, 0.250794, 0.292290, 0.332370, 0.371882,
    0.410998, 0.450680, 0.490760, 0.532426, 0.575907,
    0.620918, 0.667290, 0.713152, 0.760601, 0.807029,
    0.851927, 0.889739, 0.925283, 0.961735, 1.000000,
]
const FIG2_PAPER_CURVE_UNCERTAINTY = fill(0.02, 25)
const FIG2_REQUIRED_CANDIDATE_FILES = (
    "state.h5",
    "summary.toml",
    "convergence.tsv",
    "expansion.tsv",
    "density.tsv",
    "entanglement_spectrum.tsv",
    "schmidt_sectors.tsv",
    "progress.toml",
    "momentum_entanglement_spectrum.tsv",
    "mixed_fidelity.tsv",
    "candidate.toml",
)
const FIG2_ARTIFACT_HEADERS = Dict(
    "convergence.tsv" =>
        "stage\titeration\tmaxlinkdim\tenergy_left\tenergy_right\tenergy_mismatch\tdelta_energy\teps_left\teps_right\tprecision_error\telapsed_seconds\tconverged",
    "expansion.tsv" =>
        "stage\ttarget\tbefore\tafter\tprogressed\telapsed_seconds",
    "density.tsv" => "site\tx\ty\tdensity\tvalid",
    "entanglement_spectrum.tsv" =>
        "cut_x\tbond\tlevel\tsingular_value\tprobability\tentanglement_energy\tqn\traw_charge\tphysical_charge\tvalid",
    "schmidt_sectors.tsv" =>
        "cut_x\tbond\tqn\traw_charge\tphysical_charge\tweight\tvalid",
    "momentum_entanglement_spectrum.tsv" =>
        "level\tprobability\tentanglement_energy\traw_charge\tky_index\tky\teigenphase\tphase_residual\tunitary_residual\tpower_residual\tcommutator_residual\ttransfer_residual\tvalidation_valid\tvalidation_reason\tmomentum_valid\tmomentum_reason\tmomentum_convention\tsector_ky_offsets",
    "mixed_fidelity.tsv" =>
        "fidelity_to_previous\tmixed_value_real\tmixed_value_imag\tself_value1_real\tself_value1_imag\tself_value2_real\tself_value2_imag\tmixed_residual\tself_residual1\tself_residual2\tconverged\tvalid\treason\treference_present\treference_dimension\treference_point\treference_candidate_id\treference_directory\treference_state_sha256",
)
const FIG2_LEGALLY_EMPTY_TABLES = Set(("expansion.tsv",))

struct Fig2BenchmarkSpec
    source_path::String
    source_bytes::Vector{UInt8}
    sha256::String
    data::Dict{String,Any}
    config::InfiniteCylinderConfig
    model::CylinderModelParams
    acceptance_rows::Vector{String}
end

Base.@kwdef struct Fig2RunnerOverride
    cutoff::Float64
    vumps_tol::Float64
    energy_tol::Float64
    energy_mismatch_tol::Float64
    transfer_tol::Float64
    momentum_residual_tol::Float64
    max_iterations::Int
    stable_iterations::Int
end

struct Fig2MomentumCounting
    counts::Vector{Int}
    delta_ky::Vector{Int}
    cutoffs::Vector{Float64}
    reference_raw_charge::Union{Nothing,Int}
    reference_level::Union{Nothing,Int}
    reference_ky_index::Union{Nothing,Int}
    minimum_cutoff_margin::Float64
    valid::Bool
    reason::String
end

struct Fig2MixedReference
    present::Bool
    dimension::Int
    point::Int
    candidate_id::String
    directory::String
    state_sha256::String
end

Base.@kwdef struct Fig2CandidateEvidence
    state::Any
    requested_maxdim::Int
    achieved_maxlinkdim::Int
    checkpoint_maxlinkdim::Int=achieved_maxlinkdim
    converged::Bool
    valid::Bool
    energy_per_site::Float64
    raw_schmidt_polarization::Float64
    sector_weights::Dict{Int,Float64}
    fidelity_to_previous::Float64=NaN
    fidelity_valid::Bool=false
    mixed_fidelity::Union{Nothing,MixedTransferFidelity}=nothing
    momentum::MomentumEntanglementData
    momentum_counting::Vector{Int}=Int[]
    momentum_counting_evidence::Union{Nothing,Fig2MomentumCounting}=nothing
    restart_valid::Bool=false
    reason::String=""
end

function _default_fig2_candidate_ids end
function _default_fig2_run_candidate end
function _default_fig2_load_state end
function _default_fig2_provenance end
function _default_fig2_persisted_checkpoint_audit end
function _default_fig2_progress_audit end

Base.@kwdef struct Fig2BenchmarkOperations
    candidate_ids::Function=_default_fig2_candidate_ids
    run_candidate::Function=_default_fig2_run_candidate
    load_state::Function=_default_fig2_load_state
    provenance::Function=_default_fig2_provenance
    checkpoint_audit::Function=_default_fig2_persisted_checkpoint_audit
    progress_audit::Function=_default_fig2_progress_audit
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
    momentum_validation::MomentumValidation
    momentum_valid::Bool
    momentum_reason::String
    momentum_convention::Symbol
    momentum_sector_ky_offsets::Dict{Int,Int}
    momentum_counting::Vector{Int}
    momentum_counting_evidence::Fig2MomentumCounting
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

function _fig2_file_sha256(path::AbstractString)
    return open(path, "r") do io
        bytes2hex(SHA.sha256(io))
    end
end

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

function _fig2_number_vector(data, key::String)
    value = _fig2_required(data, key)
    value isa AbstractVector || throw(
        ArgumentError("Fig. 2 manifest $key must be an array")
    )
    result = Float64[]
    for entry in value
        entry isa Real && !(entry isa Bool) || throw(
            ArgumentError("Fig. 2 manifest $key must contain only real numbers")
        )
        converted = Float64(entry)
        isfinite(converted) || throw(
            ArgumentError("Fig. 2 manifest $key must contain only finite numbers")
        )
        push!(result, converted)
    end
    return result
end

function _fig2_integer_vector(data, key::String)
    value = _fig2_required(data, key)
    value isa AbstractVector || throw(
        ArgumentError("Fig. 2 manifest $key must be an array")
    )
    result = Int[]
    for entry in value
        entry isa Integer && !(entry isa Bool) || throw(
            ArgumentError("Fig. 2 manifest $key must contain only integers")
        )
        push!(result, Int(entry))
    end
    return result
end

function _fig2_positive_number(data, key::String, description)
    value = _fig2_number(data, key)
    value > 0 || throw(ArgumentError("$description must be positive"))
    return value
end

function _fig2_positive_integer(data, key::String, description)
    value = _fig2_integer(data, key)
    value > 0 || throw(ArgumentError("$description must be positive"))
    return value
end

function _fig2_strictly_increasing(values, description)
    length(values) >= 2 || throw(
        ArgumentError("$description must contain at least two values")
    )
    all(index -> values[index] < values[index + 1], 1:(length(values) - 1)) ||
        throw(ArgumentError("$description must be strictly increasing"))
    return values
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
    _fig2_number(data, "expected_endpoint_pump") == 1.0 || throw(
        ArgumentError("Fig. 2 expected endpoint pump is not fixed")
    )
    _fig2_integer(data, "expected_sector_shift") == 1 || throw(
        ArgumentError("Fig. 2 expected sector shift is not fixed")
    )
    _fig2_positive_number(
        data, "pump_tolerance", "Fig. 2 pump tolerance"
    ) == 0.05 || throw(ArgumentError("Fig. 2 pump tolerance is not fixed"))
    _fig2_required(data, "momentum_convention") ==
        "sector_relative_canonical_cyclic_orbit" || throw(
        ArgumentError("Fig. 2 momentum convention is not fixed")
    )
    _fig2_required(data, "raw_pump_policy") ==
        "no_offset_no_unwrap_no_sign_flip_no_branch_translation" || throw(
        ArgumentError("Fig. 2 manifest must preserve raw pump data")
    )
    counting = _fig2_integer_vector(data, "counting")
    counting == [1, 1, 2, 3, 5] || throw(
        ArgumentError("Fig. 2 benchmark counting target is not fixed")
    )
    counting_rule = _fig2_required(data, "momentum_counting_rule")
    _fig2_required(counting_rule, "format") ==
        "fqahc_fig2_momentum_counting_v1" || throw(
        ArgumentError("unsupported Fig. 2 momentum-counting rule format")
    )
    _fig2_required(counting_rule, "reference_sector") ==
        "largest_eigenvalue_raw_charge_sector" || throw(
        ArgumentError("Fig. 2 momentum counting must use the largest-eigenvalue charge sector")
    )
    _fig2_integer_vector(counting_rule, "delta_ky") == [0, -1, -2, -3, -4] || throw(
        ArgumentError("Fig. 2 momentum-counting sectors are not the fixed contract")
    )
    _fig2_number(counting_rule, "cut_intercept") == 2.15 || throw(
        ArgumentError("Fig. 2 momentum-counting intercept is not fixed")
    )
    _fig2_number(counting_rule, "cut_slope") == -1.58 || throw(
        ArgumentError("Fig. 2 momentum-counting slope is not fixed")
    )
    _fig2_number(counting_rule, "digitization_uncertainty") == 0.20 || throw(
        ArgumentError("Fig. 2 momentum-counting uncertainty is not fixed")
    )
    _fig2_required(counting_rule, "source_pdf_sha256") ==
        "23f0691570b906ae7a9394cef2e3dd6eebd4e36d79e882a4bbd2f17824cc89a0" || throw(
        ArgumentError("Fig. 2 momentum-counting source checksum is not fixed")
    )
    _fig2_integer(counting_rule, "source_page") == 3 || throw(
        ArgumentError("Fig. 2 momentum-counting source page is not fixed")
    )
    _fig2_required(counting_rule, "source_panel") == "f" || throw(
        ArgumentError("Fig. 2 momentum-counting source panel is not fixed")
    )
    paper_curve_rule = _fig2_required(data, "paper_curve_rule")
    _fig2_required(paper_curve_rule, "format") ==
        "fqahc_fig2_paper_curve_comparison_v3" || throw(
        ArgumentError("unsupported Fig. 2 paper-curve rule format")
    )
    _fig2_required(paper_curve_rule, "source_pdf_sha256") ==
        "23f0691570b906ae7a9394cef2e3dd6eebd4e36d79e882a4bbd2f17824cc89a0" || throw(
        ArgumentError("Fig. 2 paper-curve source checksum is not fixed")
    )
    _fig2_integer(paper_curve_rule, "source_page") == 3 || throw(
        ArgumentError("Fig. 2 paper-curve source page is not fixed")
    )
    _fig2_required(paper_curve_rule, "source_panel") == "d" || throw(
        ArgumentError("Fig. 2 paper-curve source panel is not fixed")
    )
    _fig2_required(paper_curve_rule, "source_vector_path_sha256") ==
        "54e5c57f88a735e5efe0c2a24de455bbaa987bbceb5a9affa44f85fab58a892c" || throw(
        ArgumentError("Fig. 2 paper-curve SVG vector-element checksum is not fixed")
    )
    _fig2_required(paper_curve_rule, "target") ==
        "nu_star=1/3,V1=10,V2=V3=2" || throw(
        ArgumentError("Fig. 2 paper-curve target is not fixed")
    )
    _fig2_required(paper_curve_rule, "extraction") ==
        "pdftocairo_svg_vector_path" || throw(
        ArgumentError("Fig. 2 paper-curve extraction method is not fixed")
    )
    _fig2_number_vector(paper_curve_rule, "flux_units_2pi") ==
        FIG2_PAPER_CURVE_FLUX_UNITS_2PI || throw(
        ArgumentError("Fig. 2 paper-curve flux grid is not fixed")
    )
    _fig2_number_vector(paper_curve_rule, "delta_q") ==
        FIG2_PAPER_CURVE_DELTA_Q || throw(
        ArgumentError("Fig. 2 paper-curve digitization is not fixed")
    )
    _fig2_number_vector(paper_curve_rule, "uncertainty") ==
        FIG2_PAPER_CURVE_UNCERTAINTY || throw(
        ArgumentError("Fig. 2 paper-curve uncertainty is not fixed")
    )

    optimization = _fig2_required(data, "optimization")
    _fig2_required(optimization, "multisite_update_alg") == "sequential" || throw(
        ArgumentError("Fig. 2 multisite update algorithm must be sequential")
    )
    optimization_targets = (
        "cutoff" => 1.0e-9,
        "vumps_tol" => 1.0e-6,
        "energy_tol" => 1.0e-6,
        "energy_mismatch_tol" => 1.0e-6,
        "transfer_tol" => 1.0e-8,
        "momentum_residual_tol" => 1.0e-6,
    )
    for (key, target) in optimization_targets
        _fig2_positive_number(
            optimization, key, "Fig. 2 optimization $key"
        ) == target || throw(
            ArgumentError("Fig. 2 optimization $key is not fixed")
        )
    end
    max_iterations = _fig2_positive_integer(
        optimization, "max_iterations", "Fig. 2 optimization max_iterations"
    )
    max_iterations == 50 || throw(
        ArgumentError("Fig. 2 optimization max_iterations is not fixed")
    )
    stable_iterations = _fig2_positive_integer(
        optimization, "stable_iterations", "Fig. 2 optimization stable_iterations"
    )
    stable_iterations == 2 || throw(
        ArgumentError("Fig. 2 optimization stable_iterations is not fixed")
    )
    progress_generations_to_keep = _fig2_positive_integer(
        optimization,
        "progress_generations_to_keep",
        "Fig. 2 optimization progress_generations_to_keep",
    )
    progress_generations_to_keep == 2 || throw(
        ArgumentError(
            "Fig. 2 optimization progress_generations_to_keep is not fixed"
        )
    )
    stable_iterations <= max_iterations || throw(
        ArgumentError("Fig. 2 stable_iterations cannot exceed max_iterations")
    )

    pilot = _fig2_required(data, "pilot")
    pilot_maxdims = _fig2_integer_vector(pilot, "maxdims")
    all(value -> value > 0, pilot_maxdims) || throw(
        ArgumentError("Fig. 2 pilot maxdims must be positive")
    )
    _fig2_strictly_increasing(pilot_maxdims, "Fig. 2 pilot maxdims")
    pilot_maxdims == [32, 64, 128] || throw(
        ArgumentError("Fig. 2 pilot maxdims are not fixed")
    )
    pilot_fluxes = _fig2_number_vector(pilot, "flux_units_2pi")
    _fig2_strictly_increasing(pilot_fluxes, "Fig. 2 pilot flux grid")
    pilot_fluxes == collect(0.0:0.5:3.0) || throw(
        ArgumentError("Fig. 2 pilot flux grid is not fixed")
    )
    _fig2_positive_integer(
        pilot, "cold_candidates", "Fig. 2 pilot cold_candidates"
    ) == 2 || throw(ArgumentError("Fig. 2 pilot cold_candidates is not fixed"))

    convergence = _fig2_required(data, "convergence")
    convergence_maxdims = _fig2_integer_vector(convergence, "maxdims")
    all(value -> value > 0, convergence_maxdims) || throw(
        ArgumentError("Fig. 2 convergence maxdims must be positive")
    )
    _fig2_strictly_increasing(
        convergence_maxdims, "Fig. 2 convergence maxdims"
    )
    convergence_maxdims == [256, 512, 1000, 2000, 3000] || throw(
        ArgumentError("Fig. 2 convergence maxdims progression is not fixed")
    )
    convergence_targets = (
        "energy_per_site_change_tol" => 1.0e-6,
        "endpoint_pump_change_tol" => 0.02,
        "sector_residual_tol" => 1.0e-6,
        "momentum_validation_residual_tol" => 1.0e-6,
    )
    for (key, target) in convergence_targets
        _fig2_positive_number(
            convergence, key, "Fig. 2 convergence $key"
        ) == target || throw(
            ArgumentError("Fig. 2 convergence $key is not fixed")
        )
    end
    early_plateau = _fig2_required(data, "early_plateau")
    _fig2_required(early_plateau, "policy") ==
        "allowed_if_all_gates_pass" || throw(
        ArgumentError("Fig. 2 early_plateau policy is not fixed")
    )
    early_plateau_dimensions = (
        "lower_dimension" => 1000,
        "upper_dimension" => 2000,
        "omitted_dimension" => 3000,
    )
    for (key, target) in early_plateau_dimensions
        _fig2_integer(early_plateau, key) == target || throw(
            ArgumentError("Fig. 2 early_plateau $key is not fixed")
        )
    end
    _fig2_required(early_plateau, "evidence_filename") ==
        "early_plateau.toml" || throw(
        ArgumentError("Fig. 2 early_plateau evidence filename is not fixed")
    )
    restart_gate = _fig2_required(data, "restart_gate")
    _fig2_required(restart_gate, "format") ==
        "fqahc_fig2_restart_gate_v1" || throw(
        ArgumentError("Fig. 2 restart-gate format is not fixed")
    )
    _fig2_required(restart_gate, "evidence_filename") ==
        "restart_gate.toml" || throw(
        ArgumentError("Fig. 2 restart-gate evidence filename is not fixed")
    )
    _fig2_required(restart_gate, "backend_commit") ==
        ITENSOR_INFINITE_MPS_COMMIT || throw(
        ArgumentError("Fig. 2 restart-gate backend commit is not fixed")
    )
    for (key, target) in (
        "energy_tolerance" => 1.0e-10,
        "density_tolerance" => 1.0e-10,
    )
        _fig2_positive_number(
            restart_gate, key, "Fig. 2 restart-gate $key"
        ) == target || throw(ArgumentError(
            "Fig. 2 restart-gate $key is not fixed"
        ))
    end
    required_restart_checks = _fig2_required(restart_gate, "required_checks")
    required_restart_checks isa AbstractVector && all(
        check -> check isa AbstractString, required_restart_checks
    ) || throw(ArgumentError(
        "Fig. 2 restart-gate required checks must be strings"
    ))
    String.(required_restart_checks) == [
        "load_completed",
        "no_expansion_iteration_completed",
        "controlled_expansion_completed",
        "energy_reproduced",
        "density_reproduced",
        "schmidt_sectors_reproduced",
        "entanglement_spectrum_reproduced",
        "site_indices_reproduced",
        "next_flux_iteration_completed",
    ] || throw(ArgumentError(
        "Fig. 2 restart-gate required checks are not the fixed contract"
    ))
    acceptance = _fig2_required(data, "acceptance")
    row_data = _fig2_required(acceptance, "rows")
    row_data isa AbstractVector && all(row -> row isa AbstractString, row_data) ||
        throw(ArgumentError("Fig. 2 acceptance rows must be strings"))
    rows = String.(row_data)
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

function _fig2_config_from_data(data)
    return InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=_fig2_integer(data, "Ny"),
        x_period=_fig2_integer(data, "x_period"),
        filling_num=_fig2_integer(data, "filling_num"),
        filling_den=_fig2_integer(data, "filling_den"),
        phi_y=_fig2_number(data, "flux_start"),
    )
end

function _fig2_model_from_data(data)
    return CylinderModelParams(;
        t1=_fig2_number(data, "t1"),
        t3=_fig2_number(data, "t3"),
        V1=_fig2_number(data, "V1"),
        V2=_fig2_number(data, "V2"),
        V3=_fig2_number(data, "V3"),
    )
end

function _fig2_argument_error(
    operation::Function,
    context;
    passthrough_argument_error=false,
)
    try
        return operation()
    catch error
        error isa InterruptException && rethrow()
        passthrough_argument_error && error isa ArgumentError && rethrow()
        throw(ArgumentError("$context: $(sprint(showerror, error))"))
    end
end

function _fig2_validated_snapshot(
    spec::Fig2BenchmarkSpec;
    parse_snapshot=TOML.parse,
    validate_snapshot=_validate_fig2_manifest,
)
    _fig2_sha256(spec.source_bytes) == spec.sha256 || throw(ArgumentError(
        "Fig. 2 manifest in-memory drift: source bytes no longer match the recorded checksum"
    ))
    snapshot = try
        Dict{String,Any}(parse_snapshot(String(copy(spec.source_bytes))))
    catch error
        error isa InterruptException && rethrow()
        throw(ArgumentError(
            "Fig. 2 manifest in-memory drift: source bytes are not valid TOML ($(sprint(showerror, error)))"
        ))
    end
    rows = try
        validate_snapshot(snapshot)
    catch error
        error isa InterruptException && rethrow()
        throw(ArgumentError(
            "Fig. 2 manifest in-memory drift: source snapshot is invalid ($(sprint(showerror, error)))"
        ))
    end
    isequal(spec.data, snapshot) || throw(ArgumentError(
        "Fig. 2 manifest in-memory drift: mutable data differs from the validated source snapshot"
    ))
    isequal(spec.config, _fig2_config_from_data(snapshot)) || throw(ArgumentError(
        "Fig. 2 manifest in-memory drift: cylinder configuration differs from the source snapshot"
    ))
    isequal(spec.model, _fig2_model_from_data(snapshot)) || throw(ArgumentError(
        "Fig. 2 manifest in-memory drift: model parameters differ from the source snapshot"
    ))
    spec.acceptance_rows == rows || throw(ArgumentError(
        "Fig. 2 manifest in-memory drift: acceptance rows differ from the source snapshot"
    ))
    return snapshot
end

function _fig2_runner_optimization(
    spec::Fig2BenchmarkSpec,
    override::Union{Nothing,Fig2RunnerOverride}=nothing,
)
    snapshot = _fig2_validated_snapshot(spec)
    if isnothing(override)
        optimization = snapshot["optimization"]
        return (;
            multisite_update_alg=Symbol(optimization["multisite_update_alg"]),
            cutoff=Float64(optimization["cutoff"]),
            vumps_tol=Float64(optimization["vumps_tol"]),
            energy_tol=Float64(optimization["energy_tol"]),
            energy_mismatch_tol=Float64(optimization["energy_mismatch_tol"]),
            transfer_tol=Float64(optimization["transfer_tol"]),
            momentum_residual_tol=Float64(
                optimization["momentum_residual_tol"]
            ),
            max_iterations=Int(optimization["max_iterations"]),
            stable_iterations=Int(optimization["stable_iterations"]),
        )
    end
    get(ENV, "IDMRG_FIG2_REAL_SMOKE", "0") == "1" || throw(ArgumentError(
        "Fig. 2 runner overrides are restricted to the opt-in real smoke"
    ))
    for (name, value) in (
        :cutoff => override.cutoff,
        :vumps_tol => override.vumps_tol,
        :energy_tol => override.energy_tol,
        :energy_mismatch_tol => override.energy_mismatch_tol,
        :transfer_tol => override.transfer_tol,
        :momentum_residual_tol => override.momentum_residual_tol,
    )
        isfinite(value) && value > 0 || throw(ArgumentError(
            "Fig. 2 smoke runner override $name must be finite and positive"
        ))
    end
    override.max_iterations > 0 && override.stable_iterations > 0 &&
        override.stable_iterations <= override.max_iterations || throw(
        ArgumentError("Fig. 2 smoke runner iteration overrides are invalid")
    )
    return (;
        multisite_update_alg=Symbol(snapshot["optimization"]["multisite_update_alg"]),
        cutoff=override.cutoff,
        vumps_tol=override.vumps_tol,
        energy_tol=override.energy_tol,
        energy_mismatch_tol=override.energy_mismatch_tol,
        transfer_tol=override.transfer_tol,
        momentum_residual_tol=override.momentum_residual_tol,
        max_iterations=override.max_iterations,
        stable_iterations=override.stable_iterations,
    )
end

function _invalid_fig2_momentum_counting(
    delta_ky,
    cutoffs,
    reason;
    reference_raw_charge=nothing,
    reference_level=nothing,
    reference_ky_index=nothing,
    minimum_cutoff_margin=NaN,
)
    return Fig2MomentumCounting(
        Int[],
        Int.(delta_ky),
        Float64.(cutoffs),
        reference_raw_charge,
        reference_level,
        reference_ky_index,
        Float64(minimum_cutoff_margin),
        false,
        String(reason),
    )
end

function _fig2_momentum_counting(
    spec::Fig2BenchmarkSpec,
    momentum::MomentumEntanglementData,
)
    _fig2_validated_snapshot(spec)
    rule = spec.data["momentum_counting_rule"]
    delta_ky = Int.(rule["delta_ky"])
    intercept = Float64(rule["cut_intercept"])
    slope = Float64(rule["cut_slope"])
    uncertainty = Float64(rule["digitization_uncertainty"])
    cutoffs = intercept .+ slope .* delta_ky
    if !(momentum.valid && momentum.validation.valid)
        return _invalid_fig2_momentum_counting(
            delta_ky,
            cutoffs,
            "momentum labels are invalid: $(momentum.reason)",
        )
    end
    string(momentum.momentum_convention) == spec.data["momentum_convention"] ||
        return _invalid_fig2_momentum_counting(
            delta_ky,
            cutoffs,
            "momentum convention does not match the immutable Fig. 2 contract",
        )
    isempty(momentum.levels) && return _invalid_fig2_momentum_counting(
        delta_ky, cutoffs, "momentum spectrum is empty"
    )
    all(level ->
        isfinite(level.probability) && level.probability >= 0 &&
        (isfinite(level.entanglement_energy) || level.entanglement_energy == Inf) &&
        !ismissing(level.ky_index), momentum.levels) ||
        return _invalid_fig2_momentum_counting(
            delta_ky,
            cutoffs,
            "momentum spectrum contains nonfinite probabilities, invalid energies, or missing labels",
        )

    reference_index = first(sortperm(eachindex(momentum.levels); by=index -> (
        -momentum.levels[index].probability,
        momentum.levels[index].level,
    )))
    reference = momentum.levels[reference_index]
    reference.probability > 0 && isfinite(reference.entanglement_energy) ||
        return _invalid_fig2_momentum_counting(
            delta_ky, cutoffs, "largest Schmidt eigenvalue is not positive and finite"
        )
    reference_ky = Int(reference.ky_index)
    Ny = spec.config.Ny
    reference_charge = reference.raw_charge
    relative_levels = NamedTuple[]
    for level in momentum.levels
        level.raw_charge == reference_charge || continue
        difference = mod(Int(level.ky_index) - reference_ky, Ny)
        relative_ky = iszero(difference) ? 0 : difference - Ny
        relative_ky in delta_ky || continue
        push!(relative_levels, (;
            delta_ky=relative_ky,
            energy=level.entanglement_energy - reference.entanglement_energy,
        ))
    end
    margins = Float64[
        abs(level.energy - cutoffs[findfirst(==(level.delta_ky), delta_ky)]) for
        level in relative_levels
    ]
    minimum_margin = isempty(margins) ? Inf : minimum(margins)
    minimum_margin > uncertainty || return _invalid_fig2_momentum_counting(
        delta_ky,
        cutoffs,
        "a low-lying level lies within the Fig. 2 digitization uncertainty of the counting cut",
        reference_raw_charge=reference_charge,
        reference_level=reference.level,
        reference_ky_index=reference_ky,
        minimum_cutoff_margin=minimum_margin,
    )
    counts = Int[
        count(level ->
            level.delta_ky == sector &&
            level.energy < cutoffs[index], relative_levels) for
        (index, sector) in enumerate(delta_ky)
    ]
    return Fig2MomentumCounting(
        counts,
        delta_ky,
        cutoffs,
        reference_charge,
        reference.level,
        reference_ky,
        minimum_margin,
        true,
        "valid gauge-invariant within-charge-sector Fig. 2 counting",
    )
end

function load_fig2_benchmark(path::AbstractString)
    source_path = abspath(path)
    isfile(source_path) || throw(ArgumentError("Fig. 2 manifest does not exist"))
    source_bytes = read(source_path)
    parsed = TOML.parse(String(copy(source_bytes)))
    data = Dict{String,Any}(parsed)
    rows = _validate_fig2_manifest(data)
    config = _fig2_config_from_data(data)
    model = _fig2_model_from_data(data)
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

function _fig2_progress_identity(
    spec,
    dimension,
    point,
    phi_y,
    candidate_id,
)
    snapshot = _fig2_validated_snapshot(spec)
    dimension isa Integer && !(dimension isa Bool) && dimension > 0 || throw(
        ArgumentError("Fig. 2 progress dimension must be a positive integer")
    )
    point isa Integer && !(point isa Bool) && point > 0 || throw(
        ArgumentError("Fig. 2 progress point must be a positive integer")
    )
    phi_y isa Real && !(phi_y isa Bool) && isfinite(phi_y) || throw(
        ArgumentError("Fig. 2 progress flux must be finite")
    )
    candidate_id isa AbstractString &&
        occursin(r"^[A-Za-z0-9_.-]+$", candidate_id) || throw(
        ArgumentError("Fig. 2 progress candidate ID is invalid")
    )
    config = with_flux(spec.config, Float64(phi_y))
    return Dict{String,Any}(
        "manifest_sha256" => spec.sha256,
        "configuration_signature" => configuration_signature(config),
        "dimension" => Int(dimension),
        "point" => Int(point),
        "phi_y" => Float64(phi_y),
        "candidate_id" => String(candidate_id),
        "multisite_update_alg" =>
            String(snapshot["optimization"]["multisite_update_alg"]),
        "progress_generations_to_keep" => Int(
            snapshot["optimization"]["progress_generations_to_keep"]
        ),
        "convergence_state_policy" =>
            FIG2_PROGRESS_CONVERGENCE_STATE_POLICY,
        "git_commit" => _fig2_repository_commit(),
        "benchmark_source_sha256" => _fig2_file_sha256(@__FILE__),
    )
end

const FIG2_PROGRESS_IDENTITY_KEYS = (
    "manifest_sha256",
    "configuration_signature",
    "dimension",
    "point",
    "phi_y",
    "candidate_id",
    "multisite_update_alg",
    "progress_generations_to_keep",
    "convergence_state_policy",
    "git_commit",
    "benchmark_source_sha256",
)

function _fig2_validate_progress_identity(data, expected)
    data isa AbstractDict || throw(
        ArgumentError("Fig. 2 progress metadata must be a table")
    )
    for key in FIG2_PROGRESS_IDENTITY_KEYS
        get(data, key, nothing) == expected[key] || throw(ArgumentError(
            "Fig. 2 progress $key disagrees with the requested candidate identity"
        ))
    end
    return nothing
end

function _fig2_progress_schedule(schedule, dimension)
    schedule isa AbstractVector || throw(
        ArgumentError("Fig. 2 progress maxdim schedule must be an array")
    )
    all(value -> value isa Integer && !(value isa Bool), schedule) || throw(
        ArgumentError("Fig. 2 progress maxdim schedule must contain integers")
    )
    targets = Int.(schedule)
    isempty(targets) && throw(
        ArgumentError("Fig. 2 progress maxdim schedule must not be empty")
    )
    full = _fig2_maxdim_schedule(dimension)
    first_stage = findfirst(==(first(targets)), full)
    !isnothing(first_stage) && full[first_stage:end] == targets || throw(
        ArgumentError("Fig. 2 progress maxdim schedule is not a canonical suffix")
    )
    return targets
end

function _fig2_progress_path(root, value, description)
    value isa AbstractString || throw(
        ArgumentError("Fig. 2 progress $description path must be a string")
    )
    relative = normpath(String(value))
    !isabspath(relative) && relative != "." || throw(
        ArgumentError("Fig. 2 progress $description path must be relative")
    )
    parts = splitpath(relative)
    all(part -> !isempty(part) && part != "." && part != "..", parts) || throw(
        ArgumentError("Fig. 2 progress $description path is unsafe")
    )
    resolved = abspath(joinpath(root, relative))
    relative_check = relpath(resolved, abspath(root))
    relative_check == relative || throw(
        ArgumentError("Fig. 2 progress $description path escapes the candidate directory")
    )
    return resolved, relative
end

function _fig2_progress_positive_integer(data, key)
    value = get(data, key, nothing)
    value isa Integer && !(value isa Bool) && value > 0 || throw(
        ArgumentError("Fig. 2 progress $key must be a positive integer")
    )
    return Int(value)
end

function _fig2_progress_nonnegative_integer(data, key)
    value = get(data, key, nothing)
    value isa Integer && !(value isa Bool) && value >= 0 || throw(
        ArgumentError("Fig. 2 progress $key must be a nonnegative integer")
    )
    return Int(value)
end

function _fig2_progress_event_payload(event::VUMPSProgressEvent)
    if event.kind == :expansion
        expansion = something(event.expansion)
        return Dict{String,Any}(
            "expansion" => Dict(
                "stage" => expansion.stage,
                "target" => expansion.target,
                "before" => expansion.before,
                "after" => expansion.after,
                "progressed" => expansion.progressed,
                "elapsed_seconds" => expansion.elapsed_seconds,
            ),
        )
    end
    record = something(event.record)
    return Dict{String,Any}(
        "iteration_record" => Dict(
            "stage" => record.stage,
            "iteration" => record.iteration,
            "maxlinkdim" => record.maxlinkdim,
            "energy_left" => record.energy_left,
            "energy_right" => record.energy_right,
            "energy_mismatch" => record.energy_mismatch,
            "delta_energy" => isfinite(record.delta_energy) ?
                record.delta_energy : "missing",
            "eps_left" => record.eps_left,
            "eps_right" => record.eps_right,
            "precision_error" => record.precision_error,
            "elapsed_seconds" => record.elapsed_seconds,
            "converged" => record.converged,
        ),
    )
end

function _fig2_read_progress_pointer(root, identity)
    pointer_path = joinpath(root, ".progress", "latest.toml")
    isfile(pointer_path) || return nothing
    pointer = Dict{String,Any}(TOML.parsefile(pointer_path))
    get(pointer, "format", "") == FIG2_PROGRESS_POINTER_FORMAT || throw(
        ArgumentError("unsupported Fig. 2 progress pointer format")
    )
    _fig2_validate_progress_identity(pointer, identity)
    _fig2_progress_positive_integer(pointer, "event_sequence")
    _fig2_progress_positive_integer(pointer, "runner_event_sequence")
    _fig2_progress_nonnegative_integer(pointer, "resume_count")
    _fig2_progress_positive_integer(pointer, "maxlinkdim")
    _fig2_progress_positive_integer(pointer, "progress_generations_to_keep")
    return pointer
end

function _fig2_prune_progress_generations!(
    root,
    state_paths,
    generations_to_keep::Integer,
)
    generations_to_keep > 0 || throw(
        ArgumentError("Fig. 2 progress retention must be positive")
    )
    generation_root = abspath(joinpath(root, ".progress", "generations"))
    isdir(generation_root) || return nothing
    retained = Set{String}()
    for relative in Iterators.take(state_paths, generations_to_keep)
        absolute, _ = _fig2_progress_path(root, relative, "state")
        dirname(absolute) == generation_root || throw(
            ArgumentError(
                "Fig. 2 progress generation is outside its state directory"
            )
        )
        push!(retained, absolute)
    end
    for name in readdir(generation_root)
        occursin(r"^state_[0-9]{8}_[0-9]+_[0-9]+\.h5$", name) || continue
        path = abspath(joinpath(generation_root, name))
        dirname(path) == generation_root || error(
            "Fig. 2 progress generation path escaped its directory"
        )
        path in retained || rm(path)
    end
    return nothing
end

function _fig2_persist_progress_event!(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    maxdim_schedule,
    event::VUMPSProgressEvent;
    resume_count::Integer,
    canonicalize_state=nothing,
    save_state=save_checkpoint,
)
    resume_count >= 0 || throw(
        ArgumentError("Fig. 2 progress resume count must be nonnegative")
    )
    root = abspath(candidate_directory)
    mkpath(root)
    identity = _fig2_progress_identity(
        spec, dimension, point, phi_y, candidate_id
    )
    schedule = _fig2_progress_schedule(maxdim_schedule, dimension)
    event.stage <= length(schedule) && schedule[event.stage] == event.target || throw(
        ArgumentError("Fig. 2 progress event target disagrees with its schedule")
    )
    previous = _fig2_read_progress_pointer(root, identity)
    event_sequence = isnothing(previous) ? 1 :
        _fig2_progress_positive_integer(previous, "event_sequence") + 1
    if isnothing(previous)
        resume_count == 0 && event.sequence == 1 || throw(ArgumentError(
            "first Fig. 2 progress event must start sequence one without a resume"
        ))
    else
        previous_resume = _fig2_progress_nonnegative_integer(
            previous, "resume_count"
        )
        previous_runner = _fig2_progress_positive_integer(
            previous, "runner_event_sequence"
        )
        same_segment = resume_count == previous_resume &&
            event.sequence == previous_runner + 1
        resumed_segment = resume_count == previous_resume + 1 &&
            event.sequence == 1
        same_segment || resumed_segment || throw(ArgumentError(
            "Fig. 2 progress runner sequence or resume count is discontinuous"
        ))
    end

    config = with_flux(spec.config, Float64(phi_y))
    canonical_seed = Int(mod(
        _derived_seed(
            _fig2_candidate_seed(dimension, point, candidate_id),
            :progress,
            event_sequence,
        ),
        UInt64(typemax(Int)),
    ))
    canonical = isnothing(canonicalize_state) ?
        _canonicalize_vumps_state(event.psi; rng_seed=canonical_seed) :
        canonicalize_state(event.psi)
    canonical isa InfiniteCanonicalMPS || throw(
        ArgumentError("Fig. 2 progress canonicalizer returned the wrong state type")
    )
    maxlinkdim = maximum(link_dimensions(canonical))
    maxlinkdim > 0 || error("Fig. 2 progress state has no positive link dimension")

    generation_root = joinpath(root, ".progress", "generations")
    event_root = joinpath(root, ".progress", "events")
    mkpath(generation_root)
    mkpath(event_root)
    stem = string(
        lpad(string(event_sequence), 8, '0'),
        "_",
        getpid(),
        "_",
        time_ns(),
    )
    state_path = joinpath(generation_root, "state_$stem.h5")
    event_path = joinpath(event_root, "event_$stem.toml")
    !ispath(state_path) && !ispath(event_path) || error(
        "Fig. 2 progress generation path collision"
    )
    save_state(state_path, canonical, config)
    isfile(state_path) && filesize(state_path) > 0 || throw(
        ArgumentError("Fig. 2 progress state writer produced no checkpoint")
    )
    state_relative = relpath(state_path, root)
    event_relative = relpath(event_path, root)
    state_sha256 = _fig2_file_sha256(state_path)
    previous_event_path = isnothing(previous) ? "" :
        String(previous["event_path"])
    previous_event_sha256 = isnothing(previous) ? "" :
        String(previous["event_sha256"])
    event_data = Dict{String,Any}(
        "format" => FIG2_PROGRESS_EVENT_FORMAT,
        identity...,
        "event_sequence" => event_sequence,
        "runner_event_sequence" => event.sequence,
        "resume_count" => Int(resume_count),
        "event_kind" => string(event.kind),
        "stage" => event.stage,
        "iteration" => event.iteration,
        "target" => event.target,
        "maxdim_schedule" => schedule,
        "maxlinkdim" => maxlinkdim,
        "state_path" => state_relative,
        "state_sha256" => state_sha256,
        "previous_event_path" => previous_event_path,
        "previous_event_sha256" => previous_event_sha256,
        "pbs_job_id" => get(ENV, "PBS_JOBID", ""),
    )
    merge!(event_data, _fig2_progress_event_payload(event))
    _write_fig2_toml(event_path, event_data)
    event_sha256 = _fig2_file_sha256(event_path)
    pointer = Dict{String,Any}(
        "format" => FIG2_PROGRESS_POINTER_FORMAT,
        identity...,
        "event_sequence" => event_sequence,
        "runner_event_sequence" => event.sequence,
        "resume_count" => Int(resume_count),
        "event_kind" => string(event.kind),
        "stage" => event.stage,
        "iteration" => event.iteration,
        "target" => event.target,
        "maxdim_schedule" => schedule,
        "maxlinkdim" => maxlinkdim,
        "state_path" => state_relative,
        "state_sha256" => state_sha256,
        "event_path" => event_relative,
        "event_sha256" => event_sha256,
    )
    _write_fig2_toml(joinpath(root, ".progress", "latest.toml"), pointer)
    chain = _fig2_validate_progress_chain(root, identity, pointer)
    _fig2_prune_progress_generations!(
        root,
        chain.state_paths,
        Int(identity["progress_generations_to_keep"]),
    )
    return (;
        event_sequence,
        runner_event_sequence=event.sequence,
        resume_count=Int(resume_count),
        verified_event_count=chain.verified_event_count,
        maxlinkdim,
        state_path=state_relative,
        event_path=event_relative,
    )
end

function _fig2_validate_progress_chain(root, identity, pointer)
    expected_sequence = _fig2_progress_positive_integer(
        pointer, "event_sequence"
    )
    event_path = String(get(pointer, "event_path", ""))
    event_sha256 = String(get(pointer, "event_sha256", ""))
    isempty(event_path) && throw(
        ArgumentError("Fig. 2 progress pointer event path is missing")
    )
    occursin(r"^[0-9a-f]{64}$", event_sha256) || throw(
        ArgumentError("Fig. 2 progress pointer event checksum is invalid")
    )
    latest = nothing
    verified = 0
    state_paths = String[]
    while true
        absolute, relative = _fig2_progress_path(
            root, event_path, "event"
        )
        isfile(absolute) || throw(
            ArgumentError("Fig. 2 progress event file is missing")
        )
        _fig2_file_sha256(absolute) == event_sha256 || throw(
            ArgumentError("Fig. 2 progress event checksum mismatch")
        )
        event = Dict{String,Any}(TOML.parsefile(absolute))
        get(event, "format", "") == FIG2_PROGRESS_EVENT_FORMAT || throw(
            ArgumentError("unsupported Fig. 2 progress event format")
        )
        _fig2_validate_progress_identity(event, identity)
        _fig2_progress_positive_integer(event, "event_sequence") ==
            expected_sequence || throw(ArgumentError(
            "Fig. 2 progress event sequence is discontinuous"
        ))
        _fig2_progress_positive_integer(event, "runner_event_sequence")
        _fig2_progress_nonnegative_integer(event, "resume_count")
        _fig2_progress_positive_integer(event, "maxlinkdim")
        _fig2_progress_schedule(event["maxdim_schedule"], identity["dimension"])
        _, state_relative = _fig2_progress_path(
            root, get(event, "state_path", nothing), "state"
        )
        push!(state_paths, state_relative)
        isnothing(latest) && (latest = event)
        verified += 1
        previous_path = get(event, "previous_event_path", nothing)
        previous_sha256 = get(event, "previous_event_sha256", nothing)
        previous_path isa AbstractString && previous_sha256 isa AbstractString ||
            throw(ArgumentError(
                "Fig. 2 progress previous-event reference is invalid"
            ))
        if expected_sequence == 1
            isempty(previous_path) && isempty(previous_sha256) || throw(
                ArgumentError("root Fig. 2 progress event has a predecessor")
            )
            break
        end
        !isempty(previous_path) && occursin(r"^[0-9a-f]{64}$", previous_sha256) ||
            throw(ArgumentError(
                "Fig. 2 progress previous-event checksum is invalid"
            ))
        event_path = String(previous_path)
        event_sha256 = String(previous_sha256)
        expected_sequence -= 1
    end
    return (;
        latest=something(latest),
        verified_event_count=verified,
        state_paths,
    )
end

function _fig2_load_progress(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id;
    load_state=load_checkpoint,
)
    root = abspath(candidate_directory)
    identity = _fig2_progress_identity(
        spec, dimension, point, phi_y, candidate_id
    )
    pointer = _fig2_read_progress_pointer(root, identity)
    isnothing(pointer) && return nothing
    chain = _fig2_validate_progress_chain(root, identity, pointer)
    latest = chain.latest
    for key in (
        "event_sequence",
        "runner_event_sequence",
        "resume_count",
        "event_kind",
        "stage",
        "iteration",
        "target",
        "maxdim_schedule",
        "maxlinkdim",
        "state_path",
        "state_sha256",
    )
        get(pointer, key, nothing) == get(latest, key, nothing) || throw(
            ArgumentError("Fig. 2 progress pointer disagrees with its latest event")
        )
    end
    state_path, state_relative = _fig2_progress_path(
        root, pointer["state_path"], "state"
    )
    isfile(state_path) && filesize(state_path) > 0 || throw(
        ArgumentError("Fig. 2 progress state checkpoint is missing")
    )
    state_sha256 = get(pointer, "state_sha256", nothing)
    state_sha256 isa AbstractString &&
        occursin(r"^[0-9a-f]{64}$", state_sha256) || throw(
        ArgumentError("Fig. 2 progress state checksum is invalid")
    )
    _fig2_file_sha256(state_path) == state_sha256 || throw(
        ArgumentError("Fig. 2 progress state checksum mismatch")
    )
    config = with_flux(spec.config, Float64(phi_y))
    state = load_state(state_path, config)
    state isa InfiniteCanonicalMPS || throw(
        ArgumentError("Fig. 2 progress loader returned the wrong state type")
    )
    maxlinkdim = maximum(link_dimensions(state))
    maxlinkdim == _fig2_progress_positive_integer(pointer, "maxlinkdim") || throw(
        ArgumentError("Fig. 2 progress state maxlinkdim disagrees with its pointer")
    )
    event_sequence = _fig2_progress_positive_integer(
        pointer, "event_sequence"
    )
    resume_count = _fig2_progress_nonnegative_integer(
        pointer, "resume_count"
    )
    return (;
        state,
        event_sequence,
        runner_event_sequence=_fig2_progress_positive_integer(
            pointer, "runner_event_sequence"
        ),
        resume_count,
        next_resume_count=resume_count + 1,
        verified_event_count=chain.verified_event_count,
        maxlinkdim,
        state_path=state_relative,
        event_path=String(pointer["event_path"]),
        latest_pointer_sha256=_fig2_file_sha256(
            joinpath(root, ".progress", "latest.toml")
        ),
        pointer,
    )
end

function _fig2_finalize_progress!(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    final_checkpoint,
)
    root = abspath(candidate_directory)
    progress = _fig2_load_progress(
        spec, root, dimension, point, phi_y, candidate_id
    )
    isnothing(progress) && throw(
        ArgumentError("Fig. 2 candidate has no restartable progress generation")
    )
    final_path = abspath(final_checkpoint)
    relpath(final_path, root) == "state.h5" || throw(
        ArgumentError("Fig. 2 final checkpoint is outside the candidate contract")
    )
    isfile(final_path) && filesize(final_path) > 0 || throw(
        ArgumentError("Fig. 2 final checkpoint is missing")
    )
    config = with_flux(spec.config, Float64(phi_y))
    final_state = load_checkpoint(final_path, config)
    final_maxlinkdim = maximum(link_dimensions(final_state))
    final_maxlinkdim == progress.maxlinkdim || throw(
        ArgumentError("Fig. 2 final and progress checkpoint maxlinkdim disagree")
    )
    final_state_sha256 = _fig2_file_sha256(final_path)
    pointer_path = joinpath(root, ".progress", "latest.toml")
    summary = Dict{String,Any}(
        "format" => FIG2_PROGRESS_SUMMARY_FORMAT,
        _fig2_progress_identity(
            spec, dimension, point, phi_y, candidate_id
        )...,
        "complete" => true,
        "event_count" => progress.verified_event_count,
        "resume_count" => progress.resume_count,
        "latest_maxlinkdim" => progress.maxlinkdim,
        "latest_pointer_path" => relpath(pointer_path, root),
        "latest_pointer_sha256" => progress.latest_pointer_sha256,
        "latest_event_path" => progress.event_path,
        "latest_event_sha256" => String(progress.pointer["event_sha256"]),
        "latest_state_path" => progress.state_path,
        "latest_state_sha256" => String(progress.pointer["state_sha256"]),
        "final_state_path" => relpath(final_path, root),
        "final_state_sha256" => final_state_sha256,
        "final_maxlinkdim" => final_maxlinkdim,
    )
    _write_fig2_toml(joinpath(root, "progress.toml"), summary)
    return _fig2_validate_progress_artifact(
        spec,
        root,
        dimension,
        point,
        phi_y,
        candidate_id,
        final_path,
    )
end

function _fig2_validate_progress_artifact(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    final_checkpoint,
)
    root = abspath(candidate_directory)
    path = joinpath(root, "progress.toml")
    isfile(path) || throw(
        ArgumentError("Fig. 2 candidate progress.toml is missing")
    )
    summary = Dict{String,Any}(TOML.parsefile(path))
    get(summary, "format", "") == FIG2_PROGRESS_SUMMARY_FORMAT || throw(
        ArgumentError("unsupported Fig. 2 progress summary format")
    )
    get(summary, "complete", nothing) === true || throw(
        ArgumentError("Fig. 2 progress summary is not complete")
    )
    identity = _fig2_progress_identity(
        spec, dimension, point, phi_y, candidate_id
    )
    _fig2_validate_progress_identity(summary, identity)
    progress = _fig2_load_progress(
        spec, root, dimension, point, phi_y, candidate_id
    )
    isnothing(progress) && throw(
        ArgumentError("Fig. 2 progress summary has no latest pointer")
    )
    get(summary, "event_count", nothing) == progress.verified_event_count ||
        throw(ArgumentError("Fig. 2 progress summary event count disagrees"))
    get(summary, "resume_count", nothing) == progress.resume_count || throw(
        ArgumentError("Fig. 2 progress summary resume count disagrees")
    )
    get(summary, "latest_maxlinkdim", nothing) == progress.maxlinkdim || throw(
        ArgumentError("Fig. 2 progress summary maxlinkdim disagrees")
    )
    get(summary, "latest_pointer_sha256", nothing) ==
        progress.latest_pointer_sha256 || throw(ArgumentError(
        "Fig. 2 progress summary pointer checksum disagrees"
    ))
    final_path = abspath(final_checkpoint)
    final_state_sha256 = _fig2_file_sha256(final_path)
    get(summary, "final_state_sha256", nothing) == final_state_sha256 || throw(
        ArgumentError("Fig. 2 progress final-state checksum disagrees")
    )
    config = with_flux(spec.config, Float64(phi_y))
    final_state = load_checkpoint(final_path, config)
    final_maxlinkdim = maximum(link_dimensions(final_state))
    get(summary, "final_maxlinkdim", nothing) == final_maxlinkdim ==
        progress.maxlinkdim || throw(ArgumentError(
        "Fig. 2 progress final-state maxlinkdim disagrees"
    ))
    return (;
        complete=true,
        event_count=progress.verified_event_count,
        resume_count=progress.resume_count,
        latest_maxlinkdim=progress.maxlinkdim,
        final_state_sha256,
        progress_sha256=_fig2_file_sha256(path),
    )
end

function _default_fig2_progress_audit(
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    final_checkpoint,
)
    return _fig2_validate_progress_artifact(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        final_checkpoint,
    )
end

function _fig2_candidate_progress_audit(
    progress_audit,
    spec,
    candidate_directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    final_checkpoint,
    expected_maxlinkdim,
)
    result = progress_audit(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        final_checkpoint,
    )
    result isa NamedTuple || throw(
        ArgumentError("Fig. 2 candidate progress audit result is invalid")
    )
    required = (
        :complete,
        :event_count,
        :resume_count,
        :latest_maxlinkdim,
        :final_state_sha256,
        :progress_sha256,
    )
    all(name -> hasproperty(result, name), required) || throw(
        ArgumentError("Fig. 2 candidate progress audit fields are incomplete")
    )
    result.complete === true || throw(
        ArgumentError("Fig. 2 candidate progress audit is not complete")
    )
    event_count = result.event_count
    event_count isa Integer && !(event_count isa Bool) && event_count > 0 || throw(
        ArgumentError("Fig. 2 candidate progress event count is invalid")
    )
    resume_count = result.resume_count
    resume_count isa Integer && !(resume_count isa Bool) && resume_count >= 0 || throw(
        ArgumentError("Fig. 2 candidate progress resume count is invalid")
    )
    latest_maxlinkdim = result.latest_maxlinkdim
    latest_maxlinkdim isa Integer && !(latest_maxlinkdim isa Bool) &&
        latest_maxlinkdim > 0 || throw(
        ArgumentError("Fig. 2 candidate progress maxlinkdim is invalid")
    )
    expected_maxlinkdim isa Integer && !(expected_maxlinkdim isa Bool) &&
        expected_maxlinkdim > 0 || throw(
        ArgumentError("Fig. 2 candidate expected progress maxlinkdim is invalid")
    )
    Int(latest_maxlinkdim) == Int(expected_maxlinkdim) || throw(
        ArgumentError(
            "Fig. 2 candidate progress maxlinkdim disagrees with achieved maxlinkdim"
        )
    )
    final_state_sha256 = result.final_state_sha256
    progress_sha256 = result.progress_sha256
    all(
        value -> value isa AbstractString && occursin(r"^[0-9a-f]{64}$", value),
        (final_state_sha256, progress_sha256),
    ) || throw(
        ArgumentError("Fig. 2 candidate progress checksum is invalid")
    )
    final_path = abspath(final_checkpoint)
    progress_path = joinpath(abspath(candidate_directory), "progress.toml")
    isfile(final_path) && filesize(final_path) > 0 || throw(
        ArgumentError("Fig. 2 candidate final checkpoint is missing")
    )
    isfile(progress_path) && filesize(progress_path) > 0 || throw(
        ArgumentError("Fig. 2 candidate progress.toml is missing")
    )
    _fig2_file_sha256(final_path) == final_state_sha256 || throw(
        ArgumentError("Fig. 2 candidate progress final-state checksum disagrees")
    )
    _fig2_file_sha256(progress_path) == progress_sha256 || throw(
        ArgumentError("Fig. 2 candidate progress.toml checksum disagrees")
    )
    return (;
        complete=true,
        event_count=Int(event_count),
        resume_count=Int(resume_count),
        latest_maxlinkdim=Int(latest_maxlinkdim),
        final_state_sha256=String(final_state_sha256),
        progress_sha256=String(progress_sha256),
    )
end

function _ensure_fig2_manifest_copy(spec::Fig2BenchmarkSpec, output)
    _fig2_validated_snapshot(spec)
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

function _fig2_git_directories(repository::AbstractString)
    dotgit = joinpath(abspath(repository), ".git")
    git_directory = if isdir(dotgit)
        normpath(dotgit)
    elseif isfile(dotgit)
        declaration = strip(read(dotgit, String))
        startswith(declaration, "gitdir: ") || throw(ArgumentError(
            "Fig. 2 repository .git file has an invalid gitdir declaration"
        ))
        target = strip(declaration[9:end])
        isempty(target) && throw(ArgumentError(
            "Fig. 2 repository .git file has an empty gitdir target"
        ))
        normpath(isabspath(target) ? target : joinpath(repository, target))
    else
        throw(ArgumentError("Fig. 2 repository has no .git directory or file"))
    end
    isdir(git_directory) || throw(ArgumentError(
        "Fig. 2 repository git directory does not exist"
    ))

    commondir_path = joinpath(git_directory, "commondir")
    common_directory = if isfile(commondir_path)
        target = strip(read(commondir_path, String))
        isempty(target) && throw(ArgumentError(
            "Fig. 2 repository commondir target is empty"
        ))
        normpath(isabspath(target) ? target : joinpath(git_directory, target))
    else
        git_directory
    end
    isdir(common_directory) || throw(ArgumentError(
        "Fig. 2 repository common git directory does not exist"
    ))
    return git_directory, common_directory
end

function _fig2_git_reference_parts(reference::AbstractString)
    parts = split(String(reference), '/'; keepempty=true)
    length(parts) >= 3 && first(parts) == "refs" &&
        all(part -> !isempty(part) && part != "." && part != "..", parts) ||
        throw(ArgumentError("Fig. 2 repository HEAD reference is invalid"))
    return parts
end

function _fig2_packed_reference(common_directory, reference)
    path = joinpath(common_directory, "packed-refs")
    isfile(path) || throw(ArgumentError(
        "Fig. 2 repository HEAD reference is missing from loose and packed refs"
    ))
    commits = String[]
    for line in eachline(path)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, '#') ||
            startswith(stripped, '^')) && continue
        fields = split(stripped)
        length(fields) == 2 && fields[2] == reference || continue
        push!(commits, fields[1])
    end
    length(commits) == 1 || throw(ArgumentError(
        "Fig. 2 repository packed HEAD reference is missing or ambiguous"
    ))
    return only(commits)
end

function _fig2_repository_commit(repository::AbstractString)
    git_directory, common_directory = _fig2_git_directories(repository)
    head_path = joinpath(git_directory, "HEAD")
    isfile(head_path) || throw(ArgumentError(
        "Fig. 2 repository git directory has no HEAD"
    ))
    head = strip(read(head_path, String))
    commit = if startswith(head, "ref: ")
        reference = strip(head[6:end])
        parts = _fig2_git_reference_parts(reference)
        loose_paths = (
            joinpath(git_directory, parts...),
            joinpath(common_directory, parts...),
        )
        loose_index = findfirst(isfile, loose_paths)
        isnothing(loose_index) ?
            _fig2_packed_reference(common_directory, reference) :
            strip(read(loose_paths[loose_index], String))
    else
        head
    end
    occursin(r"^[0-9a-f]{40}$", commit) || throw(
        ArgumentError("Fig. 2 benchmark git commit is not a full 40-hex object ID")
    )
    return commit
end

_fig2_repository_commit() = _fig2_repository_commit(
    normpath(joinpath(@__DIR__, "..", "..", ".."))
)

function _fig2_project_manifest()
    active_project = Base.active_project()
    isnothing(active_project) && throw(
        ArgumentError("Fig. 2 benchmark requires an active Julia project")
    )
    project_manifest = joinpath(dirname(active_project), "Manifest.toml")
    isfile(project_manifest) || throw(
        ArgumentError("Fig. 2 benchmark active project has no Manifest.toml")
    )
    return abspath(active_project), abspath(project_manifest)
end

function _fig2_parse_checkout_audit(contents::AbstractString)
    expected_keys = Set((
        "head_commit",
        "origin_dmrg",
        "tracked_dirty",
        "staged_dirty",
    ))
    fields = Dict{String,String}()
    for line in split(chomp(contents), '\n')
        pieces = split(line, '='; limit=2)
        length(pieces) == 2 || throw(WorkflowValidationError(
            "checkout audit contains a malformed output line"
        ))
        key, value = pieces
        key in expected_keys || throw(WorkflowValidationError(
            "checkout audit contains an unknown output key"
        ))
        haskey(fields, key) && throw(WorkflowValidationError(
            "checkout audit contains a duplicate output key"
        ))
        fields[key] = value
    end
    Set(keys(fields)) == expected_keys || throw(WorkflowValidationError(
        "checkout audit output is incomplete"
    ))

    parse_dirty(key) = if fields[key] == "true"
        true
    elseif fields[key] == "false"
        false
    else
        throw(WorkflowValidationError(
            "checkout audit dirty flags must be true or false"
        ))
    end
    head_commit = fields["head_commit"]
    origin_commit = fields["origin_dmrg"]
    occursin(r"^[0-9a-f]{40}$", head_commit) || throw(
        WorkflowValidationError("checkout audit HEAD commit is invalid")
    )
    occursin(r"^[0-9a-f]{40}$", origin_commit) || throw(
        WorkflowValidationError("checkout audit origin/DMRG commit is invalid")
    )
    return (
        head_commit,
        origin_commit,
        tracked_dirty=parse_dirty("tracked_dirty"),
        staged_dirty=parse_dirty("staged_dirty"),
    )
end

function _fig2_checkout_audit(repository::AbstractString)
    auditor = normpath(joinpath(
        @__DIR__,
        "..",
        "bin",
        "audit_production_checkout.jl",
    ))
    isfile(auditor) || throw(WorkflowValidationError(
        "Fig. 2 checkout auditor is missing"
    ))
    stdout = IOBuffer()
    stderr = IOBuffer()
    command = `$(Base.julia_cmd()) --startup-file=no $auditor $repository`
    process = run(
        pipeline(
            ignorestatus(command);
            stdout=stdout,
            stderr=stderr,
        ),
    )
    audit_output = String(take!(stdout))
    audit_errors = strip(String(take!(stderr)))
    process.exitcode == 0 || throw(WorkflowValidationError(
        isempty(audit_errors) ?
            "Fig. 2 checkout audit failed closed" :
            "Fig. 2 checkout audit failed closed: $audit_errors"
    ))
    audit = _fig2_parse_checkout_audit(audit_output)
    !audit.tracked_dirty || throw(WorkflowValidationError(
        "default Fig. 2 provenance refuses tracked worktree changes"
    ))
    !audit.staged_dirty || throw(WorkflowValidationError(
        "default Fig. 2 provenance refuses staged index changes"
    ))
    audit.head_commit == audit.origin_commit || throw(
        WorkflowValidationError(
            "default Fig. 2 provenance requires HEAD == origin/DMRG"
        )
    )
    return audit
end

function _default_fig2_provenance(spec, output, runtime_seconds)
    repository = normpath(joinpath(@__DIR__, "..", "..", ".."))
    checkout_audit = _fig2_checkout_audit(repository)
    active_project, project_manifest = _fig2_project_manifest()
    return Dict{String,Any}(
        "format" => "fqahc_fig2_provenance_v2",
        "manifest_sha256" => spec.sha256,
        "git_commit" => checkout_audit.head_commit,
        "git_tree_clean" => true,
        "julia_version" => string(VERSION),
        "pbs_job_id" => get(ENV, "PBS_JOBID", ""),
        "threads" => Threads.nthreads(),
        "blas_threads" => LinearAlgebra.BLAS.get_num_threads(),
        "strided_threads" => ITensors.Strided.get_num_threads(),
        "blocksparse_threaded" => ITensors.using_threaded_blocksparse(),
        "runtime_seconds" => runtime_seconds,
        "active_project" => active_project,
        "project_manifest" => project_manifest,
        "project_manifest_sha256" => _fig2_file_sha256(project_manifest),
        "benchmark_source" => abspath(@__FILE__),
        "benchmark_source_sha256" => _fig2_file_sha256(@__FILE__),
    )
end

const FIG2_GENERATION_PROVENANCE_KEYS = (
    "manifest_sha256",
    "git_commit",
    "git_tree_clean",
    "julia_version",
    "pbs_job_id",
    "threads",
    "blas_threads",
    "strided_threads",
    "blocksparse_threaded",
    "active_project",
    "project_manifest",
    "project_manifest_sha256",
    "benchmark_source",
    "benchmark_source_sha256",
)

function _fig2_threading_provenance_valid(provenance)
    threads = get(provenance, "threads", nothing)
    blas_threads = get(provenance, "blas_threads", nothing)
    strided_threads = get(provenance, "strided_threads", nothing)
    blocksparse_threaded = get(provenance, "blocksparse_threaded", nothing)
    return threads isa Integer && !(threads isa Bool) &&
        Int(threads) == Threads.nthreads() &&
        blas_threads isa Integer && !(blas_threads isa Bool) &&
        Int(blas_threads) == 1 &&
        strided_threads isa Integer && !(strided_threads isa Bool) &&
        Int(strided_threads) == 1 &&
        blocksparse_threaded isa Bool &&
        blocksparse_threaded == (Int(threads) > 1)
end

function _fig2_generation_provenance(spec, provenance)
    provenance isa AbstractDict || throw(
        ArgumentError("Fig. 2 provenance supplier must return a dictionary")
    )
    active_project, project_manifest = _fig2_project_manifest()
    benchmark_source = abspath(@__FILE__)
    get(provenance, "format", "") == "fqahc_fig2_provenance_v2" || throw(
        ArgumentError("unsupported Fig. 2 execution provenance format")
    )
    get(provenance, "manifest_sha256", "") == spec.sha256 || throw(
        ArgumentError("execution provenance manifest checksum does not match")
    )
    commit = get(provenance, "git_commit", "")
    commit isa AbstractString && occursin(r"^[0-9a-f]{40}$", commit) || throw(
        ArgumentError("execution provenance git commit is invalid")
    )
    commit == _fig2_repository_commit() || throw(ArgumentError(
        "execution provenance git commit is not the actual repository HEAD"
    ))
    get(provenance, "git_tree_clean", false) === true || throw(
        ArgumentError("execution provenance does not certify a clean tracked tree")
    )
    get(provenance, "julia_version", "") == "1.12.5" &&
        get(provenance, "julia_version", "") == string(VERSION) || throw(
        ArgumentError("execution provenance Julia version is invalid")
    )
    pbs_job_id = get(provenance, "pbs_job_id", "")
    pbs_job_id isa AbstractString &&
        occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", pbs_job_id) || throw(
        ArgumentError("execution provenance PBS job ID is invalid")
    )
    _fig2_threading_provenance_valid(provenance) || throw(
        ArgumentError("execution provenance measured threading state is invalid")
    )
    get(provenance, "active_project", "") == active_project &&
        get(provenance, "project_manifest", "") == project_manifest &&
        get(provenance, "project_manifest_sha256", "") ==
            _fig2_file_sha256(project_manifest) || throw(
        ArgumentError("execution provenance project manifest is invalid")
    )
    get(provenance, "benchmark_source", "") == benchmark_source &&
        get(provenance, "benchmark_source_sha256", "") ==
            _fig2_file_sha256(benchmark_source) || throw(
        ArgumentError("execution provenance benchmark source is invalid")
    )
    generation = Dict{String,Any}(
        key => deepcopy(provenance[key]) for key in FIG2_GENERATION_PROVENANCE_KEYS
    )
    generation["format"] = "fqahc_fig2_candidate_provenance_v1"
    return generation
end

function _fig2_generation_provenance_sha256(generation)
    return _fig2_sha256(collect(codeunits(_render_summary(generation))))
end

function _fig2_validate_generation_compatibility(stored, current)
    get(stored, "format", "") == "fqahc_fig2_candidate_provenance_v1" || throw(
        ArgumentError("unsupported candidate generation provenance format")
    )
    immutable_keys = (
        "manifest_sha256",
        "git_commit",
        "git_tree_clean",
        "julia_version",
        "threads",
        "blas_threads",
        "strided_threads",
        "blocksparse_threaded",
        "active_project",
        "project_manifest",
        "project_manifest_sha256",
        "benchmark_source",
        "benchmark_source_sha256",
    )
    all(key -> haskey(stored, key) && stored[key] == current[key], immutable_keys) ||
        throw(ArgumentError(
            "persisted candidate provenance is incompatible with the current commit, source, or project manifest"
        ))
    actual_commit = _fig2_repository_commit()
    get(stored, "git_commit", "") == actual_commit &&
        get(current, "git_commit", "") == actual_commit || throw(ArgumentError(
        "persisted candidate provenance git commit is not the actual repository HEAD"
    ))
    pbs_job_id = get(stored, "pbs_job_id", "")
    pbs_job_id isa AbstractString &&
        occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", pbs_job_id) || throw(
        ArgumentError("persisted candidate provenance PBS job ID is invalid")
    )
    return true
end

function _fig2_candidate_pbs_job_ids(root, ledger)
    jobs = String[]
    for row in ledger["candidate"]
        metadata = TOML.parsefile(joinpath(root, String(row["directory"]), "candidate.toml"))
        generation = get(metadata, "generation_provenance", Dict{String,Any}())
        job = get(generation, "pbs_job_id", "")
        isempty(job) || push!(jobs, String(job))
    end
    return sort!(unique(jobs))
end

function _write_fig2_provenance!(
    spec,
    output,
    started_ns,
    provenance,
    ledger,
)
    runtime_seconds = max(0.0, (time_ns() - started_ns) / 1.0e9)
    finalized = Dict{String,Any}(deepcopy(provenance))
    get(finalized, "git_commit", "") == _fig2_repository_commit() || throw(
        ArgumentError(
            "refusing to write Fig. 2 provenance for a commit other than the actual repository HEAD"
        )
    )
    finalized["runtime_seconds"] = runtime_seconds
    finalized["candidate_pbs_job_ids"] = _fig2_candidate_pbs_job_ids(output, ledger)
    _write_fig2_toml(joinpath(output, "provenance.toml"), finalized)
    return finalized
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
    safe_id = String(candidate_id)
    occursin(r"^[A-Za-z0-9_.-]+$", safe_id) || throw(
        ArgumentError("candidate ID must match [A-Za-z0-9_.-]+ exactly")
    )
    return joinpath(
        "D_$(lpad(dimension, 4, '0'))",
        "phi_$(lpad(point - 1, 4, '0'))",
        "candidate_$safe_id",
    )
end

_fig2_same_float(left, right) = isequal(Float64(left), Float64(right))
_fig2_same_complex(left, right) =
    _fig2_same_float(real(left), real(right)) &&
    _fig2_same_float(imag(left), imag(right))

function _fig2_tsv_text(value)
    return replace(String(value), '\t' => ' ', '\n' => ' ', '\r' => ' ')
end

function _fig2_missing_mixed_fidelity()
    missing_value = ComplexF64(NaN, NaN)
    return MixedTransferFidelity(
        NaN,
        missing_value,
        missing_value,
        missing_value,
        NaN,
        NaN,
        NaN,
        false,
        false,
        "no previous selected state at the reference flux",
    )
end

function _fig2_mixed_audit(evidence::Fig2CandidateEvidence, point)
    if point == 1
        isnothing(evidence.mixed_fidelity) || throw(
            ArgumentError("zero-flux candidate must not carry mixed-fidelity evidence")
        )
        evidence.fidelity_valid && throw(
            ArgumentError("zero-flux candidate must not claim previous-state fidelity")
        )
        isnan(evidence.fidelity_to_previous) || throw(
            ArgumentError("zero-flux candidate fidelity must be NaN")
        )
        return _fig2_missing_mixed_fidelity()
    end
    fidelity = evidence.mixed_fidelity
    fidelity isa MixedTransferFidelity || throw(
        ArgumentError("later Fig. 2 candidate is missing mixed-transfer evidence")
    )
    evidence.fidelity_valid == fidelity.valid || throw(
        ArgumentError("mixed-transfer validity disagrees with candidate fidelity validity")
    )
    _fig2_same_float(evidence.fidelity_to_previous, fidelity.fidelity_cell) || throw(
        ArgumentError("mixed-transfer fidelity disagrees with candidate fidelity value")
    )
    return fidelity
end

function _fig2_same_counting(left::Fig2MomentumCounting, right::Fig2MomentumCounting)
    return left.counts == right.counts &&
        left.delta_ky == right.delta_ky &&
        length(left.cutoffs) == length(right.cutoffs) &&
        all(_fig2_same_float.(left.cutoffs, right.cutoffs)) &&
        left.reference_raw_charge == right.reference_raw_charge &&
        left.reference_level == right.reference_level &&
        left.reference_ky_index == right.reference_ky_index &&
        _fig2_same_float(left.minimum_cutoff_margin, right.minimum_cutoff_margin) &&
        left.valid == right.valid && left.reason == right.reason
end

function _fig2_counting_fields(counting::Fig2MomentumCounting)
    references = (
        counting.reference_raw_charge,
        counting.reference_level,
        counting.reference_ky_index,
    )
    has_reference = all(value -> !isnothing(value), references)
    (has_reference || all(isnothing, references)) || throw(
        ArgumentError("momentum counting reference fields must be all present or all absent")
    )
    return Dict{String,Any}(
        "momentum_counting" => copy(counting.counts),
        "momentum_counting_valid" => counting.valid,
        "momentum_counting_reason" => counting.reason,
        "momentum_counting_delta_ky" => copy(counting.delta_ky),
        "momentum_counting_cutoffs" => copy(counting.cutoffs),
        "momentum_counting_has_reference" => has_reference,
        "momentum_counting_reference_raw_charge" =>
            has_reference ? Int(counting.reference_raw_charge) : 0,
        "momentum_counting_reference_level" =>
            has_reference ? Int(counting.reference_level) : 0,
        "momentum_counting_reference_ky_index" =>
            has_reference ? Int(counting.reference_ky_index) : 0,
        "momentum_counting_minimum_cutoff_margin" =>
            counting.minimum_cutoff_margin,
    )
end

function _fig2_counting_from_data(data)
    has_reference = Bool(data["momentum_counting_has_reference"])
    return Fig2MomentumCounting(
        Int.(data["momentum_counting"]),
        Int.(data["momentum_counting_delta_ky"]),
        Float64.(data["momentum_counting_cutoffs"]),
        has_reference ? Int(data["momentum_counting_reference_raw_charge"]) : nothing,
        has_reference ? Int(data["momentum_counting_reference_level"]) : nothing,
        has_reference ? Int(data["momentum_counting_reference_ky_index"]) : nothing,
        Float64(data["momentum_counting_minimum_cutoff_margin"]),
        Bool(data["momentum_counting_valid"]),
        String(data["momentum_counting_reason"]),
    )
end

function _fig2_momentum_fields(momentum)
    validation = momentum.validation
    return Dict{String,Any}(
        "momentum_unitary_residual" => validation.unitary_residual,
        "momentum_power_residual" => validation.power_residual,
        "momentum_commutator_residual" => validation.commutator_residual,
        "momentum_transfer_residual" => validation.transfer_residual,
        "momentum_validation_valid" => validation.valid,
        "momentum_validation_reason" => validation.reason,
        "momentum_valid" => momentum.valid,
        "momentum_reason" => momentum.reason,
        "momentum_convention" => string(momentum.momentum_convention),
        "momentum_sector_ky_offsets" => Dict(
            string(charge) => offset for
            (charge, offset) in momentum.sector_ky_offsets
        ),
    )
end

function _fig2_momentum_validation_from_data(data)
    return MomentumValidation(
        Float64(data["momentum_unitary_residual"]),
        Float64(data["momentum_power_residual"]),
        Float64(data["momentum_commutator_residual"]),
        Float64(data["momentum_transfer_residual"]),
        Bool(data["momentum_validation_valid"]),
        String(data["momentum_validation_reason"]),
    )
end

function _fig2_sector_offsets_from_data(data)
    return Dict(
        parse(Int, charge) => Int(offset) for
        (charge, offset) in data["momentum_sector_ky_offsets"]
    )
end

function _fig2_sector_offsets_text(offsets)
    return join(
        ("$charge:$(offsets[charge])" for charge in sort!(collect(keys(offsets)))),
        ",",
    )
end

function _fig2_mixed_fields(fidelity::MixedTransferFidelity; present::Bool)
    return Dict{String,Any}(
        "mixed_present" => present,
        "mixed_value_real" => real(fidelity.mixed_value),
        "mixed_value_imag" => imag(fidelity.mixed_value),
        "self_value1_real" => real(fidelity.self_value1),
        "self_value1_imag" => imag(fidelity.self_value1),
        "self_value2_real" => real(fidelity.self_value2),
        "self_value2_imag" => imag(fidelity.self_value2),
        "mixed_residual" => fidelity.mixed_residual,
        "self_residual1" => fidelity.self_residual1,
        "self_residual2" => fidelity.self_residual2,
        "mixed_converged" => fidelity.converged,
        "mixed_valid" => fidelity.valid,
        "mixed_reason" => fidelity.reason,
    )
end

function _fig2_absent_mixed_reference()
    return Fig2MixedReference(
        false,
        0,
        0,
        FIG2_MIXED_REFERENCE_ABSENT,
        FIG2_MIXED_REFERENCE_ABSENT,
        FIG2_MIXED_REFERENCE_ABSENT,
    )
end

function _fig2_same_mixed_reference(left, right)
    return left.present == right.present &&
        left.dimension == right.dimension &&
        left.point == right.point &&
        left.candidate_id == right.candidate_id &&
        left.directory == right.directory &&
        left.state_sha256 == right.state_sha256
end

function _fig2_mixed_reference_fields(reference::Fig2MixedReference)
    return Dict{String,Any}(
        "mixed_reference_present" => reference.present,
        "mixed_reference_dimension" => reference.dimension,
        "mixed_reference_point" => reference.point,
        "mixed_reference_candidate_id" => reference.candidate_id,
        "mixed_reference_directory" => reference.directory,
        "mixed_reference_state_sha256" => reference.state_sha256,
    )
end

function _fig2_mixed_reference_from_data(data)
    present = get(data, "mixed_reference_present", nothing)
    present isa Bool || throw(ArgumentError(
        "candidate mixed reference present flag is missing or invalid"
    ))
    dimension = get(data, "mixed_reference_dimension", nothing)
    dimension isa Integer && !(dimension isa Bool) || throw(ArgumentError(
        "candidate mixed reference dimension is missing or invalid"
    ))
    point = get(data, "mixed_reference_point", nothing)
    point isa Integer && !(point isa Bool) || throw(ArgumentError(
        "candidate mixed reference point is missing or invalid"
    ))
    candidate_id = get(data, "mixed_reference_candidate_id", nothing)
    candidate_id isa AbstractString || throw(ArgumentError(
        "candidate mixed reference candidate_id is missing or invalid"
    ))
    directory = get(data, "mixed_reference_directory", nothing)
    directory isa AbstractString || throw(ArgumentError(
        "candidate mixed reference directory is missing or invalid"
    ))
    state_sha256 = get(data, "mixed_reference_state_sha256", nothing)
    state_sha256 isa AbstractString || throw(ArgumentError(
        "candidate mixed reference state_sha256 is missing or invalid"
    ))
    reference = Fig2MixedReference(
        present,
        Int(dimension),
        Int(point),
        String(candidate_id),
        String(directory),
        String(state_sha256),
    )
    if !reference.present
        _fig2_same_mixed_reference(
            reference, _fig2_absent_mixed_reference()
        ) || throw(ArgumentError(
            "candidate absent mixed reference does not use the strict sentinel"
        ))
        return reference
    end
    reference.dimension > 0 && reference.point > 0 || throw(ArgumentError(
        "candidate present mixed reference dimension and point must be positive"
    ))
    reference.directory == _candidate_relative_directory(
        reference.dimension, reference.point, reference.candidate_id
    ) || throw(ArgumentError(
        "candidate mixed reference directory disagrees with its identity"
    ))
    occursin(r"^[0-9a-f]{64}$", reference.state_sha256) || throw(
        ArgumentError("candidate mixed reference state_sha256 is invalid")
    )
    return reference
end

function _fig2_validate_mixed_reference_point(reference, dimension, point)
    if point == 1
        _fig2_same_mixed_reference(
            reference, _fig2_absent_mixed_reference()
        ) || throw(ArgumentError(
            "reference-flux candidate must use the strict absent mixed reference"
        ))
        return nothing
    end
    reference.present || throw(ArgumentError(
        "later candidate is missing its previous selected mixed reference"
    ))
    reference.dimension == dimension && reference.point == point - 1 || throw(
        ArgumentError(
            "candidate mixed reference dimension/point is not the prior flux point"
        )
    )
    return nothing
end

function _fig2_mixed_audit_from_data(data)
    return MixedTransferFidelity(
        Float64(data["fidelity_to_previous"]),
        ComplexF64(
            Float64(data["mixed_value_real"]),
            Float64(data["mixed_value_imag"]),
        ),
        ComplexF64(
            Float64(data["self_value1_real"]),
            Float64(data["self_value1_imag"]),
        ),
        ComplexF64(
            Float64(data["self_value2_real"]),
            Float64(data["self_value2_imag"]),
        ),
        Float64(data["mixed_residual"]),
        Float64(data["self_residual1"]),
        Float64(data["self_residual2"]),
        Bool(data["mixed_converged"]),
        Bool(data["mixed_valid"]),
        String(data["mixed_reason"]),
    )
end

function _validate_fig2_evidence(
    spec::Fig2BenchmarkSpec,
    evidence::Fig2CandidateEvidence,
    point,
    dimension,
)
    evidence.requested_maxdim > 0 || throw(
        ArgumentError("candidate requested maxdim must be positive")
    )
    evidence.requested_maxdim == dimension || throw(
        ArgumentError("candidate requested maxdim disagrees with its ledger dimension")
    )
    evidence.achieved_maxlinkdim > 0 || throw(
        ArgumentError("candidate achieved maxlinkdim must be positive")
    )
    evidence.achieved_maxlinkdim <= evidence.requested_maxdim || throw(
        ArgumentError("candidate achieved maxlinkdim exceeds its requested maxdim")
    )
    evidence.checkpoint_maxlinkdim >= 0 || throw(
        ArgumentError("candidate checkpoint maxlinkdim must be nonnegative")
    )
    if evidence.restart_valid
        evidence.checkpoint_maxlinkdim == evidence.achieved_maxlinkdim || throw(
            ArgumentError("restartable candidate checkpoint maxlinkdim disagrees with achieved maxlinkdim")
        )
        evidence.checkpoint_maxlinkdim > 0 || throw(
            ArgumentError("restartable candidate checkpoint maxlinkdim must be positive")
        )
    end
    isfinite(evidence.energy_per_site) || throw(
        ArgumentError("candidate energy per site must be finite")
    )
    isfinite(evidence.raw_schmidt_polarization) || throw(
        ArgumentError("candidate raw Schmidt polarization must be finite")
    )
    normalized = _normalized_sector_weights(evidence.sector_weights)
    _fig2_mixed_audit(evidence, point)
    if point > 1 && evidence.fidelity_valid
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
    counting = evidence.momentum_counting_evidence
    counting isa Fig2MomentumCounting || throw(
        ArgumentError("candidate is missing momentum-counting audit evidence")
    )
    recomputed = _fig2_momentum_counting(spec, evidence.momentum)
    _fig2_same_counting(counting, recomputed) || throw(
        ArgumentError("candidate momentum-counting evidence does not match its spectrum")
    )
    evidence.momentum_counting == counting.counts || throw(
        ArgumentError("candidate momentum-counting array disagrees with its audit evidence")
    )
    return normalized
end

function _render_fig2_momentum(momentum)
    header = FIG2_ARTIFACT_HEADERS["momentum_entanglement_spectrum.tsv"]
    validation = momentum.validation
    offsets = _fig2_sector_offsets_text(momentum.sector_ky_offsets)
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
            validation.unitary_residual,
            validation.power_residual,
            validation.commutator_residual,
            validation.transfer_residual,
            validation.valid,
            _fig2_tsv_text(validation.reason),
            momentum.valid,
            _fig2_tsv_text(momentum.reason),
            string(momentum.momentum_convention),
            offsets,
        ) for level in momentum.levels
    )
    return _render_tsv(header, rows)
end

function _render_fig2_mixed(
    fidelity::MixedTransferFidelity,
    reference::Fig2MixedReference,
)
    return _render_tsv(
        FIG2_ARTIFACT_HEADERS["mixed_fidelity.tsv"],
        [(
            fidelity.fidelity_cell,
            real(fidelity.mixed_value),
            imag(fidelity.mixed_value),
            real(fidelity.self_value1),
            imag(fidelity.self_value1),
            real(fidelity.self_value2),
            imag(fidelity.self_value2),
            fidelity.mixed_residual,
            fidelity.self_residual1,
            fidelity.self_residual2,
            fidelity.converged,
            fidelity.valid,
            _fig2_tsv_text(fidelity.reason),
            reference.present,
            reference.dimension,
            reference.point,
            reference.candidate_id,
            reference.directory,
            reference.state_sha256,
        )],
    )
end

function _validate_fig2_tsv_structure(path, filename)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("$filename is empty"))
    expected_header = FIG2_ARTIFACT_HEADERS[filename]
    first(lines) == expected_header || throw(
        ArgumentError("$filename header does not match the production format")
    )
    rows = lines[2:end]
    filename in FIG2_LEGALLY_EMPTY_TABLES || !isempty(rows) || throw(
        ArgumentError("$filename has no data rows")
    )
    expected_columns = length(split(expected_header, '\t'))
    for (row, line) in enumerate(rows)
        length(split(line, '\t'; keepempty=true)) == expected_columns || throw(
            ArgumentError("$filename row $row has the wrong column count")
        )
    end
    return nothing
end

function _fig2_artifact_integer(text, filename, column, row; positive=false)
    value = tryparse(Int, text)
    isnothing(value) && throw(ArgumentError(
        "$filename $column at row $row is not a strict integer"
    ))
    positive && value <= 0 && throw(ArgumentError(
        "$filename $column at row $row must be positive"
    ))
    return value
end

function _fig2_artifact_float(
    text,
    filename,
    column,
    row;
    nonnegative=false,
)
    value = tryparse(Float64, text)
    isnothing(value) && throw(ArgumentError(
        "$filename $column at row $row is not a floating-point number"
    ))
    isfinite(value) || throw(ArgumentError(
        "$filename $column at row $row must be finite"
    ))
    nonnegative && value < 0 && throw(ArgumentError(
        "$filename $column at row $row must be nonnegative"
    ))
    return value
end

function _fig2_artifact_boolean(text, filename, column, row)
    text == "true" && return true
    text == "false" && return false
    throw(ArgumentError(
        "$filename $column at row $row is not a strict boolean"
    ))
end

function _fig2_validate_sector_distribution(actual, expected, filename)
    tolerance = 1e-10
    isapprox(sum(values(actual)), 1.0; atol=tolerance, rtol=0) || throw(
        ArgumentError("$filename fixed cut/bond sector weights are not normalized")
    )
    Set(keys(actual)) == Set(keys(expected)) || throw(ArgumentError(
        "$filename fixed cut/bond raw-charge sectors disagree with candidate metadata"
    ))
    all(
        isapprox(actual[charge], expected[charge]; atol=tolerance, rtol=tolerance)
        for charge in keys(expected)
    ) || throw(ArgumentError(
        "$filename fixed cut/bond sector weights disagree with candidate metadata"
    ))
    return nothing
end

function _validate_fig2_sector_artifacts(spec, directory, metadata)
    expected = _fig2_argument_error(
        "candidate sector metadata is incomplete or invalid"
    ) do
        _normalized_sector_weights(Dict(
            parse(Int, charge) => Float64(weight) for
            (charge, weight) in metadata["sector_weights"]
        ))
    end
    scale = charge_scale(spec.config)
    expected_pair = (
        1,
        Base.Checked.checked_mul(1, spec.config.Ly),
    )

    sector_filename = "schmidt_sectors.tsv"
    sector_lines = readlines(joinpath(directory, sector_filename))[2:end]
    sector_rows = NamedTuple[]
    for (row, line) in enumerate(sector_lines)
        fields = split(line, '\t'; keepempty=true)
        cut_x = _fig2_artifact_integer(
            fields[1], sector_filename, "cut_x", row; positive=true
        )
        bond = _fig2_artifact_integer(
            fields[2], sector_filename, "bond", row; positive=true
        )
        (cut_x, bond) == expected_pair || continue
        isempty(fields[3]) && throw(ArgumentError(
            "$sector_filename qn at row $row is empty"
        ))
        raw_charge = _fig2_artifact_integer(
            fields[4], sector_filename, "raw_charge", row
        )
        physical_charge = _fig2_artifact_float(
            fields[5], sector_filename, "physical_charge", row
        )
        _fig2_same_float(physical_charge, raw_charge / scale) || throw(
            ArgumentError(
                "$sector_filename physical_charge disagrees with raw_charge/charge_scale at row $row"
            )
        )
        weight = _fig2_artifact_float(
            fields[6], sector_filename, "weight", row; nonnegative=true
        )
        valid = _fig2_artifact_boolean(
            fields[7], sector_filename, "valid", row
        )
        push!(sector_rows, (;
            cut_x, bond, raw_charge, physical_charge, weight, valid
        ))
    end

    spectrum_filename = "entanglement_spectrum.tsv"
    spectrum_lines = readlines(joinpath(directory, spectrum_filename))[2:end]
    spectrum_rows = NamedTuple[]
    for (row, line) in enumerate(spectrum_lines)
        fields = split(line, '\t'; keepempty=true)
        cut_x = _fig2_artifact_integer(
            fields[1], spectrum_filename, "cut_x", row; positive=true
        )
        bond = _fig2_artifact_integer(
            fields[2], spectrum_filename, "bond", row; positive=true
        )
        (cut_x, bond) == expected_pair || continue
        level = _fig2_artifact_integer(
            fields[3], spectrum_filename, "level", row; positive=true
        )
        singular_value = _fig2_artifact_float(
            fields[4], spectrum_filename, "singular_value", row;
            nonnegative=true,
        )
        probability = _fig2_artifact_float(
            fields[5], spectrum_filename, "probability", row;
            nonnegative=true,
        )
        entanglement_energy = tryparse(Float64, fields[6])
        isnothing(entanglement_energy) && throw(ArgumentError(
            "$spectrum_filename entanglement_energy at row $row " *
            "is not a floating-point number"
        ))
        if iszero(probability)
            entanglement_energy == Inf || throw(ArgumentError(
                "$spectrum_filename zero-probability level at row $row " *
                "must have infinite entanglement energy"
            ))
        else
            isfinite(entanglement_energy) || throw(ArgumentError(
                "$spectrum_filename entanglement_energy at row $row must be finite"
            ))
            isapprox(
                entanglement_energy,
                -log(probability);
                atol=1e-12,
                rtol=1e-10,
            ) || throw(ArgumentError(
                "$spectrum_filename entanglement_energy disagrees with " *
                "probability at row $row"
            ))
        end
        isempty(fields[7]) && throw(ArgumentError(
            "$spectrum_filename qn at row $row is empty"
        ))
        raw_charge = _fig2_artifact_integer(
            fields[8], spectrum_filename, "raw_charge", row
        )
        physical_charge = _fig2_artifact_float(
            fields[9], spectrum_filename, "physical_charge", row
        )
        _fig2_same_float(physical_charge, raw_charge / scale) || throw(
            ArgumentError(
                "$spectrum_filename physical_charge disagrees with raw_charge/charge_scale at row $row"
            )
        )
        valid = _fig2_artifact_boolean(
            fields[10], spectrum_filename, "valid", row
        )
        isapprox(
            probability,
            singular_value^2;
            atol=1e-12,
            rtol=1e-10,
        ) || throw(ArgumentError(
            "$spectrum_filename probability disagrees with singular_value^2 " *
            "at row $row"
        ))
        push!(
            spectrum_rows,
            (;
                cut_x,
                bond,
                level,
                raw_charge,
                physical_charge,
                probability,
                valid,
            ),
        )
    end

    isempty(sector_rows) && throw(ArgumentError(
        "schmidt_sectors.tsv is missing the fixed cut/bond $(expected_pair)"
    ))
    isempty(spectrum_rows) && throw(ArgumentError(
        "entanglement_spectrum.tsv is missing the fixed cut/bond $(expected_pair)"
    ))
    fixed_cut_valid = all(row -> row.valid, sector_rows) &&
        all(row -> row.valid, spectrum_rows)
    levels = [row.level for row in spectrum_rows]
    length(unique(levels)) == length(levels) || throw(ArgumentError(
        "entanglement_spectrum.tsv first cut/bond has duplicate levels"
    ))

    sector_weights = Dict{Int,Float64}()
    for row in sector_rows
        sector_weights[row.raw_charge] =
            get(sector_weights, row.raw_charge, 0.0) + row.weight
    end
    spectrum_weights = Dict{Int,Float64}()
    for row in spectrum_rows
        spectrum_weights[row.raw_charge] =
            get(spectrum_weights, row.raw_charge, 0.0) + row.probability
    end
    _fig2_validate_sector_distribution(
        sector_weights, expected, sector_filename
    )
    _fig2_validate_sector_distribution(
        spectrum_weights, expected, spectrum_filename
    )
    raw_schmidt_polarization = sum(
        row.probability * row.physical_charge for row in spectrum_rows
    )
    metadata_polarization = get(metadata, "raw_schmidt_polarization", nothing)
    metadata_polarization isa Real && !(metadata_polarization isa Bool) &&
        isfinite(Float64(metadata_polarization)) || throw(ArgumentError(
        "candidate metadata raw_schmidt_polarization must be finite"
    ))
    isapprox(
        raw_schmidt_polarization,
        Float64(metadata_polarization);
        atol=1e-12,
        rtol=1e-10,
    ) || throw(ArgumentError(
        "candidate raw Schmidt polarization disagrees with entanglement_spectrum.tsv"
    ))
    return (;
        sector_pair=expected_pair,
        fixed_cut_valid,
        sector_weights,
        spectrum_weights,
        raw_schmidt_polarization,
    )
end

function _fig2_convergence_integer(text, column, row)
    value = tryparse(Int, text)
    isnothing(value) && throw(
        ArgumentError("convergence.tsv $column at row $row is not a strict integer")
    )
    value > 0 || throw(
        ArgumentError("convergence.tsv $column at row $row must be positive")
    )
    return value
end

function _fig2_convergence_float(text, column, row)
    value = tryparse(Float64, text)
    isnothing(value) && throw(
        ArgumentError("convergence.tsv $column at row $row is not a floating-point number")
    )
    isfinite(value) || throw(
        ArgumentError("convergence.tsv $column at row $row must be finite")
    )
    return value
end


function _parse_fig2_convergence_lines(lines)
    isempty(lines) && throw(ArgumentError("convergence.tsv is empty"))
    header = FIG2_ARTIFACT_HEADERS["convergence.tsv"]
    first(lines) == header || throw(
        ArgumentError("convergence.tsv header does not match the production format")
    )
    rows = lines[2:end]
    isempty(rows) && throw(ArgumentError("convergence.tsv has no data rows"))
    columns = split(header, '\t')
    parsed = NamedTuple[]
    for (row, line) in enumerate(rows)
        fields = split(line, '\t'; keepempty=true)
        length(fields) == length(columns) || throw(
            ArgumentError("convergence.tsv row $row has the wrong column count")
        )
        stage = _fig2_convergence_integer(fields[1], columns[1], row)
        iteration = _fig2_convergence_integer(fields[2], columns[2], row)
        maxlinkdim = _fig2_convergence_integer(fields[3], columns[3], row)
        energy_left = _fig2_convergence_float(fields[4], columns[4], row)
        energy_right = _fig2_convergence_float(fields[5], columns[5], row)
        energy_mismatch = _fig2_convergence_float(fields[6], columns[6], row)
        delta_energy = if fields[7] == "missing"
            missing
        else
            _fig2_convergence_float(fields[7], columns[7], row)
        end
        eps_left = _fig2_convergence_float(fields[8], columns[8], row)
        eps_right = _fig2_convergence_float(fields[9], columns[9], row)
        precision_error = _fig2_convergence_float(fields[10], columns[10], row)
        elapsed_seconds = _fig2_convergence_float(fields[11], columns[11], row)
        elapsed_seconds < 0 && throw(ArgumentError(
            "convergence.tsv elapsed_seconds at row $row must be nonnegative"
        ))
        converged = if fields[12] == "true"
            true
        elseif fields[12] == "false"
            false
        else
            throw(ArgumentError(
                "convergence.tsv converged at row $row is not a strict boolean"
            ))
        end
        push!(parsed, (;
            stage,
            iteration,
            maxlinkdim,
            energy_left,
            energy_right,
            energy_mismatch,
            delta_energy,
            eps_left,
            eps_right,
            precision_error,
            elapsed_seconds,
            converged,
        ))
    end
    return parsed
end

function _parse_fig2_convergence_tsv(path)
    return _parse_fig2_convergence_lines(readlines(path))
end

function _normalize_fig2_convergence_sentinel!(path)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("convergence.tsv is empty"))
    header = FIG2_ARTIFACT_HEADERS["convergence.tsv"]
    first(lines) == header || throw(
        ArgumentError("convergence.tsv header does not match the production format")
    )
    length(lines) >= 2 || throw(ArgumentError("convergence.tsv has no data rows"))
    columns = split(header, '\t')
    first_fields = split(lines[2], '\t'; keepempty=true)
    length(first_fields) == length(columns) || throw(
        ArgumentError("convergence.tsv row 1 has the wrong column count")
    )
    first_fields[7] == "Inf" || throw(ArgumentError(
        "convergence.tsv producer sentinel must be literal Inf at row 1 delta_energy"
    ))

    normalized = copy(lines)
    for line_index in 2:length(normalized)
        fields = split(normalized[line_index], '\t'; keepempty=true)
        if length(fields) >= 7 && fields[7] == "Inf"
            fields[7] = "missing"
            normalized[line_index] = join(fields, '\t')
        end
    end
    _parse_fig2_convergence_lines(normalized)
    _write_fig2_bytes(path, join(normalized, '\n') * '\n')
    return String(path)
end

function _validate_fig2_convergence_tsv(path)
    return last(_parse_fig2_convergence_tsv(path))
end

function _fig2_convergence_optimization(optimization)
    values = Dict{Symbol,Any}()
    for name in (:vumps_tol, :energy_tol, :energy_mismatch_tol)
        hasproperty(optimization, name) || throw(ArgumentError(
            "candidate convergence optimization is missing $name"
        ))
        value = getproperty(optimization, name)
        value isa Real && !(value isa Bool) || throw(ArgumentError(
            "candidate convergence optimization $name must be a real number"
        ))
        converted = Float64(value)
        isfinite(converted) && converted > 0 || throw(ArgumentError(
            "candidate convergence optimization $name must be finite and positive"
        ))
        values[name] = converted
    end
    hasproperty(optimization, :stable_iterations) || throw(ArgumentError(
        "candidate convergence optimization is missing stable_iterations"
    ))
    stable_iterations = getproperty(optimization, :stable_iterations)
    stable_iterations isa Integer && !(stable_iterations isa Bool) &&
        stable_iterations > 0 || throw(ArgumentError(
            "candidate convergence optimization stable_iterations must be a positive integer"
        ))
    return (;
        vumps_tol=values[:vumps_tol],
        energy_tol=values[:energy_tol],
        energy_mismatch_tol=values[:energy_mismatch_tol],
        stable_iterations=Int(stable_iterations),
    )
end

function _validate_fig2_convergence_tsv(
    path,
    summary_schedule;
    optimization,
    energy_normalization_sites,
)
    summary_schedule isa AbstractVector && !isempty(summary_schedule) || throw(
        ArgumentError("candidate summary maxdim schedule must be a nonempty array")
    )
    all(
        value -> value isa Integer && !(value isa Bool) && value > 0,
        summary_schedule,
    ) || throw(ArgumentError(
        "candidate summary maxdim schedule must contain only positive integers"
    ))
    energy_normalization_sites isa Integer &&
        !(energy_normalization_sites isa Bool) &&
        energy_normalization_sites > 0 || throw(ArgumentError(
            "candidate convergence energy_normalization_sites must be a positive integer"
        ))
    contract = _fig2_convergence_optimization(optimization)
    schedule = Int.(summary_schedule)
    parsed = _parse_fig2_convergence_tsv(path)

    current_stage = 0
    expected_iteration = 1
    stable_count = 0
    stage_converged = false
    previous_energy = nothing
    previous_maxlinkdim = nothing
    for (row_number, row) in enumerate(parsed)
        if row.stage == current_stage + 1
            current_stage > 0 && !stage_converged && throw(ArgumentError(
                "convergence.tsv stage $current_stage did not end with a valid converged row before stage $(row.stage)"
            ))
            current_stage = row.stage
            expected_iteration = 1
            stable_count = 0
            stage_converged = false
            previous_energy = nothing
            previous_maxlinkdim = nothing
        elseif row.stage == current_stage
            stage_converged && throw(ArgumentError(
                "convergence.tsv stage $(row.stage) contains rows after its converged row"
            ))
        else
            throw(ArgumentError(
                "convergence.tsv stages must begin at 1 and increase contiguously at row $row_number"
            ))
        end
        row.stage <= length(schedule) || throw(ArgumentError(
            "convergence.tsv stage $(row.stage) at row $row_number exceeds the summary schedule"
        ))
        row.iteration == expected_iteration || throw(ArgumentError(
            "convergence.tsv iteration $(row.iteration) at row $row_number is not the expected consecutive iteration $expected_iteration for stage $(row.stage)"
        ))
        row.maxlinkdim <= schedule[row.stage] || throw(ArgumentError(
            "convergence.tsv maxlinkdim $(row.maxlinkdim) at row $row_number exceeds stage $(row.stage) cap $(schedule[row.stage])"
        ))
        row.iteration == 1 && !ismissing(row.delta_energy) && throw(ArgumentError(
            "convergence.tsv first delta_energy for stage $(row.stage) must be the normalized missing sentinel"
        ))
        for (name, value) in (
            "energy_mismatch" => row.energy_mismatch,
            "eps_left" => row.eps_left,
            "eps_right" => row.eps_right,
            "precision_error" => row.precision_error,
        )
            value >= 0 || throw(ArgumentError(
                "convergence.tsv $name at row $row_number must be nonnegative"
            ))
        end
        if !ismissing(row.delta_energy)
            row.delta_energy >= 0 || throw(ArgumentError(
                "convergence.tsv delta_energy at row $row_number must be nonnegative"
            ))
        end
        replayed_energy_mismatch = abs(row.energy_left - row.energy_right)
        _fig2_same_float(
            row.energy_mismatch, replayed_energy_mismatch
        ) || throw(ArgumentError(
            "convergence.tsv energy_mismatch at row $row_number disagrees with abs(energy_left-energy_right)"
        ))
        replayed_precision_error = max(row.eps_left, row.eps_right)
        _fig2_same_float(
            row.precision_error, replayed_precision_error
        ) || throw(ArgumentError(
            "convergence.tsv precision_error at row $row_number disagrees with max(eps_left,eps_right)"
        ))
        energy = (row.energy_left + row.energy_right) / 2
        if ismissing(row.delta_energy)
            if row.iteration != 1 && !(
                !isnothing(previous_maxlinkdim) &&
                    previous_maxlinkdim < schedule[row.stage]
            )
                throw(ArgumentError(
                    "convergence.tsv delta_energy at row $row_number is missing after stage $(row.stage) already reached its target maxlinkdim"
                ))
            end
            stable_count = 0
        else
            !isnothing(previous_energy) || throw(ArgumentError(
                "convergence.tsv first delta_energy for stage $(row.stage) must be the normalized missing sentinel"
            ))
            replayed_delta_energy = abs(energy - previous_energy)
            _fig2_same_float(
                row.delta_energy, replayed_delta_energy
            ) || throw(ArgumentError(
                "convergence.tsv delta_energy at row $row_number disagrees with adjacent mean raw energies"
            ))
        end
        stable_now = !ismissing(row.delta_energy) &&
            row.maxlinkdim == schedule[row.stage] &&
            all(isfinite, (
                row.precision_error,
                row.delta_energy,
                row.energy_mismatch,
            )) &&
            row.precision_error < contract.vumps_tol &&
            row.delta_energy / energy_normalization_sites < contract.energy_tol &&
            row.energy_mismatch / energy_normalization_sites <
                contract.energy_mismatch_tol
        stable_count = stable_now ? stable_count + 1 : 0
        expected_converged = stable_now &&
            stable_count >= contract.stable_iterations
        row.converged && row.maxlinkdim != schedule[row.stage] && throw(
            ArgumentError(
                "convergence.tsv final converged row $row_number for stage $(row.stage) has the wrong maxlinkdim: achieved=$(row.maxlinkdim), target=$(schedule[row.stage])"
            )
        )
        row.converged == expected_converged || throw(ArgumentError(
            "convergence.tsv converged flag at row $row_number disagrees with replayed production tolerances and stable-iteration count"
        ))
        stage_converged = row.converged
        expected_iteration += 1
        previous_energy = energy
        previous_maxlinkdim = row.maxlinkdim
    end

    observed_stages = unique(row.stage for row in parsed)
    expected_stages = collect(eachindex(schedule))
    observed_stages == expected_stages || throw(ArgumentError(
        "convergence.tsv stages $(collect(observed_stages)) do not exactly match summary stages $expected_stages"
    ))
    final = last(parsed)
    final.stage == last(expected_stages) || throw(ArgumentError(
        "convergence.tsv final row does not belong to the final summary stage"
    ))
    return final
end

function _fig2_candidate_positive_integer(data, key)
    value = get(data, key, nothing)
    value isa Integer && !(value isa Bool) || throw(
        ArgumentError("candidate metadata $key must be an integer")
    )
    converted = Int(value)
    converted > 0 || throw(
        ArgumentError("candidate metadata $key must be positive")
    )
    return converted
end

function _fig2_candidate_nonnegative_integer(data, key)
    value = get(data, key, nothing)
    value isa Integer && !(value isa Bool) || throw(
        ArgumentError("candidate metadata $key must be an integer")
    )
    converted = Int(value)
    converted >= 0 || throw(
        ArgumentError("candidate metadata $key must be nonnegative")
    )
    return converted
end

function _fig2_summary_maxdim_schedule(summary)
    optimization = get(summary, "optimization", nothing)
    optimization isa AbstractDict || throw(
        ArgumentError("candidate summary optimization table is missing")
    )
    schedule = get(optimization, "maxdim_schedule", nothing)
    schedule isa AbstractVector || throw(
        ArgumentError("candidate summary maxdim schedule must be an array")
    )
    all(value -> value isa Integer && !(value isa Bool), schedule) || throw(
        ArgumentError("candidate summary maxdim schedule must contain only integers")
    )
    return Int.(schedule)
end

function _fig2_summary_convergence_optimization(summary, snapshot)
    optimization = get(summary, "optimization", nothing)
    optimization isa AbstractDict || throw(
        ArgumentError("candidate summary optimization table is missing")
    )
    expected = snapshot["optimization"]
    multisite_update_alg = get(
        optimization, "multisite_update_alg", nothing
    )
    multisite_update_alg isa AbstractString || throw(ArgumentError(
        "candidate summary multisite_update_alg must be a string"
    ))
    String(multisite_update_alg) == String(expected["multisite_update_alg"]) ||
        throw(ArgumentError(
            "candidate summary multisite_update_alg disagrees with the immutable manifest"
        ))
    values = Dict{Symbol,Any}()
    for key in ("vumps_tol", "energy_tol", "energy_mismatch_tol")
        value = get(optimization, key, nothing)
        value isa Real && !(value isa Bool) || throw(ArgumentError(
            "candidate summary optimization $key must be a real number"
        ))
        converted = Float64(value)
        converted == Float64(expected[key]) || throw(ArgumentError(
            "candidate summary optimization $key disagrees with the immutable manifest"
        ))
        values[Symbol(key)] = converted
    end
    stable_iterations = get(optimization, "stable_iterations", nothing)
    stable_iterations isa Integer && !(stable_iterations isa Bool) || throw(
        ArgumentError(
            "candidate summary optimization stable_iterations must be an integer"
        )
    )
    Int(stable_iterations) == Int(expected["stable_iterations"]) || throw(
        ArgumentError(
            "candidate summary optimization stable_iterations disagrees with the immutable manifest"
        )
    )
    return (;
        multisite_update_alg=Symbol(multisite_update_alg),
        vumps_tol=values[:vumps_tol],
        energy_tol=values[:energy_tol],
        energy_mismatch_tol=values[:energy_mismatch_tol],
        stable_iterations=Int(stable_iterations),
    )
end

function _fig2_same_mixed_audit(left, right)
    return _fig2_same_float(left.fidelity_cell, right.fidelity_cell) &&
        _fig2_same_complex(left.mixed_value, right.mixed_value) &&
        _fig2_same_complex(left.self_value1, right.self_value1) &&
        _fig2_same_complex(left.self_value2, right.self_value2) &&
        _fig2_same_float(left.mixed_residual, right.mixed_residual) &&
        _fig2_same_float(left.self_residual1, right.self_residual1) &&
        _fig2_same_float(left.self_residual2, right.self_residual2) &&
        left.converged == right.converged &&
        left.valid == right.valid &&
        left.reason == right.reason
end

function _fig2_mixed_raw_float(text, column, point)
    value = tryparse(Float64, text)
    isnothing(value) && throw(ArgumentError(
        "mixed_fidelity.tsv $column at point $point is not a floating-point number"
    ))
    return value
end

function _validate_fig2_mixed_tsv(spec, path, metadata, point)
    _validate_fig2_tsv_structure(path, "mixed_fidelity.tsv")
    rows = readlines(path)[2:end]
    length(rows) == 1 || throw(
        ArgumentError("mixed_fidelity.tsv must contain exactly one evidence row")
    )
    fields = split(only(rows), '\t'; keepempty=true)
    values = Float64[
        _fig2_mixed_raw_float(fields[index], index, point) for index in 1:10
    ]
    converged = _fig2_artifact_boolean(
        fields[11], "mixed_fidelity.tsv", "converged", 1
    )
    valid = _fig2_artifact_boolean(
        fields[12], "mixed_fidelity.tsv", "valid", 1
    )
    audit = if point == 1
        all(isnan, values) || throw(ArgumentError(
            "mixed_fidelity.tsv reference point must use only NaN missing-value sentinels"
        ))
        _fig2_missing_mixed_fidelity()
    else
        snapshot = _fig2_validated_snapshot(spec)
        _mixed_transfer_result(
            ComplexF64(values[2], values[3]),
            ComplexF64(values[4], values[5]),
            ComplexF64(values[6], values[7]),
            values[8],
            values[9],
            values[10],
            converged;
            residual_tolerance=Float64(
                snapshot["optimization"]["transfer_tol"]
            ),
        )
    end
    _fig2_same_float(values[1], audit.fidelity_cell) || throw(ArgumentError(
        "mixed_fidelity.tsv fidelity disagrees with replayed raw transfer values"
    ))
    converged == audit.converged || throw(ArgumentError(
        "mixed_fidelity.tsv convergence disagrees with replayed raw transfer values"
    ))
    valid == audit.valid || throw(ArgumentError(
        "mixed_fidelity.tsv validity disagrees with replayed raw transfer values"
    ))
    fields[13] == _fig2_tsv_text(audit.reason) || throw(ArgumentError(
        "mixed_fidelity.tsv reason disagrees with replayed raw transfer values"
    ))
    metadata_audit = _fig2_mixed_audit_from_data(metadata)
    _fig2_same_mixed_audit(metadata_audit, audit) || throw(ArgumentError(
        "candidate mixed-fidelity metadata disagrees with replayed raw transfer values"
    ))
    tsv_reference = _fig2_mixed_reference_from_data(Dict{String,Any}(
        "mixed_reference_present" => _fig2_artifact_boolean(
            fields[14], "mixed_fidelity.tsv", "reference_present", 1
        ),
        "mixed_reference_dimension" => _fig2_artifact_integer(
            fields[15], "mixed_fidelity.tsv", "reference_dimension", 1
        ),
        "mixed_reference_point" => _fig2_artifact_integer(
            fields[16], "mixed_fidelity.tsv", "reference_point", 1
        ),
        "mixed_reference_candidate_id" => fields[17],
        "mixed_reference_directory" => fields[18],
        "mixed_reference_state_sha256" => fields[19],
    ))
    metadata_reference = _fig2_mixed_reference_from_data(metadata)
    _fig2_same_mixed_reference(tsv_reference, metadata_reference) || throw(
        ArgumentError(
            "mixed_fidelity.tsv reference disagrees with candidate metadata"
        )
    )
    _fig2_validate_mixed_reference_point(
        tsv_reference,
        Int(metadata["dimension"]),
        point,
    )
    return audit
end

function _fig2_momentum_residual(text, column)
    value = tryparse(Float64, text)
    isnothing(value) && throw(ArgumentError(
        "momentum_entanglement_spectrum.tsv $column is not a floating-point number"
    ))
    !isnan(value) && value >= 0 || throw(ArgumentError(
        "momentum_entanglement_spectrum.tsv $column must be nonnegative and not NaN"
    ))
    return value
end

function _fig2_momentum_offsets(text)
    isempty(text) && return Dict{Int,Int}()
    offsets = Dict{Int,Int}()
    for entry in split(text, ',')
        fields = split(entry, ':'; keepempty=true)
        length(fields) == 2 || throw(ArgumentError(
            "momentum_entanglement_spectrum.tsv sector_ky_offsets is malformed"
        ))
        charge = tryparse(Int, fields[1])
        offset = tryparse(Int, fields[2])
        (isnothing(charge) || isnothing(offset)) && throw(ArgumentError(
            "momentum_entanglement_spectrum.tsv sector_ky_offsets must contain strict integers"
        ))
        haskey(offsets, charge) && throw(ArgumentError(
            "momentum_entanglement_spectrum.tsv sector_ky_offsets contains a duplicate charge"
        ))
        offsets[charge] = offset
    end
    _fig2_sector_offsets_text(offsets) == text || throw(ArgumentError(
        "momentum_entanglement_spectrum.tsv sector_ky_offsets is not canonical"
    ))
    return offsets
end

function _validate_fig2_momentum_tsv(spec, path, metadata)
    filename = "momentum_entanglement_spectrum.tsv"
    snapshot = _fig2_validated_snapshot(spec)
    momentum_residual_tolerance = _fig2_number(
        snapshot["optimization"], "momentum_residual_tol"
    )
    transfer_tolerance = _fig2_number(
        snapshot["optimization"], "transfer_tol"
    )
    _validate_fig2_tsv_structure(path, filename)
    rows = readlines(path)[2:end]
    first_fields = split(first(rows), '\t'; keepempty=true)
    residuals = ntuple(index -> _fig2_momentum_residual(
        first_fields[index + 8], split(FIG2_ARTIFACT_HEADERS[filename], '\t')[index + 8]
    ), 4)
    validation_valid = _fig2_artifact_boolean(
        first_fields[13], filename, "validation_valid", 1
    )
    validation_reason = first_fields[14]
    momentum_valid = _fig2_artifact_boolean(
        first_fields[15], filename, "momentum_valid", 1
    )
    momentum_reason = first_fields[16]
    convention = Symbol(first_fields[17])
    offsets = _fig2_momentum_offsets(first_fields[18])
    validation = MomentumValidation(
        residuals..., validation_valid, validation_reason
    )
    if validation_valid
        all(
            residual -> isfinite(residual) &&
                residual <= momentum_residual_tolerance,
            residuals[1:3],
        ) || throw(ArgumentError(
            "momentum TSV validation residual exceeds the immutable manifest tolerance"
        ))
        isfinite(residuals[4]) && residuals[4] <= transfer_tolerance || throw(
            ArgumentError(
                "momentum TSV transfer residual exceeds the immutable manifest tolerance"
            )
        )
    end
    momentum_valid == validation_valid || throw(ArgumentError(
        "momentum TSV momentum and validation validity flags disagree"
    ))
    momentum_reason == validation_reason || throw(ArgumentError(
        "momentum TSV momentum and validation reasons disagree"
    ))
    if momentum_valid
        String(convention) == spec.data["momentum_convention"] || throw(
            ArgumentError("momentum TSV convention disagrees with the immutable manifest")
        )
    else
        convention == :invalid || throw(ArgumentError(
            "invalid momentum TSV must use the invalid convention"
        ))
        isempty(offsets) || throw(ArgumentError(
            "invalid momentum TSV must not report sector offsets"
        ))
    end

    levels = MomentumEntanglementLevel[]
    tail = first_fields[9:18]
    for (row, line) in enumerate(rows)
        fields = split(line, '\t'; keepempty=true)
        fields[9:18] == tail || throw(ArgumentError(
            "momentum TSV validation diagnostics are not identical at row $row"
        ))
        level = _fig2_artifact_integer(
            fields[1], filename, "level", row; positive=true
        )
        level == row || throw(ArgumentError(
            "momentum TSV levels must be unique, ordered, and consecutive"
        ))
        probability = _fig2_artifact_float(
            fields[2], filename, "probability", row; nonnegative=true
        )
        probability > 0 || throw(ArgumentError(
            "momentum TSV probability at row $row must be positive"
        ))
        entanglement_energy = _fig2_artifact_float(
            fields[3], filename, "entanglement_energy", row
        )
        isapprox(
            entanglement_energy,
            -log(probability);
            atol=1e-12,
            rtol=1e-10,
        ) || throw(ArgumentError(
            "momentum TSV entanglement energy disagrees with probability at row $row"
        ))
        raw_charge = _fig2_artifact_integer(
            fields[4], filename, "raw_charge", row
        )
        labels_missing = fields[5] == fields[6] == fields[7] == "missing"
        if labels_missing
            momentum_valid && throw(ArgumentError(
                "valid momentum TSV contains a missing ky label at row $row"
            ))
            fields[8] == "NaN" || throw(ArgumentError(
                "missing momentum labels require a NaN phase residual at row $row"
            ))
            ky_index = missing
            ky = missing
            eigenphase = missing
            phase_residual = NaN
        else
            momentum_valid || throw(ArgumentError(
                "invalid momentum TSV contains a fabricated ky label at row $row"
            ))
            any(field -> field == "missing", fields[5:7]) && throw(ArgumentError(
                "momentum TSV ky labels must be all present or all missing at row $row"
            ))
            ky_index = _fig2_artifact_integer(
                fields[5], filename, "ky_index", row
            )
            0 <= ky_index < spec.config.Ny || throw(ArgumentError(
                "momentum TSV ky_index at row $row is outside [0, Ny)"
            ))
            ky = _fig2_artifact_float(fields[6], filename, "ky", row)
            eigenphase = _fig2_artifact_float(
                fields[7], filename, "eigenphase", row
            )
            phase_residual = _fig2_artifact_float(
                fields[8], filename, "phase_residual", row; nonnegative=true
            )
            isapprox(
                ky,
                2pi * ky_index / spec.config.Ny;
                atol=1e-12,
                rtol=1e-10,
            ) || throw(ArgumentError(
                "momentum TSV ky disagrees with ky_index at row $row"
            ))
            haskey(offsets, raw_charge) || throw(ArgumentError(
                "momentum TSV is missing the sector offset for raw charge $raw_charge"
            ))
            unshifted_index = mod(
                ky_index - offsets[raw_charge], spec.config.Ny
            )
            expected_residual = abs(
                cis(eigenphase) - cis(2pi * unshifted_index / spec.config.Ny)
            )
            isapprox(
                phase_residual,
                expected_residual;
                atol=1e-12,
                rtol=1e-10,
            ) || throw(ArgumentError(
                "momentum TSV phase residual is inconsistent at row $row"
            ))
            phase_residual <= momentum_residual_tolerance || throw(ArgumentError(
                "momentum TSV phase residual exceeds the immutable manifest tolerance at row $row"
            ))
        end
        push!(levels, MomentumEntanglementLevel(
            level,
            probability,
            entanglement_energy,
            raw_charge,
            ky_index,
            ky,
            eigenphase,
            phase_residual,
        ))
    end
    isapprox(
        sum(level.probability for level in levels),
        1.0;
        atol=1e-10,
        rtol=1e-10,
    ) || throw(ArgumentError("momentum TSV probabilities are not normalized"))

    metadata_validation = _fig2_momentum_validation_from_data(metadata)
    raw_validation_matches =
        _fig2_same_float(
            validation.unitary_residual, metadata_validation.unitary_residual
        ) &&
        _fig2_same_float(
            validation.power_residual, metadata_validation.power_residual
        ) &&
        _fig2_same_float(
            validation.commutator_residual,
            metadata_validation.commutator_residual,
        ) &&
        _fig2_same_float(
            validation.transfer_residual, metadata_validation.transfer_residual
        ) &&
        validation.valid == metadata_validation.valid &&
        validation.reason == _fig2_tsv_text(metadata_validation.reason)
    raw_validation_matches || throw(ArgumentError(
        "momentum TSV validation diagnostics disagree with candidate metadata"
    ))
    momentum_valid == Bool(metadata["momentum_valid"]) &&
        momentum_reason == _fig2_tsv_text(metadata["momentum_reason"]) &&
        convention == Symbol(String(metadata["momentum_convention"])) &&
        offsets == _fig2_sector_offsets_from_data(metadata) || throw(
        ArgumentError("momentum TSV labels disagree with candidate metadata")
    )
    data = MomentumEntanglementData(
        levels,
        validation,
        momentum_valid,
        momentum_reason,
        convention,
        offsets,
    )
    counting = _fig2_momentum_counting(spec, data)
    metadata_counting = _fig2_counting_from_data(metadata)
    _fig2_same_counting(counting, metadata_counting) &&
        Int.(metadata["momentum_counting"]) == counting.counts || throw(
        ArgumentError(
            "candidate momentum counting disagrees with the replayed raw spectrum"
        )
    )
    return (; data, validation, counting)
end

function _validate_fig2_candidate_artifacts(
    spec,
    directory,
    dimension,
    point,
    phi_y,
    candidate_id,
    current_generation_provenance,
    ;
    progress_audit=_default_fig2_progress_audit,
)
    snapshot = _fig2_validated_snapshot(spec)
    state_path = joinpath(directory, "state.h5")
    isfile(state_path) && filesize(state_path) > 0 || throw(
        ArgumentError("candidate output is missing required file state.h5")
    )
    summary_path = joinpath(directory, "summary.toml")
    isfile(summary_path) || throw(
        ArgumentError("candidate output is missing required file summary.toml")
    )
    summary = TOML.parsefile(summary_path)
    get(summary, "format", "") == TEXT_OUTPUT_FORMAT || throw(
        ArgumentError("candidate summary format is invalid")
    )
    get(summary, "algorithm", "") == "VUMPS" || throw(
        ArgumentError("candidate summary algorithm is invalid")
    )
    get(summary, "valid", nothing) isa Bool || throw(
        ArgumentError("candidate summary valid flag is missing")
    )
    get(summary, "converged", nothing) isa Bool || throw(
        ArgumentError("candidate summary converged flag is missing")
    )
    summary_schedule = _fig2_summary_maxdim_schedule(summary)
    full_schedule = _fig2_maxdim_schedule(dimension)
    suffix_start = length(full_schedule) - length(summary_schedule) + 1
    is_bounded_suffix = !isempty(summary_schedule) && suffix_start >= 1 &&
        summary_schedule == full_schedule[suffix_start:end]
    is_bounded_suffix || throw(
        ArgumentError("candidate summary maxdim schedule disagrees with requested maxdim")
    )
    convergence_optimization = _fig2_summary_convergence_optimization(
        summary, snapshot
    )
    summary_energy = get(summary, "energy", nothing)
    summary_energy isa AbstractDict || throw(
        ArgumentError("candidate summary energy table is missing")
    )
    summary_energy_per_site = get(summary_energy, "per_site", nothing)
    summary_energy_per_site isa Real &&
        !(summary_energy_per_site isa Bool) &&
        isfinite(Float64(summary_energy_per_site)) || throw(ArgumentError(
        "candidate summary energy per_site must be finite"
    ))
    summary_energy_per_site = Float64(summary_energy_per_site)

    for filename in keys(FIG2_ARTIFACT_HEADERS)
        path = joinpath(directory, filename)
        isfile(path) || throw(
            ArgumentError("candidate output is missing required file $filename")
        )
        _validate_fig2_tsv_structure(path, filename)
    end
    final_convergence = _validate_fig2_convergence_tsv(
        joinpath(directory, "convergence.tsv"),
        summary_schedule;
        optimization=convergence_optimization,
        energy_normalization_sites=sites_per_cell(spec.config),
    )
    replayed_energy_per_site = (
        final_convergence.energy_left + final_convergence.energy_right
    ) / (2 * sites_per_cell(spec.config))
    energy_tolerance = Float64(snapshot["optimization"]["energy_tol"])
    isapprox(
        summary_energy_per_site,
        replayed_energy_per_site;
        atol=energy_tolerance,
        rtol=8eps(Float64),
    ) || throw(ArgumentError(
        "candidate summary energy per_site disagrees with final convergence.tsv raw energy beyond the immutable tolerance"
    ))

    metadata_path = joinpath(directory, "candidate.toml")
    isfile(metadata_path) || throw(
        ArgumentError("candidate output is missing required file candidate.toml")
    )
    metadata = TOML.parsefile(metadata_path)
    get(metadata, "format", "") == FIG2_CANDIDATE_FORMAT || throw(
        ArgumentError("candidate metadata format is invalid")
    )
    get(metadata, "manifest_sha256", "") == spec.sha256 || throw(
        ArgumentError("candidate metadata manifest checksum is invalid")
    )
    identity = _fig2_argument_error("candidate identity fields are missing") do
        (
            Int(metadata["dimension"]),
            Int(metadata["point"]),
            Float64(metadata["phi_y"]),
            String(metadata["candidate_id"]),
        )
    end
    identity == (dimension, point, phi_y, String(candidate_id)) || throw(
        ArgumentError("candidate identity does not match its ledger row")
    )
    requested_maxdim = _fig2_candidate_positive_integer(
        metadata, "requested_maxdim"
    )
    achieved_maxlinkdim = _fig2_candidate_positive_integer(
        metadata, "achieved_maxlinkdim"
    )
    checkpoint_maxlinkdim = _fig2_candidate_nonnegative_integer(
        metadata, "checkpoint_maxlinkdim"
    )
    requested_maxdim == dimension || throw(
        ArgumentError("candidate metadata requested maxdim disagrees with its ledger dimension")
    )
    achieved_maxlinkdim <= requested_maxdim || throw(
        ArgumentError("candidate metadata achieved maxlinkdim exceeds requested maxdim")
    )
    progress = _fig2_candidate_progress_audit(
        progress_audit,
        spec,
        directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        state_path,
        achieved_maxlinkdim,
    )
    get(metadata, "progress_complete", nothing) === progress.complete || throw(
        ArgumentError("candidate metadata progress completion flag disagrees")
    )
    for (key, measured) in (
        "progress_event_count" => progress.event_count,
        "progress_resume_count" => progress.resume_count,
        "progress_latest_maxlinkdim" => progress.latest_maxlinkdim,
    )
        get(metadata, key, nothing) == measured || throw(
            ArgumentError("candidate metadata $key disagrees with progress audit")
        )
    end
    for (key, measured) in (
        "progress_final_state_sha256" => progress.final_state_sha256,
        "progress_sha256" => progress.progress_sha256,
    )
        get(metadata, key, nothing) == measured || throw(
            ArgumentError("candidate metadata $key disagrees with progress audit")
        )
    end
    sector_audit = _validate_fig2_sector_artifacts(
        spec, directory, metadata
    )
    final_convergence.maxlinkdim == achieved_maxlinkdim || throw(
        ArgumentError("final convergence.tsv maxlinkdim disagrees with achieved maxlinkdim")
    )
    restart_valid = get(metadata, "restart_valid", nothing)
    restart_valid isa Bool || throw(
        ArgumentError("candidate metadata restart flag is missing or invalid")
    )
    if restart_valid
        checkpoint_maxlinkdim == achieved_maxlinkdim || throw(ArgumentError(
            "candidate checkpoint maxlinkdim disagrees with achieved maxlinkdim"
        ))
        checkpoint_maxlinkdim > 0 || throw(ArgumentError(
            "restartable candidate checkpoint maxlinkdim must be positive"
        ))
    end
    generation = get(metadata, "generation_provenance", nothing)
    generation isa AbstractDict || throw(
        ArgumentError("candidate metadata is missing generation provenance")
    )
    _fig2_validate_generation_compatibility(
        generation, current_generation_provenance
    )
    get(metadata, "provenance_valid", false) === true || throw(
        ArgumentError("candidate metadata provenance is invalid")
    )
    metadata_energy_per_site = get(metadata, "energy_per_site", nothing)
    metadata_energy_per_site isa Real &&
        !(metadata_energy_per_site isa Bool) &&
        isfinite(Float64(metadata_energy_per_site)) || throw(ArgumentError(
        "candidate metadata energy_per_site must be finite"
    ))
    _fig2_same_float(summary_energy_per_site, metadata_energy_per_site) || throw(
        ArgumentError(
            "candidate summary energy per_site disagrees with candidate metadata"
        )
    )
    mixed_audit, momentum_audit = _fig2_argument_error(
        "candidate audit metadata is incomplete or invalid";
        passthrough_argument_error=true,
    ) do
        (
            _validate_fig2_mixed_tsv(
                spec,
                joinpath(directory, "mixed_fidelity.tsv"),
                metadata,
                point,
            ),
            _validate_fig2_momentum_tsv(
                spec,
                joinpath(directory, "momentum_entanglement_spectrum.tsv"),
                metadata,
            ),
        )
    end
    momentum_validation = momentum_audit.validation
    counting = momentum_audit.counting
    Bool(metadata["mixed_present"]) == (point > 1) || throw(
        ArgumentError("candidate mixed-evidence presence does not match its flux point")
    )
    Bool(metadata["fidelity_valid"]) == mixed_audit.valid || throw(
        ArgumentError("candidate mixed validity fields disagree")
    )
    Bool(metadata["momentum_validation_valid"]) == momentum_validation.valid || throw(
        ArgumentError("candidate momentum validation fields disagree")
    )
    Int.(metadata["momentum_counting"]) == counting.counts || throw(
        ArgumentError("candidate momentum counting fields disagree")
    )
    (!counting.valid || Bool(metadata["momentum_valid"])) || throw(
        ArgumentError("valid momentum counting requires valid momentum diagnostics")
    )
    get(metadata, "converged", nothing) === summary["converged"] || throw(
        ArgumentError("candidate summary and metadata convergence flags disagree")
    )
    summary["converged"] == final_convergence.converged || throw(
        ArgumentError("candidate summary and final convergence.tsv flags disagree")
    )
    metadata_valid = get(metadata, "valid", nothing)
    metadata_valid isa Bool || throw(
        ArgumentError("candidate metadata selection-core validity is missing")
    )
    replayed_selection_core_validity =
        final_convergence.converged &&
        isfinite(replayed_energy_per_site) &&
        sector_audit.sector_pair == (1, spec.config.Ly) &&
        sector_audit.fixed_cut_valid
    metadata_valid == replayed_selection_core_validity || throw(
        ArgumentError(
            "candidate metadata selection-core validity disagrees with " *
            "replayed raw convergence and fixed-cut artifacts"
        )
    )
    return (;
        metadata,
        replayed_energy_per_site,
        final_convergence,
        progress,
    )
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
    mixed_reference,
    relative_directory,
    generation_provenance,
    ;
    progress_audit=_default_fig2_progress_audit,
)
    directory = joinpath(output, relative_directory)
    for filename in FIG2_REQUIRED_CANDIDATE_FILES[1:7]
        path = joinpath(directory, filename)
        isfile(path) && filesize(path) > 0 || throw(
            ArgumentError("candidate output is missing required file $filename")
        )
    end
    normalized_sectors = _validate_fig2_evidence(
        spec, evidence, point, dimension
    )
    mixed_fidelity = _fig2_mixed_audit(evidence, point)
    mixed_reference isa Fig2MixedReference || throw(ArgumentError(
        "candidate mixed reference must be a Fig2MixedReference"
    ))
    _fig2_validate_mixed_reference_point(
        mixed_reference, dimension, point
    )
    progress = _fig2_candidate_progress_audit(
        progress_audit,
        spec,
        directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        joinpath(directory, "state.h5"),
        evidence.achieved_maxlinkdim,
    )
    _write_fig2_bytes(
        joinpath(directory, "momentum_entanglement_spectrum.tsv"),
        _render_fig2_momentum(evidence.momentum),
    )
    _write_fig2_bytes(
        joinpath(directory, "mixed_fidelity.tsv"),
        _render_fig2_mixed(mixed_fidelity, mixed_reference),
    )
    candidate_metadata = Dict{String,Any}(
        "format" => FIG2_CANDIDATE_FORMAT,
        "dimension" => dimension,
        "requested_maxdim" => evidence.requested_maxdim,
        "achieved_maxlinkdim" => evidence.achieved_maxlinkdim,
        "checkpoint_maxlinkdim" => evidence.checkpoint_maxlinkdim,
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
        "restart_valid" => evidence.restart_valid,
        "progress_complete" => progress.complete,
        "progress_event_count" => progress.event_count,
        "progress_resume_count" => progress.resume_count,
        "progress_latest_maxlinkdim" => progress.latest_maxlinkdim,
        "progress_final_state_sha256" => progress.final_state_sha256,
        "progress_sha256" => progress.progress_sha256,
        "provenance_valid" => true,
        "generation_provenance" => deepcopy(generation_provenance),
        "reason" => evidence.reason,
        "manifest_sha256" => spec.sha256,
    )
    merge!(
        candidate_metadata,
        _fig2_mixed_fields(mixed_fidelity; present=!isnothing(evidence.mixed_fidelity)),
        _fig2_mixed_reference_fields(mixed_reference),
        _fig2_momentum_fields(evidence.momentum),
        _fig2_counting_fields(evidence.momentum_counting_evidence),
    )
    _write_fig2_toml(joinpath(directory, "candidate.toml"), candidate_metadata)
    artifact_audit = _validate_fig2_candidate_artifacts(
        spec,
        directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        generation_provenance,
        ;
        progress_audit,
    )

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
        "requested_maxdim" => evidence.requested_maxdim,
        "achieved_maxlinkdim" => evidence.achieved_maxlinkdim,
        "checkpoint_maxlinkdim" => evidence.checkpoint_maxlinkdim,
        "point" => point,
        "phi_y" => phi_y,
        "candidate_id" => String(candidate_id),
        "directory" => relative_directory,
        "state_sha256" => checksums["state.h5"],
        "progress_complete" => artifact_audit.progress.complete,
        "progress_event_count" => artifact_audit.progress.event_count,
        "progress_resume_count" => artifact_audit.progress.resume_count,
        "progress_latest_maxlinkdim" =>
            artifact_audit.progress.latest_maxlinkdim,
        "progress_final_state_sha256" =>
            artifact_audit.progress.final_state_sha256,
        "progress_sha256" => checksums["progress.toml"],
        "generation_provenance_sha256" =>
            _fig2_generation_provenance_sha256(generation_provenance),
        "checksums" => checksums,
    )
    push!(ledger["candidate"], row)
    _write_fig2_toml(joinpath(output, "ledger.toml"), ledger)
    return (
        sectors=normalized_sectors,
        provenance_valid=true,
        authoritative_energy_per_site=
            artifact_audit.replayed_energy_per_site,
    )
end

function _fig2_authoritative_candidate_energy(candidate)
    hasproperty(candidate, :authoritative_energy_per_site) || throw(
        ArgumentError(
            "Fig. 2 candidate is missing authoritative convergence energy"
        )
    )
    energy = candidate.authoritative_energy_per_site
    energy isa Real && !(energy isa Bool) && isfinite(energy) || throw(
        ArgumentError(
            "Fig. 2 candidate authoritative convergence energy must be finite"
        )
    )
    return Float64(energy)
end

function _fig2_candidate_selection_eligible(candidate)
    return candidate.evidence.converged && candidate.evidence.valid &&
        candidate.evidence.restart_valid && candidate.provenance_valid &&
        candidate.evidence.achieved_maxlinkdim ==
            candidate.evidence.requested_maxdim &&
        candidate.evidence.checkpoint_maxlinkdim ==
            candidate.evidence.achieved_maxlinkdim &&
        isfinite(_fig2_authoritative_candidate_energy(candidate))
end

function _select_fig2_candidate(
    candidates,
    point;
    previous_sector_weights=nothing,
    sector_tolerance::Real=0.0,
)
    isfinite(sector_tolerance) && sector_tolerance >= 0 || throw(
        ArgumentError("Fig. 2 sector selection tolerance must be finite and nonnegative")
    )
    eligible = findall(_fig2_candidate_selection_eligible, candidates)
    isempty(eligible) && throw(WorkflowValidationError(
        "no converged valid restartable Fig. 2 candidate is selectable"
    ))
    if point == 1
        index = first(sort(eligible; by=i -> (
            _fig2_authoritative_candidate_energy(candidates[i]),
            candidates[i].candidate_id,
        )))
        return index, "minimum_zero_flux_energy"
    end
    isnothing(previous_sector_weights) && throw(
        ArgumentError("later Fig. 2 branch selection needs previous sector weights")
    )
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
        best_sector_shift(
            previous_sector_weights, candidates[i].sectors
        ).distance,
        _fig2_authoritative_candidate_energy(candidates[i]),
        candidates[i].candidate_id,
    )))
    return index, "maximum_mixed_fidelity"
end

function _fig2_expected_candidate_ids(
    candidate_ids_provider,
    spec,
    dimension,
    point,
)
    candidate_ids = String.(candidate_ids_provider(
        spec, dimension, point, nothing
    ))
    isempty(candidate_ids) && throw(
        WorkflowValidationError("candidate list is empty at point $point")
    )
    length(unique(candidate_ids)) == length(candidate_ids) || throw(
        ArgumentError("candidate IDs must be unique at each point")
    )
    return candidate_ids
end

function _validate_fig2_candidate_id_set(
    candidate_ids_provider,
    spec,
    dimension,
    point,
    actual_candidate_ids;
    require_complete,
)
    expected_ids = _fig2_expected_candidate_ids(
        candidate_ids_provider, spec, dimension, point
    )
    expected = Set(expected_ids)
    actual = Set(String(candidate_id) for candidate_id in actual_candidate_ids)
    valid = require_complete ? actual == expected : issubset(actual, expected)
    valid || throw(ArgumentError(
        "candidate ID set at dimension $dimension point $point does not " *
        (require_complete ? "exactly match" : "belong to") *
        " the immutable generation contract: actual=$(sort!(collect(actual))) " *
        "expected=$(sort!(collect(expected)))"
    ))
    return expected_ids
end

function _selection_dict(selection::Fig2Selection)
    row = Dict{String,Any}(
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
        "restart_valid" => selection.restart_valid,
        "provenance_valid" => selection.provenance_valid,
    )
    momentum = (;
        validation=selection.momentum_validation,
        valid=selection.momentum_valid,
        reason=selection.momentum_reason,
        momentum_convention=selection.momentum_convention,
        sector_ky_offsets=selection.momentum_sector_ky_offsets,
    )
    merge!(
        row,
        _fig2_momentum_fields(momentum),
        _fig2_counting_fields(selection.momentum_counting_evidence),
    )
    return row
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
        _fig2_momentum_validation_from_data(row),
        Bool(row["momentum_valid"]),
        String(row["momentum_reason"]),
        Symbol(String(row["momentum_convention"])),
        _fig2_sector_offsets_from_data(row),
        Int.(row["momentum_counting"]),
        _fig2_counting_from_data(row),
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

function _fig2_replayed_candidate_energy_per_site(spec, directory)
    final = last(_parse_fig2_convergence_tsv(joinpath(
        directory, "convergence.tsv"
    )))
    return (final.energy_left + final.energy_right) /
        (2 * sites_per_cell(spec.config))
end

function _fig2_mixed_reference_from_candidate_row(root, row)
    dimension, point, candidate_id = _fig2_candidate_key(row)
    directory = get(row, "directory", nothing)
    directory isa AbstractString || throw(ArgumentError(
        "mixed reference selected candidate directory is missing"
    ))
    expected_directory = _candidate_relative_directory(
        dimension, point, candidate_id
    )
    String(directory) == expected_directory || throw(ArgumentError(
        "mixed reference selected candidate directory disagrees with its identity"
    ))
    state_sha256 = get(row, "state_sha256", nothing)
    state_sha256 isa AbstractString &&
        occursin(r"^[0-9a-f]{64}$", state_sha256) || throw(ArgumentError(
        "mixed reference selected candidate state_sha256 is invalid"
    ))
    checksums = get(row, "checksums", nothing)
    checksums isa AbstractDict &&
        get(checksums, "state.h5", nothing) == state_sha256 || throw(
        ArgumentError(
            "mixed reference selected candidate state checksum fields disagree"
        )
    )
    state_path = joinpath(root, expected_directory, "state.h5")
    isfile(state_path) && filesize(state_path) > 0 || throw(ArgumentError(
        "mixed reference selected checkpoint is missing"
    ))
    _fig2_file_sha256(state_path) == state_sha256 || throw(ArgumentError(
        "mixed reference selected checkpoint checksum mismatch"
    ))
    return Fig2MixedReference(
        true,
        dimension,
        point,
        candidate_id,
        expected_directory,
        String(state_sha256),
    )
end

function _fig2_mixed_reference_for_selection(root, selection, candidate_row)
    expected = _fig2_mixed_reference_from_candidate_row(root, candidate_row)
    identity = if selection isa Fig2Selection
        (
            selection.dimension,
            selection.point,
            selection.candidate_id,
            selection.directory,
        )
    else
        (
            Int(selection["dimension"]),
            Int(selection["point"]),
            String(selection["candidate_id"]),
            String(selection["directory"]),
        )
    end
    identity == (
        expected.dimension,
        expected.point,
        expected.candidate_id,
        expected.directory,
    ) || throw(ArgumentError(
        "mixed reference does not identify the prior selected candidate"
    ))
    return expected
end

function _validate_fig2_mixed_reference_chain(
    root,
    candidate_rows,
    selection_rows,
)
    for row in values(candidate_rows)
        dimension, point, _ = _fig2_candidate_key(row)
        metadata = TOML.parsefile(joinpath(
            root, String(row["directory"]), "candidate.toml"
        ))
        actual = _fig2_mixed_reference_from_data(metadata)
        _fig2_validate_mixed_reference_point(actual, dimension, point)
        point == 1 && continue
        previous_key = (dimension, point - 1)
        haskey(selection_rows, previous_key) || throw(ArgumentError(
            "mixed reference has no prior selected ledger row"
        ))
        selected = selection_rows[previous_key]
        selected_key = (
            dimension,
            point - 1,
            String(selected["candidate_id"]),
        )
        haskey(candidate_rows, selected_key) || throw(ArgumentError(
            "mixed reference prior selected candidate row is missing"
        ))
        expected = _fig2_mixed_reference_for_selection(
            root, selected, candidate_rows[selected_key]
        )
        _fig2_same_mixed_reference(actual, expected) || throw(ArgumentError(
            "candidate mixed reference does not match the exact prior selected checkpoint"
        ))
    end
    return nothing
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

function _validate_persisted_fig2_selection(spec, root, row, candidate_row)
    selection = _fig2_argument_error(
        "persisted selection audit fields are incomplete"
    ) do
        _selection_from_dict(row)
    end
    candidate_identity = (
        Int(candidate_row["dimension"]),
        Int(candidate_row["point"]),
        Float64(candidate_row["phi_y"]),
        String(candidate_row["candidate_id"]),
    )
    selection_identity = (
        selection.dimension,
        selection.point,
        selection.phi_y,
        selection.candidate_id,
    )
    selection_identity == candidate_identity || throw(
        ArgumentError("persisted selection identity disagrees with its candidate")
    )
    selection.directory == String(candidate_row["directory"]) || throw(
        ArgumentError("persisted selection directory disagrees with its candidate")
    )
    metadata = TOML.parsefile(joinpath(root, selection.directory, "candidate.toml"))
    authoritative_energy_per_site =
        _fig2_replayed_candidate_energy_per_site(
            spec, joinpath(root, selection.directory)
        )
    scalar_matches =
        _fig2_same_float(
            selection.energy_per_site, authoritative_energy_per_site
        ) &&
        _fig2_same_float(
            selection.raw_schmidt_polarization,
            metadata["raw_schmidt_polarization"],
        ) &&
        _fig2_same_float(
            selection.fidelity_to_previous,
            metadata["fidelity_to_previous"],
        ) &&
        selection.fidelity_valid == Bool(metadata["fidelity_valid"]) &&
        selection.restart_valid == Bool(metadata["restart_valid"]) &&
        selection.provenance_valid == Bool(metadata["provenance_valid"])
    scalar_matches || throw(
        ArgumentError("persisted selection scalar evidence disagrees with candidate metadata")
    )
    metadata_sectors = _normalized_sector_weights(Dict(
        parse(Int, charge) => Float64(weight) for
        (charge, weight) in metadata["sector_weights"]
    ))
    selection.sector_weights == metadata_sectors || throw(
        ArgumentError("persisted selection sectors disagree with candidate metadata")
    )

    validation = _fig2_momentum_validation_from_data(metadata)
    momentum_matches =
        _fig2_same_float(
            selection.momentum_validation.unitary_residual,
            validation.unitary_residual,
        ) &&
        _fig2_same_float(
            selection.momentum_validation.power_residual,
            validation.power_residual,
        ) &&
        _fig2_same_float(
            selection.momentum_validation.commutator_residual,
            validation.commutator_residual,
        ) &&
        _fig2_same_float(
            selection.momentum_validation.transfer_residual,
            validation.transfer_residual,
        ) &&
        selection.momentum_validation.valid == validation.valid &&
        selection.momentum_validation.reason == validation.reason &&
        selection.momentum_valid == Bool(metadata["momentum_valid"]) &&
        selection.momentum_reason == String(metadata["momentum_reason"]) &&
        selection.momentum_convention ==
            Symbol(String(metadata["momentum_convention"])) &&
        selection.momentum_sector_ky_offsets ==
            _fig2_sector_offsets_from_data(metadata)
    momentum_matches || throw(
        ArgumentError("persisted selection momentum evidence disagrees with candidate metadata")
    )
    counting = _fig2_counting_from_data(metadata)
    selection.momentum_counting == counting.counts &&
        _fig2_same_counting(selection.momentum_counting_evidence, counting) || throw(
        ArgumentError("persisted selection momentum counting disagrees with candidate metadata")
    )
    return selection
end

function _validate_persisted_fig2_candidate_files(
    spec,
    root,
    row,
    current_generation_provenance,
    ;
    checkpoint_audit=_default_fig2_persisted_checkpoint_audit,
    progress_audit=_default_fig2_progress_audit,
)
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
    artifact_audit = _validate_fig2_candidate_artifacts(
        spec,
        directory,
        dimension,
        point,
        Float64(row["phi_y"]),
        candidate_id,
        current_generation_provenance,
        ;
        progress_audit,
    )
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
    get(row, "progress_sha256", "") == String(checksums["progress.toml"]) || throw(
        ArgumentError("persisted candidate progress checksum is inconsistent")
    )
    metadata = TOML.parsefile(joinpath(directory, "candidate.toml"))
    progress = artifact_audit.progress
    for (key, measured) in (
        "progress_complete" => progress.complete,
        "progress_event_count" => progress.event_count,
        "progress_resume_count" => progress.resume_count,
        "progress_latest_maxlinkdim" => progress.latest_maxlinkdim,
        "progress_final_state_sha256" => progress.final_state_sha256,
        "progress_sha256" => progress.progress_sha256,
    )
        get(row, key, nothing) == measured || throw(
            ArgumentError("persisted candidate ledger $key disagrees with progress audit")
        )
        get(metadata, key, nothing) == measured || throw(
            ArgumentError("persisted candidate metadata $key disagrees with progress audit")
        )
    end
    checkpoint_path = joinpath(directory, "state.h5")
    checkpoint_result = checkpoint_audit(
        spec, checkpoint_path, Float64(row["phi_y"])
    )
    checkpoint_result isa NamedTuple &&
        hasproperty(checkpoint_result, :restart_valid) &&
        hasproperty(checkpoint_result, :checkpoint_maxlinkdim) || throw(
        ArgumentError("persisted candidate checkpoint audit result is invalid")
    )
    checkpoint_result.restart_valid isa Bool || throw(
        ArgumentError("persisted candidate checkpoint restart audit is invalid")
    )
    checkpoint_result.checkpoint_maxlinkdim isa Integer &&
        !(checkpoint_result.checkpoint_maxlinkdim isa Bool) &&
        checkpoint_result.checkpoint_maxlinkdim >= 0 || throw(
        ArgumentError("persisted candidate checkpoint maxlinkdim is invalid")
    )
    measured_checkpoint_maxlinkdim =
        Int(checkpoint_result.checkpoint_maxlinkdim)
    checkpoint_result.restart_valid ==
        (measured_checkpoint_maxlinkdim > 0) || throw(
        ArgumentError(
            "persisted candidate checkpoint reload and maxlinkdim audits disagree",
        )
    )
    for key in (
        "requested_maxdim",
        "achieved_maxlinkdim",
        "checkpoint_maxlinkdim",
    )
        value = get(row, key, nothing)
        value isa Integer && !(value isa Bool) || throw(
            ArgumentError("persisted candidate ledger $key must be an integer")
        )
        Int(value) == Int(metadata[key]) || throw(
            ArgumentError("persisted candidate ledger $key disagrees with candidate metadata")
        )
    end
    stored_checkpoint_maxlinkdim = Int(metadata["checkpoint_maxlinkdim"])
    measured_checkpoint_maxlinkdim == stored_checkpoint_maxlinkdim || throw(
        ArgumentError(
            "persisted candidate checkpoint maxlinkdim disagrees with the independently reloaded state",
        )
    )
    achieved_maxlinkdim = Int(metadata["achieved_maxlinkdim"])
    measured_restart_valid = checkpoint_result.restart_valid &&
        measured_checkpoint_maxlinkdim == achieved_maxlinkdim
    Bool(metadata["restart_valid"]) == measured_restart_valid || throw(
        ArgumentError(
            "persisted candidate restart flag disagrees with the independently reloaded checkpoint",
        )
    )
    generation = metadata["generation_provenance"]
    get(row, "generation_provenance_sha256", "") ==
        _fig2_generation_provenance_sha256(generation) || throw(
        ArgumentError("persisted candidate generation provenance checksum mismatch")
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
    get(metadata, "format", "") == FIG2_CANDIDATE_FORMAT || throw(
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
    authoritative_energy_per_site =
        _fig2_replayed_candidate_energy_per_site(spec, directory)
    evidence = (
        state=nothing,
        requested_maxdim=Int(metadata["requested_maxdim"]),
        achieved_maxlinkdim=Int(metadata["achieved_maxlinkdim"]),
        checkpoint_maxlinkdim=Int(metadata["checkpoint_maxlinkdim"]),
        converged=Bool(metadata["converged"]),
        valid=Bool(metadata["valid"]),
        energy_per_site=authoritative_energy_per_site,
        raw_schmidt_polarization=Float64(metadata["raw_schmidt_polarization"]),
        sector_weights=sectors,
        fidelity_to_previous=Float64(metadata["fidelity_to_previous"]),
        fidelity_valid=Bool(metadata["fidelity_valid"]),
        mixed_fidelity=Bool(metadata["mixed_present"]) ?
            _fig2_mixed_audit_from_data(metadata) : nothing,
        momentum=(
            validation=_fig2_momentum_validation_from_data(metadata),
            valid=Bool(metadata["momentum_valid"]),
            reason=String(metadata["momentum_reason"]),
            momentum_convention=Symbol(String(metadata["momentum_convention"])),
            sector_ky_offsets=_fig2_sector_offsets_from_data(metadata),
        ),
        momentum_counting=Int.(metadata["momentum_counting"]),
        momentum_counting_evidence=_fig2_counting_from_data(metadata),
        restart_valid=Bool(metadata["restart_valid"]),
    )
    return (;
        candidate_id,
        relative_directory,
        evidence,
        sectors,
        authoritative_energy_per_site,
        provenance_valid=Bool(metadata["provenance_valid"]),
    )
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
        _fig2_authoritative_candidate_energy(candidate),
        evidence.raw_schmidt_polarization,
        evidence.fidelity_to_previous,
        evidence.fidelity_valid,
        candidate.sectors,
        evidence.momentum.validation,
        evidence.momentum.valid,
        evidence.momentum.reason,
        Symbol(string(evidence.momentum.momentum_convention)),
        Dict(evidence.momentum.sector_ky_offsets),
        copy(evidence.momentum_counting),
        evidence.momentum_counting_evidence,
        evidence.restart_valid,
        candidate.provenance_valid,
    )
end

function _write_fig2_pump_tables(spec, output, selections)
    ordered = sort(selections; by=row -> (row.dimension, row.point))
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
                isfinite(row.fidelity_to_previous)
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
                    "mixed-fidelity branch is valid; adjacent sector residual is diagnostic"
                elseif !(row.fidelity_valid && isfinite(row.fidelity_to_previous))
                    "invalid because mixed fidelity is unavailable"
                else
                    "invalid because an earlier mixed-fidelity branch step was invalid"
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
    _fig2_validated_snapshot(spec)
    started_ns = time_ns()
    isempty(strip(stage)) && throw(ArgumentError("benchmark stage must not be empty"))
    dims, phis = _validate_fig2_schedule(dimensions, fluxes)
    root = abspath(output)
    mkpath(root)
    _ensure_fig2_manifest_copy(spec, root)
    execution_provenance = operations.provenance(spec, root, 0.0)
    generation_provenance = _fig2_generation_provenance(
        spec, execution_provenance
    )
    ledger = _load_fig2_ledger(spec, root, stage)
    persisted_candidates = _fig2_unique_rows(
        ledger["candidate"], _fig2_candidate_key, "candidate"
    )
    foreach(
        row -> _validate_persisted_fig2_candidate_files(
            spec,
            root,
            row,
            generation_provenance;
            checkpoint_audit=operations.checkpoint_audit,
            progress_audit=operations.progress_audit,
        ),
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
    validated_selections = Dict{Tuple{Int,Int},Fig2Selection}()
    for row in values(persisted_selections)
        candidate_key = (
            Int(row["dimension"]),
            Int(row["point"]),
            String(row["candidate_id"]),
        )
        haskey(persisted_candidates, candidate_key) || throw(
            ArgumentError("persisted selection has no complete candidate row")
        )
        selection = _validate_persisted_fig2_selection(
            spec, root, row, persisted_candidates[candidate_key]
        )
        validated_selections[(selection.dimension, selection.point)] = selection
    end
    selected_candidate_rows = Dict(
        key => row for (key, row) in persisted_candidates
        if (key[1], key[2]) in keys(persisted_selections)
    )
    _validate_fig2_selection_replay(
        spec,
        root,
        selected_candidate_rows,
        persisted_selections,
        validated_selections,
        ;
        candidate_ids_provider=operations.candidate_ids,
    )
    _validate_fig2_mixed_reference_chain(
        root, persisted_candidates, persisted_selections
    )
    selections = Fig2Selection[]
    for dimension in dims
        previous_state = nothing
        previous_selection = nothing
        for (point, phi_y) in enumerate(phis)
            selection_key = (dimension, point)
            persisted_candidate_ids = String[
                candidate_id for (candidate_dimension, candidate_point, candidate_id)
                in keys(persisted_candidates)
                if candidate_dimension == dimension && candidate_point == point
            ]
            candidate_ids = _validate_fig2_candidate_id_set(
                operations.candidate_ids,
                spec,
                dimension,
                point,
                persisted_candidate_ids;
                require_complete=haskey(persisted_selections, selection_key),
            )
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
            mixed_reference = if point == 1
                _fig2_absent_mixed_reference()
            else
                isnothing(previous_selection) && throw(ArgumentError(
                    "later Fig. 2 point has no prior selected state reference"
                ))
                previous_key = (
                    dimension,
                    point - 1,
                    previous_selection.candidate_id,
                )
                haskey(persisted_candidates, previous_key) || throw(
                    ArgumentError(
                        "later Fig. 2 point has no prior selected candidate row"
                    )
                )
                _fig2_mixed_reference_for_selection(
                    root,
                    previous_selection,
                    persisted_candidates[previous_key],
                )
            end
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
                completion = _complete_fig2_candidate!(
                    spec,
                    root,
                    ledger,
                    dimension,
                    point,
                    phi_y,
                    candidate_id,
                    evidence,
                    mixed_reference,
                    relative_directory,
                    generation_provenance,
                    ;
                    progress_audit=operations.progress_audit,
                )
                push!(candidates, (;
                    candidate_id,
                    relative_directory,
                    evidence,
                    sectors=completion.sectors,
                    authoritative_energy_per_site=
                        completion.authoritative_energy_per_site,
                    provenance_valid=completion.provenance_valid,
                ))
                persisted_candidates[candidate_key] = last(ledger["candidate"])
            end
            selected_index, selection_basis = _select_fig2_candidate(
                candidates,
                point;
                previous_sector_weights=isnothing(previous_selection) ? nothing :
                    previous_selection.sector_weights,
                sector_tolerance=_fig2_number(
                    spec.data["convergence"], "sector_residual_tol"
                ),
            )
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
    end
    _write_fig2_pump_tables(spec, root, selections)
    _write_fig2_provenance!(
        spec,
        root,
        started_ns,
        execution_provenance,
        ledger,
    )
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

function _latest_dimension_selections(rows::AbstractVector{Fig2Selection})
    isempty(rows) && return Fig2Selection[]
    dimension = maximum(row.dimension for row in rows)
    return sort!(filter(row -> row.dimension == dimension, rows); by=row -> row.point)
end

function _validate_fig2_selection_replay_grid(spec, grouped_candidates)
    _fig2_validated_snapshot(spec)
    by_dimension = Dict{Int,Vector{Tuple{Int,Float64}}}()
    for (key, rows) in grouped_candidates
        dimension, point = key
        group_fluxes = unique(Float64(row["phi_y"]) for row in rows)
        length(group_fluxes) == 1 || throw(ArgumentError(
            "selection replay candidate group at dimension $dimension " *
            "point $point has inconsistent flux values"
        ))
        phi_y = only(group_fluxes)
        isfinite(phi_y) || throw(
            ArgumentError(
                "selection replay candidate flux at dimension $dimension " *
                "point $point is not finite",
            )
        )
        push!(
            get!(by_dimension, dimension, Tuple{Int,Float64}[]),
            (point, phi_y),
        )
    end
    for (dimension, grid) in by_dimension
        sort!(grid; by=first)
        points = first.(grid)
        phis = last.(grid)
        points == collect(eachindex(points)) || throw(ArgumentError(
            "selection replay candidate points for dimension $dimension " *
            "are not a complete contiguous grid"
        ))
        issorted(phis) && length(unique(phis)) == length(phis) || throw(
            ArgumentError(
                "selection replay candidate fluxes for dimension $dimension " *
                "must be strictly increasing",
            )
        )
    end
    return nothing
end

function _validate_fig2_selection_replay(
    spec,
    root,
    candidate_rows,
    selection_rows,
    validated_selections,
    ;
    candidate_ids_provider=_default_fig2_candidate_ids,
)
    grouped_candidates = Dict{Tuple{Int,Int},Vector{Any}}()
    for row in values(candidate_rows)
        key = (Int(row["dimension"]), Int(row["point"]))
        push!(get!(grouped_candidates, key, Any[]), row)
    end
    for rows in values(grouped_candidates)
        sort!(rows; by=row -> String(row["candidate_id"]))
    end

    _validate_fig2_selection_replay_grid(spec, grouped_candidates)
    candidate_group_keys = Set{Tuple{Int,Int}}(keys(grouped_candidates))
    selection_keys = Set{Tuple{Int,Int}}(keys(selection_rows))
    validated_selection_keys =
        Set{Tuple{Int,Int}}(keys(validated_selections))
    candidate_group_keys == selection_keys == validated_selection_keys || throw(
        ArgumentError(
            "selection replay candidate-group, selection, and validated " *
            "selection keys do not exactly match",
        )
    )
    _validate_fig2_mixed_reference_chain(
        root, candidate_rows, selection_rows
    )

    candidate_groups = Dict{Tuple{Int,Int},Vector{Any}}()
    for (key, rows) in grouped_candidates
        dimension, point = key
        _validate_fig2_candidate_id_set(
            candidate_ids_provider,
            spec,
            dimension,
            point,
            (String(row["candidate_id"]) for row in rows);
            require_complete=true,
        )
        candidates = [
            _persisted_fig2_candidate(
                spec,
                root,
                row,
                dimension,
                point,
                Float64(row["phi_y"]),
                String(row["candidate_id"]),
            ) for row in rows
        ]
        candidate_groups[key] = candidates
    end

    previous_sector_weights = Dict{Int,Dict{Int,Float64}}()
    sector_tolerance = _fig2_number(
        spec.data["convergence"], "sector_residual_tol"
    )
    for key in sort!(collect(keys(selection_rows)))
        dimension, point = key
        ledger_selection = validated_selections[key]
        candidates = get(candidate_groups, key, nothing)
        isnothing(candidates) && throw(ArgumentError(
            "selection replay has no complete candidates at dimension " *
            "$dimension point $point"
        ))
        previous = point == 1 ? nothing :
            get(previous_sector_weights, dimension, nothing)
        selected_index, replay_basis = _fig2_argument_error(
            "selection replay failed at dimension $dimension point $point"
        ) do
            _select_fig2_candidate(
                candidates,
                point;
                previous_sector_weights=previous,
                sector_tolerance,
            )
        end
        replay_winner = candidates[selected_index]
        ledger_selection.candidate_id == replay_winner.candidate_id || throw(
            ArgumentError(
                "selection replay winner mismatch at dimension $dimension " *
                "point $point: ledger selected " *
                "$(ledger_selection.candidate_id), replay selected " *
                "$(replay_winner.candidate_id)"
            )
        )
        ledger_selection.selection_basis == replay_basis || throw(
            ArgumentError(
                "selection replay basis mismatch at dimension $dimension " *
                "point $point: ledger recorded " *
                "$(ledger_selection.selection_basis), replay requires " *
                "$replay_basis"
            )
        )
        previous_sector_weights[dimension] = replay_winner.sectors
    end
    return nothing
end

function _fig2_acceptance_integrity_inputs(
    spec,
    root;
    checkpoint_audit=_default_fig2_persisted_checkpoint_audit,
    progress_audit=_default_fig2_progress_audit,
    candidate_ids_provider=_default_fig2_candidate_ids,
)
    _fig2_validated_snapshot(spec)
    _ensure_fig2_manifest_copy(spec, root)
    ledger_path = joinpath(root, "ledger.toml")
    isfile(ledger_path) || throw(ArgumentError("acceptance ledger.toml is missing"))
    raw_ledger = TOML.parsefile(ledger_path)
    stage = get(raw_ledger, "stage", nothing)
    stage isa AbstractString && !isempty(strip(stage)) || throw(
        ArgumentError("acceptance ledger stage is missing")
    )
    ledger = _load_fig2_ledger(spec, root, stage)

    provenance_path = joinpath(root, "provenance.toml")
    isfile(provenance_path) || throw(
        ArgumentError("acceptance root provenance.toml is missing")
    )
    execution_provenance = TOML.parsefile(provenance_path)
    provenance_result = _fig2_provenance_result(spec, root)
    provenance_result.valid || throw(ArgumentError(
        "acceptance root provenance is invalid: $(provenance_result.reason)"
    ))
    generation_provenance = _fig2_generation_provenance(
        spec, execution_provenance
    )

    candidates = _fig2_unique_rows(
        ledger["candidate"], _fig2_candidate_key, "candidate"
    )
    foreach(
        row -> _validate_persisted_fig2_candidate_files(
            spec,
            root,
            row,
            generation_provenance;
            checkpoint_audit,
            progress_audit,
        ),
        values(candidates),
    )
    expected_jobs = _fig2_candidate_pbs_job_ids(root, ledger)
    recorded_jobs = get(execution_provenance, "candidate_pbs_job_ids", nothing)
    recorded_jobs isa AbstractVector && String.(recorded_jobs) == expected_jobs ||
        throw(ArgumentError(
            "acceptance root provenance candidate PBS job IDs are invalid"
        ))

    selection_rows = _fig2_unique_rows(
        ledger["selection"], _fig2_selection_key, "selection"
    )
    selections = Fig2Selection[]
    validated_selections = Dict{Tuple{Int,Int},Fig2Selection}()
    for row in values(selection_rows)
        candidate_key = (
            Int(row["dimension"]),
            Int(row["point"]),
            String(row["candidate_id"]),
        )
        haskey(candidates, candidate_key) || throw(
            ArgumentError("persisted selection has no complete candidate row")
        )
        selection = _validate_persisted_fig2_selection(
            spec, root, row, candidates[candidate_key]
        )
        push!(selections, selection)
        validated_selections[(selection.dimension, selection.point)] = selection
    end
    _validate_fig2_selection_replay(
        spec,
        root,
        candidates,
        selection_rows,
        validated_selections,
        ;
        candidate_ids_provider,
    )
    sort!(selections; by=row -> (row.dimension, row.point))
    _write_fig2_pump_tables(spec, root, selections)
    return (; ledger, selections, execution_provenance)
end

function _fig2_integrity_digest(path; sha256=_fig2_file_sha256)
    isfile(path) || return "missing"
    return try
        sha256(path)
    catch error
        error isa InterruptException && rethrow()
        "unreadable"
    end
end

function _write_fig2_integrity_failure!(spec, root, error)
    evidence_path = joinpath(root, "acceptance_integrity_failure.toml")
    reason = "acceptance integrity validation failed ($(typeof(error))): $(sprint(showerror, error))"
    evidence = Dict{String,Any}(
        "format" => "fqahc_fig2_acceptance_integrity_failure_v1",
        "pass" => false,
        "reason" => reason,
        "exception_type" => string(typeof(error)),
        "exception_message" => sprint(showerror, error),
        "manifest_sha256" => spec.sha256,
        "manifest_copy_sha256" =>
            _fig2_integrity_digest(joinpath(root, "manifest.toml")),
        "ledger_sha256" =>
            _fig2_integrity_digest(joinpath(root, "ledger.toml")),
        "provenance_sha256" =>
            _fig2_integrity_digest(joinpath(root, "provenance.toml")),
    )
    _write_fig2_toml(evidence_path, evidence)
    report = Dict{String,Any}(
        row => _acceptance_row(
            false,
            "integrity_failure",
            "valid candidate/selection/provenance checksums",
            "acceptance_integrity_failure.toml",
            reason,
        ) for row in spec.acceptance_rows
    )
    _write_fig2_toml(joinpath(root, "acceptance.toml"), report)
    return report
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

function _fig2_table_integer(fields, index, table, line_number, name)
    value = tryparse(Int, fields[index])
    isnothing(value) && throw(ArgumentError(
        "$table line $line_number has invalid integer field $name"
    ))
    return value
end

function _fig2_table_float(fields, index, table, line_number, name)
    value = tryparse(Float64, fields[index])
    isnothing(value) && throw(ArgumentError(
        "$table line $line_number has invalid floating-point field $name"
    ))
    isfinite(value) || throw(ArgumentError(
        "$table line $line_number has nonfinite field $name"
    ))
    return value
end

function _fig2_table_bool(fields, index, table, line_number, name)
    value = fields[index]
    value in ("true", "false") || throw(ArgumentError(
        "$table line $line_number has invalid Boolean field $name"
    ))
    return value == "true"
end

function _validate_fig2_table_grid(rows, table)
    by_dimension = Dict{Int,Vector{eltype(rows)}}()
    for row in rows
        push!(get!(by_dimension, row.dimension, eltype(rows)[]), row)
    end
    for (dimension, dimension_rows) in by_dimension
        points = [row.point for row in dimension_rows]
        phis = [row.phi_y for row in dimension_rows]
        length(unique(points)) == length(points) || throw(ArgumentError(
            "$table has duplicate point rows for D=$dimension"
        ))
        length(unique(phis)) == length(phis) || throw(ArgumentError(
            "$table has duplicate flux rows for D=$dimension"
        ))
        issorted(points) && issorted(phis) || throw(ArgumentError(
            "$table rows for D=$dimension must be sorted by point and flux"
        ))
    end
    return nothing
end

function _read_fig2_raw_pump(path, dimension)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("pump_raw.tsv is empty"))
    expected_header =
        "dimension\tpoint\tphi_y\tcandidate_id\traw_schmidt_polarization\tdelta_raw_schmidt_polarization"
    first(lines) == expected_header || throw(
        ArgumentError("pump_raw.tsv header does not match the immutable raw format")
    )
    rows = NamedTuple[]
    for (offset, line) in enumerate(lines[2:end])
        line_number = offset + 1
        fields = split(line, '\t')
        length(fields) == 6 || throw(ArgumentError(
            "pump_raw.tsv line $line_number does not have six fields"
        ))
        parsed_dimension = _fig2_table_integer(
            fields, 1, "pump_raw.tsv", line_number, "dimension"
        )
        parsed_dimension > 0 || throw(ArgumentError(
            "pump_raw.tsv line $line_number has nonpositive dimension"
        ))
        point = _fig2_table_integer(
            fields, 2, "pump_raw.tsv", line_number, "point"
        )
        point > 0 || throw(ArgumentError(
            "pump_raw.tsv line $line_number has nonpositive point"
        ))
        phi_y = _fig2_table_float(
            fields, 3, "pump_raw.tsv", line_number, "phi_y"
        )
        candidate_id = String(fields[4])
        isempty(candidate_id) && throw(ArgumentError(
            "pump_raw.tsv line $line_number has empty candidate_id"
        ))
        raw = _fig2_table_float(
            fields, 5, "pump_raw.tsv", line_number,
            "raw_schmidt_polarization",
        )
        delta = _fig2_table_float(
            fields, 6, "pump_raw.tsv", line_number,
            "delta_raw_schmidt_polarization",
        )
        push!(rows, (;
            dimension=parsed_dimension,
            point,
            phi_y,
            candidate_id,
            raw_schmidt_polarization=raw,
            delta_raw_schmidt_polarization=delta,
        ))
    end
    _validate_fig2_table_grid(rows, "pump_raw.tsv")
    by_dimension = Dict{Int,Vector{eltype(rows)}}()
    for row in rows
        push!(get!(by_dimension, row.dimension, eltype(rows)[]), row)
    end
    for (parsed_dimension, dimension_rows) in by_dimension
        baseline = first(dimension_rows).raw_schmidt_polarization
        for row in dimension_rows
            expected_delta = row.raw_schmidt_polarization - baseline
            _fig2_same_float(
                row.delta_raw_schmidt_polarization, expected_delta
            ) || throw(ArgumentError(
                "pump_raw.tsv delta disagrees with the D=$parsed_dimension baseline at point $(row.point)"
            ))
        end
    end
    selected = filter(row -> row.dimension == dimension, rows)
    isempty(selected) && throw(ArgumentError(
        "pump_raw.tsv has no selected rows for D=$dimension"
    ))
    return selected
end

function _read_fig2_tracked_pump(path, dimension)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("pump_sector_tracked.tsv is empty"))
    expected_header =
        "dimension\tpoint\tphi_y\tcandidate_id\talignment_shift_to_previous_raw\tcumulative_sector_shift\tsector_distance_raw\tvalid\treason"
    first(lines) == expected_header || throw(ArgumentError(
        "pump_sector_tracked.tsv header does not match the immutable tracked format"
    ))
    rows = NamedTuple[]
    for (offset, line) in enumerate(lines[2:end])
        line_number = offset + 1
        fields = split(line, '\t')
        length(fields) == 9 || throw(ArgumentError(
            "pump_sector_tracked.tsv line $line_number does not have nine fields"
        ))
        parsed_dimension = _fig2_table_integer(
            fields, 1, "pump_sector_tracked.tsv", line_number, "dimension"
        )
        parsed_dimension > 0 || throw(ArgumentError(
            "pump_sector_tracked.tsv line $line_number has nonpositive dimension"
        ))
        point = _fig2_table_integer(
            fields, 2, "pump_sector_tracked.tsv", line_number, "point"
        )
        point > 0 || throw(ArgumentError(
            "pump_sector_tracked.tsv line $line_number has nonpositive point"
        ))
        phi_y = _fig2_table_float(
            fields, 3, "pump_sector_tracked.tsv", line_number, "phi_y"
        )
        candidate_id = String(fields[4])
        isempty(candidate_id) && throw(ArgumentError(
            "pump_sector_tracked.tsv line $line_number has empty candidate_id"
        ))
        alignment_shift = _fig2_table_integer(
            fields, 5, "pump_sector_tracked.tsv", line_number,
            "alignment_shift_to_previous_raw",
        )
        cumulative_shift = _fig2_table_integer(
            fields, 6, "pump_sector_tracked.tsv", line_number,
            "cumulative_sector_shift",
        )
        sector_distance = _fig2_table_float(
            fields, 7, "pump_sector_tracked.tsv", line_number,
            "sector_distance_raw",
        )
        sector_distance >= 0 || throw(ArgumentError(
            "pump_sector_tracked.tsv line $line_number has negative sector residual"
        ))
        valid = _fig2_table_bool(
            fields, 8, "pump_sector_tracked.tsv", line_number, "valid"
        )
        reason = String(fields[9])
        isempty(strip(reason)) && throw(ArgumentError(
            "pump_sector_tracked.tsv line $line_number has empty reason"
        ))
        push!(rows, (;
            dimension=parsed_dimension,
            point,
            phi_y,
            candidate_id,
            alignment_shift_to_previous_raw=alignment_shift,
            cumulative_sector_shift=cumulative_shift,
            sector_distance_raw=sector_distance,
            valid,
            reason,
        ))
    end
    _validate_fig2_table_grid(rows, "pump_sector_tracked.tsv")
    selected = filter(row -> row.dimension == dimension, rows)
    isempty(selected) && throw(ArgumentError(
        "pump_sector_tracked.tsv has no selected rows for D=$dimension"
    ))
    return selected
end

function _fig2_aligned_curve_rows(
    spec,
    dimension,
    raw_rows,
    tracked_rows,
    validated_selections,
)
    selections = sort!(filter(
        row -> row.dimension == dimension, collect(validated_selections)
    ); by=row -> row.point)
    isempty(selections) && throw(ArgumentError(
        "paper curve has no integrity-validated selections for D=$dimension"
    ))
    length(raw_rows) == length(selections) || throw(ArgumentError(
        "pump_raw.tsv row count disagrees with integrity-validated selections for D=$dimension"
    ))
    length(tracked_rows) == length(selections) || throw(ArgumentError(
        "pump_sector_tracked.tsv row count disagrees with integrity-validated selections for D=$dimension"
    ))

    cumulative_shift = 0
    path_valid = true
    aligned = NamedTuple[]
    for index in eachindex(selections)
        selection = selections[index]
        raw = raw_rows[index]
        tracked = tracked_rows[index]
        expected_identity = (
            dimension,
            selection.point,
            selection.phi_y,
            selection.candidate_id,
        )
        raw_identity = (
            raw.dimension, raw.point, raw.phi_y, raw.candidate_id
        )
        raw_identity == expected_identity || throw(ArgumentError(
            "pump_raw.tsv dimension/point/phi/candidate identity disagrees with integrity-validated selection at D=$dimension point $(selection.point)"
        ))
        _fig2_same_float(
            raw.raw_schmidt_polarization,
            selection.raw_schmidt_polarization,
        ) || throw(ArgumentError(
            "pump_raw.tsv raw polarization disagrees with integrity-validated candidate at D=$dimension point $(selection.point)"
        ))
        tracked_identity = (
            tracked.dimension,
            tracked.point,
            tracked.phi_y,
            tracked.candidate_id,
        )
        tracked_identity == expected_identity || throw(ArgumentError(
            "pump_sector_tracked.tsv dimension/point/phi/candidate identity disagrees with integrity-validated selection at D=$dimension point $(selection.point)"
        ))

        expected_alignment = 0
        expected_distance = 0.0
        expected_valid = true
        if index > 1
            previous = selections[index - 1]
            alignment = best_sector_shift(
                previous.sector_weights, selection.sector_weights
            )
            expected_alignment = alignment.shift
            expected_distance = alignment.distance
            step_valid = selection.fidelity_valid &&
                isfinite(selection.fidelity_to_previous)
            expected_valid = path_valid && step_valid
            expected_valid && (cumulative_shift -= expected_alignment)
        end
        path_valid = expected_valid
        tracked.alignment_shift_to_previous_raw == expected_alignment &&
            tracked.cumulative_sector_shift == cumulative_shift &&
            _fig2_same_float(tracked.sector_distance_raw, expected_distance) &&
            tracked.valid == expected_valid || throw(ArgumentError(
                "pump_sector_tracked.tsv branch evidence disagrees with integrity-validated mixed-fidelity sector correspondence at D=$dimension point $(selection.point); refinement required"
            ))
        push!(aligned, (; raw, tracked, selection))
    end
    return aligned
end

function _fig2_interpolation_bracket(rows, target_phi, sector_tolerance)
    phis = [row.raw.phi_y for row in rows]
    target_phi < first(phis) - 1e-10 && throw(ArgumentError(
        "paper-curve target requires left extrapolation; refinement covering the target interval is required"
    ))
    target_phi > last(phis) + 1e-10 && throw(ArgumentError(
        "paper-curve target requires right extrapolation; refinement covering the target interval is required"
    ))
    exact = findfirst(phi -> isapprox(
        phi, target_phi; atol=1e-10, rtol=0
    ), phis)
    if !isnothing(exact)
        left_index = exact
        right_index = exact
        right_weight = 0.0
        sampling = "exact"
    else
        right_index = searchsortedfirst(phis, target_phi)
        1 < right_index <= length(rows) || throw(ArgumentError(
            "paper-curve interpolation bracket is missing; refinement is required"
        ))
        left_index = right_index - 1
        left_phi = phis[left_index]
        right_phi = phis[right_index]
        right_weight = (target_phi - left_phi) / (right_phi - left_phi)
        0 < right_weight < 1 || throw(ArgumentError(
            "paper-curve interpolation weight is invalid; refinement is required"
        ))
        sampling = "linear_interpolation"
    end
    left = rows[left_index]
    right = rows[right_index]
    branch_valid = left.tracked.valid && right.tracked.valid
    branch_valid || throw(ArgumentError(
        "paper-curve interpolation crosses an invalid mixed-fidelity charge-sector branch; refinement is required"
    ))
    left_weight = 1 - right_weight
    interpolated_raw =
        left_weight * left.raw.raw_schmidt_polarization +
        right_weight * right.raw.raw_schmidt_polarization
    return (;
        left,
        right,
        left_weight,
        right_weight,
        interpolated_raw,
        sampling,
        exact_match=!isnothing(exact),
    )
end

function _fig2_paper_curve_result(
    spec,
    root,
    dimension,
    validated_selections;
    comparison_filename="paper_curve_comparison.toml",
    raw_reader=_read_fig2_raw_pump,
)
    _fig2_validated_snapshot(spec)
    rule = spec.data["paper_curve_rule"]
    raw_path = joinpath(root, "pump_raw.tsv")
    tracked_path = joinpath(root, "pump_sector_tracked.tsv")
    comparison_path = joinpath(root, comparison_filename)
    uncertainty = _fig2_number_vector(rule, "uncertainty")
    comparison = Dict{String,Any}(
        "format" => rule["format"],
        "manifest_sha256" => spec.sha256,
        "source_pdf_sha256" => rule["source_pdf_sha256"],
        "source_page" => rule["source_page"],
        "source_panel" => rule["source_panel"],
        "source_vector_path_sha256" => rule["source_vector_path_sha256"],
        "target" => rule["target"],
        "extraction" => rule["extraction"],
        "raw_pump_policy" => spec.data["raw_pump_policy"],
        "sampling_policy" =>
            "linear_interpolation_on_validated_raw_grid_without_extrapolation",
        "dimension" => dimension,
        "point" => Dict{String,Any}[],
    )
    try
        isfile(raw_path) || throw(ArgumentError("pump_raw.tsv is missing"))
        isfile(tracked_path) || throw(ArgumentError(
            "pump_sector_tracked.tsv is missing"
        ))
        comparison["pump_raw_sha256"] = _fig2_file_sha256(raw_path)
        comparison["pump_sector_tracked_sha256"] =
            _fig2_file_sha256(tracked_path)
        raw_rows = raw_reader(raw_path, dimension)
        tracked_rows = _read_fig2_tracked_pump(tracked_path, dimension)
        aligned_rows = _fig2_aligned_curve_rows(
            spec,
            dimension,
            raw_rows,
            tracked_rows,
            validated_selections,
        )
        reference_fluxes = _fig2_number_vector(rule, "flux_units_2pi")
        reference_values = _fig2_number_vector(rule, "delta_q")
        sector_tolerance = _fig2_number(
            spec.data["convergence"], "sector_residual_tol"
        )
        baseline = _fig2_interpolation_bracket(
            aligned_rows,
            2pi * first(reference_fluxes),
            sector_tolerance,
        )
        points = Dict{String,Any}[]
        for index in eachindex(reference_fluxes)
            flux_units = reference_fluxes[index]
            target_phi = 2pi * flux_units
            bracket = _fig2_interpolation_bracket(
                aligned_rows, target_phi, sector_tolerance
            )
            measured = bracket.interpolated_raw - baseline.interpolated_raw
            residual = abs(measured - reference_values[index])
            point_pass = residual <= uncertainty[index]
            left = bracket.left
            right = bracket.right
            push!(points, Dict{String,Any}(
                "flux_units_2pi" => flux_units,
                "phi_y" => target_phi,
                "sampling" => bracket.sampling,
                "exact_match" => bracket.exact_match,
                "left_point" => left.raw.point,
                "right_point" => right.raw.point,
                "left_candidate_id" => left.raw.candidate_id,
                "right_candidate_id" => right.raw.candidate_id,
                "left_phi_y" => left.raw.phi_y,
                "right_phi_y" => right.raw.phi_y,
                "left_raw_schmidt_polarization" =>
                    left.raw.raw_schmidt_polarization,
                "right_raw_schmidt_polarization" =>
                    right.raw.raw_schmidt_polarization,
                "left_weight" => bracket.left_weight,
                "right_weight" => bracket.right_weight,
                "left_branch_valid" => left.tracked.valid,
                "right_branch_valid" => right.tracked.valid,
                "left_sector_residual" => left.tracked.sector_distance_raw,
                "right_sector_residual" => right.tracked.sector_distance_raw,
                "left_alignment_shift_to_previous_raw" =>
                    left.tracked.alignment_shift_to_previous_raw,
                "right_alignment_shift_to_previous_raw" =>
                    right.tracked.alignment_shift_to_previous_raw,
                "left_cumulative_sector_shift" =>
                    left.tracked.cumulative_sector_shift,
                "right_cumulative_sector_shift" =>
                    right.tracked.cumulative_sector_shift,
                "raw_schmidt_polarization" => bracket.interpolated_raw,
                "baseline_raw_schmidt_polarization" =>
                    baseline.interpolated_raw,
                "measured" => measured,
                "reference" => reference_values[index],
                "uncertainty" => uncertainty[index],
                "absolute_residual" => residual,
                "pass" => point_pass,
            ))
        end
        maximum_residual = maximum(point["absolute_residual"] for point in points)
        pass = all(point["pass"] for point in points)
        reason = pass ?
            "all linearly interpolated raw Fig. 2(d) residuals satisfy their digitization uncertainties on valid tracked branches" :
            "at least one linearly interpolated raw Fig. 2(d) residual exceeds its digitization uncertainty"
        comparison["point"] = points
        comparison["maximum_residual"] = maximum_residual
        comparison["tolerance"] = maximum(uncertainty)
        comparison["pass"] = pass
        comparison["reason"] = reason
        _write_fig2_toml(comparison_path, comparison)
        return (;
            pass,
            measured=maximum_residual,
            reason,
            evidence=String(comparison_filename),
            sha256=_fig2_file_sha256(comparison_path),
        )
    catch error
        error isa InterruptException && rethrow()
        reason = "paper curve recomputation from validated raw/tracked pump interpolation failed: $(sprint(showerror, error))"
        comparison["maximum_residual"] = "missing"
        comparison["tolerance"] = maximum(uncertainty)
        comparison["pass"] = false
        comparison["reason"] = reason
        _write_fig2_toml(comparison_path, comparison)
        return (
            pass=false,
            measured="missing",
            reason,
            evidence=String(comparison_filename),
            sha256=_fig2_file_sha256(comparison_path),
        )
    end
end

function _fig2_provenance_result(
    spec,
    root;
    toml_parsefile=TOML.parsefile,
)
    path = joinpath(root, "provenance.toml")
    isfile(path) || return (
        valid=false,
        measured="missing",
        reason="root provenance.toml is missing",
    )
    try
        provenance = toml_parsefile(path)
        active_project, project_manifest = _fig2_project_manifest()
        benchmark_source = abspath(@__FILE__)
        commit = get(provenance, "git_commit", "")
        pbs_job_id = get(provenance, "pbs_job_id", "")
        runtime_seconds = get(provenance, "runtime_seconds", nothing)
        valid =
            get(provenance, "format", "") == "fqahc_fig2_provenance_v2" &&
            get(provenance, "manifest_sha256", "") == spec.sha256 &&
            commit isa AbstractString && occursin(r"^[0-9a-f]{40}$", commit) &&
            commit == _fig2_repository_commit() &&
            get(provenance, "git_tree_clean", false) === true &&
            get(provenance, "julia_version", "") == "1.12.5" &&
            get(provenance, "julia_version", "") == string(VERSION) &&
            pbs_job_id isa AbstractString &&
            occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", pbs_job_id) &&
            _fig2_threading_provenance_valid(provenance) &&
            runtime_seconds isa Real && !(runtime_seconds isa Bool) &&
            isfinite(Float64(runtime_seconds)) && Float64(runtime_seconds) >= 0 &&
            get(provenance, "active_project", "") == active_project &&
            get(provenance, "project_manifest", "") == project_manifest &&
            get(provenance, "project_manifest_sha256", "") ==
                _fig2_file_sha256(project_manifest) &&
            get(toml_parsefile(project_manifest), "julia_version", "") == "1.12.5" &&
            get(provenance, "benchmark_source", "") == benchmark_source &&
            get(provenance, "benchmark_source_sha256", "") ==
                _fig2_file_sha256(benchmark_source)
        return (;
            valid,
            measured=valid ? String(commit) : "invalid",
            reason=valid ?
                "root provenance contract and source checksums are valid" :
                "root provenance commit, contract, or source checksums are invalid",
        )
    catch error
        error isa InterruptException && rethrow()
        return (
            valid=false,
            measured="invalid",
            reason="root provenance validation failed: $(sprint(showerror, error))",
        )
    end
end

function _fig2_pinned_backend_revision(active_project, project_manifest)
    project = TOML.parsefile(active_project)
    source = get(get(project, "sources", Dict{String,Any}()),
        "ITensorInfiniteMPS", nothing)
    source isa AbstractDict || throw(ArgumentError(
        "active project has no pinned ITensorInfiniteMPS backend source"
    ))
    project_revision = get(source, "rev", "")

    manifest = TOML.parsefile(project_manifest)
    dependencies = get(manifest, "deps", nothing)
    dependencies isa AbstractDict || throw(ArgumentError(
        "project manifest has no dependency table for the pinned backend"
    ))
    entries = get(dependencies, "ITensorInfiniteMPS", nothing)
    entries isa AbstractVector && length(entries) == 1 &&
        only(entries) isa AbstractDict || throw(ArgumentError(
        "project manifest ITensorInfiniteMPS backend entry is missing or ambiguous"
    ))
    manifest_revision = get(only(entries), "repo-rev", "")
    project_revision == ITENSOR_INFINITE_MPS_COMMIT &&
        manifest_revision == ITENSOR_INFINITE_MPS_COMMIT || throw(ArgumentError(
        "project dependency provenance does not use the pinned ITensorInfiniteMPS backend commit"
    ))
    return ITENSOR_INFINITE_MPS_COMMIT
end

function _fig2_restart_gate_result(
    spec,
    root,
    root_provenance;
    parse_gate=TOML.parsefile,
)
    path = joinpath(root, "restart_gate.toml")
    isfile(path) || return (
        valid=false,
        measured="missing",
        reason="global two-process restart gate restart_gate.toml is missing",
    )
    try
        snapshot = _fig2_validated_snapshot(spec)
        contract = snapshot["restart_gate"]
        gate = parse_gate(path)
        root_provenance isa AbstractDict || throw(ArgumentError(
            "root provenance is unavailable for restart-gate comparison"
        ))

        get(gate, "format", "") == contract["format"] || throw(
            ArgumentError("restart-gate format does not match the fixed contract")
        )
        get(gate, "pass", nothing) === true || throw(
            ArgumentError("restart-gate pass flag is not exactly true")
        )
        get(gate, "manifest_sha256", "") == spec.sha256 || throw(
            ArgumentError("restart-gate manifest checksum does not match")
        )

        commit = get(gate, "git_commit", "")
        commit isa AbstractString && occursin(r"^[0-9a-f]{40}$", commit) &&
            commit == get(root_provenance, "git_commit", "") &&
            commit == _fig2_repository_commit() || throw(ArgumentError(
            "restart-gate commit does not match root provenance and actual HEAD"
        ))

        julia_version = get(gate, "julia_version", "")
        julia_version == "1.12.5" &&
            julia_version == string(VERSION) &&
            julia_version == get(root_provenance, "julia_version", "") || throw(
            ArgumentError("restart-gate Julia version is invalid")
        )

        active_project, project_manifest = _fig2_project_manifest()
        project_manifest_sha256 = _fig2_file_sha256(project_manifest)
        get(gate, "project_manifest_sha256", "") ==
            project_manifest_sha256 &&
            get(gate, "project_manifest_sha256", "") ==
                get(root_provenance, "project_manifest_sha256", "") || throw(
            ArgumentError(
                "restart-gate project manifest dependency provenance is invalid"
            )
        )
        backend_commit = _fig2_pinned_backend_revision(
            active_project, project_manifest
        )
        get(gate, "backend_commit", "") == backend_commit &&
            get(gate, "backend_commit", "") == contract["backend_commit"] ||
            throw(ArgumentError(
                "restart-gate backend commit is not the pinned dependency"
            ))

        function process_exit_code(key, description)
            value = get(gate, key, nothing)
            value isa Integer && !(value isa Bool) && Int(value) == 0 || throw(
                ArgumentError("restart-gate $description process exit code is not zero")
            )
            return Int(value)
        end
        process_exit_code("save_process_exit_code", "save")
        process_exit_code("resume_process_exit_code", "resume")

        function process_id(key, description)
            value = get(gate, key, nothing)
            value isa AbstractString &&
                occursin(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$", value) || throw(
                ArgumentError("restart-gate $description process ID is invalid")
            )
            return String(value)
        end
        save_process_id = process_id("save_process_id", "save")
        resume_process_id = process_id("resume_process_id", "resume")
        save_process_id != resume_process_id || throw(ArgumentError(
            "restart-gate save and resume process IDs must be distinct"
        ))

        required_checks = String.(contract["required_checks"])
        checks = Dict{String,Bool}()
        for check in required_checks
            get(gate, check, nothing) === true || throw(ArgumentError(
                "restart-gate required check $check is not exactly true"
            ))
            checks[check] = true
        end

        function positive_integer(key)
            value = get(gate, key, nothing)
            value isa Integer && !(value isa Bool) && Int(value) > 0 || throw(
                ArgumentError("restart-gate $key must be a positive integer")
            )
            return Int(value)
        end
        pre_maxlinkdim = positive_integer("pre_maxlinkdim")
        post_maxlinkdim = positive_integer("post_maxlinkdim")
        post_maxlinkdim > pre_maxlinkdim || throw(ArgumentError(
            "restart-gate post_maxlinkdim must exceed pre_maxlinkdim"
        ))

        quantitative = Dict{String,Float64}()
        for (error_key, tolerance_key, description) in (
            ("energy_absolute_error", "energy_tolerance", "energy"),
            ("density_max_absolute_error", "density_tolerance", "density"),
        )
            has_error = haskey(gate, error_key)
            has_tolerance = haskey(gate, tolerance_key)
            has_error && has_tolerance || throw(ArgumentError(
                "restart-gate $description error and tolerance are mandatory"
            ))
            error_value = gate[error_key]
            tolerance = gate[tolerance_key]
            error_value isa Real && !(error_value isa Bool) &&
                isfinite(Float64(error_value)) && Float64(error_value) >= 0 ||
                throw(ArgumentError(
                    "restart-gate $description error must be finite and nonnegative"
                ))
            tolerance isa Real && !(tolerance isa Bool) &&
                isfinite(Float64(tolerance)) && Float64(tolerance) >= 0 ||
                throw(ArgumentError(
                    "restart-gate $description tolerance must be finite and nonnegative"
                ))
            manifest_tolerance = _fig2_positive_number(
                contract,
                tolerance_key,
                "restart-gate immutable $description tolerance",
            )
            Float64(tolerance) == manifest_tolerance || throw(ArgumentError(
                "restart-gate $description tolerance does not match the immutable manifest"
            ))
            Float64(error_value) <= manifest_tolerance || throw(ArgumentError(
                "restart-gate $description reproduction error exceeds tolerance"
            ))
            quantitative[error_key] = Float64(error_value)
            quantitative[tolerance_key] = Float64(tolerance)
        end

        measured = Dict{String,Any}(
            "restart_gate_sha256" => _fig2_file_sha256(path),
            "save_process_id" => save_process_id,
            "resume_process_id" => resume_process_id,
            "pre_maxlinkdim" => pre_maxlinkdim,
            "post_maxlinkdim" => post_maxlinkdim,
            "required_checks" => checks,
        )
        merge!(measured, quantitative)
        return (
            valid=true,
            measured,
            reason="validated global two-process load, no-expansion iteration, controlled expansion, observable reproduction, site identity, and next-flux evidence",
        )
    catch error
        error isa InterruptException && rethrow()
        return (
            valid=false,
            measured="invalid",
            reason="global two-process restart-gate validation failed: $(sprint(showerror, error))",
        )
    end
end

function _fig2_sector_path_result(spec, rows)
    tolerance = _fig2_number(spec.data["convergence"], "sector_residual_tol")
    for index in 2:length(rows)
        current = rows[index]
        step_valid = current.fidelity_valid &&
            isfinite(current.fidelity_to_previous)
        step_valid || return (
            valid=false,
            shift=0,
            raw_shift=0,
            maximum_distance=Inf,
        )
    end
    isempty(rows) && return (
        valid=false,
        shift=0,
        raw_shift=0,
        maximum_distance=Inf,
    )
    alignment = best_sector_shift(
        first(rows).sector_weights, last(rows).sector_weights
    )
    raw_shift = -alignment.shift
    scale = charge_scale(spec.config)
    divisible = rem(raw_shift, scale) == 0
    shift = div(raw_shift, scale)
    valid = divisible && isfinite(alignment.distance) &&
        alignment.distance <= tolerance
    return (; valid, shift, raw_shift, maximum_distance=alignment.distance)
end

function _fig2_momentum_endpoint_valid(spec, row)
    snapshot = _fig2_validated_snapshot(spec)
    momentum_residual_tolerance = _fig2_number(
        snapshot["optimization"], "momentum_residual_tol"
    )
    transfer_tolerance = _fig2_number(
        snapshot["optimization"], "transfer_tol"
    )
    validation = row.momentum_validation
    momentum_residuals = (
        validation.unitary_residual,
        validation.power_residual,
        validation.commutator_residual,
    )
    return row.momentum_valid && validation.valid &&
        all(
            residual -> isfinite(residual) &&
                residual <= momentum_residual_tolerance,
            momentum_residuals,
        ) &&
        isfinite(validation.transfer_residual) &&
        validation.transfer_residual <= transfer_tolerance &&
        string(row.momentum_convention) == spec.data["momentum_convention"] &&
        row.momentum_counting_evidence.valid &&
        row.momentum_counting == Int.(spec.data["counting"])
end

function _fig2_high_dimension_gate(
    spec,
    previous,
    current,
    previous_curve,
    current_curve,
)
    energy_tolerance = _fig2_number(
        spec.data["convergence"], "energy_per_site_change_tol"
    )
    pump_change_tolerance = _fig2_number(
        spec.data["convergence"], "endpoint_pump_change_tol"
    )
    pump_target_tolerance = _fig2_number(spec.data, "pump_tolerance")
    expected_pump = _fig2_number(spec.data, "expected_endpoint_pump")
    expected_shift = _fig2_integer(spec.data, "expected_sector_shift")
    energy_change = abs(
        last(current).energy_per_site - last(previous).energy_per_site
    )
    previous_pump = last(previous).raw_schmidt_polarization -
        first(previous).raw_schmidt_polarization
    current_pump = last(current).raw_schmidt_polarization -
        first(current).raw_schmidt_polarization
    pump_change = abs(current_pump - previous_pump)
    measured = Dict{String,Any}(
        "lower_dimension" => first(previous).dimension,
        "upper_dimension" => first(current).dimension,
        "energy_per_site_change" => energy_change,
        "lower_endpoint_pump" => previous_pump,
        "upper_endpoint_pump" => current_pump,
        "endpoint_pump_change" => pump_change,
    )
    energy_change <= energy_tolerance || return (
        pass=false,
        measured,
        reason="required high-D endpoint energy change exceeds tolerance",
    )
    abs(previous_pump - expected_pump) <= pump_target_tolerance &&
        abs(current_pump - expected_pump) <= pump_target_tolerance &&
        pump_change <= pump_change_tolerance || return (
        pass=false,
        measured,
        reason="required high-D endpoint pump target or pump convergence failed",
    )
    previous_sector = _fig2_sector_path_result(spec, previous)
    current_sector = _fig2_sector_path_result(spec, current)
    merge!(measured, Dict(
        "lower_sector_shift" => previous_sector.shift,
        "upper_sector_shift" => current_sector.shift,
        "lower_sector_residual" => previous_sector.maximum_distance,
        "upper_sector_residual" => current_sector.maximum_distance,
    ))
    previous_sector.valid && current_sector.valid &&
        previous_sector.shift == expected_shift &&
        current_sector.shift == expected_shift || return (
        pass=false,
        measured,
        reason="required high-D sector shift or sector residual gate failed",
    )
    _fig2_momentum_endpoint_valid(spec, last(previous)) &&
        _fig2_momentum_endpoint_valid(spec, last(current)) || return (
        pass=false,
        measured,
        reason="required high-D momentum residual, convention, or counting gate failed",
    )
    all(row -> row.restart_valid && row.provenance_valid, (previous..., current...)) ||
        return (
            pass=false,
            measured,
            reason="required high-D restart or provenance gate failed",
        )
    merge!(measured, Dict(
        "lower_paper_curve_maximum_residual" => previous_curve.measured,
        "upper_paper_curve_maximum_residual" => current_curve.measured,
        "lower_paper_curve_evidence" => previous_curve.evidence,
        "upper_paper_curve_evidence" => current_curve.evidence,
        "lower_paper_curve_sha256" => previous_curve.sha256,
        "upper_paper_curve_sha256" => current_curve.sha256,
    ))
    previous_curve.pass && current_curve.pass || return (
        pass=false,
        measured,
        reason="required high-D paper Fig. 2(d) curve gate failed for one or both dimensions",
    )
    return (
        pass=true,
        measured,
        reason="required high-D energy, pump, sector, momentum, curve, restart, and provenance gates pass",
    )
end

function _fig2_early_plateau_evidence(
    spec,
    root;
    toml_parsefile=TOML.parsefile,
)
    snapshot = _fig2_validated_snapshot(spec)
    contract = snapshot["early_plateau"]
    filename = String(contract["evidence_filename"])
    lower_dimension = Int(contract["lower_dimension"])
    upper_dimension = Int(contract["upper_dimension"])
    omitted_dimension = Int(contract["omitted_dimension"])
    path = joinpath(root, filename)
    isfile(path) || return (
        valid=false,
        reason="$filename is missing",
        evidence=filename,
        sha256="missing",
    )
    try
        evidence = toml_parsefile(path)
        get(evidence, "format", "") == "fqahc_fig2_early_plateau_v1" ||
            return (
                valid=false,
                reason="early_plateau format is invalid",
                evidence=filename,
                sha256=_fig2_file_sha256(path),
            )
        get(evidence, "manifest_sha256", "") == spec.sha256 || return (
            valid=false,
            reason="early_plateau manifest checksum is invalid",
            evidence=filename,
            sha256=_fig2_file_sha256(path),
        )
        dimensions_valid =
            _fig2_integer(evidence, "omit_dimension") == omitted_dimension &&
            _fig2_integer(evidence, "lower_dimension") == lower_dimension &&
            _fig2_integer(evidence, "upper_dimension") == upper_dimension
        dimensions_valid || return (
            valid=false,
            reason="early_plateau dimension contract is invalid",
            evidence=filename,
            sha256=_fig2_file_sha256(path),
        )
        reason = get(evidence, "reason", nothing)
        reason isa AbstractString && !isempty(strip(reason)) || return (
            valid=false,
            reason="early_plateau reason is missing",
            evidence=filename,
            sha256=_fig2_file_sha256(path),
        )
        occursin("D=$lower_dimension", reason) &&
            occursin("D=$upper_dimension", reason) || return (
            valid=false,
            reason="early_plateau reason must explain the " *
                "D=$lower_dimension/D=$upper_dimension plateau",
            evidence=filename,
            sha256=_fig2_file_sha256(path),
        )
        return (
            valid=true,
            reason=String(reason),
            evidence=filename,
            sha256=_fig2_file_sha256(path),
        )
    catch error
        error isa InterruptException && rethrow()
        return (
            valid=false,
            reason="early_plateau validation failed: $(sprint(showerror, error))",
            evidence=filename,
            sha256=_fig2_file_sha256(path),
        )
    end
end

function _fig2_bond_dimension_pair(spec, root, dimensions)
    snapshot = _fig2_validated_snapshot(spec)
    declared = _fig2_integer_vector(snapshot["convergence"], "maxdims")
    early_plateau = snapshot["early_plateau"]
    required_progression = declared[1:(end - 1)]
    missing = filter(dimension -> !(dimension in dimensions), required_progression)
    isempty(missing) || return (
        valid=false,
        pair=nothing,
        mode="missing_progression",
        evidence="ledger.toml",
        evidence_sha256=_fig2_file_sha256(joinpath(root, "ledger.toml")),
        reason="bond convergence is missing required declared progression D=$(join(missing, ","))",
    )
    final_dimension = Int(early_plateau["omitted_dimension"])
    if final_dimension in dimensions
        return (
            valid=true,
            pair=(declared[end - 1], final_dimension),
            mode="full_D3000_progression",
            evidence="ledger.toml",
            evidence_sha256=_fig2_file_sha256(joinpath(root, "ledger.toml")),
            reason="full declared progression includes D=$final_dimension",
        )
    end
    early = _fig2_early_plateau_evidence(spec, root)
    early.valid || return (
        valid=false,
        pair=nothing,
        mode="missing_D3000",
        evidence=early.evidence,
        evidence_sha256=early.sha256,
        reason="D=$final_dimension is required unless explicit early_plateau evidence is valid: $(early.reason)",
    )
    return (
        valid=true,
        pair=(
            Int(early_plateau["lower_dimension"]),
            Int(early_plateau["upper_dimension"]),
        ),
        mode="audited_early_plateau",
        evidence=early.evidence,
        evidence_sha256=early.sha256,
        reason=early.reason,
    )
end

function write_fig2_acceptance_report!(
    spec::Fig2BenchmarkSpec,
    output::AbstractString,
    ;
    checkpoint_audit=_default_fig2_persisted_checkpoint_audit,
    progress_audit=_default_fig2_progress_audit,
    candidate_ids_provider=_default_fig2_candidate_ids,
)
    _fig2_validated_snapshot(spec)
    root = abspath(output)
    mkpath(root)
    integrity = try
        _fig2_acceptance_integrity_inputs(
            spec,
            root;
            checkpoint_audit,
            progress_audit,
            candidate_ids_provider,
        )
    catch error
        error isa InterruptException && rethrow()
        return _write_fig2_integrity_failure!(spec, root, error)
    end
    ledger = integrity.ledger
    selection_rows = integrity.selections
    selected = _latest_dimension_selections(selection_rows)
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

    sector_path = isempty(selected) ? nothing :
        _fig2_sector_path_result(spec, selected)
    expected_shift = _fig2_integer(spec.data, "expected_sector_shift")
    sector_tolerance = _fig2_number(spec.data["convergence"], "sector_residual_tol")
    sector_pass = !isnothing(sector_path) && sector_path.valid &&
        sector_path.shift == expected_shift
    sector_reason = if isnothing(sector_path)
        "insufficient mixed-fidelity endpoint entanglement-spectrum evidence"
    elseif sector_pass
        "mixed-fidelity path and endpoint entanglement-spectrum raw shift $(sector_path.raw_shift)/charge scale $(charge_scale(spec.config)) match the physical charge target"
    else
        "mixed-fidelity path or endpoint entanglement-spectrum raw shift $(sector_path.raw_shift)/charge scale $(charge_scale(spec.config)) is invalid"
    end
    report["sector_shift"] = _acceptance_row(
        sector_pass,
        isnothing(sector_path) ? "missing" : sector_path.shift,
        sector_tolerance,
        "pump_sector_tracked.tsv",
        sector_reason,
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

    curve_result = if isempty(selected)
        _fig2_paper_curve_result(spec, root, 0, selection_rows)
    else
        _fig2_paper_curve_result(
            spec, root, last(selected).dimension, selection_rows
        )
    end
    curve_tolerance = maximum(_fig2_number_vector(
        spec.data["paper_curve_rule"], "uncertainty"
    ))
    report["paper_curve"] = _acceptance_row(
        curve_result.pass, curve_result.measured, curve_tolerance,
        "paper_curve_comparison.toml", curve_result.reason,
    )

    dimensions = sort!(unique(row.dimension for row in selection_rows))
    schedule = _fig2_bond_dimension_pair(spec, root, dimensions)
    convergence_pass = schedule.valid
    convergence_evidence = schedule.evidence
    convergence_measured = Dict{String,Any}(
        "available_dimensions" => dimensions,
        "mode" => schedule.mode,
        "schedule_evidence_sha256" => schedule.evidence_sha256,
        "ledger_sha256" => _fig2_file_sha256(joinpath(root, "ledger.toml")),
    )
    convergence_reason = schedule.reason
    if schedule.valid
        high_dimensions = schedule.pair
        previous = sort!(filter(
            row -> row.dimension == first(high_dimensions), selection_rows
        ); by=row -> row.point)
        current = sort!(filter(
            row -> row.dimension == last(high_dimensions), selection_rows
        ); by=row -> row.point)
        start_flux = _fig2_number(spec.data, "flux_start")
        stop_flux = _fig2_number(spec.data, "flux_stop")
        complete(rows) = length(rows) >= 2 &&
            first(rows).phi_y == start_flux && last(rows).phi_y == stop_flux
        if !(complete(previous) && complete(current))
            convergence_pass = false
            convergence_reason =
                "required high-D start/stop endpoint evidence is incomplete"
        else
            previous_curve = _fig2_paper_curve_result(
                spec,
                root,
                first(high_dimensions),
                selection_rows;
                comparison_filename=
                    "paper_curve_comparison_D$(first(high_dimensions)).toml",
            )
            current_curve = _fig2_paper_curve_result(
                spec,
                root,
                last(high_dimensions),
                selection_rows;
                comparison_filename=
                    "paper_curve_comparison_D$(last(high_dimensions)).toml",
            )
            convergence = _fig2_high_dimension_gate(
                spec, previous, current, previous_curve, current_curve
            )
            convergence_pass = convergence.pass
            merge!(convergence_measured, convergence.measured)
            if schedule.mode == "audited_early_plateau"
                convergence_measured["early_plateau_sha256"] =
                    schedule.evidence_sha256
                convergence_measured["early_plateau_reason"] = schedule.reason
            end
            convergence_reason = convergence.reason
        end
    end
    report["bond_convergence"] = _acceptance_row(
        convergence_pass,
        convergence_measured,
        _fig2_number(spec.data["convergence"], "energy_per_site_change_tol"),
        convergence_evidence,
        convergence_reason,
    )

    root_provenance = _fig2_provenance_result(spec, root)
    checkpoint_restart_pass =
        !isempty(selected) && all(row -> row.restart_valid, selected)
    restart_gate = _fig2_restart_gate_result(
        spec, root, integrity.execution_provenance
    )
    restart_pass = checkpoint_restart_pass && restart_gate.valid
    restart_reason = if !checkpoint_restart_pass
        "at least one selected checkpoint failed its reload/integrity audit"
    elseif !restart_gate.valid
        restart_gate.reason
    else
        "every selected checkpoint reload/integrity audit and the global two-process restart gate passed"
    end
    report["restart"] = _acceptance_row(
        restart_pass,
        Dict{String,Any}(
            "selected_checkpoint_count" => length(selected),
            "reload_integrity_audited_checkpoint_count" =>
                count(row -> row.restart_valid, selected),
            "restart_gate" => restart_gate.measured,
        ),
        "all selected checkpoint reload/integrity audits plus a valid global two-process gate",
        "restart_gate.toml",
        restart_reason,
    )
    provenance_pass = root_provenance.valid && !isempty(selected) &&
        all(row -> row.provenance_valid, selected) &&
        _fig2_file_sha256(joinpath(root, "manifest.toml")) == spec.sha256
    report["provenance"] = _acceptance_row(
        provenance_pass,
        root_provenance.measured,
        "valid fqahc_fig2_provenance_v2",
        "provenance.toml",
        provenance_pass ?
            "root provenance, manifest checksum, and candidate provenance are valid" :
            "provenance validation failed: $(root_provenance.reason)",
    )

    sort!(collect(keys(report))) == sort!(copy(spec.acceptance_rows)) || error(
        "acceptance report rows drifted from the immutable manifest"
    )
    _write_fig2_toml(joinpath(root, "acceptance.toml"), report)
    return report
end

function _default_fig2_candidate_ids(spec, dimension, point, previous_state)
    snapshot = _fig2_validated_snapshot(spec)
    candidates = fig2_initial_candidates(spec.config)
    point == 1 && return [candidate.id for candidate in candidates]
    cold_count = Int(snapshot["pilot"]["cold_candidates"])
    bounded = first(candidates, min(cold_count, length(candidates)))
    return ["warm"; ["cold_$(candidate.id)" for candidate in bounded]]
end

function _default_fig2_load_state(spec, dimension, point, selection, checkpoint)
    _fig2_validated_snapshot(spec)
    config = with_flux(spec.config, selection.phi_y)
    return load_checkpoint(checkpoint, config)
end

function _fig2_checkpoint_restart_valid(
    checkpoint,
    config;
    loader=load_checkpoint,
)
    try
        loader(checkpoint, config)
        return true
    catch error
        error isa InterruptException && rethrow()
        return false
    end
end

function _fig2_checkpoint_audit(
    checkpoint,
    config;
    loader=load_checkpoint,
    maxlinkdim_of=state -> maximum(link_dimensions(state)),
)
    try
        reloaded = loader(checkpoint, config)
        checkpoint_maxlinkdim = maxlinkdim_of(reloaded)
        checkpoint_maxlinkdim isa Integer &&
            !(checkpoint_maxlinkdim isa Bool) && checkpoint_maxlinkdim > 0 ||
            error("checkpoint maxlinkdim must be a positive integer")
        return (;
            restart_valid=true,
            checkpoint_maxlinkdim=Int(checkpoint_maxlinkdim),
        )
    catch error
        error isa InterruptException && rethrow()
        return (; restart_valid=false, checkpoint_maxlinkdim=0)
    end
end

function _default_fig2_persisted_checkpoint_audit(spec, checkpoint, phi_y)
    _fig2_validated_snapshot(spec)
    config = with_flux(spec.config, phi_y)
    return _fig2_checkpoint_audit(checkpoint, config)
end

function _fig2_candidate_pattern(spec, candidate_id)
    id = startswith(candidate_id, "cold_") ? candidate_id[6:end] : candidate_id
    matches = filter(candidate -> candidate.id == id, fig2_initial_candidates(spec.config))
    length(matches) == 1 || throw(ArgumentError("unknown Fig. 2 candidate ID $candidate_id"))
    return only(matches).occupied_sites
end

function _fig2_product_state_on_sites(config, sites, occupied_sites)
    length(sites) == sites_per_cell(config) || throw(
        ArgumentError("reused Fig. 2 site count does not match the reference cell")
    )
    validate_centered_site_charges(config, sites)
    pattern = collect(occupied_sites)
    length(pattern) == particles_per_cell(config) || throw(
        ArgumentError("Fig. 2 product-state particle count mismatch")
    )
    all(site -> site isa Integer && !(site isa Bool), pattern) || throw(
        ArgumentError("Fig. 2 occupied site indices must be integers")
    )
    all(site -> 1 <= site <= sites_per_cell(config), pattern) || throw(
        ArgumentError("Fig. 2 occupied site lies outside the reference cell")
    )
    occupied = Set(Int.(pattern))
    length(occupied) == length(pattern) || throw(
        ArgumentError("Fig. 2 initial pattern contains duplicate occupied sites")
    )
    initstate(n::Integer) = mod1(n, sites_per_cell(config)) in occupied ? "Occ" : "Emp"
    return InfMPS(sites, initstate)
end

function _prepare_fig2_candidate_state(
    spec,
    config,
    candidate_id,
    previous_state,
)
    _fig2_validated_snapshot(spec)
    if candidate_id == "warm"
        isnothing(previous_state) && throw(
            ArgumentError("warm candidate needs previous state")
        )
        return (;
            sites=siteinds(only, previous_state.AL),
            psi=previous_state,
            occupied_sites=default_occupied_sites(config),
        )
    end

    occupied_sites = _fig2_candidate_pattern(spec, candidate_id)
    if isnothing(previous_state)
        sites, _, psi = initial_infinite_mps(config; occupied_sites)
        return (; sites, psi, occupied_sites=collect(occupied_sites))
    end

    sites = siteinds(only, previous_state.AL)
    psi = _fig2_product_state_on_sites(config, sites, occupied_sites)
    siteinds(only, psi.AL) == sites || error(
        "cold Fig. 2 product state did not preserve exact site-index identities"
    )
    return (; sites, psi, occupied_sites=collect(occupied_sites))
end

function _fig2_selection_core_validity(result, expected_pair)
    result.optimization.converged || return false
    energy_per_site = result.energy.per_site
    energy_per_site isa Real && !(energy_per_site isa Bool) &&
        isfinite(energy_per_site) || return false
    fixed_cuts = filter(
        data -> (data.cut_x, data.bond) == expected_pair,
        result.entanglements,
    )
    length(fixed_cuts) == 1 || return false
    fixed_cut = only(fixed_cuts)
    isfinite(fixed_cut.entropy) &&
        isfinite(fixed_cut.raw_schmidt_polarization) || return false
    isempty(fixed_cut.levels) && return false
    isempty(fixed_cut.sectors) && return false
    all(level ->
        (level.cut_x, level.bond) == expected_pair &&
        isfinite(level.singular_value) && level.singular_value >= 0 &&
        isfinite(level.probability) && level.probability >= 0 &&
        (isfinite(level.entanglement_energy) ||
            (iszero(level.probability) && level.entanglement_energy == Inf)) &&
        level.raw_charge isa Integer && !(level.raw_charge isa Bool) &&
        isfinite(level.physical_charge), fixed_cut.levels) || return false
    all(sector ->
        (sector.cut_x, sector.bond) == expected_pair &&
        sector.raw_charge isa Integer && !(sector.raw_charge isa Bool) &&
        isfinite(sector.physical_charge) &&
        isfinite(sector.weight) && sector.weight >= 0,
        fixed_cut.sectors,
    ) || return false
    return true
end

function _fig2_candidate_seed(
    dimension::Integer,
    point::Integer,
    candidate_id::AbstractString,
)
    !(dimension isa Bool) && dimension > 0 || throw(
        ArgumentError("Fig. 2 candidate dimension must be positive")
    )
    !(point isa Bool) && point > 0 || throw(
        ArgumentError("Fig. 2 candidate point must be positive")
    )
    isempty(candidate_id) && throw(
        ArgumentError("Fig. 2 candidate ID must not be empty")
    )
    seed = _derived_seed(0, :fig2_canonicalization, dimension, point, candidate_id)
    return Int(mod(seed, UInt64(typemax(Int))))
end

function _fig2_maxdim_schedule(target::Integer)
    !(target isa Bool) && target > 0 || throw(
        ArgumentError("Fig. 2 maxdim target must be a positive integer")
    )
    target <= typemax(Int) || throw(
        ArgumentError("Fig. 2 maxdim target exceeds the supported Int range")
    )
    final_target = Int(target)
    final_target <= 4 && return [final_target]

    schedule = Int[4]
    while last(schedule) <= fld(final_target, 2)
        push!(schedule, 2 * last(schedule))
    end
    last(schedule) == final_target || push!(schedule, final_target)
    return schedule
end

function _fig2_maxdim_schedule(initial::Integer, target::Integer)
    !(initial isa Bool) && initial > 0 || throw(
        ArgumentError("Fig. 2 initial maxdim must be a positive integer")
    )
    initial <= typemax(Int) || throw(
        ArgumentError("Fig. 2 initial maxdim exceeds the supported Int range")
    )
    schedule = _fig2_maxdim_schedule(target)
    initial_maxdim = Int(initial)
    initial_maxdim <= last(schedule) || throw(
        ArgumentError(
            "Fig. 2 initial maxdim=$initial_maxdim exceeds target maxdim=$(last(schedule))"
        )
    )
    first_stage = findfirst(stage -> stage >= initial_maxdim, schedule)
    isnothing(first_stage) && error("Fig. 2 maxdim schedule has no compatible stage")
    return schedule[first_stage:end]
end

function _default_fig2_run_candidate(
    spec,
    dimension,
    point,
    phi_y,
    candidate_id,
    previous_state,
    candidate_directory;
    runner_override::Union{Nothing,Fig2RunnerOverride}=nothing,
    build_hamiltonian=build_infinite_mpo,
    run_prepared_point=_run_prepared_point,
)
    optimization = _fig2_runner_optimization(spec, runner_override)
    config = with_flux(spec.config, phi_y)
    prepared = _prepare_fig2_candidate_state(
        spec, config, candidate_id, previous_state
    )
    progress = _fig2_load_progress(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
    )
    psi = isnothing(progress) ? prepared.psi : progress.state
    sites = isnothing(progress) ? prepared.sites : siteinds(only, psi.AL)
    occupied_sites = prepared.occupied_sites
    initial_maxdim = maximum(link_dimensions(psi))
    settings = SinglePointSettings(
        config,
        spec.model,
        _fig2_maxdim_schedule(initial_maxdim, dimension),
        optimization.cutoff,
        optimization.vumps_tol,
        optimization.energy_tol,
        optimization.energy_mismatch_tol,
        1e-12,
        optimization.transfer_tol,
        4,
        optimization.max_iterations,
        optimization.stable_iterations,
        optimization.multisite_update_alg,
        Base.Threads.nthreads(),
        candidate_directory,
        joinpath(candidate_directory, "state.h5"),
        nothing,
        collect(occupied_sites),
        _fig2_candidate_seed(dimension, point, candidate_id),
        true,
    )
    H = build_hamiltonian(config, spec.model, sites)
    progress_resume_count = isnothing(progress) ? 0 :
        progress.next_resume_count
    progress_callback = event -> _fig2_persist_progress_event!(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        settings.maxdim_schedule,
        event;
        resume_count=progress_resume_count,
    )
    candidate_operations = SinglePointOperations(
        optimize=(candidate_H, candidate_psi, candidate_settings) ->
            _default_optimize(
                candidate_H,
                candidate_psi,
                candidate_settings;
                progress_callback,
            ),
    )
    result = run_prepared_point(
        settings,
        H,
        psi,
        candidate_operations;
        transfer_rng=Random.Xoshiro(_derived_seed(0, :fig2, dimension, point, candidate_id)),
    )
    _fig2_finalize_progress!(
        spec,
        candidate_directory,
        dimension,
        point,
        phi_y,
        candidate_id,
        joinpath(candidate_directory, "state.h5"),
    )
    _normalize_fig2_convergence_sentinel!(
        joinpath(candidate_directory, "convergence.tsv")
    )
    entanglement = first(result.entanglements)
    sectors = Dict(sector.raw_charge => sector.weight for sector in entanglement.sectors)
    momentum = momentum_entanglement_data(
        result.optimization.psi,
        config;
        cut_x=1,
        residual_tol=optimization.momentum_residual_tol,
        transfer_tol=optimization.transfer_tol,
        rng_seed=_derived_seed(0, :momentum, dimension, point, candidate_id),
    )
    momentum_counting = _fig2_momentum_counting(spec, momentum)
    fidelity = if isnothing(previous_state)
        nothing
    else
        mixed_transfer_fidelity(
            previous_state,
            result.optimization.psi,
            config;
            tol=optimization.transfer_tol,
            rng=Random.Xoshiro(_derived_seed(0, :fidelity, dimension, point, candidate_id)),
        )
    end
    checkpoint_audit = _fig2_checkpoint_audit(
        joinpath(candidate_directory, "state.h5"), config
    )
    achieved_maxlinkdim = maximum(link_dimensions(result.optimization.psi))
    restart_valid = checkpoint_audit.restart_valid &&
        checkpoint_audit.checkpoint_maxlinkdim == achieved_maxlinkdim
    reasons = [result.optimization.reason; result.observable_reasons]
    !isnothing(fidelity) && !fidelity.valid && push!(reasons, fidelity.reason)
    checkpoint_audit.restart_valid && !restart_valid && push!(
        reasons,
        "checkpoint maxlinkdim $(checkpoint_audit.checkpoint_maxlinkdim) " *
        "disagrees with achieved maxlinkdim $achieved_maxlinkdim",
    )
    return Fig2CandidateEvidence(;
        state=result.optimization.psi,
        requested_maxdim=dimension,
        achieved_maxlinkdim,
        checkpoint_maxlinkdim=checkpoint_audit.checkpoint_maxlinkdim,
        converged=result.optimization.converged,
        valid=_fig2_selection_core_validity(
            result, (1, spec.config.Ly)
        ),
        energy_per_site=result.energy.per_site,
        raw_schmidt_polarization=entanglement.raw_schmidt_polarization,
        sector_weights=sectors,
        fidelity_to_previous=isnothing(fidelity) ? NaN : fidelity.fidelity_cell,
        fidelity_valid=!isnothing(fidelity) && fidelity.valid,
        mixed_fidelity=fidelity,
        momentum,
        momentum_counting=momentum_counting.valid ? momentum_counting.counts : Int[],
        momentum_counting_evidence=momentum_counting,
        restart_valid,
        reason=join(reasons, "; "),
    )
end
