const VUMPS_SOLVER_TOLERANCE_POLICY =
    "max_previous_precision_error_or_vumps_tol_over_100"

function _configuration_signature(c::InfiniteCylinderConfig)
    geometry_code = c.geometry == :legacy_sheared ? "IC$(c.Ly)" : "IP$(c.Ny)"
    signature = "$(geometry_code)x$(c.x_period)n$(c.filling_num)d$(c.filling_den)"
    ncodeunits(signature) <= 16 || throw(
        ArgumentError(
            "cylinder configuration signature exceeds the ITensor 16-character tag limit"
        )
    )
    signature_tags = TagSet(signature)
    length(signature_tags) == 1 || throw(
        ArgumentError("cylinder configuration signature must form exactly one ITensor tag")
    )
    return signature
end

configuration_signature(c::InfiniteCylinderConfig) = _configuration_signature(c)

function _tag_configuration_sites(sites, c::InfiniteCylinderConfig)
    signature = _configuration_signature(c)
    tagged = Index[]
    sizehint!(tagged, length(sites))
    for site in sites
        hastags(site, "Site") ||
            throw(ArgumentError("centered site index is missing the generic Site tag"))
        configured_site = replacetags(site, "Site" => signature)
        hastags(configured_site, signature) && !hastags(configured_site, "Site") ||
            throw(ArgumentError("failed to bind cylinder configuration to site index"))
        push!(tagged, configured_site)
    end
    return CelledVector(tagged, translator(sites))
end

function initial_infinite_mps(
    c::InfiniteCylinderConfig; occupied_sites=default_occupied_sites(c)
)
    nsites = sites_per_cell(c)
    pattern = collect(occupied_sites)
    length(pattern) == particles_per_cell(c) ||
        throw(ArgumentError("initial pattern particle count mismatch"))
    all(site -> site isa Integer && !(site isa Bool), pattern) ||
        throw(ArgumentError("occupied site indices must be integers"))
    all(site -> 1 <= site <= nsites, pattern) ||
        throw(ArgumentError("occupied site lies outside the reference cell"))

    occupied = Set(Int.(pattern))
    length(occupied) == length(pattern) ||
        throw(ArgumentError("initial pattern contains duplicate occupied sites"))
    initstate(n::Integer) = mod1(n, nsites) in occupied ? "Occ" : "Emp"
    centered_sites = infsiteinds("Fermion", nsites; initstate, conserve_qns=true)
    sites = _tag_configuration_sites(centered_sites, c)
    validate_centered_site_charges(c, sites)
    return sites, initstate, InfMPS(sites, initstate)
end

function centered_site_charges(c::InfiniteCylinderConfig)
    nsites = sites_per_cell(c)
    nparticles = particles_per_cell(c)
    common = gcd(nsites, nparticles)
    empty = -(nparticles ÷ common)
    occupied = nsites ÷ common + empty
    return (; empty, occupied)
end

function validate_centered_site_charges(c::InfiniteCylinderConfig, sites)
    length(sites) == sites_per_cell(c) ||
        throw(ArgumentError("site-index count does not match the reference cell"))
    charges = centered_site_charges(c)
    expected_empty = QN("Nf", charges.empty, -1)
    expected_occupied = QN("Nf", charges.occupied, -1)
    for (n, site) in enumerate(sites)
        hasqns(site) || throw(ArgumentError("site index $n does not conserve QNs"))
        flux(state(site, "Emp")) == expected_empty ||
            throw(ArgumentError("site index $n has the wrong centered empty charge"))
        flux(state(site, "Occ")) == expected_occupied ||
            throw(ArgumentError("site index $n has the wrong centered occupied charge"))
    end
    return nothing
end

function link_dimensions(psi::InfiniteCanonicalMPS)
    links = linkinds(only, psi.AL)
    return [dim(links[n]) for n in 1:nsites(psi)]
end

function _validate_convergence_parameters(
    vumps_tol::Real,
    energy_tol::Real,
    energy_mismatch_tol::Real,
    stable_iterations::Integer,
)
    isfinite(vumps_tol) && vumps_tol > 0 ||
        throw(ArgumentError("vumps_tol must be finite and positive"))
    isfinite(energy_tol) && energy_tol > 0 ||
        throw(ArgumentError("energy_tol must be finite and positive"))
    isfinite(energy_mismatch_tol) && energy_mismatch_tol > 0 ||
        throw(ArgumentError("energy_mismatch_tol must be finite and positive"))
    !(stable_iterations isa Bool) && stable_iterations > 0 ||
        throw(ArgumentError("stable_iterations must be positive"))
    return nothing
end

function vumps_converged(
    residual::Real,
    delta_energy::Real,
    energy_mismatch::Real,
    consecutive_iterations::Integer;
    vumps_tol::Real,
    energy_tol::Real,
    energy_mismatch_tol::Real,
    stable_iterations::Integer=2,
    energy_normalization_sites::Integer=1,
)
    _validate_convergence_parameters(
        vumps_tol, energy_tol, energy_mismatch_tol, stable_iterations
    )
    !(energy_normalization_sites isa Bool) && energy_normalization_sites > 0 || throw(
        ArgumentError("energy_normalization_sites must be positive")
    )
    diagnostics = (residual, delta_energy, energy_mismatch)
    return all(value -> isfinite(value) && value >= 0, diagnostics) &&
           residual < vumps_tol &&
           delta_energy / energy_normalization_sites < energy_tol &&
           energy_mismatch / energy_normalization_sites < energy_mismatch_tol &&
           consecutive_iterations >= stable_iterations
end

# Convenience form for callers that have one conservative energy diagnostic.
# It is checked against both energy thresholds.
function vumps_converged(
    residual::Real,
    energy_error::Real,
    consecutive_iterations::Integer;
    kwargs...,
)
    return vumps_converged(
        residual, energy_error, energy_error, consecutive_iterations; kwargs...
    )
end

struct VUMPSRecord
    stage::Int
    iteration::Int
    maxlinkdim::Int
    energy_left::Float64
    energy_right::Float64
    energy_mismatch::Float64
    delta_energy::Float64
    eps_left::Float64
    eps_right::Float64
    precision_error::Float64
    elapsed_seconds::Float64
    converged::Bool
end

struct SubspaceExpansionRecord
    stage::Int
    target::Int
    before::Vector{Int}
    after::Vector{Int}
    progressed::Bool
    elapsed_seconds::Float64

    function SubspaceExpansionRecord(
        stage::Integer,
        target::Integer,
        before::Vector{Int},
        after::Vector{Int},
        progressed::Bool,
        elapsed_seconds::Real,
    )
        return new(
            Int(stage),
            Int(target),
            copy(before),
            copy(after),
            progressed,
            Float64(elapsed_seconds),
        )
    end
end

struct VUMPSProgressEvent
    sequence::Int
    kind::Symbol
    stage::Int
    iteration::Int
    target::Int
    psi::InfiniteCanonicalMPS
    expansion::Union{Nothing,SubspaceExpansionRecord}
    record::Union{Nothing,VUMPSRecord}

    function VUMPSProgressEvent(
        sequence::Integer,
        kind::Symbol,
        stage::Integer,
        iteration::Integer,
        target::Integer,
        psi::InfiniteCanonicalMPS,
        expansion::Union{Nothing,SubspaceExpansionRecord},
        record::Union{Nothing,VUMPSRecord},
    )
        sequence > 0 || throw(ArgumentError("progress event sequence must be positive"))
        stage > 0 || throw(ArgumentError("progress event stage must be positive"))
        iteration > 0 || throw(ArgumentError("progress event iteration must be positive"))
        target > 0 || throw(ArgumentError("progress event target must be positive"))
        kind in (:expansion, :iteration) || throw(
            ArgumentError("progress event kind must be :expansion or :iteration")
        )
        if kind == :expansion
            !isnothing(expansion) && isnothing(record) || throw(
                ArgumentError("expansion progress event has inconsistent records")
            )
        else
            isnothing(expansion) && !isnothing(record) || throw(
                ArgumentError("iteration progress event has inconsistent records")
            )
        end
        return new(
            Int(sequence),
            kind,
            Int(stage),
            Int(iteration),
            Int(target),
            psi,
            expansion,
            record,
        )
    end
end

struct VUMPSProgressContinuation
    psi::InfiniteCanonicalMPS
end

struct VUMPSResult
    psi::InfiniteCanonicalMPS
    records::Vector{VUMPSRecord}
    expansions::Vector{SubspaceExpansionRecord}
    converged::Bool
    reason::String
end

VUMPSResult(psi, records, converged, reason) =
    VUMPSResult(psi, records, SubspaceExpansionRecord[], converged, reason)

function _validate_imaginary_tolerance(imaginary_tol::Real)
    tolerance = Float64(imaginary_tol)
    isfinite(tolerance) && tolerance >= 0 ||
        throw(ArgumentError("imaginary_tol must be finite and nonnegative"))
    return tolerance
end

function _numerically_real_energy(energy::Number, imaginary_tol::Float64)
    isfinite(real(energy)) && isfinite(imag(energy)) ||
        throw(ArgumentError("energy estimate must be finite"))
    abs(imag(energy)) <= imaginary_tol ||
        throw(ArgumentError("energy estimate has imaginary part larger than imaginary_tol"))
    real_energy = Float64(real(energy))
    isfinite(real_energy) ||
        throw(ArgumentError("energy estimate does not fit finite Float64 output"))
    return real_energy
end

const SUPPORTED_MULTISITE_UPDATE_ALGS = (:sequential, :parallel)

function _validate_multisite_update_alg(multisite_update_alg)
    multisite_update_alg isa Symbol &&
        multisite_update_alg in SUPPORTED_MULTISITE_UPDATE_ALGS || throw(
        ArgumentError(
            "multisite_update_alg must be :sequential or :parallel"
        )
    )
    return multisite_update_alg
end

function unit_cell_energy(
    energy::Number,
    nsites::Integer;
    imaginary_tol::Real=1e-12,
)
    nsites > 0 || throw(ArgumentError("site count must be positive"))
    tolerance = _validate_imaginary_tolerance(imaginary_tol)
    cell_energy = _numerically_real_energy(energy, tolerance) * nsites
    isfinite(cell_energy) ||
        throw(ArgumentError("unit-cell energy does not fit finite Float64 output"))
    return cell_energy
end

function unit_cell_energy(
    energies::AbstractVector{<:Number},
    nsites::Integer;
    imaginary_tol::Real=1e-12,
)
    nsites > 0 || throw(ArgumentError("site count must be positive"))
    length(energies) == nsites ||
        throw(ArgumentError("energy vector must contain one estimate per MPS site"))
    tolerance = _validate_imaginary_tolerance(imaginary_tol)
    cell_energy = sum(energy -> _numerically_real_energy(energy, tolerance), energies)
    isfinite(cell_energy) ||
        throw(ArgumentError("unit-cell energy does not fit finite Float64 output"))
    return cell_energy
end

function vumps_iteration(
    H,
    psi::InfiniteCanonicalMPS;
    vumps_tol::Real,
    imaginary_tol::Real=1e-12,
    multisite_update_alg=:sequential,
    solver_tol=(x -> x / 100),
    eager::Bool=true,
)
    isfinite(vumps_tol) && vumps_tol > 0 ||
        throw(ArgumentError("vumps_tol must be finite and positive"))
    tolerance = _validate_imaginary_tolerance(imaginary_tol)
    update_alg = _validate_multisite_update_alg(multisite_update_alg)
    eps_left = fill(Float64(vumps_tol), nsites(psi))
    eps_right = fill(Float64(vumps_tol), nsites(psi))
    result = nothing
    elapsed_seconds = @elapsed result = ITensorInfiniteMPS.tdvp_iteration(
        ITensorInfiniteMPS.vumps_solver,
        H,
        psi;
        (ϵᴸ!)=eps_left,
        (ϵᴿ!)=eps_right,
        time_step=-Inf,
        multisite_update_alg=string(update_alg),
        solver_tol,
        eager,
    )
    newpsi, (energy_left, energy_right) = result
    n = nsites(newpsi)
    return (
        psi=newpsi,
        energy_left=unit_cell_energy(energy_left, n; imaginary_tol=tolerance),
        energy_right=unit_cell_energy(energy_right, n; imaginary_tol=tolerance),
        eps_left,
        eps_right,
        elapsed_seconds,
    )
end

function expand_subspace(
    psi::InfiniteCanonicalMPS,
    H,
    target::Integer;
    cutoff::Real,
)
    target > 0 || throw(ArgumentError("target maxdim must be positive"))
    isfinite(cutoff) && cutoff > 0 ||
        throw(ArgumentError("cutoff must be finite and positive"))
    before = link_dimensions(psi)
    maximum(before) < target || return psi

    expanded = subspace_expansion(psi, H; maxdim=target, cutoff)
    after = link_dimensions(expanded)
    maximum(after) <= target || error(
        "subspace expansion exceeded target maxdim=$target: after=$after"
    )
    _subspace_expansion_progressed(before, after) || error(
        "subspace expansion stalled at achieved maxdim=$(maximum(after)) " *
        "before target maxdim=$target; before=$before after=$after"
    )
    return expanded
end

function _subspace_expansion_progressed(before, after)
    length(before) == length(after) || throw(
        DimensionMismatch("subspace expansion changed the number of MPS links")
    )
    return all(after[index] >= before[index] for index in eachindex(before, after)) &&
        any(after[index] > before[index] for index in eachindex(before, after))
end

function _validate_vumps_schedule(
    maxdim_schedule,
    cutoff::Real,
    max_iterations::Integer,
    vumps_tol::Real,
    energy_tol::Real,
    energy_mismatch_tol::Real,
    stable_iterations::Integer,
    imaginary_tol::Real,
)
    targets = collect(maxdim_schedule)
    isempty(targets) && throw(ArgumentError("maxdim_schedule must not be empty"))
    all(
        target -> target isa Integer && !(target isa Bool) && target > 0,
        targets,
    ) ||
        throw(ArgumentError("maxdim_schedule targets must be positive integers"))
    all(targets[index] < targets[index + 1] for index in 1:(length(targets) - 1)) ||
        throw(ArgumentError("maxdim_schedule must be strictly increasing"))
    isfinite(cutoff) && cutoff > 0 ||
        throw(ArgumentError("cutoff must be finite and positive"))
    !(max_iterations isa Bool) && max_iterations > 0 ||
        throw(ArgumentError("max_iterations must be positive"))
    _validate_convergence_parameters(
        vumps_tol, energy_tol, energy_mismatch_tol, stable_iterations
    )
    _validate_imaginary_tolerance(imaginary_tol)
    return Int.(targets)
end

function _identity_transfer_initial(psi::InfiniteMPS)
    transfer = ITensorInfiniteMPS.TransferMatrix(psi)
    initial = ITensor(dag(ITensorInfiniteMPS.input_inds(transfer)))
    initial_indices = inds(initial)
    length(initial_indices) == 2 || error(
        "canonical transfer fixed point must have two indices"
    )
    dim(initial_indices[1]) == dim(initial_indices[2]) || error(
        "canonical transfer fixed-point indices must have equal dimensions"
    )
    for value in 1:dim(initial_indices[1])
        initial[
            initial_indices[1] => value,
            initial_indices[2] => value,
        ] = 1.0
    end
    normalize!(initial)
    return transfer, initial
end

function _arnoldi_transfer_fixed_point(
    psi::InfiniteMPS;
    tol::Real=1e-12,
    max_iterations::Integer=100,
)
    isfinite(tol) && tol > 0 || throw(
        ArgumentError("fixed-point tolerance must be finite and positive")
    )
    !(max_iterations isa Bool) && max_iterations > 0 || throw(
        ArgumentError("fixed-point maximum iterations must be positive")
    )

    transfer, initial = _identity_transfer_initial(psi)
    _, vectors, values, _ = KrylovKit.schursolve(
        transfer,
        initial,
        1,
        :LM,
        KrylovKit.Arnoldi(;
            tol,
            maxiter=max_iterations,
            eager=false,
        ),
    )
    isempty(values) && error("canonical transfer Arnoldi solve returned no eigenvalue")
    isempty(vectors) && error("canonical transfer Arnoldi solve returned no eigenvector")
    eigenvalue = values[1]
    fixed_point = vectors[1]
    isfinite(eigenvalue) && !iszero(eigenvalue) || error(
        "canonical transfer Arnoldi solve returned an invalid eigenvalue"
    )
    applied = transfer(fixed_point)
    eigen_residual = norm(applied - eigenvalue * fixed_point) /
        max(norm(applied), abs(eigenvalue) * norm(fixed_point), 1.0)
    isfinite(eigen_residual) && eigen_residual <= tol || error(
        "canonical transfer Arnoldi eigenpair failed its residual check"
    )

    overlap = inner(initial, fixed_point)
    isfinite(overlap) && !iszero(overlap) || error(
        "canonical transfer Arnoldi fixed point lost identity overlap"
    )
    fixed_point *= conj(sign(overlap))

    adjoint_fixed_point = swapinds(
        dag(fixed_point), reverse(Pair(inds(fixed_point)...))
    )
    hermitian_residual = norm(fixed_point - adjoint_fixed_point) / norm(fixed_point)
    hermitian_residual <= max(100 * tol, 1e-10) || error(
        "canonical transfer fixed point is not Hermitian"
    )
    fixed_point = (fixed_point + adjoint_fixed_point) / 2
    normalize!(fixed_point)

    checked = transfer(fixed_point)
    eigenvalue = inner(fixed_point, checked) / inner(fixed_point, fixed_point)
    isfinite(eigenvalue) && real(eigenvalue) > 0 || error(
        "canonical transfer fixed-point check produced a nonpositive eigenvalue"
    )
    checked_residual = norm(checked - eigenvalue * fixed_point) /
        max(norm(checked), abs(eigenvalue) * norm(fixed_point), 1.0)
    isfinite(checked_residual) && checked_residual <= tol || error(
        "canonical transfer fixed point failed its eigen-residual check"
    )
    return fixed_point, eigenvalue
end

function _positive_transfer_sqrt(
    fixed_point;
    tol::Real=1e-12,
)
    tol isa Bool && throw(
        ArgumentError("canonical transfer square-root tolerance must be real")
    )
    isfinite(tol) && tol > 0 || throw(
        ArgumentError(
            "canonical transfer square-root tolerance must be finite and positive"
        )
    )

    spectrum, eigenvectors = eigen(fixed_point; ishermitian=true)
    eigenvalues = [spectrum[n, n] for n in 1:mindim(spectrum)]
    isempty(eigenvalues) && error(
        "canonical transfer fixed-point eigendecomposition returned no eigenvalue"
    )
    all(isfinite, eigenvalues) || error(
        "canonical transfer fixed-point eigendecomposition returned a nonfinite eigenvalue"
    )

    spectrum_scale = max(maximum(abs, eigenvalues), 1.0)
    negative_threshold = tol * spectrum_scale
    isfinite(negative_threshold) && negative_threshold > 0 || throw(
        ArgumentError(
            "canonical transfer square-root tolerance produced an invalid negative-eigenvalue threshold"
        )
    )
    minimum_eigenvalue = minimum(eigenvalues)
    minimum_eigenvalue >= -negative_threshold || error(
        "canonical transfer fixed point has a negative eigenvalue: " *
        "minimum=$minimum_eigenvalue, threshold=-$negative_threshold"
    )

    sqrt_spectrum = copy(spectrum)
    for n in eachindex(eigenvalues)
        eigenvalue = eigenvalues[n]
        sqrt_spectrum[n, n] = sqrt(max(eigenvalue, zero(eigenvalue)))
    end
    positive_sqrt = eigenvectors' * sqrt_spectrum * dag(eigenvectors)

    input_indices = inds(fixed_point)
    hassameinds(positive_sqrt, input_indices) || error(
        "canonical transfer square root changed the fixed-point indices"
    )
    positive_sqrt = ITensors.permute(positive_sqrt, input_indices...)
    output_indices = inds(positive_sqrt)
    for (output_index, input_index) in zip(output_indices, input_indices)
        id(output_index) == id(input_index) &&
            tags(output_index) == tags(input_index) &&
            plev(output_index) == plev(input_index) &&
            dir(output_index) == dir(input_index) &&
            space(output_index) == space(input_index) || error(
            "canonical transfer square root changed fixed-point index metadata"
        )
    end
    return positive_sqrt
end

function _canonical_from_left_isometries(psi::InfiniteMPS)
    fixed_point, eigenvalue = _arnoldi_transfer_fixed_point(psi)
    center = _positive_transfer_sqrt(fixed_point)
    center = replacetags(center, ts"" => ts"Right"; plev=1)
    center = noprime(center, ts"Right")
    normalize!(center)
    centers, right, polar_eigenvalue =
        ITensorInfiniteMPS.right_orthogonalize_polar(
            psi,
            center;
            left_tags=ts"",
            right_tags=ts"Right",
        )
    isapprox(polar_eigenvalue, sqrt(real(eigenvalue))) || error(
        "canonical transfer eigenvalue mismatch: $polar_eigenvalue versus $eigenvalue"
    )
    return InfiniteCanonicalMPS(psi, centers, right)
end

function _index_qn_dimensions(index)
    return Dict(qn(index, block) => blockdim(index, block) for block in 1:nblocks(index))
end

function _dual_qn_index(index)
    blocks = [
        -qn(index, block) => blockdim(index, block) for block in 1:nblocks(index)
    ]
    return Index(
        blocks;
        tags=tags(index),
        plev=plev(index),
        dir=dir(dag(index)),
    )
end

function _replace_with_dual_index(tensor, old_index, new_index)
    isomorphism = ITensor(dag(old_index), new_index)
    for value in 1:dim(old_index)
        isomorphism[dag(old_index) => value, new_index => value] = 1.0
    end
    isomorphism_flux = flux(isomorphism)
    (isnothing(isomorphism_flux) || isomorphism_flux == QN()) || error(
        "right-link remapping isomorphism must have neutral QN flux"
    )
    return tensor * isomorphism
end

function _right_link_occurrence(new_base, old_base, occurrence)
    matched = settags(new_base, tags(occurrence))
    matched = setprime(matched, plev(occurrence))
    return dir(occurrence) == dir(old_base) ? matched : dag(matched)
end

function _remap_right_canonical_links(
    psi::InfiniteCanonicalMPS,
    new_index::Function,
)
    old_links = [
        only(commoninds(psi.C[site], psi.AR[site])) for
        site in 1:nsites(psi)
    ]
    new_links = new_index.(old_links)
    centers = copy(psi.C)
    right = copy(psi.AR)

    for site in 1:nsites(psi)
        previous = mod1(site - 1, nsites(psi))
        centers[site] = _replace_with_dual_index(
            centers[site], old_links[site], new_links[site]
        )

        old_left = only(commoninds(psi.AR[site], psi.AR[site - 1]))
        new_left = _right_link_occurrence(
            new_links[previous], old_links[previous], old_left
        )
        new_right = _right_link_occurrence(
            new_links[site], old_links[site], old_links[site]
        )
        tensor = _replace_with_dual_index(
            right[site], old_left, new_left
        )
        right[site] = _replace_with_dual_index(
            tensor, old_links[site], new_right
        )
    end
    return InfiniteCanonicalMPS(psi.AL, centers, right)
end

_dualize_right_canonical_links(psi::InfiniteCanonicalMPS) =
    _remap_right_canonical_links(psi, _dual_qn_index)

_freshen_right_canonical_links(psi::InfiniteCanonicalMPS) =
    _remap_right_canonical_links(psi, sim)

function _normalize_right_link_convention(psi::InfiniteCanonicalMPS)
    center_links = Index[]
    for site in 1:nsites(psi)
        push!(center_links, only(commoninds(psi.C[site], psi.AL[site])))
        push!(center_links, only(commoninds(psi.C[site], psi.AR[site])))
    end
    all(index -> !hasqns(index), center_links) && return (
        _freshen_right_canonical_links(psi)
    )
    all(hasqns, center_links) || error(
        "canonical center links mix QN-conserving and dense index spaces"
    )

    same = Bool[]
    dual = Bool[]
    for site in 1:nsites(psi)
        left_index = only(commoninds(psi.C[site], psi.AL[site]))
        right_index = only(commoninds(psi.C[site], psi.AR[site]))
        left_space = _index_qn_dimensions(left_index)
        right_space = _index_qn_dimensions(right_index)
        push!(same, right_space == left_space)
        push!(
            dual,
            right_space == Dict(
                -charge => dimension for (charge, dimension) in left_space
            ),
        )
    end
    all(dual) && return _freshen_right_canonical_links(psi)
    all(same) || error(
        "canonical right-link QN spaces are neither uniformly equal nor dual to left links"
    )
    return _dualize_right_canonical_links(psi)
end

function _isometry_residual(tensor, retained_index)
    gram = tensor * dag(prime(tensor, retained_index))
    identity = ITensors.denseblocks(delta(inds(gram)...))
    return norm(gram - identity)
end

function _is_strictly_canonical_vumps_state(
    psi::InfiniteCanonicalMPS;
    tolerance::Real=1e-10,
)
    isfinite(tolerance) && tolerance > 0 || throw(
        ArgumentError("canonical-state tolerance must be finite and positive")
    )
    site_count = nsites(psi)
    all(component -> nsites(component) == site_count, (psi.AL, psi.C, psi.AR)) ||
        return false
    all(
        component -> translator(component) === translator(psi.AL),
        (psi.C, psi.AR),
    ) || return false
    siteinds(only, psi.AL) == siteinds(only, psi.AR) || return false
    left_ids = Set(id.(linkinds(only, psi.AL)))
    right_ids = Set(id.(linkinds(only, psi.AR)))
    isempty(intersect(left_ids, right_ids)) || return false

    for site in 1:site_count
        order(psi.AL[site]) == 3 || return false
        order(psi.AR[site]) == 3 || return false
        order(psi.C[site]) == 2 || return false
        left_common = commoninds(psi.AL[site], psi.AL[site + 1])
        right_common = commoninds(psi.AR[site - 1], psi.AR[site])
        length(left_common) == 1 || return false
        length(right_common) == 1 || return false
        left_retained = only(left_common)
        right_retained = only(right_common)
        left_residual = _isometry_residual(psi.AL[site], left_retained)
        right_residual = _isometry_residual(psi.AR[site], right_retained)
        all(isfinite, (left_residual, right_residual)) || return false
        max(left_residual, right_residual) <= tolerance || return false

        left_center = psi.AL[site] * psi.C[site]
        right_center = psi.C[site - 1] * psi.AR[site]
        center_norms = (norm(left_center), norm(right_center), norm(psi.C[site]))
        all(value -> isfinite(value) && value > 0, center_norms) || return false
        center_scale = max(center_norms[1], center_norms[2], 1.0)
        center_residual = norm(left_center - right_center) / center_scale
        isfinite(center_residual) && center_residual <= tolerance || return false
    end
    return true
end

function _canonical_rng_seed(seed::Integer)
    seed isa Bool && throw(ArgumentError("canonical RNG seed must be an integer"))
    0 <= seed <= typemax(UInt64) || throw(
        ArgumentError("canonical RNG seed must be in the UInt64 range")
    )
    return UInt64(seed)
end

function _with_canonical_rng(f::Function, seed::Integer)
    canonical_seed = _canonical_rng_seed(seed)
    default_rng = Random.default_rng()
    index_rng = ITensors.index_id_rng()
    saved_default_rng = copy(default_rng)
    saved_index_rng = copy(index_rng)
    try
        Random.seed!(default_rng, canonical_seed)
        Random.seed!(index_rng, xor(canonical_seed, UInt64(0x49444d5053494458)))
        return f()
    finally
        copy!(default_rng, saved_default_rng)
        copy!(index_rng, saved_index_rng)
    end
end

function _canonicalize_vumps_state(
    psi::InfiniteCanonicalMPS;
    rng_seed::Integer=0,
)
    expected_sites = siteinds(only, psi.AL)
    expected_dims = link_dimensions(psi)
    expected_flux = flux(psi.AL)
    _is_strictly_canonical_vumps_state(psi) && return psi
    canonical = _with_canonical_rng(rng_seed) do
        _normalize_right_link_convention(
            _canonical_from_left_isometries(psi.AL)
        )
    end

    canonical isa InfiniteCanonicalMPS || error(
        "VUMPS canonicalization did not return an InfiniteCanonicalMPS"
    )
    nsites(canonical) == nsites(psi) || error(
        "VUMPS canonicalization changed the reference-cell size"
    )
    siteinds(only, canonical.AL) == expected_sites || error(
        "VUMPS canonicalization changed the site index identities"
    )
    link_dimensions(canonical) == expected_dims || error(
        "VUMPS canonicalization changed the link dimensions"
    )
    flux(canonical.AL) == expected_flux || error(
        "VUMPS canonicalization changed the conserved QN flux"
    )
    left_ids = Set(id.(linkinds(only, canonical.AL)))
    right_ids = Set(id.(linkinds(only, canonical.AR)))
    isempty(intersect(left_ids, right_ids)) || error(
        "VUMPS canonicalization reused a left link ID on the right-canonical chain"
    )
    for site in 1:nsites(canonical)
        left_retained = only(commoninds(
            canonical.AL[site], canonical.AL[site + 1]
        ))
        right_retained = only(commoninds(
            canonical.AR[site - 1], canonical.AR[site]
        ))
        left_residual = _isometry_residual(
            canonical.AL[site], left_retained
        )
        right_residual = _isometry_residual(
            canonical.AR[site], right_retained
        )
        isfinite(left_residual) && left_residual <= 1e-10 || error(
            "VUMPS canonicalization left a nonisometric AL tensor at site $site"
        )
        isfinite(right_residual) && right_residual <= 1e-10 || error(
            "VUMPS canonicalization left a nonisometric AR tensor at site $site"
        )
        left_center = canonical.AL[site] * canonical.C[site]
        right_center = canonical.C[site - 1] * canonical.AR[site]
        scale = max(norm(left_center), norm(right_center), 1.0)
        residual = norm(left_center - right_center) / scale
        isfinite(residual) && residual <= 1e-10 || error(
            "VUMPS canonicalization left a noncanonical center at site $site"
        )
    end
    return canonical
end

function _apply_vumps_progress_continuation(
    current::InfiniteCanonicalMPS,
    response,
    target::Integer,
)
    !(target isa Bool) && target > 0 || throw(
        ArgumentError("progress continuation target must be positive")
    )
    response isa VUMPSProgressContinuation || return current
    replacement = response.psi
    nsites(replacement) == nsites(current) || throw(
        ArgumentError("progress continuation changed the reference-cell size")
    )
    translator(replacement) === translator(current) || throw(
        ArgumentError("progress continuation changed the cell translator")
    )
    siteinds(only, replacement.AL) == siteinds(only, current.AL) || throw(
        ArgumentError("progress continuation changed the site index identities")
    )
    linkinds(only, replacement.AL) == linkinds(only, current.AL) || throw(
        ArgumentError("progress continuation changed the left link identities")
    )
    link_dimensions(replacement) == link_dimensions(current) || throw(
        ArgumentError("progress continuation changed the left link dimensions")
    )
    maximum(link_dimensions(replacement)) <= target || throw(
        ArgumentError("progress continuation exceeds the active maxdim target")
    )
    flux(replacement.AL) == flux(current.AL) || throw(
        ArgumentError("progress continuation changed the conserved QN flux")
    )
    _is_strictly_canonical_vumps_state(replacement) || throw(
        ArgumentError("progress continuation state is not strictly canonical")
    )
    return replacement
end

function run_vumps(
    H,
    psi::InfiniteCanonicalMPS;
    maxdim_schedule=[maximum(link_dimensions(psi))],
    cutoff::Real=1e-8,
    max_iterations::Integer=50,
    vumps_tol::Real=1e-5,
    energy_tol::Real=10 * vumps_tol,
    energy_mismatch_tol::Real=10 * vumps_tol,
    stable_iterations::Integer=2,
    imaginary_tol::Real=1e-12,
    multisite_update_alg=:sequential,
    solver_tol=(x -> x / 100),
    eager::Bool=true,
    canonical_seed::Integer=0,
    progress_callback=nothing,
)
    update_alg = _validate_multisite_update_alg(multisite_update_alg)
    targets = _validate_vumps_schedule(
        maxdim_schedule,
        cutoff,
        max_iterations,
        vumps_tol,
        energy_tol,
        energy_mismatch_tol,
        stable_iterations,
        imaginary_tol,
    )
    _canonical_rng_seed(canonical_seed)
    initial_maxdim = maximum(link_dimensions(psi))
    initial_maxdim <= targets[1] || throw(
        ArgumentError(
            "initial/achieved maxdim=$initial_maxdim exceeds first target=$(targets[1]); " *
            "exact-target stages cannot lower bond dimension"
        )
    )
    records = VUMPSRecord[]
    expansions = SubspaceExpansionRecord[]
    current = psi
    event_sequence = 0

    for (stage, target) in enumerate(targets)
        previous_energy = nothing
        solver_tolerance_seed = Float64(vumps_tol)
        stable_count = 0
        stage_converged = false

        for iteration in 1:max_iterations
            before = link_dimensions(current)
            if maximum(before) < target
                elapsed_seconds = @elapsed current = expand_subspace(
                    current, H, target; cutoff
                )
                after = link_dimensions(current)
                expansion_record = SubspaceExpansionRecord(
                        stage,
                        target,
                        before,
                        after,
                        true,
                        elapsed_seconds,
                    )
                push!(expansions, expansion_record)
                event_sequence += 1
                if !isnothing(progress_callback)
                    response = progress_callback(VUMPSProgressEvent(
                        event_sequence,
                        :expansion,
                        stage,
                        iteration,
                        target,
                        current,
                        expansion_record,
                        nothing,
                    ))
                    current = _apply_vumps_progress_continuation(
                        current, response, target
                    )
                end
                previous_energy = nothing
                solver_tolerance_seed = Float64(vumps_tol)
                stable_count = 0
            end

            effective_solver_tol = _ -> solver_tol(solver_tolerance_seed)
            step = vumps_iteration(
                H,
                current;
                vumps_tol,
                imaginary_tol,
                multisite_update_alg=update_alg,
                solver_tol=effective_solver_tol,
                eager,
            )
            current = step.psi
            energy_left = step.energy_left
            energy_right = step.energy_right
            energy = (energy_left + energy_right) / 2
            energy_mismatch = abs(energy_left - energy_right)
            delta_energy =
                isnothing(previous_energy) ? Inf : abs(energy - previous_energy)
            eps_left = maximum(step.eps_left)
            eps_right = maximum(step.eps_right)
            precision_error = max(eps_left, eps_right)
            energy_normalization_sites = nsites(current)
            current_maxdim = maximum(link_dimensions(current))

            stable_now =
                current_maxdim == target &&
                all(isfinite, (precision_error, delta_energy, energy_mismatch)) &&
                precision_error < vumps_tol &&
                delta_energy / energy_normalization_sites < energy_tol &&
                energy_mismatch / energy_normalization_sites < energy_mismatch_tol
            stable_count = stable_now ? stable_count + 1 : 0
            converged =
                stable_now &&
                vumps_converged(
                    precision_error,
                    delta_energy,
                    energy_mismatch,
                    stable_count;
                    vumps_tol,
                    energy_tol,
                    energy_mismatch_tol,
                    stable_iterations,
                    energy_normalization_sites,
                )
            iteration_record = VUMPSRecord(
                    stage,
                    iteration,
                    current_maxdim,
                    energy_left,
                    energy_right,
                    energy_mismatch,
                    delta_energy,
                    eps_left,
                    eps_right,
                    precision_error,
                    step.elapsed_seconds,
                    converged,
                )
            push!(records, iteration_record)
            event_sequence += 1
            if !isnothing(progress_callback)
                response = progress_callback(VUMPSProgressEvent(
                    event_sequence,
                    :iteration,
                    stage,
                    iteration,
                    target,
                    current,
                    nothing,
                    iteration_record,
                ))
                current = _apply_vumps_progress_continuation(
                    current, response, target
                )
            end
            previous_energy = energy
            solver_tolerance_seed = max(precision_error, Float64(vumps_tol))
            if converged
                stage_converged = true
                break
            end
        end

        if !stage_converged
            achieved = maximum(link_dimensions(current))
            reason = achieved < target ?
                "stage $stage reached maximum iterations ($max_iterations) before target maxdim=$target (achieved maxdim=$achieved)" :
                "stage $stage reached maximum iterations ($max_iterations) without convergence"
            return VUMPSResult(
                _canonicalize_vumps_state(current; rng_seed=canonical_seed),
                records,
                expansions,
                false,
                reason,
            )
        end
    end

    count = length(targets)
    suffix = count == 1 ? "stage" : "stages"
    return VUMPSResult(
        _canonicalize_vumps_state(current; rng_seed=canonical_seed),
        records,
        expansions,
        true,
        "converged after $count $suffix",
    )
end
