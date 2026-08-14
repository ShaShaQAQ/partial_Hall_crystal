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
    solver_tol=(x -> x / 100),
    eager::Bool=true,
)
    isfinite(vumps_tol) && vumps_tol > 0 ||
        throw(ArgumentError("vumps_tol must be finite and positive"))
    tolerance = _validate_imaginary_tolerance(imaginary_tol)
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
        multisite_update_alg="sequential",
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
    _subspace_expansion_progressed(before, after) ||
        error("subspace expansion made no progress toward target maxdim=$target")
    return expanded
end

function _subspace_expansion_progressed(before, after)
    length(before) == length(after) || throw(
        DimensionMismatch("subspace expansion changed the number of MPS links")
    )
    return any(after[index] > before[index] for index in eachindex(before, after))
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

function _mixed_canonical_from_left_links(psi::InfiniteMPS)
    _, right, _ = ITensorInfiniteMPS.right_orthogonalize(
        psi;
        left_tags=ts"Left",
        right_tags=ts"Right",
    )
    left, centers, eigenvalue = ITensorInfiniteMPS.left_orthogonalize(
        right;
        left_tags=ts"Left",
        right_tags=ts"Right",
    )
    isapprox(eigenvalue, one(eigenvalue)) || error(
        "canonical transfer eigenvalue must be one, got $eigenvalue"
    )
    return InfiniteCanonicalMPS(left, centers, right)
end

function _index_qn_dimensions(index)
    return Dict(qn(index, block) => blockdim(index, block) for block in 1:nblocks(index))
end

function _dual_qn_index(index)
    blocks = [
        -qn(index, block) => blockdim(index, block) for block in 1:nblocks(index)
    ]
    return Index(blocks; tags=tags(index), dir=dir(dag(index)))
end

function _replace_with_dual_index(tensor, old_index, new_index)
    isomorphism = ITensor(dag(old_index), new_index)
    for value in 1:dim(old_index)
        isomorphism[dag(old_index) => value, new_index => value] = 1.0
    end
    flux(isomorphism) == QN() || error(
        "right-link dualization isomorphism must have neutral QN flux"
    )
    return tensor * isomorphism
end

function _right_link_occurrence(new_base, old_base, occurrence)
    tagged = settags(new_base, tags(occurrence))
    return dir(occurrence) == dir(old_base) ? tagged : dag(tagged)
end

function _dualize_right_canonical_links(psi::InfiniteCanonicalMPS)
    centers = copy(psi.C)
    right = copy(psi.AR)
    replacements = Dict{Any,Tuple{Index,Index}}()

    for site in 1:nsites(psi)
        old_index = only(filterinds(centers[site]; tags="Right"))
        replacements[id(old_index)] = (old_index, _dual_qn_index(old_index))
    end

    for site in 1:nsites(psi)
        tensor = centers[site]
        for old_index in filterinds(tensor; tags="Right")
            old_base, new_base = replacements[id(old_index)]
            new_index = _right_link_occurrence(new_base, old_base, old_index)
            tensor = _replace_with_dual_index(tensor, old_index, new_index)
        end
        centers[site] = tensor

        tensor = right[site]
        for old_index in filterinds(tensor; tags="Right")
            old_base, new_base = replacements[id(old_index)]
            new_index = _right_link_occurrence(new_base, old_base, old_index)
            tensor = _replace_with_dual_index(tensor, old_index, new_index)
        end
        right[site] = tensor
    end
    return InfiniteCanonicalMPS(psi.AL, centers, right)
end

function _normalize_right_link_convention(psi::InfiniteCanonicalMPS)
    same = Bool[]
    dual = Bool[]
    for site in 1:nsites(psi)
        left_index = only(filterinds(psi.C[site]; tags="Left"))
        right_index = only(filterinds(psi.C[site]; tags="Right"))
        left_space = _index_qn_dimensions(left_index)
        right_space = _index_qn_dimensions(right_index)
        push!(same, right_space == left_space)
        push!(dual, right_space == Dict(-charge => dimension for (charge, dimension) in left_space))
    end
    all(dual) && return psi
    all(same) || error(
        "canonical right-link QN spaces are neither uniformly equal nor dual to left links"
    )
    return _dualize_right_canonical_links(psi)
end

function _canonicalize_vumps_state(psi::InfiniteCanonicalMPS)
    expected_sites = siteinds(only, psi.AL)
    expected_dims = link_dimensions(psi)
    expected_flux = flux(psi.AL)
    left_links = linkinds(only, psi.AL)
    canonical = if all(index -> hastags(index, "Left"), left_links)
        _mixed_canonical_from_left_links(psi.AL)
    else
        ITensorMPS.orthogonalize(psi.AL, :)
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
    canonical = _normalize_right_link_convention(canonical)

    for site in 1:nsites(canonical)
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
    solver_tol=(x -> x / 100),
    eager::Bool=true,
)
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
    records = VUMPSRecord[]
    expansions = SubspaceExpansionRecord[]
    current = psi

    for (stage, target) in enumerate(targets)
        before = link_dimensions(current)
        if maximum(before) < target
            elapsed_seconds = @elapsed current = expand_subspace(
                current, H, target; cutoff
            )
            after = link_dimensions(current)
            push!(
                expansions,
                SubspaceExpansionRecord(
                    stage,
                    target,
                    before,
                    after,
                    _subspace_expansion_progressed(before, after),
                    elapsed_seconds,
                ),
            )
        end
        previous_energy = nothing
        stable_count = 0
        stage_converged = false

        for iteration in 1:max_iterations
            step = vumps_iteration(
                H,
                current;
                vumps_tol,
                imaginary_tol,
                solver_tol,
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

            stable_now =
                all(isfinite, (precision_error, delta_energy, energy_mismatch)) &&
                precision_error < vumps_tol &&
                delta_energy / energy_normalization_sites < energy_tol &&
                energy_mismatch / energy_normalization_sites < energy_mismatch_tol
            stable_count = stable_now ? stable_count + 1 : 0
            converged = vumps_converged(
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
            push!(
                records,
                VUMPSRecord(
                    stage,
                    iteration,
                    maximum(link_dimensions(current)),
                    energy_left,
                    energy_right,
                    energy_mismatch,
                    delta_energy,
                    eps_left,
                    eps_right,
                    precision_error,
                    step.elapsed_seconds,
                    converged,
                ),
            )
            previous_energy = energy
            if converged
                stage_converged = true
                break
            end
        end

        if !stage_converged
            reason =
                "stage $stage reached maximum iterations ($max_iterations) without convergence"
            return VUMPSResult(
                _canonicalize_vumps_state(current), records, expansions, false, reason
            )
        end
    end

    count = length(targets)
    suffix = count == 1 ? "stage" : "stages"
    return VUMPSResult(
        _canonicalize_vumps_state(current),
        records,
        expansions,
        true,
        "converged after $count $suffix",
    )
end
