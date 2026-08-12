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
    sites = infsiteinds("Fermion", nsites; initstate, conserve_qns=true)
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
)
    _validate_convergence_parameters(
        vumps_tol, energy_tol, energy_mismatch_tol, stable_iterations
    )
    diagnostics = (residual, delta_energy, energy_mismatch)
    return all(value -> isfinite(value) && value >= 0, diagnostics) &&
           residual < vumps_tol &&
           delta_energy < energy_tol &&
           energy_mismatch < energy_mismatch_tol &&
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

struct VUMPSResult
    psi::InfiniteCanonicalMPS
    records::Vector{VUMPSRecord}
    converged::Bool
    reason::String
end

function unit_cell_energy(energy::Real, nsites::Integer)
    nsites > 0 || throw(ArgumentError("site count must be positive"))
    isfinite(energy) || throw(ArgumentError("energy estimate must be finite"))
    return Float64(energy) * nsites
end

function unit_cell_energy(energies::AbstractVector{<:Real}, nsites::Integer)
    nsites > 0 || throw(ArgumentError("site count must be positive"))
    length(energies) == nsites ||
        throw(ArgumentError("energy vector must contain one estimate per MPS site"))
    all(isfinite, energies) ||
        throw(ArgumentError("energy estimates must be finite"))
    return sum(Float64, energies)
end

function vumps_iteration(
    H,
    psi::InfiniteCanonicalMPS;
    vumps_tol::Real,
    solver_tol=(x -> x / 100),
    eager::Bool=true,
)
    isfinite(vumps_tol) && vumps_tol > 0 ||
        throw(ArgumentError("vumps_tol must be finite and positive"))
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
        energy_left=unit_cell_energy(energy_left, n),
        energy_right=unit_cell_energy(energy_right, n),
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
    maximum(after) > maximum(before) ||
        error("subspace expansion made no progress toward target maxdim=$target")
    return expanded
end

function _validate_vumps_schedule(
    maxdim_schedule,
    cutoff::Real,
    max_iterations::Integer,
    vumps_tol::Real,
    energy_tol::Real,
    energy_mismatch_tol::Real,
    stable_iterations::Integer,
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
    return Int.(targets)
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
    )
    records = VUMPSRecord[]
    current = psi

    for (stage, target) in enumerate(targets)
        current = expand_subspace(current, H, target; cutoff)
        previous_energy = nothing
        stable_count = 0
        stage_converged = false

        for iteration in 1:max_iterations
            step = vumps_iteration(
                H,
                current;
                vumps_tol,
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

            stable_now =
                all(isfinite, (precision_error, delta_energy, energy_mismatch)) &&
                precision_error < vumps_tol &&
                delta_energy < energy_tol &&
                energy_mismatch < energy_mismatch_tol
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
            return VUMPSResult(current, records, false, reason)
        end
    end

    count = length(targets)
    suffix = count == 1 ? "stage" : "stages"
    return VUMPSResult(current, records, true, "converged after $count $suffix")
end
