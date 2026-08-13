function _normalized_sector_weights(weights::AbstractDict)
    normalized = Dict{Int,Float64}()
    isempty(weights) && throw(ArgumentError("sector weights must not be empty"))
    for (raw_charge, weight) in weights
        raw_charge isa Integer && !(raw_charge isa Bool) || throw(
            ArgumentError("sector keys must be integer raw charges")
        )
        charge = try
            Int(raw_charge)
        catch
            throw(ArgumentError("raw sector charge does not fit Int"))
        end
        weight isa Real || throw(ArgumentError("sector weights must be real"))
        value = try
            Float64(weight)
        catch
            throw(ArgumentError("sector weight does not fit Float64"))
        end
        isfinite(value) && value >= 0 || throw(
            ArgumentError("sector weights must be finite and nonnegative")
        )
        normalized[charge] = get(normalized, charge, 0.0) + value
        isfinite(normalized[charge]) || throw(
            ArgumentError("summed sector weight must be finite")
        )
    end
    total = sum(values(normalized))
    isfinite(total) && total > 0 || throw(
        ArgumentError("sector weights must have finite nonzero total weight")
    )
    for charge in keys(normalized)
        normalized[charge] /= total
    end
    return normalized
end

function _bhattacharyya_distance_normalized(
    first_weights::Dict{Int,Float64},
    second_weights::Dict{Int,Float64},
)
    coefficient = sum(
        sqrt(get(first_weights, charge, 0.0) * get(second_weights, charge, 0.0)) for
        charge in union(keys(first_weights), keys(second_weights))
    )
    distance = 1.0 - coefficient
    tolerance = 64 * eps(Float64)
    -tolerance <= distance <= 1.0 + tolerance || error(
        "normalized Bhattacharyya distance lies outside [0, 1] beyond roundoff"
    )
    return clamp(distance, 0.0, 1.0)
end

"""
    bhattacharyya_distance(first_weights, second_weights)

Normalize two nonnegative raw-charge weight tables and return
`1 - sum(sqrt(p[q] * r[q]))` over the union of their keys.
"""
function bhattacharyya_distance(
    first_weights::AbstractDict,
    second_weights::AbstractDict,
)
    first_normalized = _normalized_sector_weights(first_weights)
    second_normalized = _normalized_sector_weights(second_weights)
    return _bhattacharyya_distance_normalized(first_normalized, second_normalized)
end

"""
    best_sector_shift(first_weights, second_weights)

Find the integer `shift` that best aligns sector weights. The reported convention is
diagnostic only: `shift` is added to the keys of `second_weights` to align them with
`first_weights`. It never modifies either input or any stored raw observable. Equal
distances prefer smaller `abs(shift)`, then the smaller signed shift.
"""
function best_sector_shift(first_weights::AbstractDict, second_weights::AbstractDict)
    first_normalized = _normalized_sector_weights(first_weights)
    second_normalized = _normalized_sector_weights(second_weights)
    shifts = Int[]
    for first_charge in keys(first_normalized), second_charge in keys(second_normalized)
        difference = BigInt(first_charge) - BigInt(second_charge)
        typemin(Int) <= difference <= typemax(Int) || continue
        push!(shifts, Int(difference))
    end
    isempty(shifts) && throw(ArgumentError("no candidate sector shift fits Int"))
    unique!(shifts)
    sort!(shifts; by=shift -> (abs(BigInt(shift)), BigInt(shift)))

    best_shift = first(shifts)
    best_distance = Inf
    for shift in shifts
        shifted = Dict{Int,Float64}()
        valid_shift = true
        for (raw_charge, weight) in second_normalized
            shifted_charge = BigInt(raw_charge) + shift
            if !(typemin(Int) <= shifted_charge <= typemax(Int))
                valid_shift = false
                break
            end
            shifted[Int(shifted_charge)] = weight
        end
        valid_shift || continue
        distance = _bhattacharyya_distance_normalized(first_normalized, shifted)
        if distance < best_distance
            best_distance = distance
            best_shift = shift
        end
    end
    isfinite(best_distance) || throw(ArgumentError("no valid sector alignment shift"))
    return (shift=best_shift, distance=best_distance)
end

function _validated_spectrum_levels(levels)
    validated = NamedTuple{(:raw_charge, :entanglement_energy, :probability),Tuple{Int,Float64,Float64}}[]
    for level in levels
        hasproperty(level, :raw_charge) &&
            hasproperty(level, :entanglement_energy) &&
            hasproperty(level, :probability) || throw(
                ArgumentError(
                    "spectrum levels need raw_charge, entanglement_energy, and probability fields",
                )
            )
        raw_charge = getproperty(level, :raw_charge)
        raw_charge isa Integer && !(raw_charge isa Bool) || throw(
            ArgumentError("spectrum raw charges must be integers")
        )
        charge = try
            Int(raw_charge)
        catch
            throw(ArgumentError("spectrum raw charge does not fit Int"))
        end
        energy = try
            Float64(getproperty(level, :entanglement_energy))
        catch
            throw(ArgumentError("entanglement energy does not fit Float64"))
        end
        probability = try
            Float64(getproperty(level, :probability))
        catch
            throw(ArgumentError("spectrum probability does not fit Float64"))
        end
        isfinite(energy) || throw(
            ArgumentError("entanglement energies must be finite")
        )
        isfinite(probability) && probability >= 0 || throw(
            ArgumentError("spectrum probabilities must be finite and nonnegative")
        )
        push!(
            validated,
            (; raw_charge=charge, entanglement_energy=energy, probability),
        )
    end
    sort!(validated; by=level -> (level.entanglement_energy, level.raw_charge))
    return validated
end

"""
    charge_matched_spectrum_distance(first_levels, second_levels; nlevels=20)

Compare the lowest `nlevels` entanglement energies from each raw spectrum,
matching levels by their integer raw charge and their energy rank within that charge.
The distance is the largest matched absolute energy change. Fewer matches than
requested remain valid and are reported through `matched` and `requested`; zero
matches return an explicit invalid result.
"""
function charge_matched_spectrum_distance(
    first_levels,
    second_levels;
    nlevels::Integer=20,
)
    !(nlevels isa Bool) && nlevels > 0 || throw(
        ArgumentError("nlevels must be a positive integer")
    )
    requested = try
        Int(nlevels)
    catch
        throw(ArgumentError("nlevels does not fit Int"))
    end
    first_lowest = Iterators.take(_validated_spectrum_levels(first_levels), requested)
    second_lowest = Iterators.take(_validated_spectrum_levels(second_levels), requested)
    first_by_charge = Dict{Int,Vector{Float64}}()
    second_by_charge = Dict{Int,Vector{Float64}}()
    for level in first_lowest
        push!(get!(first_by_charge, level.raw_charge, Float64[]), level.entanglement_energy)
    end
    for level in second_lowest
        push!(get!(second_by_charge, level.raw_charge, Float64[]), level.entanglement_energy)
    end

    differences = Float64[]
    for charge in sort!(collect(intersect(keys(first_by_charge), keys(second_by_charge))))
        first_energies = first_by_charge[charge]
        second_energies = second_by_charge[charge]
        count = min(length(first_energies), length(second_energies))
        append!(
            differences,
            abs.(first_energies[1:count] .- second_energies[1:count]),
        )
    end
    matched = length(differences)
    if iszero(matched)
        return (
            distance=NaN,
            matched,
            requested,
            valid=false,
            reason="no charge-matched levels among the lowest $requested levels",
        )
    end
    return (
        distance=maximum(differences),
        matched,
        requested,
        valid=true,
        reason="matched $matched of $requested requested levels by raw charge",
    )
end

function branch_flags(
    ;
    fidelity,
    charge_step,
    sector_distance,
    spectrum_distance,
    energy_crossing,
    fidelity_drop_tol=1e-3,
    charge_jump_tol=0.1,
    sector_tol=0.1,
    spectrum_tol=0.5,
)
    values = Float64[]
    for (name, value) in (
        (:fidelity, fidelity),
        (:charge_step, charge_step),
        (:sector_distance, sector_distance),
        (:spectrum_distance, spectrum_distance),
        (:fidelity_drop_tol, fidelity_drop_tol),
        (:charge_jump_tol, charge_jump_tol),
        (:sector_tol, sector_tol),
        (:spectrum_tol, spectrum_tol),
    )
        value isa Real || throw(ArgumentError("$name must be real"))
        converted = try
            Float64(value)
        catch
            throw(ArgumentError("$name does not fit Float64"))
        end
        isfinite(converted) || throw(ArgumentError("$name must be finite"))
        push!(values, converted)
    end
    fidelity_value,
    charge_value,
    sector_value,
    spectrum_value,
    fidelity_tolerance,
    charge_tolerance,
    sector_tolerance,
    spectrum_tolerance = values
    0 <= fidelity_value <= 1 || throw(ArgumentError("fidelity must lie in [0, 1]"))
    sector_value >= 0 || throw(ArgumentError("sector_distance must be nonnegative"))
    spectrum_value >= 0 || throw(
        ArgumentError("spectrum_distance must be nonnegative")
    )
    all(>=(0.0), values[5:8]) || throw(
        ArgumentError("branch flag tolerances must be nonnegative")
    )
    energy_crossing isa Bool || throw(ArgumentError("energy_crossing must be Bool"))
    return (
        fidelity=(1 - fidelity_value) > fidelity_tolerance,
        charge=abs(charge_value) > charge_tolerance,
        sector=sector_value > sector_tolerance,
        spectrum=spectrum_value > spectrum_tolerance,
        energy_crossing,
    )
end

struct MixedTransferFidelity
    fidelity_cell::Float64
    mixed_value::ComplexF64
    self_value1::ComplexF64
    self_value2::ComplexF64
    mixed_residual::Float64
    self_residual1::Float64
    self_residual2::Float64
    converged::Bool
    valid::Bool
    reason::String
end

function _invalid_mixed_transfer(
    mixed_value,
    self_value1,
    self_value2,
    mixed_residual,
    self_residual1,
    self_residual2,
    reason::AbstractString,
    ;
    converged::Bool=false,
)
    return MixedTransferFidelity(
        NaN,
        ComplexF64(mixed_value),
        ComplexF64(self_value1),
        ComplexF64(self_value2),
        Float64(mixed_residual),
        Float64(self_residual1),
        Float64(self_residual2),
        converged,
        false,
        String(reason),
    )
end

function _mixed_transfer_result(
    mixed_value,
    self_value1,
    self_value2,
    mixed_residual,
    self_residual1,
    self_residual2,
    converged::Bool;
    residual_tolerance::Real,
    clamp_tolerance::Real=64 * eps(Float64),
)
    isfinite(residual_tolerance) && residual_tolerance > 0 || throw(
        ArgumentError("residual_tolerance must be finite and positive")
    )
    isfinite(clamp_tolerance) && clamp_tolerance >= 0 || throw(
        ArgumentError("clamp_tolerance must be finite and nonnegative")
    )
    values = ComplexF64[mixed_value, self_value1, self_value2]
    residuals = Float64[mixed_residual, self_residual1, self_residual2]
    all(isfinite, values) || return _invalid_mixed_transfer(
        values...,
        residuals...,
        "mixed and self transfer dominant values must be finite",
        converged=converged,
    )
    all(residual -> isfinite(residual) && residual >= 0, residuals) || return _invalid_mixed_transfer(
        values...,
        residuals...,
        "mixed and self transfer residuals must be finite and nonnegative",
        converged=converged,
    )
    converged || return _invalid_mixed_transfer(
        values...,
        residuals...,
        "at least one dominant transfer solve did not converge",
    )
    all(<=(Float64(residual_tolerance)), residuals) || return _invalid_mixed_transfer(
        values...,
        residuals...,
        "at least one dominant transfer residual exceeds tolerance",
        converged=converged,
    )
    denominator = sqrt(abs(values[2] * values[3]))
    isfinite(denominator) && denominator > 0 || return _invalid_mixed_transfer(
        values...,
        residuals...,
        "self transfer normalization must be finite and positive",
        converged=converged,
    )
    fidelity = abs(values[1]) / denominator
    isfinite(fidelity) || return _invalid_mixed_transfer(
        values...,
        residuals...,
        "normalized mixed-transfer fidelity must be finite",
        converged=converged,
    )
    tolerance = Float64(clamp_tolerance)
    -tolerance <= fidelity <= 1 + tolerance || return _invalid_mixed_transfer(
        values...,
        residuals...,
        "normalized mixed-transfer fidelity lies outside [0, 1] beyond roundoff",
        converged=converged,
    )
    fidelity = clamp(fidelity, 0.0, 1.0)
    return MixedTransferFidelity(
        fidelity,
        values...,
        residuals...,
        true,
        true,
        "valid normalized mixed-transfer fidelity",
    )
end

function _validate_mixed_transfer_states(
    psi1::InfiniteCanonicalMPS,
    psi2::InfiniteCanonicalMPS,
    c::InfiniteCylinderConfig,
)
    _validate_observable_cell(psi1, c)
    _validate_observable_cell(psi2, c)
    first_sites = siteinds(only, psi1.AL)
    second_sites = siteinds(only, psi2.AL)
    length(first_sites) == length(second_sites) || throw(
        ArgumentError("mixed-transfer states have different cell geometry")
    )
    all(first_sites .== second_sites) || throw(
        ArgumentError(
            "mixed-transfer states must reuse the exact same site indices and configuration signature",
        )
    )
    return nothing
end

function _mixed_transfer_map(psi1::InfiniteCanonicalMPS, psi2::InfiniteCanonicalMPS)
    n = nsites(psi1)
    bra2 = prime(linkinds, dag(psi2.AL))
    ket_cell = psi1.AL[Cell(1)]
    bra_cell = bra2[Cell(1)]
    right_boundary = unioninds(
        linkinds(psi1.AL, n => n + 1),
        linkinds(bra2, n => n + 1),
    )
    left_boundary = unioninds(
        linkinds(psi1.AL, 1 => 0),
        linkinds(bra2, 1 => 0),
    )
    return ITensorMap(
        ket_cell,
        bra_cell;
        input_inds=right_boundary,
        output_inds=left_boundary,
    )
end

function _dominant_transfer_value(
    transfer;
    tol::Float64,
    rng::AbstractRNG=Random.default_rng(),
    kwargs...,
)
    initial = random_itensor(rng, QN(), dag(input_inds(transfer)))
    flux(initial) == QN() || error("neutral transfer initial tensor has nonzero flux")
    values, vectors, info = KrylovKit.eigsolve(
        transfer,
        initial,
        1,
        :LM;
        tol,
        verbosity=0,
        kwargs...,
    )
    isempty(values) && error("dominant transfer solve returned no eigenvalue")
    isempty(vectors) && error("dominant transfer solve returned no eigenvector")
    value = ComplexF64(first(values))
    residual = Float64(norm(transfer(first(vectors)) - value * first(vectors)))
    return (; value, residual, converged=info.converged >= 1)
end

"""
    mixed_transfer_fidelity(psi1, psi2, config; tol=1e-10, kwargs...)

Compute normalized fidelity per complete infinite MPS cell from the dominant
mixed-transfer eigenvalue and both self-transfer normalizations. States must
have identical geometry, configuration signatures, and exact reused site indices.
Eigensolver failures and nonfinite diagnostics return `valid=false`.
"""
function mixed_transfer_fidelity(
    psi1::InfiniteCanonicalMPS,
    psi2::InfiniteCanonicalMPS,
    c::InfiniteCylinderConfig;
    tol::Real=1e-10,
    rng::AbstractRNG=Random.default_rng(),
    kwargs...,
)
    isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be finite and positive"))
    _validate_mixed_transfer_states(psi1, psi2, c)
    tolerance = Float64(tol)
    try
        mixed = _dominant_transfer_value(
            _mixed_transfer_map(psi1, psi2); tol=tolerance, rng, kwargs...
        )
        self1 = _dominant_transfer_value(
            TransferMatrix(psi1.AL); tol=tolerance, rng, kwargs...
        )
        self2 = _dominant_transfer_value(
            TransferMatrix(psi2.AL); tol=tolerance, rng, kwargs...
        )
        return _mixed_transfer_result(
            mixed.value,
            self1.value,
            self2.value,
            mixed.residual,
            self1.residual,
            self2.residual,
            mixed.converged && self1.converged && self2.converged;
            residual_tolerance=tolerance,
        )
    catch exception
        reason = sprint(showerror, exception)
        return _invalid_mixed_transfer(
            ComplexF64(NaN),
            ComplexF64(NaN),
            ComplexF64(NaN),
            NaN,
            NaN,
            NaN,
            "invalid mixed-transfer solve: $reason",
        )
    end
end

struct BranchCandidate
    label::String
    psi::InfiniteCanonicalMPS
    energy::EnergyData
    entanglement::Vector{EntanglementData}
    converged::Bool
    valid::Bool
    reason::String

    function BranchCandidate(
        label::AbstractString,
        psi::InfiniteCanonicalMPS,
        energy::EnergyData,
        entanglement::AbstractVector{<:EntanglementData},
        converged::Bool,
        valid::Bool,
        reason::AbstractString,
    )
        isempty(strip(label)) && throw(ArgumentError("candidate label must not be empty"))
        all(
            isfinite,
            (energy.per_cell, energy.per_x, energy.per_unit_cell, energy.per_site),
        ) || throw(ArgumentError("candidate energy must be finite"))
        isempty(strip(reason)) && throw(ArgumentError("candidate reason must not be empty"))
        return new(
            String(label),
            psi,
            energy,
            entanglement isa Vector{EntanglementData} ?
            entanglement : Vector{EntanglementData}(entanglement),
            converged,
            valid,
            String(reason),
        )
    end
end

struct CandidateSelection
    candidates::Vector{BranchCandidate}
    selected_index::Int
    mode::Symbol
    reason::String
    fidelities::Vector{Union{Nothing,MixedTransferFidelity}}
end

function _selectable_candidate(candidate::BranchCandidate)
    return candidate.converged && candidate.valid && isfinite(candidate.energy.per_site)
end

function _ground_candidate_selection(candidates::Vector{BranchCandidate})
    valid_indices = findall(_selectable_candidate, candidates)
    isempty(valid_indices) && throw(ArgumentError("no converged valid candidate is available"))
    selected_index = first(valid_indices)
    selected_energy = candidates[selected_index].energy.per_site
    for index in Iterators.drop(valid_indices, 1)
        energy = candidates[index].energy.per_site
        if energy < selected_energy
            selected_index = index
            selected_energy = energy
        end
    end
    return CandidateSelection(
        candidates,
        selected_index,
        :ground,
        "selected 1-based candidate $selected_index by minimum valid converged energy per site",
        fill(nothing, length(candidates)),
    )
end

function _adiabatic_candidate_selection(
    candidates::Vector{BranchCandidate},
    previous_state::InfiniteCanonicalMPS,
    c::InfiniteCylinderConfig;
    tol::Real,
    kwargs...,
)
    fidelities = Vector{Union{Nothing,MixedTransferFidelity}}(
        nothing,
        length(candidates),
    )
    selected_index = 0
    selected_fidelity = -Inf
    for (index, candidate) in pairs(candidates)
        _selectable_candidate(candidate) || continue
        fidelity = mixed_transfer_fidelity(
            previous_state,
            candidate.psi,
            c;
            tol,
            kwargs...,
        )
        fidelities[index] = fidelity
        fidelity.valid || continue
        if fidelity.fidelity_cell > selected_fidelity
            selected_index = index
            selected_fidelity = fidelity.fidelity_cell
        end
    end
    iszero(selected_index) && throw(
        ArgumentError("no converged valid candidate has a valid mixed-transfer fidelity")
    )
    return CandidateSelection(
        candidates,
        selected_index,
        :adiabatic,
        "selected 1-based candidate $selected_index by maximum valid mixed-transfer fidelity",
        fidelities,
    )
end

"""
    select_candidate(candidates; mode, previous_state=nothing, config=nothing, kwargs...)

Retain all candidates and return a 1-based selected index. `:ground` chooses the
first deterministic minimum among converged, valid energy densities. `:adiabatic`
requires a previous state and configuration and chooses the first deterministic
maximum valid mixed-transfer fidelity. Candidate entanglement and raw charges are
never shifted, unwrapped, or otherwise modified.
"""
function select_candidate(
    candidates::AbstractVector{<:BranchCandidate};
    mode,
    previous_state=nothing,
    config=nothing,
    tol::Real=1e-10,
    kwargs...,
)
    isempty(candidates) && throw(ArgumentError("candidate list must not be empty"))
    mode isa Symbol || throw(ArgumentError("candidate selection mode must be a Symbol"))
    retained = Vector{BranchCandidate}(candidates)
    if mode === :ground
        return _ground_candidate_selection(retained)
    elseif mode === :adiabatic
        previous_state isa InfiniteCanonicalMPS || throw(
            ArgumentError("adiabatic selection requires a previous infinite state")
        )
        config isa InfiniteCylinderConfig || throw(
            ArgumentError("adiabatic selection requires an infinite-cylinder configuration")
        )
        return _adiabatic_candidate_selection(
            retained,
            previous_state,
            config;
            tol,
            kwargs...,
        )
    end
    throw(ArgumentError("unsupported candidate selection mode: $mode"))
end
