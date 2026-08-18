struct MomentumValidation
    unitary_residual::Float64
    power_residual::Float64
    commutator_residual::Float64
    transfer_residual::Float64
    valid::Bool
    reason::String
end

struct MomentumEntanglementLevel
    level::Int
    probability::Float64
    entanglement_energy::Float64
    raw_charge::Int
    ky_index::Union{Missing,Int}
    ky::Union{Missing,Float64}
    eigenphase::Union{Missing,Float64}
    phase_residual::Float64
end

struct MomentumEntanglementData
    levels::Vector{MomentumEntanglementLevel}
    validation::MomentumValidation
    valid::Bool
    reason::String
    momentum_convention::Symbol
    sector_ky_offsets::Dict{Int,Int}
end

# `sector_ky_offsets` fixes the independent momentum origin of each charge
# sector at one parameter point. It is gauge-dependent metadata and must not
# be compared across flux points as a physical momentum flow.

struct MomentumConstructionError <: Exception
    message::String
end

Base.showerror(io::IO, error::MomentumConstructionError) =
    print(io, error.message)

const DEFAULT_MOMENTUM_RNG_SEED::UInt64 = 0x4d4f4d454e54554d

function _resolve_momentum_rng(rng, rng_seed)
    !isnothing(rng) && !isnothing(rng_seed) && throw(
        ArgumentError("rng and rng_seed are mutually exclusive")
    )
    if !isnothing(rng)
        rng isa AbstractRNG || throw(
            ArgumentError("rng must be an AbstractRNG or nothing")
        )
        return rng
    end

    seed = isnothing(rng_seed) ? DEFAULT_MOMENTUM_RNG_SEED : rng_seed
    seed isa Integer && !(seed isa Bool) || throw(
        ArgumentError("rng_seed must be an integer in the UInt64 range")
    )
    zero(seed) <= seed <= typemax(UInt64) || throw(
        ArgumentError("rng_seed must be an integer in the UInt64 range")
    )
    return Random.Xoshiro(UInt64(seed))
end

function _momentum_parameters(Ny, residual_tol, degeneracy_tol)
    Ny isa Integer && !(Ny isa Bool) && Ny > 0 || throw(
        ArgumentError("Ny must be a positive integer")
    )
    isfinite(residual_tol) && residual_tol > 0 || throw(
        ArgumentError("residual_tol must be finite and positive")
    )
    isfinite(degeneracy_tol) && degeneracy_tol > 0 || throw(
        ArgumentError("degeneracy_tol must be finite and positive")
    )
    return Int(Ny), Float64(residual_tol), Float64(degeneracy_tol)
end

function _probability_threshold(probability_tol)
    probability_tol isa Real && !(probability_tol isa Bool) &&
        isfinite(probability_tol) && probability_tol >= 0 || throw(
        ArgumentError("probability_tol must be finite and nonnegative")
    )
    return Float64(probability_tol)
end

function _momentum_inputs(rho, Ty, raw_charges, tolerance)
    ndims(rho) == 2 && size(rho, 1) == size(rho, 2) || throw(
        DimensionMismatch("Schmidt density matrix must be square")
    )
    ndims(Ty) == 2 && size(Ty, 1) == size(Ty, 2) || throw(
        DimensionMismatch("virtual translation matrix must be square")
    )
    size(rho) == size(Ty) || throw(
        DimensionMismatch("density and translation matrices must have equal size")
    )
    dimension = size(rho, 1)
    dimension > 0 || throw(ArgumentError("Schmidt space must not be empty"))
    charges = isnothing(raw_charges) ? zeros(Int, dimension) : collect(raw_charges)
    length(charges) == dimension || throw(
        DimensionMismatch("raw charge count must match the Schmidt dimension")
    )
    all(charge -> charge isa Integer && !(charge isa Bool), charges) || throw(
        ArgumentError("raw charges must be integers")
    )
    charges = Int.(charges)

    density = Matrix{ComplexF64}(rho)
    translation = Matrix{ComplexF64}(Ty)
    all(isfinite, density) || throw(
        ArgumentError("Schmidt density matrix entries must be finite")
    )
    all(isfinite, translation) || throw(
        ArgumentError("virtual translation matrix entries must be finite")
    )
    density_scale = max(LinearAlgebra.norm(density), 1.0)
    LinearAlgebra.norm(density - density') / density_scale <= tolerance || throw(
        ArgumentError("Schmidt density matrix must be Hermitian")
    )
    density = (density + density') / 2
    trace_density = Float64(real(LinearAlgebra.tr(density)))
    isfinite(trace_density) && trace_density > 0 || throw(
        ArgumentError("Schmidt density matrix must have positive finite trace")
    )
    density ./= trace_density
    return density, translation, charges
end

function _charge_density_basis(density, charges, tolerance)
    dimension = size(density, 1)
    basis = zeros(ComplexF64, dimension, dimension)
    probabilities = Float64[]
    basis_charges = Int[]
    column = 1
    for charge in sort!(unique(charges))
        indices = findall(==(charge), charges)
        block = LinearAlgebra.Hermitian(density[indices, indices])
        decomposition = LinearAlgebra.eigen(block)
        order = sortperm(real.(decomposition.values); rev=true)
        values = Float64.(real.(decomposition.values[order]))
        minimum(values) >= -tolerance || throw(
            ArgumentError("Schmidt density matrix must be positive semidefinite")
        )
        values = max.(values, 0.0)
        count = length(indices)
        basis[indices, column:(column + count - 1)] = decomposition.vectors[:, order]
        append!(probabilities, values)
        append!(basis_charges, fill(charge, count))
        column += count
    end
    abs(sum(probabilities) - 1.0) <= 10tolerance || throw(
        ArgumentError("Schmidt probabilities are not normalized")
    )
    return probabilities, basis_charges, basis
end

function _momentum_levels(probabilities, charges; labels=nothing)
    levels = MomentumEntanglementLevel[]
    for level in eachindex(probabilities)
        probability = probabilities[level]
        if isnothing(labels)
            ky_index = missing
            ky = missing
            eigenphase = missing
            phase_residual = NaN
        else
            label = labels[level]
            ky_index = label.ky_index
            ky = label.ky
            eigenphase = label.eigenphase
            phase_residual = label.phase_residual
        end
        push!(
            levels,
            MomentumEntanglementLevel(
                level,
                probability,
                iszero(probability) ? Inf : -log(probability),
                charges[level],
                ky_index,
                ky,
                eigenphase,
                phase_residual,
            ),
        )
    end
    return levels
end

function _invalid_momentum_data(probabilities, charges, validation)
    return MomentumEntanglementData(
        _momentum_levels(probabilities, charges),
        validation,
        false,
        validation.reason,
        :invalid,
        Dict{Int,Int}(),
    )
end

function _relative_matrix_residual(value, scale)
    return Float64(LinearAlgebra.norm(value) / max(Float64(scale), eps(Float64)))
end

function _degenerate_momentum_groups(probabilities, charges, tolerance)
    groups = Vector{Vector{Int}}()
    first = 1
    while first <= length(probabilities)
        group = Int[first]
        next = first + 1
        while next <= length(probabilities) &&
                charges[next] == charges[first] &&
                abs(probabilities[next] - probabilities[first]) <= tolerance
            push!(group, next)
            next += 1
        end
        push!(groups, group)
        first = next
    end
    return groups
end

function _quantized_momentum_label(value, Ny)
    magnitude = abs(value)
    isfinite(magnitude) && magnitude > 0 || return (
        ky_index=0,
        ky=0.0,
        eigenphase=NaN,
        phase_residual=Inf,
    )
    normalized = value / magnitude
    phase = mod(angle(normalized), 2pi)
    momentum = mod(round(Int, phase * Ny / (2pi)), Ny)
    target = cis(2pi * momentum / Ny)
    return (
        ky_index=momentum,
        ky=2pi * momentum / Ny,
        eigenphase=phase,
        phase_residual=Float64(abs(normalized - target)),
    )
end

function _lexicographically_less(first, second)
    for (left, right) in zip(first, second)
        left == right || return left < right
    end
    return length(first) < length(second)
end

function _sector_momentum_signature(labels, groups, offset, Ny)
    signature = Int[]
    for group in groups
        shifted = sort!(
            [mod(labels[level].ky_index + offset, Ny) for level in group]
        )
        append!(signature, shifted)
    end
    return signature
end

function _shift_momentum_label(label, offset, Ny)
    momentum = mod(label.ky_index + offset, Ny)
    return (;
        ky_index=momentum,
        ky=2pi * momentum / Ny,
        eigenphase=label.eigenphase,
        phase_residual=label.phase_residual,
    )
end

function _canonical_sector_relative_labels!(labels, groups, charges, Ny)
    sector_offsets = Dict{Int,Int}()
    for charge in sort!(unique(charges))
        sector_groups = [
            group for group in groups if charges[first(group)] == charge
        ]
        best_offset = 0
        best_signature = _sector_momentum_signature(
            labels, sector_groups, best_offset, Ny
        )
        for offset in 1:(Ny - 1)
            signature = _sector_momentum_signature(
                labels, sector_groups, offset, Ny
            )
            if _lexicographically_less(signature, best_signature)
                best_offset = offset
                best_signature = signature
            end
        end
        sector_offsets[charge] = best_offset

        for group in sector_groups
            shifted_labels = [
                _shift_momentum_label(labels[level], best_offset, Ny) for
                level in group
            ]
            sort!(shifted_labels; by=label -> (label.ky_index, label.eigenphase))
            for (level, label) in zip(group, shifted_labels)
                labels[level] = label
            end
        end
    end
    return sector_offsets
end

function _block_momentum_inputs(
    density_blocks,
    translation_blocks,
    raw_charges,
    tolerance;
    prepared_density=false,
)
    densities = collect(density_blocks)
    translations = collect(translation_blocks)
    isempty(densities) && throw(ArgumentError("Schmidt space must not be empty"))
    length(densities) == length(translations) || throw(
        DimensionMismatch("density and translation block counts must match")
    )
    charges = isnothing(raw_charges) ? zeros(Int, length(densities)) : collect(raw_charges)
    length(charges) == length(densities) || throw(
        DimensionMismatch("raw charge count must match the QN block count")
    )
    all(charge -> charge isa Integer && !(charge isa Bool), charges) || throw(
        ArgumentError("raw charges must be integers")
    )

    density_matrices = Matrix{ComplexF64}[]
    translation_matrices = Matrix{ComplexF64}[]
    for block in eachindex(densities)
        density = densities[block]
        translation = translations[block]
        ndims(density) == 2 && size(density, 1) == size(density, 2) || throw(
            DimensionMismatch("Schmidt density block $block must be square")
        )
        ndims(translation) == 2 && size(translation, 1) == size(translation, 2) || throw(
            DimensionMismatch("virtual translation block $block must be square")
        )
        size(density) == size(translation) || throw(
            DimensionMismatch("density and translation block $block must have equal size")
        )
        size(density, 1) > 0 || throw(
            ArgumentError("Schmidt density blocks must not be empty")
        )
        density_matrix = Matrix{ComplexF64}(density)
        translation_matrix = Matrix{ComplexF64}(translation)
        all(isfinite, density_matrix) || throw(
            ArgumentError("Schmidt density block $block entries must be finite")
        )
        all(isfinite, translation_matrix) || throw(
            ArgumentError("virtual translation block $block entries must be finite")
        )
        if !prepared_density
            density_scale = max(LinearAlgebra.norm(density_matrix), 1.0)
            LinearAlgebra.norm(density_matrix - density_matrix') / density_scale <= tolerance || throw(
                ArgumentError("Schmidt density block $block must be Hermitian")
            )
            density_matrix = (density_matrix + density_matrix') / 2
        end
        push!(density_matrices, density_matrix)
        push!(translation_matrices, translation_matrix)
    end

    trace_density = Float64(sum(real(LinearAlgebra.tr(block)) for block in density_matrices))
    isfinite(trace_density) && trace_density > 0 || throw(
        ArgumentError("Schmidt density matrix must have positive finite trace")
    )
    if !prepared_density
        for block in density_matrices
            block ./= trace_density
        end
    end
    return density_matrices, translation_matrices, Int.(charges)
end

function _block_density_basis(density_blocks, block_charges, tolerance)
    block_probabilities = Vector{Vector{Float64}}(undef, length(density_blocks))
    basis_blocks = Vector{Matrix{ComplexF64}}(undef, length(density_blocks))
    for block in eachindex(density_blocks)
        decomposition = LinearAlgebra.eigen(
            LinearAlgebra.Hermitian(density_blocks[block])
        )
        order = sortperm(real.(decomposition.values); rev=true)
        values = Float64.(real.(decomposition.values[order]))
        minimum(values) >= -tolerance || throw(
            ArgumentError("Schmidt density matrix must be positive semidefinite")
        )
        block_probabilities[block] = max.(values, 0.0)
        basis_blocks[block] = Matrix{ComplexF64}(decomposition.vectors[:, order])
    end

    probabilities = Float64[]
    basis_charges = Int[]
    locations = Tuple{Int,Int}[]
    for charge in sort!(unique(block_charges))
        records = NamedTuple[]
        for block in eachindex(density_blocks)
            block_charges[block] == charge || continue
            for local_level in eachindex(block_probabilities[block])
                push!(records, (;
                    probability=block_probabilities[block][local_level],
                    block,
                    local_level,
                ))
            end
        end
        sort!(records; by=record -> record.probability, rev=true)
        for record in records
            push!(probabilities, record.probability)
            push!(basis_charges, charge)
            push!(locations, (record.block, record.local_level))
        end
    end
    abs(sum(probabilities) - 1.0) <= 10tolerance || throw(
        ArgumentError("Schmidt probabilities are not normalized")
    )
    return probabilities, basis_charges, basis_blocks, locations
end

function _block_frobenius_norm(blocks)
    squared_norm = 0.0
    for block in blocks
        squared_norm += sum(abs2, block)
    end
    return Float64(sqrt(squared_norm))
end

function _block_validation_residuals(density_blocks, translation_blocks, circumference)
    total_dimension = sum(block -> size(block, 1), density_blocks)
    unitary_squared = 0.0
    power_squared = 0.0
    commutator_squared = 0.0
    for (density, translation) in zip(density_blocks, translation_blocks)
        block_dimension = size(density, 1)
        identity_block = Matrix{ComplexF64}(
            LinearAlgebra.I, block_dimension, block_dimension
        )
        unitary_squared += sum(abs2, translation' * translation - identity_block)
        power_squared += sum(abs2, translation^circumference - identity_block)
        commutator_squared += sum(abs2, density * translation - translation * density)
    end
    identity_norm = sqrt(Float64(total_dimension))
    density_norm = _block_frobenius_norm(density_blocks)
    translation_norm = _block_frobenius_norm(translation_blocks)
    return (;
        unitary=Float64(sqrt(unitary_squared) / identity_norm),
        power=Float64(sqrt(power_squared) / identity_norm),
        density_translation=Float64(
            sqrt(commutator_squared) /
            max(density_norm * translation_norm, 1.0)
        ),
        charge_density=0.0,
        charge_translation=0.0,
    )
end

function _block_momentum_labels(
    probabilities,
    charges,
    locations,
    translation_blocks,
    circumference,
    degeneracy,
    tolerance,
)
    labels = Vector{NamedTuple}(undef, length(probabilities))
    phase_failure = false
    groups = _degenerate_momentum_groups(probabilities, charges, degeneracy)
    for group in groups
        block_values = ComplexF64[]
        contributing_blocks = unique([
            first(locations[level]) for level in group
        ])
        for block in contributing_blocks
            local_levels = [
                last(locations[level]) for level in group if
                first(locations[level]) == block
            ]
            append!(
                block_values,
                LinearAlgebra.eigen(
                    translation_blocks[block][local_levels, local_levels]
                ).values,
            )
        end
        length(block_values) == length(group) || error(
            "momentum block locations do not cover a degenerate Schmidt group"
        )
        block_labels = [
            _quantized_momentum_label(value, circumference) for
            value in block_values
        ]
        sort!(block_labels; by=label -> (label.ky_index, label.eigenphase))
        for (level, label) in zip(group, block_labels)
            labels[level] = label
            phase_failure |= label.phase_residual > tolerance
        end
    end
    return labels, groups, phase_failure
end

function _validate_momentum_block_vectors(
    density_blocks,
    translation_blocks;
    Ny,
    raw_charges=nothing,
    residual_tol::Real=1e-10,
    degeneracy_tol::Real=residual_tol,
    transfer_residual::Real=0.0,
    probability_tol::Real=0.0,
    _prepared_density::Bool=false,
    _dense_residuals=nothing,
)
    circumference, tolerance, degeneracy = _momentum_parameters(
        Ny, residual_tol, degeneracy_tol
    )
    probability_threshold = _probability_threshold(probability_tol)
    transfer_residual isa Real && !isnan(transfer_residual) && transfer_residual >= 0 || throw(
        ArgumentError("transfer_residual must be nonnegative and not NaN")
    )
    densities, translations, block_charges = _block_momentum_inputs(
        density_blocks,
        translation_blocks,
        raw_charges,
        tolerance;
        prepared_density=_prepared_density,
    )
    probabilities, basis_charges, density_bases, locations =
        _block_density_basis(densities, block_charges, tolerance)
    density_basis_translations = [
        density_bases[block]' * translations[block] * density_bases[block] for
        block in eachindex(translations)
    ]
    block_residuals = _block_validation_residuals(
        densities, translations, circumference
    )
    residuals = isnothing(_dense_residuals) ? block_residuals : _dense_residuals
    retained = findall(probability -> probability > probability_threshold, probabilities)
    isempty(retained) && throw(
        ArgumentError(
            "probability_tol removed all positive Schmidt support"
        )
    )
    physical_probabilities = probabilities[retained]
    physical_charges = basis_charges[retained]
    physical_locations = locations[retained]
    unitary_residual = Float64(residuals.unitary)
    power_residual = Float64(residuals.power)
    charge_density_commutator = Float64(residuals.charge_density)
    density_translation_commutator = Float64(residuals.density_translation)
    charge_translation_commutator = Float64(residuals.charge_translation)
    commutator_residual = max(
        charge_density_commutator,
        density_translation_commutator,
        charge_translation_commutator,
    )

    failures = String[]
    unitary_residual <= tolerance || push!(failures, "unitarity residual exceeds tolerance")
    power_residual <= tolerance || push!(failures, "Ty^Ny residual exceeds tolerance")
    charge_density_commutator <= tolerance || push!(
        failures, "charge-density commutator residual exceeds tolerance"
    )
    density_translation_commutator <= tolerance || push!(
        failures, "Schmidt-density translation commutator residual exceeds tolerance"
    )
    charge_translation_commutator <= tolerance || push!(
        failures, "charge-translation commutator residual exceeds tolerance"
    )
    transfer_residual <= tolerance || push!(failures, "twisted-transfer residual exceeds tolerance")
    if !isempty(failures)
        reason = join(failures, "; ")
        validation = MomentumValidation(
            unitary_residual,
            power_residual,
            commutator_residual,
            Float64(transfer_residual),
            false,
            reason,
        )
        return _invalid_momentum_data(
            physical_probabilities, physical_charges, validation
        )
    end

    labels, groups, phase_failure = _block_momentum_labels(
        physical_probabilities,
        physical_charges,
        physical_locations,
        density_basis_translations,
        circumference,
        degeneracy,
        tolerance,
    )
    if phase_failure
        reason = "at least one virtual-translation phase is not quantized to 2pi*m/Ny"
        validation = MomentumValidation(
            unitary_residual,
            power_residual,
            commutator_residual,
            Float64(transfer_residual),
            false,
            reason,
        )
        return _invalid_momentum_data(
            physical_probabilities, physical_charges, validation
        )
    end

    sector_ky_offsets = _canonical_sector_relative_labels!(
        labels, groups, physical_charges, circumference
    )

    validation = MomentumValidation(
        unitary_residual,
        power_residual,
        commutator_residual,
        Float64(transfer_residual),
        true,
        "valid charge-resolved Schmidt spectrum with sector-relative Δky labels",
    )
    return MomentumEntanglementData(
        _momentum_levels(physical_probabilities, physical_charges; labels),
        validation,
        true,
        validation.reason,
        :sector_relative_canonical_cyclic_orbit,
        sector_ky_offsets,
    )
end

function validate_momentum_blocks(
    rho,
    Ty;
    Ny,
    raw_charges=nothing,
    residual_tol::Real=1e-10,
    degeneracy_tol::Real=residual_tol,
    transfer_residual::Real=0.0,
    probability_tol::Real=0.0,
)
    circumference, tolerance, _ = _momentum_parameters(
        Ny, residual_tol, degeneracy_tol
    )
    _probability_threshold(probability_tol)
    transfer_residual isa Real && !isnan(transfer_residual) && transfer_residual >= 0 || throw(
        ArgumentError("transfer_residual must be nonnegative and not NaN")
    )
    density, translation, charges = _momentum_inputs(
        rho, Ty, raw_charges, tolerance
    )
    raw_charge_diagonal = LinearAlgebra.Diagonal(Float64.(charges))
    identity_matrix = Matrix{ComplexF64}(
        LinearAlgebra.I, size(density, 1), size(density, 1)
    )
    dense_residuals = (;
        unitary=_relative_matrix_residual(
            translation' * translation - identity_matrix,
            LinearAlgebra.norm(identity_matrix),
        ),
        power=_relative_matrix_residual(
            translation^circumference - identity_matrix,
            LinearAlgebra.norm(identity_matrix),
        ),
        charge_density=_relative_matrix_residual(
            raw_charge_diagonal * density - density * raw_charge_diagonal,
            max(
                LinearAlgebra.norm(raw_charge_diagonal) * LinearAlgebra.norm(density),
                1.0,
            ),
        ),
        density_translation=_relative_matrix_residual(
            density * translation - translation * density,
            max(
                LinearAlgebra.norm(density) * LinearAlgebra.norm(translation),
                1.0,
            ),
        ),
        charge_translation=_relative_matrix_residual(
            raw_charge_diagonal * translation - translation * raw_charge_diagonal,
            max(
                LinearAlgebra.norm(raw_charge_diagonal) * LinearAlgebra.norm(translation),
                1.0,
            ),
        ),
    )
    charge_values = sort!(unique(charges))
    charge_indices = [findall(==(charge), charges) for charge in charge_values]
    density_blocks = [density[indices, indices] for indices in charge_indices]
    translation_blocks = [translation[indices, indices] for indices in charge_indices]
    return _validate_momentum_block_vectors(
        density_blocks,
        translation_blocks;
        Ny=circumference,
        raw_charges=charge_values,
        residual_tol=tolerance,
        degeneracy_tol,
        transfer_residual,
        probability_tol,
        _prepared_density=true,
        _dense_residuals=dense_residuals,
    )
end

function _state_value(site_index, state_name)
    basis_state = state(site_index, state_name)
    values = findall(value -> !iszero(value), [basis_state[site_index => n] for n in 1:dim(site_index)])
    length(values) == 1 || error("$state_name is not a unique local basis state")
    return only(values)
end

function _physical_fermion_swap_phase(first_site, second_site)
    first_occupied = _state_value(first_site, "Occ")
    second_occupied = _state_value(second_site, "Occ")
    gate = ITensor(
        QN(),
        prime(first_site),
        prime(second_site),
        dag(first_site),
        dag(second_site),
    )
    for first_value in 1:dim(first_site), second_value in 1:dim(second_site)
        sign = first_value == first_occupied && second_value == second_occupied ? -1.0 : 1.0
        gate[
            prime(first_site) => first_value,
            prime(second_site) => second_value,
            dag(first_site) => first_value,
            dag(second_site) => second_value,
        ] = sign
    end
    return gate
end

function _physical_site_index(psi::MPS, site::Int, configuration_tag)
    physical = filter(index -> hastags(index, configuration_tag), inds(psi[site]))
    length(physical) == 1 || throw(ArgumentError(
        "finite-cell MPS site $site must have exactly one configuration-tagged physical index"
    ))
    return only(physical)
end

function _validate_finite_cell_topology(psi::MPS, configuration_tag)
    count = length(psi)
    count >= 2 || throw(ArgumentError(
        "finite-cell MPS topology validation requires at least two sites"
    ))
    for site in eachindex(psi)
        _physical_site_index(psi, site, configuration_tag)
    end

    internal_links = Vector{Index}(undef, count - 1)
    for bond in 1:(count - 1)
        shared = commoninds(psi[bond], psi[bond + 1])
        length(shared) == 1 || throw(ArgumentError(
            "finite-cell MPS bond $bond must have exactly one shared index"
        ))
        link = only(shared)
        !hastags(link, configuration_tag) || throw(ArgumentError(
            "finite-cell MPS bond $bond must share a nonphysical link index"
        ))
        internal_links[bond] = link
    end
    length(unique(internal_links)) == length(internal_links) || throw(
        ArgumentError("finite-cell MPS internal link indices must be distinct")
    )

    nonphysical(site) = filter(
        index -> !hastags(index, configuration_tag), inds(psi[site])
    )
    first_nonphysical = nonphysical(1)
    length(first_nonphysical) == 2 &&
        internal_links[1] in first_nonphysical || throw(ArgumentError(
            "finite-cell MPS first site must have one left boundary and one right link"
        ))
    left_boundary = only(filter(
        index -> index != internal_links[1], first_nonphysical
    ))
    all(site -> !hasind(psi[site], left_boundary), 2:count) || throw(
        ArgumentError("finite-cell MPS left boundary index must be dangling")
    )

    for site in 2:(count - 1)
        site_nonphysical = nonphysical(site)
        left_link = internal_links[site - 1]
        right_link = internal_links[site]
        length(site_nonphysical) == 2 &&
            left_link in site_nonphysical &&
            right_link in site_nonphysical || throw(ArgumentError(
                "finite-cell MPS interior site $site must have exactly its adjacent links"
            ))
    end

    last_nonphysical = nonphysical(count)
    length(last_nonphysical) == 2 &&
        internal_links[end] in last_nonphysical || throw(ArgumentError(
            "finite-cell MPS last site must have one left link and one right boundary"
        ))
    right_boundary = only(filter(
        index -> index != internal_links[end], last_nonphysical
    ))
    all(site -> !hasind(psi[site], right_boundary), 1:(count - 1)) || throw(
        ArgumentError("finite-cell MPS right boundary index must be dangling")
    )
    return nothing
end

function _left_virtual_index(psi::MPS, bond::Int, configuration_tag)
    if bond > 1
        links = linkinds(psi, bond - 1)
        length(links) == 1 || throw(ArgumentError(
            "finite-cell MPS bond $(bond - 1) must have exactly one link index"
        ))
        return only(links)
    end
    boundary = filter(
        index -> !hastags(index, configuration_tag),
        uniqueinds(psi[1], psi[2]),
    )
    length(boundary) == 1 || throw(ArgumentError(
        "finite-cell MPS must have exactly one left dangling virtual index"
    ))
    return only(boundary)
end

function _validate_fermionic_swap_cutoff(cutoff)
    cutoff isa Real && !(cutoff isa Bool) || throw(
        ArgumentError("fermionic MPS swap cutoff must be a real number")
    )
    isfinite(cutoff) && cutoff >= 0 || throw(
        ArgumentError("fermionic MPS swap cutoff must be finite and nonnegative")
    )
    return cutoff
end

function _validate_fermionic_swap(
    psi::MPS, bond::Int; configuration_tag, ortho, cutoff
)
    1 <= bond < length(psi) || throw(
        ArgumentError("fermionic MPS swap bond must lie in 1:$(length(psi) - 1)")
    )
    ortho == "left" || ortho == "right" || throw(
        ArgumentError("swap orthogonality direction must be left or right")
    )
    _validate_fermionic_swap_cutoff(cutoff)
    ITensors.using_auto_fermion() && error(
        "explicit transverse FSWAP requires ITensor auto-fermion mode to be disabled"
    )
    _physical_site_index(psi, bond, configuration_tag)
    _physical_site_index(psi, bond + 1, configuration_tag)
    _left_virtual_index(psi, bond, configuration_tag)
    links = linkinds(psi, bond)
    length(links) == 1 || throw(ArgumentError(
        "finite-cell MPS swap bond must have exactly one link index"
    ))
    return nothing
end

function _swap_fermionic_mps_sites_unchecked!(
    psi::MPS, bond::Int; configuration_tag, ortho, cutoff
)
    first_site = _physical_site_index(psi, bond, configuration_tag)
    second_site = _physical_site_index(psi, bond + 1, configuration_tag)
    orthocenter = ortho == "left" ? bond + 1 : bond
    if ITensorMPS.leftlim(psi) < bond - 1
        ITensorMPS.orthogonalize!(psi, bond)
    elseif ITensorMPS.rightlim(psi) > bond + 2
        ITensorMPS.orthogonalize!(psi, bond + 1)
    end
    combined = psi[bond] * psi[bond + 1]
    combined = noprime(
        combined * _physical_fermion_swap_phase(first_site, second_site)
    )
    left_index = _left_virtual_index(psi, bond, configuration_tag)
    link_tags = TagSet[tags(linkind(psi, bond))]
    swapped = MPS(
        combined,
        [second_site, first_site];
        leftinds=left_index,
        orthocenter=orthocenter - bond + 1,
        tags=link_tags,
        cutoff,
    )
    psi[bond:(bond + 1)] = swapped
    return psi
end

function _swap_fermionic_mps_sites!(
    psi::MPS, bond::Int; configuration_tag, ortho, cutoff
)
    _validate_fermionic_swap(
        psi, bond; configuration_tag, ortho, cutoff
    )
    return _swap_fermionic_mps_sites_unchecked!(
        psi, bond; configuration_tag, ortho, cutoff
    )
end

function _validate_fermionic_move(
    psi::MPS, source::Int, target::Int; configuration_tag, cutoff
)
    1 <= source <= length(psi) || throw(
        ArgumentError("fermionic MPS move source must lie in 1:$(length(psi))")
    )
    1 <= target <= length(psi) || throw(
        ArgumentError("fermionic MPS move target must lie in 1:$(length(psi))")
    )
    _validate_fermionic_swap_cutoff(cutoff)
    _validate_finite_cell_topology(psi, configuration_tag)
    route = if source < target
        [(bond, "left") for bond in source:(target - 1)]
    elseif source > target
        [(bond, "right") for bond in reverse(target:(source - 1))]
    else
        Tuple{Int,String}[]
    end
    for (bond, ortho) in route
        _validate_fermionic_swap(
            psi, bond; configuration_tag, ortho, cutoff
        )
    end
    return route
end

function _move_fermionic_mps_site!(
    psi::MPS, source::Int, target::Int; configuration_tag, cutoff
)
    route = _validate_fermionic_move(
        psi, source, target; configuration_tag, cutoff
    )
    for (bond, ortho) in route
        _swap_fermionic_mps_sites_unchecked!(
            psi, bond; configuration_tag, ortho, cutoff
        )
    end
    return psi
end

function _validate_left_cell_translation(
    translated::MPS,
    original::MPS,
    c;
    configuration_tag,
    cutoff,
)
    translated !== original || throw(ArgumentError(
        "translated and original finite cells must be distinct objects"
    ))
    expected_count = c.Ly * c.x_period
    length(translated) == expected_count || throw(ArgumentError(
        "translated finite cell must have length $expected_count"
    ))
    length(original) == expected_count || throw(ArgumentError(
        "original finite cell must have length $expected_count"
    ))
    _validate_fermionic_swap_cutoff(cutoff)
    for cell in (translated, original)
        _validate_finite_cell_topology(cell, configuration_tag)
    end

    route = Tuple{Int,String}[]
    for x in 0:(c.x_period - 1)
        first_site = x * c.Ly + 1
        last_site = first_site + c.Ly - 1
        append!(
            route,
            _validate_fermionic_move(
                translated,
                last_site - 1,
                first_site;
                configuration_tag,
                cutoff,
            ),
        )
        append!(
            route,
            _validate_fermionic_move(
                translated,
                last_site,
                first_site + 1;
                configuration_tag,
                cutoff,
            ),
        )
    end
    return route
end

function _translate_left_cell!(
    translated::MPS,
    original::MPS,
    c;
    configuration_tag,
    cutoff,
)
    route = _validate_left_cell_translation(
        translated,
        original,
        c;
        configuration_tag,
        cutoff,
    )
    count = length(translated)
    for (bond, ortho) in route
        _swap_fermionic_mps_sites_unchecked!(
            translated, bond; configuration_tag, ortho, cutoff
        )
    end
    for site in 1:count
        actual = _physical_site_index(translated, site, configuration_tag)
        desired = _physical_site_index(original, site, configuration_tag)
        actual == desired && continue
        translated[site] = replaceinds(translated[site], actual => desired)
    end
    return translated
end

function _translated_left_cell(psi::InfiniteCanonicalMPS, c, bond; cutoff)
    count = nsites(psi)
    configuration_tag = _configuration_signature(c)
    original = MPS([psi.AL[bond + offset] for offset in 1:count])
    translated = copy(original)
    _translate_left_cell!(
        translated,
        original,
        c;
        configuration_tag,
        cutoff,
    )
    return original, translated
end

function _dag_prime_finite_virtual_inds(psi::MPS, configuration_tag)
    bra = dag(psi)
    for site in eachindex(bra)
        virtual = filter(
            index -> !hastags(index, configuration_tag), inds(bra[site])
        )
        bra[site] = replaceinds(bra[site], virtual, prime.(virtual))
    end
    return bra
end

function _qn_block_ranges(index)
    ranges = Vector{UnitRange{Int}}(undef, nblocks(index))
    first_level = 1
    for block in eachindex(ranges)
        last_level = first_level + blockdim(index, block) - 1
        ranges[block] = first_level:last_level
        first_level = last_level + 1
    end
    first_level == dim(index) + 1 || error(
        "QN block dimensions do not match index dimension"
    )
    return ranges
end

function _same_qn_block_order(index, reference)
    hasqns(index) && hasqns(reference) || return false
    nblocks(index) == nblocks(reference) || return false
    for block in 1:nblocks(reference)
        qn(index, block) == qn(reference, block) || return false
        blockdim(index, block) == blockdim(reference, block) || return false
    end
    return true
end

function _qn_block_groups(index)
    hasqns(index) || error("QN block grouping requires a QN index")
    groups = Vector{Vector{Int}}()
    for block in 1:nblocks(index)
        group_position = findfirst(
            group -> qn(index, first(group)) == qn(index, block),
            groups,
        )
        if isnothing(group_position)
            push!(groups, Int[block])
        else
            push!(groups[group_position], block)
        end
    end
    return groups
end

function _neutral_qn_block_groups(row_index, column_index)
    _same_qn_block_order(row_index, column_index) || error(
        "neutral QN matrix indices must have identical block order"
    )
    qn_groups = _qn_block_groups(row_index)
    neutral_blocks = collect(
        ITensors.nzblocks(QN(), (row_index, column_index))
    )
    actual_blocks = [
        Tuple(Int.(Tuple(block))) for block in neutral_blocks
    ]
    expected_blocks = Tuple{Int,Int}[]
    for group in qn_groups, column_block in group, row_block in group
        push!(expected_blocks, (row_block, column_block))
    end
    sort(actual_blocks) == sort(expected_blocks) || error(
        "neutral matrix blocks do not match the complete repeated-QN block groups " *
        "(actual=$actual_blocks, expected=$expected_blocks)"
    )
    structural_groups = [
        [
            block for block in neutral_blocks if
            Int(block[1]) in group && Int(block[2]) in group
        ] for group in qn_groups
    ]
    return qn_groups, structural_groups
end

function _twisted_group_dimension(vector_indices, blocks)
    return sum(
        prod(
            blockdim(vector_indices[dimension], Int(block[dimension])) for
            dimension in eachindex(vector_indices)
        ) for block in blocks
    )
end

function _logical_qn_block_matrix(
    ordered,
    row_index,
    column_index,
    qn_group,
)
    row_dimensions = [blockdim(row_index, block) for block in qn_group]
    column_dimensions = [blockdim(column_index, block) for block in qn_group]
    row_dimensions == column_dimensions || error(
        "repeated-QN row and column block dimensions must match"
    )
    dimension = sum(row_dimensions)
    matrix = zeros(ComplexF64, dimension, dimension)
    row_first = 1
    for (row_position, row_block) in enumerate(qn_group)
        row_last = row_first + row_dimensions[row_position] - 1
        column_first = 1
        for (column_position, column_block) in enumerate(qn_group)
            column_last = column_first + column_dimensions[column_position] - 1
            block_value = ordered[ITensors.Block(row_block, column_block)]
            if !isnothing(block_value)
                matrix_block = Matrix{ComplexF64}(Array(block_value))
                size(matrix_block) == (
                    row_dimensions[row_position],
                    column_dimensions[column_position],
                ) || error(
                    "QN matrix structural block ($row_block, $column_block) " *
                    "has inconsistent dimensions $(size(matrix_block))"
                )
                matrix[row_first:row_last, column_first:column_last] = matrix_block
            end
            column_first = column_last + 1
        end
        row_first = row_last + 1
    end
    return matrix
end

function _twisted_matrix_indices(fixed_point, schmidt_index)
    vector_indices = collect(inds(fixed_point))
    length(vector_indices) == 2 || error(
        "twisted transfer fixed point must have two virtual indices " *
        "(vector=$vector_indices)"
    )
    order = sortperm(vector_indices; by=plev, rev=true)
    row_index, column_index = vector_indices[order]
    plev(row_index) > plev(column_index) || error(
        "twisted transfer fixed-point indices must have distinct prime levels"
    )
    dir(row_index) == dir(dag(schmidt_index)) || error(
        "twisted transfer row index is not the primed Schmidt-link dual"
    )
    dir(column_index) == dir(schmidt_index) || error(
        "twisted transfer column index is not the Schmidt link"
    )
    _same_qn_block_order(row_index, schmidt_index) || error(
        "twisted transfer row-index QN blocks do not match Schmidt-link order"
    )
    _same_qn_block_order(column_index, schmidt_index) || error(
        "twisted transfer column-index QN blocks do not match Schmidt-link order"
    )
    return row_index, column_index
end

function _twisted_block_matrices(fixed_point, schmidt_index, Ny::Int)
    Ny > 0 || throw(ArgumentError("Ny must be positive"))
    row_index, column_index = _twisted_matrix_indices(
        fixed_point, schmidt_index
    )
    qn_groups, _ = _neutral_qn_block_groups(row_index, column_index)

    ordered = ITensors.permute(fixed_point, row_index, column_index)
    matrices = Matrix{ComplexF64}[]
    for (sector, qn_group) in enumerate(qn_groups)
        matrix_block = _logical_qn_block_matrix(
            ordered, row_index, column_index, qn_group
        )
        block_dimension = size(matrix_block, 1)
        amplitude = Float64(LinearAlgebra.norm(matrix_block) / sqrt(block_dimension))
        isfinite(amplitude) && amplitude > eps(Float64) || throw(
            MomentumConstructionError(
                "twisted transfer charge sector $sector is missing or has invalid norm"
            )
        )
        matrix_block ./= amplitude

        anchor_position = argmax(abs.(matrix_block))
        anchor = matrix_block[anchor_position]
        isfinite(anchor) && abs(anchor) > eps(Float64) || throw(
            MomentumConstructionError(
                "twisted transfer charge sector $sector has no deterministic phase anchor"
            )
        )
        matrix_block .*= cis(-angle(anchor))

        power_phase = LinearAlgebra.tr(matrix_block^Ny) / block_dimension
        isfinite(power_phase) && abs(power_phase) > eps(Float64) || throw(
            MomentumConstructionError(
                "twisted transfer charge sector $sector has invalid Ny-power phase"
            )
        )
        matrix_block .*= cis(-angle(power_phase) / Ny)
        push!(matrices, matrix_block)
    end
    return matrices
end

function _twisted_block_matrix(fixed_point, schmidt_index, Ny::Int)
    matrices = _twisted_block_matrices(fixed_point, schmidt_index, Ny)
    dimension = sum(block -> size(block, 1), matrices)
    matrix = zeros(ComplexF64, dimension, dimension)
    first_level = 1
    for block in matrices
        last_level = first_level + size(block, 1) - 1
        matrix[first_level:last_level, first_level:last_level] = block
        first_level = last_level + 1
    end
    return matrix
end

function _random_twisted_group(rng, vector_indices, blocks)
    index_ranges = _qn_block_ranges.(vector_indices)
    initial = ITensor(ComplexF64, vector_indices...)
    for block in blocks
        selected_ranges = [
            index_ranges[dimension][Int(block[dimension])] for
            dimension in eachindex(vector_indices)
        ]
        block_dimensions = Tuple(length.(selected_ranges))
        selectors = [
            vector_indices[dimension] => selected_ranges[dimension] for
            dimension in eachindex(vector_indices)
        ]
        setindex!(
            initial,
            randn(rng, ComplexF64, block_dimensions...),
            selectors...,
        )
    end
    flux(initial) == QN() || error(
        "twisted transfer charge-group initial tensor has nonzero flux"
    )
    LinearAlgebra.norm(initial) > 0 || error(
        "twisted transfer charge-group initial tensor has zero norm"
    )
    return initial
end

function _random_twisted_block(rng, vector_indices, block)
    return _random_twisted_group(rng, vector_indices, [block])
end

function _project_twisted_group(value, vector_indices, blocks)
    ordered = ITensors.permute(value, vector_indices...)
    projected = ITensor(ComplexF64, vector_indices...)
    for block in blocks
        block_value = ordered[block]
        isnothing(block_value) || setindex!(projected, block_value, block)
    end
    return projected
end

function _project_twisted_block(value, vector_indices, block)
    return _project_twisted_group(value, vector_indices, [block])
end

function _relative_transfer_residual(image, value, vector)
    scale = max(
        Float64(LinearAlgebra.norm(image)),
        Float64(abs(value) * LinearAlgebra.norm(vector)),
        eps(Float64),
    )
    return Float64(LinearAlgebra.norm(image - value * vector) / scale)
end

function _mixed_symmetry_eigenvalue_residual(eigenvalues)
    values = ComplexF64.(collect(eigenvalues))
    isempty(values) && throw(
        ArgumentError("mixed-symmetry eigenvalues must not be empty")
    )
    all(isfinite, values) || throw(
        ArgumentError("mixed-symmetry eigenvalues must be finite")
    )
    magnitudes = abs.(values)
    any(iszero, magnitudes) && return Inf
    magnitude_residual = maximum(abs.(magnitudes .- 1.0))
    normalized_phases = values ./ magnitudes
    phase_residual = maximum(abs.(normalized_phases .- first(normalized_phases)))
    return Float64(max(magnitude_residual, phase_residual))
end

function _dominant_mixed_eigenpair(values, vectors, tol::Real)
    tolerance = Float64(tol)
    isfinite(tolerance) && tolerance >= 0 || throw(
        ArgumentError("dominant mixed-transfer gap tolerance must be finite and nonnegative")
    )
    eigenvalues = ComplexF64.(collect(values))
    eigenvectors = collect(vectors)
    isempty(eigenvalues) && throw(
        ArgumentError("dominant mixed-transfer eigenvalues must not be empty")
    )
    length(eigenvectors) == length(eigenvalues) || throw(
        DimensionMismatch("mixed-transfer eigenvalue and eigenvector counts must match")
    )
    all(isfinite, eigenvalues) || throw(
        ArgumentError("dominant mixed-transfer eigenvalues must be finite")
    )

    order = sortperm(eachindex(eigenvalues); by=index -> abs(eigenvalues[index]), rev=true)
    dominant_index = first(order)
    dominant_magnitude = abs(eigenvalues[dominant_index])
    gap = if length(order) == 1
        Inf
    else
        subleading_magnitude = abs(eigenvalues[order[2]])
        Float64(
            (dominant_magnitude - subleading_magnitude) /
            max(dominant_magnitude, eps(Float64))
        )
    end
    gap > tolerance || throw(
        MomentumConstructionError(
            "transverse-translation symmetry has an unresolved dominant " *
            "mixed-transfer eigenvalue: relative magnitude gap $gap does not " *
            "exceed tolerance $tolerance"
        )
    )
    return (;
        value=eigenvalues[dominant_index],
        vector=eigenvectors[dominant_index],
        gap,
    )
end

function _left_twisted_transfer_map(original, translated, configuration_tag)
    original_bra = _dag_prime_finite_virtual_inds(
        original, configuration_tag
    )
    return transpose(ITensorMap(translated, original_bra))
end

function _primal_twisted_fixed_point(fixed_point)
    vector_indices = collect(inds(fixed_point))
    return ITensors.setinds(fixed_point, dag.(vector_indices))
end

function _twisted_virtual_translation(
    psi::InfiniteCanonicalMPS,
    c::InfiniteCylinderConfig,
    bond::Int;
    schmidt_index,
    tol::Float64,
    rng::AbstractRNG,
    cutoff::Float64,
    kwargs...,
)
    original, translated = _translated_left_cell(psi, c, bond; cutoff)
    transfer = _left_twisted_transfer_map(
        original, translated, _configuration_signature(c)
    )
    vector_indices = collect(dag(input_inds(transfer)))
    length(vector_indices) == 2 || error(
        "twisted transfer must have two virtual input indices " *
        "(input=$(input_inds(transfer)), output=$(output_inds(transfer)))"
    )
    _, neutral_groups = _neutral_qn_block_groups(
        vector_indices[1], vector_indices[2]
    )
    isempty(neutral_groups) && throw(
        MomentumConstructionError(
            "twisted transfer has no neutral virtual-charge blocks"
        )
    )

    fixed_points = ITensor[]
    eigenvalues = ComplexF64[]
    block_residuals = Float64[]
    requested_eigenpairs = Int[]
    dominant_gaps = Float64[]
    for (sector, structural_blocks) in enumerate(neutral_groups)
        initial = _random_twisted_group(
            rng, vector_indices, structural_blocks
        )
        function block_transfer(vector)
            image = transfer(vector)
            hassameinds(image, vector) || error(
                "twisted transfer charge block $sector does not preserve its vector space " *
                "(input=$(input_inds(transfer)), output=$(output_inds(transfer)), " *
                "vector=$(collect(inds(vector))), image=$(collect(inds(image))))"
            )
            return _project_twisted_group(
                image, vector_indices, structural_blocks
            )
        end
        initial_image = block_transfer(initial)
        LinearAlgebra.norm(initial_image) > eps(Float64) || throw(
            MomentumConstructionError(
                "transverse-translation symmetry charge block $sector is " *
                "missing from the mixed transfer map"
            )
        )
        operator_dimension = _twisted_group_dimension(
            vector_indices, structural_blocks
        )
        requested_count = operator_dimension > 1 ? 2 : 1
        values, vectors, info = KrylovKit.eigsolve(
            block_transfer,
            initial,
            requested_count,
            :LM;
            tol,
            verbosity=0,
            kwargs...,
        )
        info.converged >= requested_count || throw(
            MomentumConstructionError(
                "twisted transfer charge block $sector converged $(info.converged) " *
                "of $requested_count requested eigenpairs"
            )
        )
        length(values) >= requested_count || throw(
            MomentumConstructionError(
                "twisted transfer charge block $sector returned $(length(values)) " *
                "of $requested_count requested eigenvalues"
            )
        )
        length(vectors) >= requested_count || throw(
            MomentumConstructionError(
                "twisted transfer charge block $sector returned $(length(vectors)) " *
                "of $requested_count requested eigenvectors"
            )
        )
        all(isfinite, values[1:requested_count]) || throw(
            MomentumConstructionError(
                "twisted transfer charge block $sector returned nonfinite eigenvalues"
            )
        )
        all(
            vector -> isfinite(LinearAlgebra.norm(vector)) &&
                LinearAlgebra.norm(vector) > eps(Float64),
            vectors[1:requested_count],
        ) || throw(
            MomentumConstructionError(
                "twisted transfer charge block $sector returned invalid eigenvectors"
            )
        )
        dominant = _dominant_mixed_eigenpair(values, vectors, tol)
        value = dominant.value
        vector = dominant.vector
        image = transfer(vector)
        hassameinds(image, vector) || error(
            "twisted transfer charge block $sector eigenvector left its vector space"
        )
        residual = _relative_transfer_residual(image, value, vector)
        isfinite(residual) && residual <= tol || throw(
            MomentumConstructionError(
                "twisted transfer charge block $sector residual $residual exceeds tolerance $tol"
            )
        )
        push!(fixed_points, vector)
        push!(eigenvalues, value)
        push!(block_residuals, residual)
        push!(requested_eigenpairs, requested_count)
        push!(dominant_gaps, dominant.gap)
    end

    fixed_point = copy(first(fixed_points))
    for vector in Iterators.drop(fixed_points, 1)
        fixed_point += vector
    end
    fixed_point = _primal_twisted_fixed_point(fixed_point)
    blocks = _twisted_block_matrices(fixed_point, schmidt_index, c.Ny)
    eigenpair_residual = maximum(block_residuals)
    symmetry_residual = _mixed_symmetry_eigenvalue_residual(eigenvalues)
    residual = max(eigenpair_residual, symmetry_residual)
    return (;
        blocks,
        residual,
        eigenvalues,
        eigenpair_residual,
        symmetry_residual,
        requested_eigenpairs,
        dominant_gaps,
    )
end

function _schmidt_density_blocks(center, left_index, right_index)
    hasqns(left_index) || error("Schmidt left link must carry QN blocks")
    hasqns(right_index) || error("Schmidt right link must carry QN blocks")
    hassameinds(center, ITensor(left_index, right_index)) || error(
        "Schmidt center indices do not match the requested left and right links"
    )
    left_bra_index = dag(prime(left_index))
    density_tensor = ITensors.permute(
        center * dag(prime(center, left_index)),
        left_index,
        left_bra_index,
    )
    flux(density_tensor) == QN() || error(
        "Schmidt density tensor must have neutral QN flux"
    )
    qn_groups, _ = _neutral_qn_block_groups(
        left_index, left_bra_index
    )

    density_blocks = Matrix{ComplexF64}[]
    block_charges = Int[]
    for qn_group in qn_groups
        density_block = _logical_qn_block_matrix(
            density_tensor,
            left_index,
            left_bra_index,
            qn_group,
        )
        all(isfinite, density_block) || error(
            "Schmidt density charge sector contains nonfinite entries"
        )
        push!(density_blocks, density_block)
        push!(
            block_charges,
            Int(val(qn(left_index, first(qn_group)), "Nf")),
        )
    end
    trace_density = Float64(sum(real(LinearAlgebra.tr(block)) for block in density_blocks))
    isfinite(trace_density) && trace_density > 0 || error(
        "Schmidt density blocks must have positive finite total trace"
    )
    return density_blocks, block_charges
end

function _schmidt_density_on_left_link(psi, c, bond)
    center = psi.C[bond]
    left_index = only(commoninds(center, psi.AL[bond]))
    right_index = only(uniqueinds(center, left_index))
    density_blocks, block_charges = _schmidt_density_blocks(
        center, left_index, right_index
    )
    return density_blocks, block_charges, left_index
end

"""
    momentum_entanglement_data(psi, c; cut_x, rng=nothing, rng_seed=nothing, kwargs...)

Construct the charge-resolved momentum entanglement data at a transverse cut.
When neither RNG keyword is supplied, each call creates a new `Random.Xoshiro`
with the fixed default seed `0x4d4f4d454e54554d`, independent of task-global RNG
state. Passing `rng` uses and consumes that supplied RNG object. The `rng` and
`rng_seed` keywords are mutually exclusive. An explicit `rng_seed` must be a
non-`Bool` integer in the inclusive `UInt64` range.
"""
function momentum_entanglement_data(
    psi::InfiniteCanonicalMPS,
    c::InfiniteCylinderConfig;
    cut_x::Integer,
    residual_tol::Real=1e-8,
    transfer_tol::Real=residual_tol,
    degeneracy_tol::Real=residual_tol,
    probability_tol::Real=0.0,
    translation_cutoff::Real=0.0,
    rng=nothing,
    rng_seed=nothing,
    kwargs...,
)
    resolved_rng = _resolve_momentum_rng(rng, rng_seed)
    _paper_translation_Ny(c)
    !(cut_x isa Bool) && 1 <= cut_x <= c.x_period || throw(
        ArgumentError("cut_x must identify a transverse bond in the reference cell")
    )
    isfinite(transfer_tol) && transfer_tol > 0 || throw(
        ArgumentError("transfer_tol must be finite and positive")
    )
    isfinite(translation_cutoff) && translation_cutoff >= 0 || throw(
        ArgumentError("translation_cutoff must be finite and nonnegative")
    )
    probability_threshold = _probability_threshold(probability_tol)
    _validate_observable_cell(psi, c)
    bond = Int(cut_x) * c.Ly
    density_blocks, block_charges, schmidt_index = _schmidt_density_on_left_link(
        psi, c, bond
    )
    try
        virtual = _twisted_virtual_translation(
            psi,
            c,
            bond;
            schmidt_index,
            tol=Float64(transfer_tol),
            rng=resolved_rng,
            cutoff=Float64(translation_cutoff),
            kwargs...,
        )
        return _validate_momentum_block_vectors(
            density_blocks,
            virtual.blocks;
            Ny=c.Ny,
            raw_charges=block_charges,
            residual_tol,
            degeneracy_tol,
            transfer_residual=virtual.residual,
            probability_tol=probability_threshold,
        )
    catch error
        error isa MomentumConstructionError || rethrow()
        normalized_blocks, _, normalized_charges = _block_momentum_inputs(
            density_blocks,
            [
                Matrix{ComplexF64}(LinearAlgebra.I, size(block, 1), size(block, 1)) for
                block in density_blocks
            ],
            block_charges,
            Float64(residual_tol),
        )
        probabilities, basis_charges, _, _ = _block_density_basis(
            normalized_blocks, normalized_charges, Float64(residual_tol)
        )
        retained = findall(
            probability -> probability > probability_threshold,
            probabilities,
        )
        isempty(retained) && throw(
            ArgumentError(
                "probability_tol removed all positive Schmidt support"
            )
        )
        reason = "twisted-transfer construction failed: $(sprint(showerror, error))"
        validation = MomentumValidation(Inf, Inf, Inf, Inf, false, reason)
        return _invalid_momentum_data(
            probabilities[retained], basis_charges[retained], validation
        )
    end
end
