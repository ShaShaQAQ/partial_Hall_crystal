struct MPSKitEnergyData
    per_cell::Float64
    per_x::Float64
    per_unit_cell::Float64
    per_site::Float64
    imaginary_per_cell::Float64
end

struct MPSKitDensityRow
    site::Int
    x::Int
    y::Int
    orbital::Symbol
    density::Float64
end

struct MPSKitDensityData
    rows::Vector{MPSKitDensityRow}
    densities::Vector{Float64}
    ring_sums::Vector{Float64}
    total_particles::Float64
    expected_particles::Int
    particle_error::Float64
end

struct MPSKitEntanglementLevel
    cut_x::Int
    bond::Int
    level::Int
    sector::String
    odd_parity::Bool
    raw_charge::Int
    physical_charge::Float64
    singular_value::Float64
    probability::Float64
    entanglement_energy::Float64
end

struct MPSKitEntanglementSector
    cut_x::Int
    bond::Int
    sector::String
    odd_parity::Bool
    raw_charge::Int
    physical_charge::Float64
    weight::Float64
end

struct MPSKitEntanglementData
    cut_x::Int
    bond::Int
    entropy::Float64
    mean_raw_charge::Float64
    mean_physical_charge::Float64
    levels::Vector{MPSKitEntanglementLevel}
    sectors::Vector{MPSKitEntanglementSector}
end

struct MPSKitNeutralTransferData
    eigenvalues::Vector{ComplexF64}
    phases::Vector{Float64}
    residual_norms::Vector{Float64}
    converged::Int
    valid::Bool
    reason::String
    ratio::Float64
    xi_cell::Float64
    xi_x::Float64
end

struct MPSKitSchmidtMomentumData
    spectrum::MomentumEntanglementData
    closure_residual::Float64
    valid::Bool
    reason::String
end

struct MPSKitMomentumData
    cut_x::Int
    bond::Int
    spectrum::MomentumEntanglementData
    closure_residual::Float64
    transfer_eigenvalue::ComplexF64
    transfer_residual::Float64
    transfer_converged::Bool
    dominant_gap::Float64
    translation_phase_correction::ComplexF64
    valid::Bool
    reason::String
end

function _validate_mpskit_observable_state(
    state,
    c::InfiniteCylinderConfig;
    physical_spaces::Bool,
)
    length(state) == sites_per_cell(c) || throw(
        ArgumentError("state and configuration unit-cell lengths do not match"),
    )
    if physical_spaces
        expected = mpskit_physical_spaces(c)
        all(
            site -> MPSKit.physicalspace(state, site) == expected[site],
            1:length(state),
        ) || throw(
            ArgumentError("state does not use the configuration physical spaces"),
        )
    end
    return nothing
end

function mpskit_energy_data(
    state,
    hamiltonian,
    c::InfiniteCylinderConfig;
    environments=nothing,
    imaginary_tol::Real=1e-12,
)
    isfinite(imaginary_tol) && imaginary_tol >= 0 || throw(
        ArgumentError("imaginary_tol must be finite and nonnegative"),
    )
    _validate_mpskit_observable_state(state, c; physical_spaces=true)
    length(hamiltonian) == length(state) || throw(
        ArgumentError("Hamiltonian and state unit-cell lengths do not match"),
    )
    local_environments = isnothing(environments) ?
        MPSKit.environments(state, hamiltonian, state) : environments
    energy = ComplexF64(
        MPSKit.expectation_value(state, hamiltonian, local_environments),
    )
    isfinite(real(energy)) && isfinite(imag(energy)) || throw(
        ArgumentError("MPSKit energy per cell must be finite"),
    )
    abs(imag(energy)) <= imaginary_tol || throw(
        ArgumentError(
            "MPSKit energy imaginary part $(imag(energy)) exceeds tolerance $imaginary_tol",
        ),
    )
    per_cell = real(energy)
    return MPSKitEnergyData(
        per_cell,
        per_cell / c.x_period,
        per_cell / unit_cells_per_cell(c),
        per_cell / sites_per_cell(c),
        imag(energy),
    )
end

function _mpskit_density_coordinates(c::InfiniteCylinderConfig, site::Int)
    if c.geometry == :paper_straight
        return paper_coordinates(c, site)
    end
    x, y = site_coordinates(c, site)
    return x, y, :site
end

function mpskit_density_data(
    state,
    c::InfiniteCylinderConfig;
    imaginary_tol::Real=1e-12,
)
    isfinite(imaginary_tol) && imaginary_tol >= 0 || throw(
        ArgumentError("imaginary_tol must be finite and nonnegative"),
    )
    _validate_mpskit_observable_state(state, c; physical_spaces=true)
    number_operator = _mpskit_fermion_operators(
        first(mpskit_physical_spaces(c));
        centered=true,
    ).number
    densities = Float64[]
    sizehint!(densities, length(state))
    rows = MPSKitDensityRow[]
    sizehint!(rows, length(state))
    for site in 1:length(state)
        value = ComplexF64(
            MPSKit.expectation_value(state, site => number_operator),
        )
        isfinite(real(value)) && isfinite(imag(value)) || throw(
            ArgumentError("density at site $site must be finite"),
        )
        abs(imag(value)) <= imaginary_tol || throw(
            ArgumentError(
                "density imaginary part at site $site exceeds tolerance $imaginary_tol",
            ),
        )
        density = real(value)
        x, y, orbital = _mpskit_density_coordinates(c, site)
        push!(densities, density)
        push!(rows, MPSKitDensityRow(site, x, y, orbital, density))
    end
    ring_sums = [
        sum(@view densities[(x * c.Ly + 1):((x + 1) * c.Ly)]) for
        x in 0:(c.x_period - 1)
    ]
    total_particles = sum(densities)
    expected_particles = particles_per_cell(c)
    return MPSKitDensityData(
        rows,
        densities,
        ring_sums,
        total_particles,
        expected_particles,
        total_particles - expected_particles,
    )
end

function _mpskit_entanglement_sector(sector)
    length(sector) == 2 || throw(
        ArgumentError("Schmidt sector does not use parity-and-charge labels"),
    )
    parity = sector[1]
    parity isa TensorKit.FermionParity || throw(
        ArgumentError("Schmidt sector does not carry fermion parity"),
    )
    raw = TensorKit.charge(sector[2])
    isinteger(raw) || throw(
        ArgumentError("Schmidt U(1) charge must be integer-valued"),
    )
    return Bool(getfield(parity, :isodd)), Int(raw)
end

function _mpskit_entanglement_cut_data(
    state,
    c::InfiniteCylinderConfig;
    cut_x::Integer,
)
    !(cut_x isa Bool) && 1 <= cut_x <= c.x_period || throw(
        ArgumentError("cut_x must identify a physical x cut in the unit cell"),
    )
    _validate_mpskit_observable_state(state, c; physical_spaces=true)
    cut = Int(cut_x)
    bond = Base.Checked.checked_mul(cut, c.Ly)
    spectrum = MPSKit.entanglement_spectrum(state, bond)
    normalization = sum(
        TensorKit.dim(sector) * sum(abs2, singular_values) for
        (sector, singular_values) in pairs(spectrum)
    )
    isfinite(normalization) && normalization > 0 || throw(
        ArgumentError("Schmidt spectrum has invalid or zero norm"),
    )

    levels = MPSKitEntanglementLevel[]
    sectors = MPSKitEntanglementSector[]
    entropy = 0.0
    mean_raw_charge = 0.0
    scale = charge_scale(c)
    level_index = 0
    for (sector, singular_values) in pairs(spectrum)
        odd_parity, raw_charge = _mpskit_entanglement_sector(sector)
        sector_weight = TensorKit.dim(sector) *
            sum(abs2, singular_values) / normalization
        push!(
            sectors,
            MPSKitEntanglementSector(
                cut,
                bond,
                string(sector),
                odd_parity,
                raw_charge,
                raw_charge / scale,
                Float64(sector_weight),
            ),
        )
        mean_raw_charge += sector_weight * raw_charge
        for singular_value in singular_values
            level_index += 1
            value = Float64(real(singular_value))
            isfinite(value) && value >= 0 || throw(
                ArgumentError("Schmidt singular values must be finite and nonnegative"),
            )
            probability = abs2(value) / normalization
            if !iszero(probability)
                entropy -= TensorKit.dim(sector) * probability * log(probability)
            end
            push!(
                levels,
                MPSKitEntanglementLevel(
                    cut,
                    bond,
                    level_index,
                    string(sector),
                    odd_parity,
                    raw_charge,
                    raw_charge / scale,
                    value,
                    probability,
                    iszero(probability) ? Inf : -log(probability),
                ),
            )
        end
    end
    sort!(sectors; by=sector -> (sector.raw_charge, sector.odd_parity))
    return MPSKitEntanglementData(
        cut,
        bond,
        entropy,
        mean_raw_charge,
        mean_raw_charge / scale,
        levels,
        sectors,
    )
end

function mpskit_entanglement_data(
    state,
    c::InfiniteCylinderConfig;
    cut_x::Union{Nothing,Integer}=nothing,
)
    if isnothing(cut_x)
        return [
            _mpskit_entanglement_cut_data(state, c; cut_x=value) for
            value in 1:c.x_period
        ]
    end
    return _mpskit_entanglement_cut_data(state, c; cut_x)
end

function _invalid_mpskit_neutral_transfer(
    values,
    residual_norms,
    converged::Integer,
    reason::AbstractString,
)
    eigenvalues = ComplexF64.(values)
    return MPSKitNeutralTransferData(
        eigenvalues,
        Float64[angle(value) for value in eigenvalues],
        Float64.(residual_norms),
        Int(converged),
        false,
        String(reason),
        NaN,
        NaN,
        NaN,
    )
end

function _mpskit_neutral_transfer_result(
    c::InfiniteCylinderConfig,
    values,
    residual_norms,
    converged::Integer;
    residual_tolerance::Real,
    ratio_tolerance::Real,
)
    eigenvalues = ComplexF64.(values)
    residuals = Float64.(residual_norms)
    length(eigenvalues) == length(residuals) || throw(
        ArgumentError("transfer eigenvalue and residual counts must match"),
    )
    converged >= 2 && length(eigenvalues) >= 2 || return
        _invalid_mpskit_neutral_transfer(
            eigenvalues,
            residuals,
            converged,
            "fewer than two neutral transfer eigenpairs converged",
        )
    all(isfinite, eigenvalues) || return _invalid_mpskit_neutral_transfer(
        eigenvalues,
        residuals,
        converged,
        "neutral transfer eigenvalues must be finite",
    )
    all(value -> isfinite(value) && value >= 0, residuals) || return
        _invalid_mpskit_neutral_transfer(
            eigenvalues,
            residuals,
            converged,
            "neutral transfer residuals must be finite and nonnegative",
        )
    all(<=(residual_tolerance), residuals[1:2]) || return
        _invalid_mpskit_neutral_transfer(
            eigenvalues,
            residuals,
            converged,
            "leading neutral transfer residual exceeds tolerance",
        )
    abs(eigenvalues[1]) > 0 || return _invalid_mpskit_neutral_transfer(
        eigenvalues,
        residuals,
        converged,
        "leading neutral transfer eigenvalue must be nonzero",
    )
    ratio = abs(eigenvalues[2] / eigenvalues[1])
    isfinite(ratio) && ratio > 0 || return _invalid_mpskit_neutral_transfer(
        eigenvalues,
        residuals,
        converged,
        "neutral transfer ratio must be finite and positive",
    )
    ratio < 1 - ratio_tolerance || return _invalid_mpskit_neutral_transfer(
        eigenvalues,
        residuals,
        converged,
        "neutral transfer ratio must be strictly below one",
    )
    xi_cell = -1 / log(ratio)
    xi_x = c.x_period * xi_cell
    all(isfinite, (xi_cell, xi_x)) || return _invalid_mpskit_neutral_transfer(
        eigenvalues,
        residuals,
        converged,
        "neutral transfer correlation length must be finite",
    )
    return MPSKitNeutralTransferData(
        eigenvalues,
        Float64[angle(value) for value in eigenvalues],
        residuals,
        Int(converged),
        true,
        "valid neutral transfer spectrum",
        ratio,
        xi_cell,
        xi_x,
    )
end

function _mpskit_randomize_transfer_initial!(tensor, seed::Integer)
    rng = Random.Xoshiro(seed)
    for (_, block) in TensorKit.blocks(tensor)
        randn!(rng, block)
    end
    TensorKit.normalize!(tensor)
    return tensor
end

function _mpskit_neutral_transfer_initial(state, seed::Integer)
    virtual = MPSKit.left_virtualspace(state, 1)
    space_type = TensorKit.spacetype(virtual)
    neutral = TensorKit.leftunit(state)
    auxiliary = space_type(neutral => 1)
    domain_space = TensorKit.:(⊗)(TensorKit.dual(auxiliary), virtual)
    map_space = TensorKit.:(←)(virtual, domain_space)
    initial = similar(state.AL[1], map_space)
    return _mpskit_randomize_transfer_initial!(initial, seed)
end

function mpskit_neutral_transfer_data(
    state,
    c::InfiniteCylinderConfig;
    neigs::Integer=4,
    tol::Real=1e-10,
    maxiter::Integer=500,
    krylovdim::Integer=30,
    seed::Integer=0x4d50534b,
    ratio_tolerance::Real=1e-10,
)
    !(neigs isa Bool) && neigs >= 2 || throw(
        ArgumentError("neigs must be an integer of at least two"),
    )
    isfinite(tol) && tol > 0 || throw(
        ArgumentError("tol must be finite and positive"),
    )
    !(maxiter isa Bool) && maxiter > 0 || throw(
        ArgumentError("maxiter must be positive"),
    )
    !(krylovdim isa Bool) && krylovdim > neigs || throw(
        ArgumentError("krylovdim must be greater than neigs"),
    )
    isfinite(ratio_tolerance) && ratio_tolerance >= 0 || throw(
        ArgumentError("ratio_tolerance must be finite and nonnegative"),
    )
    _validate_mpskit_observable_state(state, c; physical_spaces=false)
    initial = _mpskit_neutral_transfer_initial(state, seed)
    transfer = TensorKit.flip(MPSKit.TransferMatrix(state.AL, state.AL))
    available = TensorKit.dim(TensorKit.space(initial))
    requested = min(Int(neigs), available)
    if requested == 1
        applied = transfer(initial)
        eigenvalue = TensorKit.dot(initial, applied) /
            TensorKit.dot(initial, initial)
        residual = TensorKit.norm(applied - eigenvalue * initial) /
            TensorKit.norm(initial)
        return _invalid_mpskit_neutral_transfer(
            [eigenvalue],
            [residual],
            1,
            "fewer than two neutral transfer eigenpairs converged",
        )
    end
    algorithm = KrylovKit.Arnoldi(;
        tol=Float64(tol),
        maxiter=Int(maxiter),
        eager=true,
        krylovdim=Int(krylovdim),
        verbosity=0,
    )
    values, vectors, info = KrylovKit.eigsolve(
        transfer,
        initial,
        requested,
        :LM,
        algorithm,
    )
    order = sortperm(values; by=abs, rev=true)[1:requested]
    ordered_values = ComplexF64[values[index] for index in order]
    ordered_vectors = [vectors[index] for index in order]
    residual_norms = Float64[
        TensorKit.norm(transfer(vector) - value * vector) /
        TensorKit.norm(vector) for
        (value, vector) in zip(ordered_values, ordered_vectors)
    ]
    return _mpskit_neutral_transfer_result(
        c,
        ordered_values,
        residual_norms,
        min(info.converged, length(ordered_values));
        residual_tolerance=10Float64(tol),
        ratio_tolerance,
    )
end

function _mpskit_identity_finite_mpo(physical_spaces)
    tensors = [
        MPSKit.add_util_leg(TensorKit.id(physical_space)) for
        physical_space in physical_spaces
    ]
    return MPSKit.FiniteMPO(tensors)
end

function _mpskit_adjacent_fermionic_swap_mpo(physical_spaces, site::Int)
    1 <= site < length(physical_spaces) || throw(
        ArgumentError("adjacent swap site lies outside the finite MPO"),
    )
    physical_spaces[site] == physical_spaces[site + 1] || throw(
        ArgumentError("adjacent fermionic swap requires equal physical spaces"),
    )
    local_gate = TensorKit.TensorMap(
        TensorKit.BraidingTensor(
            physical_spaces[site],
            physical_spaces[site + 1],
        ),
    )
    local_mpo = MPSKit.FiniteMPO(local_gate)
    tensors = [
        MPSKit.add_util_leg(TensorKit.id(physical_space)) for
        physical_space in physical_spaces
    ]
    tensors[site] = local_mpo[1]
    tensors[site + 1] = local_mpo[2]
    return MPSKit.FiniteMPO(tensors)
end

function _mpskit_finite_transverse_translation_mpo(
    c::InfiniteCylinderConfig,
)
    _paper_translation_Ny(c)
    physical_spaces = mpskit_physical_spaces(c)
    translation = _mpskit_identity_finite_mpo(physical_spaces)
    for x in 0:(c.x_period - 1)
        first_site = x * c.Ly + 1
        last_site = (x + 1) * c.Ly
        for swap_site in (last_site - 2):-1:first_site
            swap_gate = _mpskit_adjacent_fermionic_swap_mpo(
                physical_spaces,
                swap_site,
            )
            translation = swap_gate * translation
        end
        for swap_site in (last_site - 1):-1:(first_site + 1)
            swap_gate = _mpskit_adjacent_fermionic_swap_mpo(
                physical_spaces,
                swap_site,
            )
            translation = swap_gate * translation
        end
    end
    return translation
end

function mpskit_transverse_translation_mpo(c::InfiniteCylinderConfig)
    finite = _mpskit_finite_transverse_translation_mpo(c)
    return MPSKit.InfiniteMPO(copy.(collect(parent(finite))))
end

function _mpskit_retained_block_residuals(
    density_blocks,
    translation_blocks,
    probability_tol::Float64,
    Ny::Int,
)
    probabilities, _, density_bases, locations = _block_density_basis(
        density_blocks,
        zeros(Int, length(density_blocks)),
        max(probability_tol, eps(Float64)),
    )
    density_basis_translations = [
        density_bases[block]' * translation_blocks[block] * density_bases[block]
        for block in eachindex(translation_blocks)
    ]
    retained_by_block = [Int[] for _ in density_blocks]
    for level in eachindex(probabilities)
        probabilities[level] > probability_tol || continue
        block, local_level = locations[level]
        push!(retained_by_block[block], local_level)
    end

    closure_squared = 0.0
    action_squared = 0.0
    retained_densities = Matrix{ComplexF64}[]
    retained_translations = Matrix{ComplexF64}[]
    for block in eachindex(density_blocks)
        retained = sort!(retained_by_block[block])
        isempty(retained) && continue
        translation = density_basis_translations[block]
        discarded = setdiff(axes(translation, 1), retained)
        if !isempty(discarded)
            closure_squared += sum(abs2, translation[discarded, retained])
        end
        action_squared += sum(abs2, translation[:, retained])
        density = density_bases[block]' * density_blocks[block] *
            density_bases[block]
        push!(retained_densities, Matrix{ComplexF64}(density[retained, retained]))
        push!(
            retained_translations,
            Matrix{ComplexF64}(translation[retained, retained]),
        )
    end
    isempty(retained_densities) && throw(
        ArgumentError("probability_tol removed all positive Schmidt support"),
    )
    closure = sqrt(closure_squared) / max(sqrt(action_squared), eps(Float64))
    residuals = _block_validation_residuals(
        retained_densities,
        retained_translations,
        Ny,
    )
    return Float64(closure), residuals
end

function _mpskit_invalid_schmidt_momentum(
    spectrum::MomentumEntanglementData,
    validation::MomentumValidation,
    closure_residual::Float64,
)
    probabilities = getproperty.(spectrum.levels, :probability)
    charges = getproperty.(spectrum.levels, :raw_charge)
    invalid = _invalid_momentum_data(probabilities, charges, validation)
    return MPSKitSchmidtMomentumData(
        invalid,
        closure_residual,
        false,
        validation.reason,
    )
end

function mpskit_schmidt_momentum_data(
    density_blocks,
    translation_blocks;
    Ny,
    raw_charges=nothing,
    residual_tol::Real=1e-10,
    degeneracy_tol::Real=residual_tol,
    transfer_residual::Real=0.0,
    probability_tol::Real=0.0,
)
    circumference, tolerance, degeneracy = _momentum_parameters(
        Ny,
        residual_tol,
        degeneracy_tol,
    )
    threshold = _probability_threshold(probability_tol)
    densities, translations, charges = _block_momentum_inputs(
        density_blocks,
        translation_blocks,
        raw_charges,
        tolerance,
    )
    spectrum = _validate_momentum_block_vectors(
        densities,
        translations;
        Ny=circumference,
        raw_charges=charges,
        residual_tol=tolerance,
        degeneracy_tol=degeneracy,
        transfer_residual,
        probability_tol=threshold,
        _prepared_density=true,
    )
    closure, retained = _mpskit_retained_block_residuals(
        densities,
        translations,
        threshold,
        circumference,
    )
    combined_unitary = max(
        spectrum.validation.unitary_residual,
        retained.unitary,
    )
    combined_power = max(
        spectrum.validation.power_residual,
        retained.power,
    )
    combined_commutator = max(
        spectrum.validation.commutator_residual,
        retained.density_translation,
    )
    failures = String[]
    closure <= tolerance || push!(
        failures,
        "retained Schmidt subspace closure residual exceeds tolerance",
    )
    retained.unitary <= tolerance || push!(
        failures,
        "retained virtual translation unitarity residual exceeds tolerance",
    )
    retained.power <= tolerance || push!(
        failures,
        "retained Ty^Ny residual exceeds tolerance",
    )
    retained.density_translation <= tolerance || push!(
        failures,
        "retained Schmidt-density translation commutator residual exceeds tolerance",
    )
    if !isempty(failures)
        reason = join(failures, "; ")
        validation = MomentumValidation(
            combined_unitary,
            combined_power,
            combined_commutator,
            spectrum.validation.transfer_residual,
            false,
            reason,
        )
        return _mpskit_invalid_schmidt_momentum(
            spectrum,
            validation,
            closure,
        )
    end
    if !spectrum.valid
        return MPSKitSchmidtMomentumData(
            spectrum,
            closure,
            false,
            spectrum.reason,
        )
    end
    return MPSKitSchmidtMomentumData(
        spectrum,
        closure,
        true,
        spectrum.reason,
    )
end

function _mpskit_validate_mixed_states(
    first_state,
    second_state,
    c::InfiniteCylinderConfig,
)
    _validate_mpskit_observable_state(first_state, c; physical_spaces=true)
    _validate_mpskit_observable_state(second_state, c; physical_spaces=true)
    all(
        site -> MPSKit.physicalspace(first_state, site) ==
            MPSKit.physicalspace(second_state, site),
        1:length(first_state),
    ) || throw(
        ArgumentError("mixed-transfer states have different physical spaces"),
    )
    return nothing
end

function _mpskit_mixed_transfer_initial(ket, bra, seed::Integer)
    map_space = TensorKit.:(←)(
        MPSKit.left_virtualspace(bra, 1),
        MPSKit.left_virtualspace(ket, 1),
    )
    initial = similar(ket.AL[1], map_space)
    return _mpskit_randomize_transfer_initial!(initial, seed)
end

function _mpskit_dominant_overlap_transfer(
    ket,
    bra;
    tol::Float64,
    maxiter::Int,
    krylovdim::Int,
    seed::Integer,
)
    initial = _mpskit_mixed_transfer_initial(ket, bra, seed)
    transfer = TensorKit.flip(MPSKit.TransferMatrix(ket.AL, bra.AL))
    available = TensorKit.dim(TensorKit.space(initial))
    if available == 1
        image = transfer(initial)
        denominator = TensorKit.dot(initial, initial)
        value = ComplexF64(TensorKit.dot(initial, image) / denominator)
        residual = _relative_transfer_residual(image, value, initial)
        return (; value, residual, converged=true)
    end
    algorithm = KrylovKit.Arnoldi(;
        tol,
        maxiter,
        eager=true,
        krylovdim,
        verbosity=0,
    )
    values, vectors, info = KrylovKit.eigsolve(
        transfer,
        initial,
        1,
        :LM,
        algorithm,
    )
    isempty(values) && error("dominant mixed-transfer solve returned no eigenvalue")
    isempty(vectors) && error("dominant mixed-transfer solve returned no eigenvector")
    value = ComplexF64(first(values))
    vector = first(vectors)
    residual = _relative_transfer_residual(
        transfer(vector),
        value,
        vector,
    )
    return (; value, residual, converged=info.converged >= 1)
end

function mpskit_mixed_transfer_fidelity(
    first_state,
    second_state,
    c::InfiniteCylinderConfig;
    tol::Real=1e-10,
    maxiter::Integer=500,
    krylovdim::Integer=30,
    seed::Integer=0x4d4958,
)
    isfinite(tol) && tol > 0 || throw(
        ArgumentError("tol must be finite and positive"),
    )
    !(maxiter isa Bool) && maxiter > 0 || throw(
        ArgumentError("maxiter must be positive"),
    )
    !(krylovdim isa Bool) && krylovdim > 1 || throw(
        ArgumentError("krylovdim must be greater than one"),
    )
    seed isa Integer && !(seed isa Bool) || throw(
        ArgumentError("seed must be an integer"),
    )
    _mpskit_validate_mixed_states(first_state, second_state, c)
    tolerance = Float64(tol)
    try
        mixed = _mpskit_dominant_overlap_transfer(
            first_state,
            second_state;
            tol=tolerance,
            maxiter=Int(maxiter),
            krylovdim=Int(krylovdim),
            seed,
        )
        self1 = _mpskit_dominant_overlap_transfer(
            first_state,
            first_state;
            tol=tolerance,
            maxiter=Int(maxiter),
            krylovdim=Int(krylovdim),
            seed=seed + 1,
        )
        self2 = _mpskit_dominant_overlap_transfer(
            second_state,
            second_state;
            tol=tolerance,
            maxiter=Int(maxiter),
            krylovdim=Int(krylovdim),
            seed=seed + 2,
        )
        return _mixed_transfer_result(
            mixed.value,
            self1.value,
            self2.value,
            mixed.residual,
            self1.residual,
            self2.residual,
            mixed.converged && self1.converged && self2.converged;
            residual_tolerance=10tolerance,
            clamp_tolerance=10tolerance,
        )
    catch exception
        return _invalid_mixed_transfer(
            ComplexF64(NaN),
            ComplexF64(NaN),
            ComplexF64(NaN),
            NaN,
            NaN,
            NaN,
            "invalid MPSKit mixed-transfer solve: $(sprint(showerror, exception))",
        )
    end
end

function _mpskit_twisted_transfer_solution(
    state,
    translation;
    tol::Float64,
    maxiter::Int,
    krylovdim::Int,
    seed::Integer,
)
    initial = MPSKit.allocate_GL(state, translation, state, 1)
    _mpskit_randomize_transfer_initial!(initial, seed)
    transfer = TensorKit.flip(
        MPSKit.TransferMatrix(state.AL, parent(translation), state.AL),
    )
    available = TensorKit.dim(TensorKit.space(initial))
    if available == 1
        image = transfer(initial)
        value = ComplexF64(
            TensorKit.dot(initial, image) / TensorKit.dot(initial, initial),
        )
        residual = _relative_transfer_residual(image, value, initial)
        return (;
            value,
            vector=initial,
            residual,
            converged=true,
            gap=Inf,
        )
    end
    requested = min(2, available)
    algorithm = KrylovKit.Arnoldi(;
        tol,
        maxiter,
        eager=true,
        krylovdim=max(krylovdim, requested + 1),
        verbosity=0,
    )
    values, vectors, info = KrylovKit.eigsolve(
        transfer,
        initial,
        requested,
        :LM,
        algorithm,
    )
    length(values) >= requested || error(
        "twisted transfer returned too few eigenvalues",
    )
    length(vectors) >= requested || error(
        "twisted transfer returned too few eigenvectors",
    )
    order = sortperm(values; by=abs, rev=true)
    first_index = order[1]
    value = ComplexF64(values[first_index])
    vector = vectors[first_index]
    residual = _relative_transfer_residual(
        transfer(vector),
        value,
        vector,
    )
    gap = requested == 1 ? Inf : Float64(
        (abs(value) - abs(values[order[2]])) / max(abs(value), eps(Float64)),
    )
    return (;
        value,
        vector,
        residual,
        converged=info.converged >= requested,
        gap,
    )
end

function _mpskit_virtual_translation_at_bond(
    state,
    translation,
    fixed_point,
    bond::Int,
)
    value = fixed_point
    for site in 1:bond
        value = value * MPSKit.TransferMatrix(
            state.AL[site],
            translation[site],
            state.AL[site],
        )
    end
    return TensorKit.removeunit(value, 2)
end

function _mpskit_schmidt_translation_blocks(
    state,
    virtual_translation,
    bond::Int,
    Ny::Int,
)
    density = state.C[bond] * state.C[bond]'
    if TensorKit.space(virtual_translation) != TensorKit.space(density)
        adjoint_translation = virtual_translation'
        TensorKit.space(adjoint_translation) == TensorKit.space(density) || throw(
            ArgumentError(
                "twisted fixed point and Schmidt density use incompatible spaces",
            ),
        )
        virtual_translation = adjoint_translation
    end
    density_by_sector = Dict(
        sector => Matrix{ComplexF64}(block) for
        (sector, block) in TensorKit.blocks(density)
    )
    translation_by_sector = Dict(
        sector => Matrix{ComplexF64}(block) for
        (sector, block) in TensorKit.blocks(virtual_translation)
    )
    sectors = collect(keys(density_by_sector))
    sort!(sectors; by=sector -> begin
        odd, raw = _mpskit_entanglement_sector(sector)
        (raw, odd)
    end)
    all(haskey(translation_by_sector, sector) for sector in sectors) || throw(
        ArgumentError("twisted fixed point is missing a Schmidt charge block"),
    )
    density_blocks = [density_by_sector[sector] for sector in sectors]
    translation_blocks = [translation_by_sector[sector] for sector in sectors]
    raw_charges = [last(_mpskit_entanglement_sector(sector)) for sector in sectors]

    total_dimension = sum(block -> size(block, 1), translation_blocks)
    translation_norm = sqrt(sum(sum(abs2, block) for block in translation_blocks))
    isfinite(translation_norm) && translation_norm > 0 || throw(
        ArgumentError("virtual translation has invalid or zero norm"),
    )
    magnitude_correction = sqrt(Float64(total_dimension)) / translation_norm
    translation_blocks = [magnitude_correction .* block for block in translation_blocks]
    cycle_average = sum(
        LinearAlgebra.tr(block^Ny) for block in translation_blocks
    ) / total_dimension
    phase_correction = iszero(cycle_average) ? 1.0 + 0.0im :
        cis(-angle(cycle_average) / Ny)
    translation_blocks = [phase_correction .* block for block in translation_blocks]
    return (;
        density_blocks,
        translation_blocks,
        raw_charges,
        phase_correction=ComplexF64(magnitude_correction * phase_correction),
    )
end

function _mpskit_invalidate_spectrum(
    spectrum::MomentumEntanglementData,
    reason::AbstractString,
    transfer_residual::Float64,
)
    validation = MomentumValidation(
        spectrum.validation.unitary_residual,
        spectrum.validation.power_residual,
        spectrum.validation.commutator_residual,
        transfer_residual,
        false,
        String(reason),
    )
    return _invalid_momentum_data(
        getproperty.(spectrum.levels, :probability),
        getproperty.(spectrum.levels, :raw_charge),
        validation,
    )
end

function mpskit_momentum_entanglement_data(
    state,
    c::InfiniteCylinderConfig;
    cut_x::Integer,
    tol::Real=1e-10,
    maxiter::Integer=500,
    krylovdim::Integer=30,
    seed::Integer=0x4d4f4d454e54554d,
    residual_tol::Real=10tol,
    degeneracy_tol::Real=residual_tol,
    probability_tol::Real=0.0,
    gap_tol::Real=residual_tol,
)
    _paper_translation_Ny(c)
    !(cut_x isa Bool) && 1 <= cut_x <= c.x_period || throw(
        ArgumentError("cut_x must identify a physical x cut in the unit cell"),
    )
    isfinite(tol) && tol > 0 || throw(
        ArgumentError("tol must be finite and positive"),
    )
    isfinite(residual_tol) && residual_tol > 0 || throw(
        ArgumentError("residual_tol must be finite and positive"),
    )
    isfinite(gap_tol) && gap_tol >= 0 || throw(
        ArgumentError("gap_tol must be finite and nonnegative"),
    )
    !(maxiter isa Bool) && maxiter > 0 || throw(
        ArgumentError("maxiter must be positive"),
    )
    !(krylovdim isa Bool) && krylovdim > 2 || throw(
        ArgumentError("krylovdim must exceed two"),
    )
    seed isa Integer && !(seed isa Bool) || throw(
        ArgumentError("seed must be an integer"),
    )
    _validate_mpskit_observable_state(state, c; physical_spaces=true)
    cut = Int(cut_x)
    bond = cut * c.Ly
    translation = mpskit_transverse_translation_mpo(c)
    solution = _mpskit_twisted_transfer_solution(
        state,
        translation;
        tol=Float64(tol),
        maxiter=Int(maxiter),
        krylovdim=Int(krylovdim),
        seed,
    )
    virtual = _mpskit_virtual_translation_at_bond(
        state,
        translation,
        solution.vector,
        bond,
    )
    blocks = _mpskit_schmidt_translation_blocks(
        state,
        virtual,
        bond,
        c.Ny,
    )
    schmidt = mpskit_schmidt_momentum_data(
        blocks.density_blocks,
        blocks.translation_blocks;
        Ny=c.Ny,
        raw_charges=blocks.raw_charges,
        residual_tol,
        degeneracy_tol,
        transfer_residual=solution.residual,
        probability_tol,
    )
    transfer_failures = String[]
    solution.converged || push!(
        transfer_failures,
        "twisted transfer did not converge",
    )
    solution.residual <= residual_tol || push!(
        transfer_failures,
        "twisted-transfer residual exceeds tolerance",
    )
    magnitude_residual = abs(abs(solution.value) - 1)
    magnitude_residual <= residual_tol || push!(
        transfer_failures,
        "twisted-transfer eigenvalue magnitude differs from one",
    )
    solution.gap > gap_tol || push!(
        transfer_failures,
        "twisted-transfer dominant eigenvalue is unresolved",
    )
    spectrum = schmidt.spectrum
    if !isempty(transfer_failures)
        reason = join(transfer_failures, "; ")
        spectrum = _mpskit_invalidate_spectrum(
            spectrum,
            reason,
            max(solution.residual, magnitude_residual),
        )
    end
    valid = isempty(transfer_failures) && schmidt.valid && spectrum.valid
    reason = valid ? spectrum.reason : join(
        filter(
            value -> !isempty(value),
            [schmidt.reason, join(transfer_failures, "; ")],
        ),
        "; ",
    )
    return MPSKitMomentumData(
        cut,
        bond,
        spectrum,
        schmidt.closure_residual,
        solution.value,
        solution.residual,
        solution.converged,
        solution.gap,
        blocks.phase_correction,
        valid,
        reason,
    )
end
