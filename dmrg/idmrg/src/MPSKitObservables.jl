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
