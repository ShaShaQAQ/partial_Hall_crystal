struct EnergyData
    per_cell::Float64
    per_x::Float64
    per_unit_cell::Float64
    per_site::Float64
end

struct DensityRow
    site::Int
    x::Int
    y::Int
    density::Float64
end

struct EntanglementLevel
    cut_x::Int
    bond::Int
    level::Int
    singular_value::Float64
    probability::Float64
    entanglement_energy::Float64
    qn::QN
    raw_charge::Int
    physical_charge::Float64
end

struct EntanglementSector
    cut_x::Int
    bond::Int
    qn::QN
    raw_charge::Int
    physical_charge::Float64
    weight::Float64
end

struct EntanglementData
    cut_x::Int
    bond::Int
    entropy::Float64
    raw_schmidt_polarization::Float64
    levels::Vector{EntanglementLevel}
    sectors::Vector{EntanglementSector}
end

struct NeutralTransferData
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

function normalize_energy(c::InfiniteCylinderConfig, energy_per_cell::Real)
    energy = Float64(energy_per_cell)
    isfinite(energy) ||
        throw(ArgumentError("energy per cell must be finite after Float64 conversion"))
    return EnergyData(
        energy,
        energy / c.x_period,
        energy / unit_cells_per_cell(c),
        energy / sites_per_cell(c),
    )
end

function _validate_observable_cell(psi, c::InfiniteCylinderConfig)
    nsites(psi) == sites_per_cell(c) || throw(
        ArgumentError("MPS cell size does not match the cylinder configuration")
    )
    expected_signature = _configuration_signature(c)
    sites = siteinds(only, psi.AL)
    all(site -> hastags(site, expected_signature), sites) || throw(
        ArgumentError("MPS site indices do not match the cylinder configuration")
    )
    return nothing
end

function _observable_energy_per_cell(
    estimates;
    imaginary_tol::Real=1e-12,
)
    tolerance = _validate_imaginary_tolerance(imaginary_tol)
    total = sum(estimates)
    return _numerically_real_energy(total, tolerance)
end

function energy_data(
    psi,
    H,
    c::InfiniteCylinderConfig;
    imaginary_tol::Real=1e-12,
)
    _validate_observable_cell(psi, c)
    return normalize_energy(
        c,
        _observable_energy_per_cell(expect(psi, H); imaginary_tol),
    )
end

function density_data(psi, c::InfiniteCylinderConfig)
    _validate_observable_cell(psi, c)
    return [
        DensityRow(site, site_coordinates(c, site)..., real(expect(psi, "N", site))) for
        site in 1:sites_per_cell(c)
    ]
end

function _normalized_schmidt_probabilities(singular_values)
    values = Float64.(singular_values)
    isempty(values) && throw(ArgumentError("Schmidt spectrum is empty"))
    all(isfinite, values) ||
        throw(ArgumentError("Schmidt singular values must be finite"))
    all(>=(0.0), values) ||
        throw(ArgumentError("Schmidt singular values must be nonnegative"))
    squared_norm = sum(abs2, values)
    isfinite(squared_norm) && squared_norm > 0 ||
        throw(ArgumentError("Schmidt spectrum has invalid or zero norm"))
    probabilities = abs2.(values) ./ squared_norm
    all(isfinite, probabilities) ||
        throw(ArgumentError("Schmidt probabilities must be finite"))
    abs(sum(probabilities) - 1.0) <= 1e-10 ||
        error("Schmidt probabilities are not normalized to tolerance 1e-10")
    return probabilities
end

function _left_schmidt_index(psi::InfiniteCanonicalMPS, bond::Int, C, S, V)
    left_center_index = only(commoninds(C, psi.AL[bond]))
    left_singular_index = only(uniqueinds(S, V))
    dir(left_singular_index) == dir(left_center_index) ||
        error("left Schmidt index orientation does not match the MPS left bond")
    return left_singular_index
end

function _schmidt_level_qns(schmidt_index)
    level_qns = QN[]
    for block in 1:nblocks(schmidt_index)
        append!(
            level_qns,
            fill(qn(schmidt_index, block), blockdim(schmidt_index, block)),
        )
    end
    length(level_qns) == dim(schmidt_index) ||
        error("Schmidt QN block dimensions do not match the index dimension")
    return level_qns
end

function _schmidt_sector_weights(level_qns, probabilities)
    length(level_qns) == length(probabilities) || throw(
        ArgumentError("Schmidt QN and probability counts must match")
    )
    weights = Dict{QN,Float64}()
    for (sector_qn, probability) in zip(level_qns, probabilities)
        weights[sector_qn] = get(weights, sector_qn, 0.0) + probability
    end
    return weights
end

function entanglement_data(
    psi::InfiniteCanonicalMPS,
    c::InfiniteCylinderConfig;
    cut_x::Integer,
)
    !(cut_x isa Bool) && 1 <= cut_x <= c.x_period || throw(
        ArgumentError("cut_x must identify a transverse bond in the reference cell")
    )
    _validate_observable_cell(psi, c)

    cut = Int(cut_x)
    bond = Base.Checked.checked_mul(cut, c.Ly)
    C = psi.C[bond]
    left_center_index = only(commoninds(C, psi.AL[bond]))
    _, S, V = svd(C, left_center_index)
    left_singular_index = _left_schmidt_index(psi, bond, C, S, V)
    singular_values = [Float64(real(S[level, level])) for level in 1:dim(left_singular_index)]
    probabilities = _normalized_schmidt_probabilities(singular_values)

    levels = EntanglementLevel[]
    level_qns = _schmidt_level_qns(left_singular_index)
    sector_weights = _schmidt_sector_weights(level_qns, probabilities)
    scale = charge_scale(c)
    for (level, (sector_qn, singular_value, probability)) in enumerate(
        zip(level_qns, singular_values, probabilities)
    )
        raw_charge = val(sector_qn, "Nf")
        physical_charge = raw_charge / scale
        push!(
            levels,
            EntanglementLevel(
                cut,
                bond,
                level,
                singular_value,
                probability,
                iszero(probability) ? Inf : -log(probability),
                sector_qn,
                raw_charge,
                physical_charge,
            ),
        )
    end

    sectors = [
        EntanglementSector(
            cut,
            bond,
            sector_qn,
            val(sector_qn, "Nf"),
            val(sector_qn, "Nf") / scale,
            weight,
        ) for (sector_qn, weight) in sector_weights
    ]
    sort!(sectors; by=sector -> sector.raw_charge)
    entropy = -sum(
        probability * log(probability) for probability in probabilities if !iszero(probability)
    )
    raw_schmidt_polarization =
        sum(level.probability * level.physical_charge for level in levels)
    return EntanglementData(
        cut,
        bond,
        entropy,
        raw_schmidt_polarization,
        levels,
        sectors,
    )
end

function _invalid_neutral_transfer(
    values,
    residual_norms,
    converged::Integer,
    reason::AbstractString,
)
    eigenvalues = ComplexF64.(values)
    phases = Float64[angle(value) for value in eigenvalues]
    return NeutralTransferData(
        eigenvalues,
        phases,
        Float64.(residual_norms),
        Int(converged),
        false,
        String(reason),
        NaN,
        NaN,
        NaN,
    )
end

function _neutral_transfer_result(
    c::InfiniteCylinderConfig,
    values,
    residual_norms,
    converged::Integer;
    ratio_tolerance::Real=1e-10,
    residual_tolerance::Real=1e-10,
)
    isfinite(ratio_tolerance) && ratio_tolerance >= 0 ||
        throw(ArgumentError("ratio_tolerance must be finite and nonnegative"))
    isfinite(residual_tolerance) && residual_tolerance > 0 ||
        throw(ArgumentError("residual_tolerance must be finite and positive"))
    eigenvalues = ComplexF64.(values)
    residuals = Float64.(residual_norms)
    !(converged isa Bool) && converged >= 0 ||
        throw(ArgumentError("converged eigenpair count must be nonnegative"))
    length(residuals) == length(eigenvalues) ||
        throw(ArgumentError("transfer residual count must match eigenvalue count"))
    converged <= length(eigenvalues) || throw(
        ArgumentError("converged eigenpair count exceeds returned eigenvalue count")
    )
    all(isfinite, eigenvalues) || return _invalid_neutral_transfer(
        eigenvalues, residuals, converged, "transfer eigenvalues must be finite"
    )
    all(residual -> isfinite(residual) && residual >= 0, residuals) ||
        return _invalid_neutral_transfer(
            eigenvalues,
            residuals,
            converged,
            "transfer residual norms must be finite and nonnegative",
        )
    converged >= 2 && length(eigenvalues) >= 2 ||
        return _invalid_neutral_transfer(
            eigenvalues,
            residuals,
            converged,
            "fewer than two neutral transfer eigenpairs converged",
        )
    all(<=(residual_tolerance), residuals[1:2]) ||
        return _invalid_neutral_transfer(
            eigenvalues,
            residuals,
            converged,
            "leading neutral transfer eigenpair residual exceeds tolerance",
        )

    ratio = abs(eigenvalues[2] / eigenvalues[1])
    isfinite(ratio) || return _invalid_neutral_transfer(
        eigenvalues, residuals, converged, "neutral transfer ratio must be finite"
    )
    ratio > 0 || return _invalid_neutral_transfer(
        eigenvalues, residuals, converged, "neutral transfer ratio must be positive"
    )
    ratio < 1 - ratio_tolerance || return _invalid_neutral_transfer(
        eigenvalues,
        residuals,
        converged,
        "neutral transfer ratio must be strictly below one outside tolerance",
    )

    xi_cell = -1 / log(ratio)
    xi_x = c.x_period * xi_cell
    all(isfinite, (xi_cell, xi_x)) || return _invalid_neutral_transfer(
        eigenvalues,
        residuals,
        converged,
        "neutral transfer correlation length must be finite",
    )
    return NeutralTransferData(
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

function neutral_transfer_data(
    psi::InfiniteCanonicalMPS,
    c::InfiniteCylinderConfig;
    neigs::Integer=4,
    tol::Real=1e-10,
    ratio_tolerance::Real=1e-10,
    rng::AbstractRNG=Random.default_rng(),
    kwargs...,
)
    !(neigs isa Bool) && neigs >= 2 ||
        throw(ArgumentError("neigs must be an integer of at least two"))
    isfinite(tol) && tol > 0 ||
        throw(ArgumentError("tol must be finite and positive"))
    _validate_observable_cell(psi, c)

    transfer = TransferMatrix(psi.AL)
    initial = random_itensor(rng, QN(), dag(input_inds(transfer)))
    flux(initial) == QN() || error("neutral transfer initial tensor has nonzero flux")
    values, vectors, info = KrylovKit.eigsolve(
        transfer,
        initial,
        Int(neigs),
        :LM;
        tol,
        verbosity=0,
        kwargs...,
    )
    explicit_residual_norms = Float64[
        norm(transfer(vector) - value * vector) for
        (value, vector) in zip(values, vectors)
    ]
    return _neutral_transfer_result(
        c,
        values,
        explicit_residual_norms,
        info.converged;
        ratio_tolerance,
        residual_tolerance=tol,
    )
end
