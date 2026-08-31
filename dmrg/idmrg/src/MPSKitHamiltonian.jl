const _MPSKIT_FERMION_OPERATORS = TensorKitTensors.FermionOperators

function _mpskit_uncentered_physical_space()
    return _MPSKIT_FERMION_OPERATORS.fermion_space(TensorKit.U1Irrep)
end

function _mpskit_centered_basis_transform(physical_space)
    uncentered_space = _mpskit_uncentered_physical_space()
    dense_centered = TensorKitTensors.desymmetrize(physical_space)
    dense_uncentered = TensorKitTensors.desymmetrize(uncentered_space)
    transform_space = TensorKit.:(←)(dense_centered, dense_uncentered)
    return TensorKit.TensorMap(Matrix{Int}(I, 2, 2), transform_space)
end

function _mpskit_fermion_operators(physical_space; centered::Bool)
    number_operator = _MPSKIT_FERMION_OPERATORS.f_num(
        ComplexF64,
        TensorKit.U1Irrep,
    )
    plus_min = _MPSKIT_FERMION_OPERATORS.f_plus_f_min(
        ComplexF64,
        TensorKit.U1Irrep,
    )
    min_plus = _MPSKIT_FERMION_OPERATORS.f_min_f_plus(
        ComplexF64,
        TensorKit.U1Irrep,
    )
    centered || return (
        number=number_operator,
        plus_min=plus_min,
        min_plus=min_plus,
    )

    transform = _mpskit_centered_basis_transform(physical_space)
    return (
        number=TensorKitTensors.symmetrize(
            number_operator,
            transform,
            physical_space,
        ),
        plus_min=TensorKitTensors.symmetrize(
            plus_min,
            transform,
            physical_space,
        ),
        min_plus=TensorKitTensors.symmetrize(
            min_plus,
            transform,
            physical_space,
        ),
    )
end

function _mpskit_validate_site(site::Int, nsites::Int, label::AbstractString)
    1 <= site <= nsites || throw(
        ArgumentError("$label site $site must lie in 1:$nsites")
    )
    return site
end

function _mpskit_local_terms(
    nsites::Int,
    hoppings::AbstractVector{<:InfiniteHoppingTerm},
    interactions::AbstractVector{<:InfiniteInteractionTerm},
    operators;
    atol::Real,
)
    onsite = Dict{Int,ComplexF64}()
    directed = Dict{Tuple{Int,Int},ComplexF64}()
    for hopping in hoppings
        target = _mpskit_validate_site(hopping.target, nsites, "hopping target")
        source = _mpskit_validate_site(hopping.source, nsites, "hopping source")
        isfinite(real(hopping.amp)) && isfinite(imag(hopping.amp)) || throw(
            ArgumentError("hopping amplitudes must be finite")
        )
        if target == source
            onsite[target] = get(onsite, target, 0.0 + 0.0im) + hopping.amp
        else
            key = (target, source)
            directed[key] = get(directed, key, 0.0 + 0.0im) + hopping.amp
        end
    end

    local_terms = Any[]
    for site in sort!(collect(keys(onsite)))
        amplitude = onsite[site]
        abs(imag(amplitude)) <= atol || throw(
            ArgumentError("onsite hopping at site $site must be real")
        )
        abs(real(amplitude)) <= atol && continue
        push!(local_terms, (site,) => real(amplitude) * operators.number)
    end

    unordered_pairs = sort!(
        unique(minmax(target, source) for (target, source) in keys(directed)),
    )
    for (first_site, second_site) in unordered_pairs
        forward_key = (first_site, second_site)
        reverse_key = (second_site, first_site)
        haskey(directed, forward_key) && haskey(directed, reverse_key) || throw(
            ArgumentError(
                "hopping pair ($first_site, $second_site) is missing its Hermitian reverse",
            ),
        )
        amplitude = directed[forward_key]
        reverse_amplitude = directed[reverse_key]
        isapprox(reverse_amplitude, conj(amplitude); atol, rtol=0) || throw(
            ArgumentError(
                "hopping pair ($first_site, $second_site) is not Hermitian",
            ),
        )
        abs(amplitude) <= atol && abs(reverse_amplitude) <= atol && continue
        hopping_operator = amplitude * operators.plus_min -
            conj(amplitude) * operators.min_plus
        push!(local_terms, (first_site, second_site) => hopping_operator)
    end

    interaction_strengths = Dict{Tuple{Int,Int},Float64}()
    for interaction in interactions
        first_site = _mpskit_validate_site(
            interaction.i,
            nsites,
            "interaction first",
        )
        second_site = _mpskit_validate_site(
            interaction.j,
            nsites,
            "interaction second",
        )
        first_site != second_site || throw(
            ArgumentError("density interactions require two distinct sites")
        )
        isfinite(interaction.V) || throw(
            ArgumentError("interaction strengths must be finite")
        )
        key = minmax(first_site, second_site)
        interaction_strengths[key] = get(interaction_strengths, key, 0.0) +
            interaction.V
    end
    density_density = TensorKit.:(⊗)(operators.number, operators.number)
    for sites in sort!(collect(keys(interaction_strengths)))
        interaction = interaction_strengths[sites]
        abs(interaction) <= atol && continue
        push!(local_terms, sites => interaction * density_density)
    end
    return local_terms
end

function _mpskit_finite_hamiltonian(
    physical_spaces,
    hoppings::AbstractVector{<:InfiniteHoppingTerm},
    interactions::AbstractVector{<:InfiniteInteractionTerm};
    centered::Bool,
    atol::Real,
)
    nsites = length(physical_spaces)
    nsites >= 2 || throw(
        ArgumentError("finite Hamiltonian requires at least two sites")
    )
    all(==(first(physical_spaces)), physical_spaces) || throw(
        ArgumentError("finite Hamiltonian requires a uniform physical space")
    )
    operators = _mpskit_fermion_operators(first(physical_spaces); centered)
    local_terms = _mpskit_local_terms(
        nsites,
        hoppings,
        interactions,
        operators;
        atol,
    )
    isempty(local_terms) && throw(
        ArgumentError("finite Hamiltonian term list cannot be empty")
    )
    return MPSKit.FiniteMPOHamiltonian(physical_spaces, local_terms)
end

function mpskit_finite_hamiltonian(
    nsites::Int,
    hoppings::AbstractVector{<:InfiniteHoppingTerm},
    interactions::AbstractVector{<:InfiniteInteractionTerm};
    charge_mode::Symbol=:uncentered,
    atol::Real=1e-12,
)
    charge_mode == :uncentered || throw(
        ArgumentError("centered finite Hamiltonians require a configuration")
    )
    physical_space = _mpskit_uncentered_physical_space()
    return _mpskit_finite_hamiltonian(
        fill(physical_space, nsites),
        hoppings,
        interactions;
        centered=false,
        atol,
    )
end

function mpskit_finite_hamiltonian(
    c::InfiniteCylinderConfig,
    hoppings::AbstractVector{<:InfiniteHoppingTerm},
    interactions::AbstractVector{<:InfiniteInteractionTerm};
    charge_mode::Symbol=:centered,
    atol::Real=1e-12,
)
    charge_mode in (:centered, :uncentered) || throw(
        ArgumentError("charge_mode must be :centered or :uncentered")
    )
    physical_spaces = if charge_mode == :centered
        mpskit_physical_spaces(c)
    else
        fill(_mpskit_uncentered_physical_space(), sites_per_cell(c))
    end
    return _mpskit_finite_hamiltonian(
        physical_spaces,
        hoppings,
        interactions;
        centered=charge_mode == :centered,
        atol,
    )
end

function mpskit_dense_matrix(hamiltonian)
    tensor = convert(TensorKit.TensorMap, hamiltonian)
    dense_tensor = TensorKitTensors.desymmetrize(tensor)
    dense_array = convert(Array, dense_tensor)
    rows = TensorKit.dim(TensorKit.codomain(dense_tensor))
    columns = TensorKit.dim(TensorKit.domain(dense_tensor))
    return reshape(dense_array, rows, columns)
end
