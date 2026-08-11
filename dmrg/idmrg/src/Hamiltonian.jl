struct InfiniteHoppingTerm
    target::Int
    source::Int
    amp::ComplexF64
end

struct InfiniteInteractionTerm
    i::Int
    j::Int
    V::Float64
    shell::Int
end

function normalize_directed(c::InfiniteCylinderConfig, target::Int, source::Int)
    n = sites_per_cell(c)
    shift = fld(min(target, source) - 1, n) * n
    return target - shift, source - shift
end

function normalize_undirected(c::InfiniteCylinderConfig, i::Int, j::Int)
    i, j = normalize_directed(c, i, j)
    return minmax(i, j)
end

function build_infinite_hoppings(
    c::InfiniteCylinderConfig,
    params::CylinderModelParams;
    atol::Real=1e-12,
)
    tR = hopping_matrices_from_get_Hk(params.t1, params.t3)
    accum = Dict{Tuple{Int,Int},ComplexF64}()

    for source in 1:sites_per_cell(c)
        source_x, source_y = site_coordinates(c, source)
        source_orb0 = mod(source_y, 2)
        source_cell_y = (source_y - source_orb0) ÷ 2

        for ((dx_uc, dy_uc), tmat) in tR
            target_x = source_x + dx_uc
            target_cell_y_unwrapped = source_cell_y + dy_uc
            winding_y = fld(target_cell_y_unwrapped, c.Ly ÷ 2)
            target_cell_y = mod(target_cell_y_unwrapped, c.Ly ÷ 2)

            for target_orb in 1:2
                amp0 = tmat[target_orb, source_orb0 + 1]
                abs(amp0) < atol && continue
                target_y = 2 * target_cell_y + target_orb - 1
                target = linear_site(c, target_x, target_y)
                key = normalize_directed(c, target, source)
                amp = amp0 * cis(winding_y * c.phi_y)
                accum[key] = get(accum, key, 0.0 + 0.0im) + amp
            end
        end
    end

    hops = InfiniteHoppingTerm[]
    for ((target, source), amp) in accum
        abs(amp) < atol && continue
        push!(hops, InfiniteHoppingTerm(target, source, amp))
    end
    sort!(hops; by=h -> (h.source, h.target, real(h.amp), imag(h.amp)))
    return hops
end

function build_infinite_interactions(
    c::InfiniteCylinderConfig,
    params::CylinderModelParams;
    atol::Real=1e-14,
)
    seen = Set{Tuple{Int,Int,Int}}()
    interactions = InfiniteInteractionTerm[]

    for (shell, V) in ((1, params.V1), (2, params.V2), (3, params.V3))
        abs(V) < atol && continue
        for source in 1:sites_per_cell(c)
            source_x, source_y = site_coordinates(c, source)
            for (dx, dy) in NB_DISPS[shell]
                target = linear_site(c, source_x + dx, mod(source_y + dy, c.Ly))
                target == source && continue
                i, j = normalize_undirected(c, target, source)
                key = (i, j, shell)
                key in seen && continue
                push!(seen, key)
                push!(interactions, InfiniteInteractionTerm(i, j, V, shell))
            end
        end
    end

    sort!(interactions; by=v -> (v.shell, v.i, v.j, v.V))
    return interactions
end

function build_infinite_model_terms(
    c::InfiniteCylinderConfig,
    params::CylinderModelParams,
)
    hops = build_infinite_hoppings(c, params)
    hamiltonian_is_hermitian(c, hops) ||
        error("translation-normalized hopping Hamiltonian is not Hermitian")
    return hops, build_infinite_interactions(c, params)
end

function hamiltonian_is_hermitian(
    c::InfiniteCylinderConfig,
    hops::AbstractVector{<:InfiniteHoppingTerm};
    atol::Real=1e-10,
)
    amplitudes = Dict((h.target, h.source) => h.amp for h in hops)
    for h in hops
        reverse_key = normalize_directed(c, h.source, h.target)
        haskey(amplitudes, reverse_key) || return false
        isapprox(amplitudes[reverse_key], conj(h.amp); atol, rtol=0) || return false
    end
    return true
end

const CanonicalTermKey = Tuple{Symbol,Int,Int,Int}

struct CanonicalTermDictionary <: AbstractDict{CanonicalTermKey,ComplexF64}
    data::Dict{CanonicalTermKey,ComplexF64}
end

Base.length(d::CanonicalTermDictionary) = length(d.data)
Base.iterate(d::CanonicalTermDictionary, state...) = iterate(d.data, state...)
Base.getindex(d::CanonicalTermDictionary, key) = getindex(d.data, key)
Base.haskey(d::CanonicalTermDictionary, key) = haskey(d.data, key)
Base.keys(d::CanonicalTermDictionary) = keys(d.data)

function Base.isapprox(
    a::CanonicalTermDictionary,
    b::CanonicalTermDictionary;
    kwargs...,
)
    keys(a.data) == keys(b.data) || return false
    return all(isapprox(a[key], b[key]; kwargs...) for key in keys(a.data))
end

function canonical_term_dict(
    hops::AbstractVector{<:InfiniteHoppingTerm},
    interactions::AbstractVector{<:InfiniteInteractionTerm},
)
    terms = Dict{CanonicalTermKey,ComplexF64}()
    for h in hops
        terms[(:hopping, h.target, h.source, 0)] = h.amp
    end
    for v in interactions
        terms[(:interaction, v.i, v.j, v.shell)] = ComplexF64(v.V)
    end
    return CanonicalTermDictionary(terms)
end

function tile_to_finite_window(
    c::InfiniteCylinderConfig,
    hops::AbstractVector{<:InfiniteHoppingTerm},
    interactions::AbstractVector{<:InfiniteInteractionTerm};
    Lx::Int,
)
    Lx > 0 || throw(ArgumentError("Lx must be positive"))
    ncell = sites_per_cell(c)
    nsites_window = Base.Checked.checked_mul(Lx, c.Ly)

    hopping_accum = Dict{Tuple{Int,Int},ComplexF64}()
    for h in hops
        first_translation = max(cld(1 - h.target, ncell), cld(1 - h.source, ncell))
        last_translation = min(
            fld(nsites_window - h.target, ncell),
            fld(nsites_window - h.source, ncell),
        )
        for cell_shift in first_translation:last_translation
            target = h.target + cell_shift * ncell
            source = h.source + cell_shift * ncell
            key = (target, source)
            hopping_accum[key] = get(hopping_accum, key, 0.0 + 0.0im) + h.amp
        end
    end

    interaction_accum = Dict{Tuple{Int,Int,Int},Float64}()
    for v in interactions
        first_translation = max(cld(1 - v.i, ncell), cld(1 - v.j, ncell))
        last_translation = min(
            fld(nsites_window - v.i, ncell),
            fld(nsites_window - v.j, ncell),
        )
        for cell_shift in first_translation:last_translation
            i = v.i + cell_shift * ncell
            j = v.j + cell_shift * ncell
            key = (i, j, v.shell)
            interaction_accum[key] = get(interaction_accum, key, 0.0) + v.V
        end
    end

    finite_hops = [
        InfiniteHoppingTerm(target, source, amp) for
        ((target, source), amp) in hopping_accum
    ]
    finite_interactions = [
        InfiniteInteractionTerm(i, j, V, shell) for
        ((i, j, shell), V) in interaction_accum
    ]
    sort!(finite_hops; by=h -> (h.source, h.target, real(h.amp), imag(h.amp)))
    sort!(finite_interactions; by=v -> (v.shell, v.i, v.j, v.V))
    return finite_hops, finite_interactions
end

const _finite_parity_model = Ref{Union{Nothing,Module}}(nothing)

function _finite_model_module()
    finite_model = _finite_parity_model[]
    finite_model === nothing || return finite_model

    sandbox = Module(:InfiniteCylinderFiniteParity)
    Base.include(sandbox, normpath(joinpath(@__DIR__, "..", "..", "DMRGFullModel.jl")))
    finite_model = Base.invokelatest(() -> getfield(sandbox, :DMRGFullModel))
    _finite_parity_model[] = finite_model
    return finite_model
end

function _dictionaries_are_approx(a::AbstractDict, b::AbstractDict; atol::Real)
    keys(a) == keys(b) || return false
    return all(isapprox(a[key], b[key]; atol, rtol=0) for key in keys(a))
end

function finite_term_parity(
    c::InfiniteCylinderConfig,
    params::CylinderModelParams;
    Lx::Int,
    atol::Real=1e-10,
)
    hops, interactions = build_infinite_model_terms(c, params)
    tiled_hops, tiled_interactions =
        tile_to_finite_window(c, hops, interactions; Lx)

    finite = _finite_model_module()
    finite_hops, finite_interactions = Base.invokelatest() do
        finite_lat = getfield(finite, :cylinder_lattice)(Lx; Ly=c.Ly)
        finite_params = getfield(finite, :CylinderModelParams)(;
            t1=params.t1,
            t3=params.t3,
            V1=params.V1,
            V2=params.V2,
            V3=params.V3,
        )
        return getfield(finite, :build_model_terms)(
            finite_lat,
            finite_params;
            phi_y=c.phi_y,
        )
    end

    tiled_hopping_dict = Dict((h.target, h.source) => h.amp for h in tiled_hops)
    finite_hopping_dict = Dict((h.target, h.source) => h.amp for h in finite_hops)
    tiled_interaction_dict =
        Dict((v.i, v.j, v.shell) => v.V for v in tiled_interactions)
    finite_interaction_dict =
        Dict((v.i, v.j, v.shell) => v.V for v in finite_interactions)

    return _dictionaries_are_approx(tiled_hopping_dict, finite_hopping_dict; atol) &&
        _dictionaries_are_approx(
            tiled_interaction_dict,
            finite_interaction_dict;
            atol,
        )
end

function infinite_opsum(c::InfiniteCylinderConfig, params::CylinderModelParams)
    hops, interactions = build_infinite_model_terms(c, params)
    os = OpSum()
    for h in hops
        if h.target == h.source
            os += h.amp, "N", h.source
        else
            os += h.amp, "Cdag", h.target, "C", h.source
        end
    end
    for v in interactions
        os += v.V, "N", v.i, "N", v.j
    end

    starts = Set(minimum(ITensors.sites(term)) for term in ITensors.terms(os))
    for j in 1:sites_per_cell(c)
        if !(j in starts)
            os += 0.0, "Id", j
        end
    end
    return os
end

build_infinite_mpo(c::InfiniteCylinderConfig, params::CylinderModelParams, sites) =
    InfiniteSum{MPO}(infinite_opsum(c, params), sites)
