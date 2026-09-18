function _fermionic_translation_shape(Ny::Integer)
    Ny isa Bool && throw(ArgumentError("Ny must be a positive integer"))
    ny = BigInt(Ny)
    ny > 0 || throw(ArgumentError("Ny must be a positive integer"))
    ny <= typemax(Int) ÷ 2 || throw(ArgumentError("Ny is too large"))
    mode_count = 2 * Int(ny)
    dimension = BigInt(1) << mode_count
    dimension <= typemax(Int) || throw(
        ArgumentError("fermionic Fock-space dimension does not fit Int")
    )
    return mode_count, Int(dimension)
end

function translated_fock_basis(bits::Integer, Ny::Integer)
    bits isa Bool && throw(ArgumentError("bits must be a nonnegative integer"))
    mode_count, dimension = _fermionic_translation_shape(Ny)
    bit_pattern = BigInt(bits)
    0 <= bit_pattern < dimension || throw(
        ArgumentError("bits lie outside the $mode_count-mode Fock basis")
    )
    source = Int(bit_pattern)
    target = 0
    translated_modes = Int[]

    for mode in 0:(mode_count - 1)
        iszero(source & (1 << mode)) && continue
        translated_mode = mod(mode + 2, mode_count)
        target |= 1 << translated_mode
        push!(translated_modes, translated_mode)
    end

    inversions = 0
    for first in 1:(length(translated_modes) - 1)
        for second in (first + 1):length(translated_modes)
            inversions += translated_modes[first] > translated_modes[second]
        end
    end
    return target, isodd(inversions) ? -1 : 1
end

function fermionic_ring_translation(Ny::Integer)
    _, dimension = _fermionic_translation_shape(Ny)
    rows = Vector{Int}(undef, dimension)
    values = Vector{ComplexF64}(undef, dimension)
    for bits in 0:(dimension - 1)
        target, sign = translated_fock_basis(bits, Ny)
        rows[bits + 1] = target + 1
        values[bits + 1] = sign
    end
    return SparseArrays.sparse(rows, 1:dimension, values, dimension, dimension)
end

function _paper_translation_Ny(c::InfiniteCylinderConfig)
    c.geometry == :paper_straight || throw(
        ArgumentError(
            "one-cell transverse translation requires paper_straight geometry"
        )
    )
    return c.Ny
end

fermionic_ring_translation(c::InfiniteCylinderConfig) =
    fermionic_ring_translation(_paper_translation_Ny(c))

function apply_fermionic_ring_translation(
    coefficients::AbstractVector,
    Ny::Integer,
)
    _, dimension = _fermionic_translation_shape(Ny)
    length(coefficients) == dimension || throw(
        DimensionMismatch(
            "coefficient vector has length $(length(coefficients)); expected $dimension"
        )
    )
    translated = similar(coefficients, dimension)
    for bits in 0:(dimension - 1)
        target, sign = translated_fock_basis(bits, Ny)
        translated[target + 1] = sign * coefficients[bits + 1]
    end
    return translated
end

apply_fermionic_ring_translation(
    coefficients::AbstractVector,
    c::InfiniteCylinderConfig,
) = apply_fermionic_ring_translation(coefficients, _paper_translation_Ny(c))
