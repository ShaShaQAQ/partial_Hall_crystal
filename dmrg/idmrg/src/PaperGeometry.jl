const PAPER_A1 = (1.0, 0.0)
const PAPER_A2 = (0.0, sqrt(3.0))
const PAPER_SUBLATTICE_POSITIONS = ((0.0, 0.0), (0.5, sqrt(3.0) / 2))

struct PaperHoppingOrbit
    dx::Int
    dy::Int
    target_orbital::Int
    source_orbital::Int
    amplitude::ComplexF64
end

function paper_triangular_coordinates(c::InfiniteCylinderConfig, site::Int)
    x, y, orbital = paper_coordinates(c, site)
    orbital_offset = orbital === :A ? 0 : 1
    return x - y, 2y + orbital_offset
end

function paper_physical_position(c::InfiniteCylinderConfig, site::Int)
    x, y, orbital = paper_coordinates(c, site)
    basis = PAPER_SUBLATTICE_POSITIONS[orbital === :A ? 1 : 2]
    return (
        x * PAPER_A1[1] + y * PAPER_A2[1] + basis[1],
        x * PAPER_A1[2] + y * PAPER_A2[2] + basis[2],
    )
end

function paper_periodic_distance_squared(
    c::InfiniteCylinderConfig,
    first_site::Int,
    second_site::Int,
)
    c.geometry == :paper_straight || throw(
        ArgumentError("paper distances require paper_straight geometry")
    )
    first_position = paper_physical_position(c, first_site)
    second_position = paper_physical_position(c, second_site)
    circumference = c.Ny * PAPER_A2[2]
    return minimum(
        (first_position[1] - second_position[1])^2 +
            (first_position[2] - second_position[2] + winding * circumference)^2 for
        winding in -1:1
    )
end

function paper_hopping_orbits(
    params::CylinderModelParams;
    atol::Real=1e-12,
)
    tR = hopping_matrices_from_get_Hk(params.t1, params.t3)
    orbits = PaperHoppingOrbit[]
    for ((legacy_dx, legacy_dy), matrix) in tR
        paper_dx = legacy_dx + legacy_dy
        paper_dy = legacy_dy
        for target_orbital in 1:2, source_orbital in 1:2
            amplitude = matrix[target_orbital, source_orbital]
            abs(amplitude) < atol && continue
            push!(
                orbits,
                PaperHoppingOrbit(
                    paper_dx,
                    paper_dy,
                    target_orbital,
                    source_orbital,
                    amplitude,
                ),
            )
        end
    end
    sort!(
        orbits;
        by=orbit -> (
            orbit.source_orbital,
            orbit.target_orbital,
            orbit.dx,
            orbit.dy,
            real(orbit.amplitude),
            imag(orbit.amplitude),
        ),
    )
    return orbits
end

function paper_orbits_are_hermitian(
    orbits::AbstractVector{<:PaperHoppingOrbit};
    atol::Real=1e-12,
)
    amplitudes = Dict(
        (orbit.dx, orbit.dy, orbit.target_orbital, orbit.source_orbital) =>
            orbit.amplitude for orbit in orbits
    )
    for orbit in orbits
        reverse_key = (
            -orbit.dx,
            -orbit.dy,
            orbit.source_orbital,
            orbit.target_orbital,
        )
        haskey(amplitudes, reverse_key) || return false
        isapprox(
            amplitudes[reverse_key],
            conj(orbit.amplitude);
            atol,
            rtol=0,
        ) || return false
    end
    return true
end

function paper_bloch_from_orbits(k, params::CylinderModelParams)
    length(k) == 2 || throw(ArgumentError("momentum must have two components"))
    momentum = (Float64(k[1]), Float64(k[2]))
    hamiltonian = zeros(ComplexF64, 2, 2)
    for orbit in paper_hopping_orbits(params)
        target_position = PAPER_SUBLATTICE_POSITIONS[orbit.target_orbital]
        source_position = PAPER_SUBLATTICE_POSITIONS[orbit.source_orbital]
        displacement = (
            orbit.dx * PAPER_A1[1] + orbit.dy * PAPER_A2[1] +
                target_position[1] - source_position[1],
            orbit.dx * PAPER_A1[2] + orbit.dy * PAPER_A2[2] +
                target_position[2] - source_position[2],
        )
        phase = cis(momentum[1] * displacement[1] + momentum[2] * displacement[2])
        hamiltonian[orbit.target_orbital, orbit.source_orbital] +=
            orbit.amplitude * phase
    end
    return hamiltonian
end

function paper_reference_Hk(k, params::CylinderModelParams)
    length(k) == 2 || throw(ArgumentError("momentum must have two components"))
    momentum = Float64[Float64(k[1]), Float64(k[2])]
    return CylinderModelCore.get_Hk(momentum, params.t1, params.t3)
end
