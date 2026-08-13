struct InfiniteCylinderConfig
    Ly::Int
    x_period::Int
    filling_num::Int
    filling_den::Int
    phi_y::Float64

    function InfiniteCylinderConfig(
        Ly::Int,
        x_period::Int,
        filling_num::Int,
        filling_den::Int,
        phi_y::Float64,
    )
        Ly > 0 && iseven(Ly) || throw(ArgumentError("Ly must be positive and even"))
        x_period > 0 || throw(ArgumentError("x_period must be positive"))
        filling_num > 0 || throw(ArgumentError("filling_num must be positive"))
        filling_den > 0 || throw(ArgumentError("filling_den must be positive"))
        isfinite(phi_y) || throw(ArgumentError("phi_y must be finite"))

        max_int = BigInt(typemax(Int))
        n_sites = BigInt(Ly) * x_period
        n_sites <= max_int || throw(ArgumentError("site count does not fit Int"))
        n_unit_cells = n_sites ÷ 2
        n_unit_cells <= max_int || throw(ArgumentError("unit-cell count does not fit Int"))
        particle_numerator = BigInt(filling_num) * n_unit_cells
        particle_numerator % filling_den == 0 ||
            throw(ArgumentError("unit cell has noninteger particle number"))
        n_particles = particle_numerator ÷ filling_den
        n_particles <= max_int || throw(ArgumentError("particle count does not fit Int"))
        n_particles <= n_sites || throw(ArgumentError("particle number exceeds site count"))

        return new(Ly, x_period, filling_num, filling_den, phi_y)
    end
end

function InfiniteCylinderConfig(;
    Ly::Int=6,
    x_period::Int=1,
    filling_num::Int=1,
    filling_den::Int=3,
    phi_y::Real=0.0,
)
    return InfiniteCylinderConfig(Ly, x_period, filling_num, filling_den, Float64(phi_y))
end

sites_per_cell(c::InfiniteCylinderConfig) = c.Ly * c.x_period
unit_cells_per_cell(c::InfiniteCylinderConfig) = sites_per_cell(c) ÷ 2
particles_per_cell(c::InfiniteCylinderConfig) =
    Int(BigInt(c.filling_num) * unit_cells_per_cell(c) ÷ c.filling_den)
physical_site_density(c::InfiniteCylinderConfig) = (c.filling_num // c.filling_den) / 2
charge_scale(c::InfiniteCylinderConfig) =
    sites_per_cell(c) ÷ gcd(sites_per_cell(c), particles_per_cell(c))

function linear_site(c::InfiniteCylinderConfig, x::Int, y::Int)
    x_offset = Base.Checked.checked_mul(x, c.Ly)
    coordinate = Base.Checked.checked_add(x_offset, mod(y, c.Ly))
    return Base.Checked.checked_add(coordinate, 1)
end

function site_coordinates(c::InfiniteCylinderConfig, n::Int)
    zero_based = Base.Checked.checked_sub(n, 1)
    x, y = fldmod(zero_based, c.Ly)
    return x, y
end

function canonical_site(c::InfiniteCylinderConfig, n::Int)
    zero_based = Base.Checked.checked_sub(n, 1)
    cell, offset = fldmod(zero_based, sites_per_cell(c))
    return offset + 1, cell
end

function default_occupied_sites(c::InfiniteCylinderConfig)
    np = particles_per_cell(c)
    n_sites = BigInt(sites_per_cell(c))
    return [
        Int(fld(BigInt(i - 1) * n_sites, BigInt(np))) + 1 for i in 1:np
    ]
end
