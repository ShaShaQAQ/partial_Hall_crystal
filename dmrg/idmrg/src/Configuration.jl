struct InfiniteCylinderConfig
    Ly::Int
    x_period::Int
    filling_num::Int
    filling_den::Int
    phi_y::Float64
end

function InfiniteCylinderConfig(;
    Ly::Int=6,
    x_period::Int=1,
    filling_num::Int=1,
    filling_den::Int=3,
    phi_y::Real=0.0,
)
    Ly > 0 && iseven(Ly) || throw(ArgumentError("Ly must be positive and even"))
    x_period > 0 || throw(ArgumentError("x_period must be positive"))
    filling_num > 0 || throw(ArgumentError("filling_num must be positive"))
    filling_den > 0 || throw(ArgumentError("filling_den must be positive"))
    c = InfiniteCylinderConfig(Ly, x_period, filling_num, filling_den, Float64(phi_y))
    numerator = filling_num * unit_cells_per_cell(c)
    numerator % filling_den == 0 ||
        throw(ArgumentError("unit cell has noninteger particle number"))
    return c
end

sites_per_cell(c::InfiniteCylinderConfig) = c.Ly * c.x_period
unit_cells_per_cell(c::InfiniteCylinderConfig) = sites_per_cell(c) ÷ 2
particles_per_cell(c::InfiniteCylinderConfig) =
    c.filling_num * unit_cells_per_cell(c) ÷ c.filling_den
physical_site_density(c::InfiniteCylinderConfig) = c.filling_num // (2 * c.filling_den)
charge_scale(c::InfiniteCylinderConfig) =
    sites_per_cell(c) ÷ gcd(sites_per_cell(c), particles_per_cell(c))

linear_site(c::InfiniteCylinderConfig, x::Int, y::Int) = x * c.Ly + mod(y, c.Ly) + 1

function site_coordinates(c::InfiniteCylinderConfig, n::Int)
    x, y = fldmod(n - 1, c.Ly)
    return x, y
end

function canonical_site(c::InfiniteCylinderConfig, n::Int)
    cell, offset = fldmod(n - 1, sites_per_cell(c))
    return offset + 1, cell
end

function default_occupied_sites(c::InfiniteCylinderConfig)
    np = particles_per_cell(c)
    return unique(round.(Int, range(1, sites_per_cell(c); length=np + 2))[2:(end - 1)])
end
