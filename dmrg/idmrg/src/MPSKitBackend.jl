const BLOCKTENSORKIT_COMMIT =
    "230cec77c9c7527817d2216b9c6de87f6d8bdda8"
const MPSKIT_BACKEND_COMMIT =
    "811ecf6c06c1f7c1bc656da61abcd679effcd428"
const TENSORKITTENSORS_COMMIT =
    "3755705a1c44a3d5e32086e7d89b2c561b268cb1"

const _MPSKIT_CHARGE_SECTOR = TensorKit.ProductSector{
    Tuple{TensorKit.FermionParity,TensorKit.U1Irrep},
}

function mpskit_backend_provenance()
    return (
        backend="mpskit_idmrg_v1",
        blocktensorkit_commit=BLOCKTENSORKIT_COMMIT,
        blocktensorkit_version=string(Base.pkgversion(BlockTensorKit)),
        mpskit_commit=MPSKIT_BACKEND_COMMIT,
        mpskit_version=string(Base.pkgversion(MPSKit)),
        tensorkit_version=string(Base.pkgversion(TensorKit)),
        tensorkittensors_commit=TENSORKITTENSORS_COMMIT,
        tensorkittensors_version=string(Base.pkgversion(TensorKitTensors)),
        jld2_version=string(Base.pkgversion(JLD2)),
    )
end

function mpskit_centered_charges(c::InfiniteCylinderConfig)
    nsites = sites_per_cell(c)
    nparticles = particles_per_cell(c)
    divisor = gcd(nsites, nparticles)
    return (
        scale=nsites ÷ divisor,
        empty=-(nparticles ÷ divisor),
        occupied=(nsites - nparticles) ÷ divisor,
    )
end

function _mpskit_charge_sector(parity::Bool, raw_charge::Integer)
    return _MPSKIT_CHARGE_SECTOR(
        TensorKit.FermionParity(parity),
        TensorKit.U1Irrep(raw_charge),
    )
end

function _mpskit_charge_space(parity::Bool, raw_charge::Integer)
    sector = _mpskit_charge_sector(parity, raw_charge)
    return TensorKit.Vect[_MPSKIT_CHARGE_SECTOR](sector => 1)
end

function _mpskit_physical_space(c::InfiniteCylinderConfig)
    charges = mpskit_centered_charges(c)
    empty = _mpskit_charge_sector(false, charges.empty)
    occupied = _mpskit_charge_sector(true, charges.occupied)
    return TensorKit.Vect[_MPSKIT_CHARGE_SECTOR](
        empty => 1,
        occupied => 1,
    )
end

function mpskit_physical_spaces(c::InfiniteCylinderConfig)
    physical_space = _mpskit_physical_space(c)
    return fill(physical_space, sites_per_cell(c))
end

function mpskit_physical_raw_charges(physical_space)
    raw_charges = Int[]
    for sector in TensorKit.sectors(physical_space)
        length(sector) == 2 || throw(
            ArgumentError("physical space does not use parity-and-charge sectors")
        )
        TensorKit.dim(physical_space, sector) == 1 || throw(
            ArgumentError("physical charge sectors must be one dimensional")
        )
        raw_charge = TensorKit.charge(sector[2])
        isinteger(raw_charge) || throw(
            ArgumentError("physical U(1) charge must be an integer")
        )
        push!(raw_charges, Int(raw_charge))
    end
    length(raw_charges) == 2 || throw(
        ArgumentError("physical space must contain exactly two charge sectors")
    )
    return sort!(raw_charges)
end

function _mpskit_occupied_sites(occupied_sites, c::InfiniteCylinderConfig)
    occupied = try
        Int[site for site in occupied_sites]
    catch error
        error isa InexactError || rethrow()
        throw(ArgumentError("occupied sites must be integer-valued"))
    end
    nsites = sites_per_cell(c)
    all(site -> 1 <= site <= nsites, occupied) || throw(
        ArgumentError("occupied sites must lie in 1:$nsites")
    )
    length(unique(occupied)) == length(occupied) || throw(
        ArgumentError("occupied sites must be distinct")
    )
    return sort!(occupied)
end

function mpskit_total_raw_charge(occupied_sites, c::InfiniteCylinderConfig)
    occupied = _mpskit_occupied_sites(occupied_sites, c)
    charges = mpskit_centered_charges(c)
    nempty = sites_per_cell(c) - length(occupied)
    return length(occupied) * charges.occupied + nempty * charges.empty
end

function mpskit_product_state(c::InfiniteCylinderConfig, occupied_sites)
    occupied = _mpskit_occupied_sites(occupied_sites, c)
    length(occupied) == particles_per_cell(c) || throw(
        ArgumentError(
            "product state requires $(particles_per_cell(c)) occupied sites",
        ),
    )
    mpskit_total_raw_charge(occupied, c) == 0 || throw(
        ArgumentError("product-state charge path must close in the zero sector")
    )

    occupied_set = Set(occupied)
    charges = mpskit_centered_charges(c)
    physical_spaces = mpskit_physical_spaces(c)
    bond_spaces = Vector{typeof(first(physical_spaces))}(
        undef,
        sites_per_cell(c),
    )
    cumulative_parity = false
    cumulative_charge = 0
    for site in eachindex(bond_spaces)
        isoccupied = site in occupied_set
        cumulative_parity = xor(cumulative_parity, isoccupied)
        cumulative_charge += isoccupied ? charges.occupied : charges.empty
        bond_spaces[site] = _mpskit_charge_space(
            cumulative_parity,
            cumulative_charge,
        )
    end
    cumulative_parity && throw(
        ArgumentError("product-state fermion-parity path must close in the even sector")
    )
    cumulative_charge == 0 || throw(
        ArgumentError("product-state charge path must close in the zero sector")
    )

    tensors = map(eachindex(bond_spaces)) do site
        left_space = bond_spaces[mod1(site - 1, length(bond_spaces))]
        codomain_space = TensorKit.:(⊗)(left_space, physical_spaces[site])
        return TensorKit.isometry(
            ComplexF64,
            codomain_space,
            bond_spaces[site],
        )
    end
    psi = MPSKit.InfiniteMPS(tensors)
    TensorKit.normalize!(psi)
    return psi
end

function _mpskit_only_raw_charge(space)
    TensorKit.dim(space) == 1 || throw(
        ArgumentError("product-state virtual spaces must have dimension one")
    )
    sectors = collect(TensorKit.sectors(space))
    length(sectors) == 1 || throw(
        ArgumentError("product-state virtual spaces must contain one sector")
    )
    sector = only(sectors)
    length(sector) == 2 || throw(
        ArgumentError("virtual space does not use parity-and-charge sectors")
    )
    raw_charge = TensorKit.charge(sector[2])
    isinteger(raw_charge) || throw(
        ArgumentError("virtual U(1) charge must be an integer")
    )
    return Int(raw_charge)
end

function mpskit_product_occupations(psi, c::InfiniteCylinderConfig)
    length(psi) == sites_per_cell(c) || throw(
        ArgumentError("state and configuration unit-cell lengths do not match")
    )
    expected_physical = _mpskit_physical_space(c)
    all(
        site -> MPSKit.physicalspace(psi, site) == expected_physical,
        1:length(psi),
    ) || throw(ArgumentError("state does not use the expected physical space"))

    charges = mpskit_centered_charges(c)
    previous_charge = _mpskit_only_raw_charge(
        MPSKit.right_virtualspace(psi, length(psi)),
    )
    occupied = Int[]
    for site in 1:length(psi)
        current_charge = _mpskit_only_raw_charge(
            MPSKit.right_virtualspace(psi, site),
        )
        local_charge = current_charge - previous_charge
        if local_charge == charges.occupied
            push!(occupied, site)
        elseif local_charge != charges.empty
            throw(
                ArgumentError(
                    "virtual charge path is not a bond-dimension-one product state",
                ),
            )
        end
        previous_charge = current_charge
    end
    return occupied
end
