using TOML
using SHA
using InfiniteCylinderDMRG

import KrylovKit
import MPSKit
import TensorKit

function mpskit_checkpoint_test_fixture()
    config = InfiniteCylinderConfig(;
        Ly=4,
        x_period=1,
        filling_num=1,
        filling_den=1,
    )
    physical_spaces = mpskit_physical_spaces(config)
    operators = InfiniteCylinderDMRG._mpskit_fermion_operators(
        first(physical_spaces);
        centered=true,
    )
    amplitude = -1.0 + 0.0im
    hopping = amplitude * operators.plus_min -
        conj(amplitude) * operators.min_plus
    hamiltonian = MPSKit.InfiniteMPOHamiltonian(
        physical_spaces,
        (1, 2) => hopping,
        (2, 3) => hopping,
        (3, 4) => hopping,
    )
    state = mpskit_product_state(config, [1, 3])
    return config, hamiltonian, state, operators.number
end

function mpskit_checkpoint_observations(state, hamiltonian, number_operator)
    environments = MPSKit.environments(state, hamiltonian, state)
    energy_per_site = real(
        MPSKit.expectation_value(state, hamiltonian, environments) /
        length(state),
    )
    galerkin_residual = MPSKit.calc_galerkin(
        state,
        hamiltonian,
        state,
        environments,
    )
    densities = [
        real(MPSKit.expectation_value(state, site => number_operator)) for
        site in 1:length(state)
    ]
    spectrum = MPSKit.entanglement_spectrum(state, 2)
    sector_weights = [
        Dict(
            "sector" => string(sector),
            "weight" => Float64(
                TensorKit.dim(sector) * sum(abs2, singular_values),
            ),
        ) for (sector, singular_values) in pairs(spectrum)
    ]
    sort!(sector_weights; by=row -> row["sector"])
    return (;
        energy_per_site,
        galerkin_residual,
        densities,
        sector_weights,
    )
end

checkpoint_sha256(path::AbstractString) = bytes2hex(SHA.sha256(read(path)))

function write_mpskit_checkpoint_toml(path::AbstractString, data)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end

function static_checkpoint_lanczos()
    return KrylovKit.Lanczos(;
        tol=1e-10,
        maxiter=200,
        eager=true,
        krylovdim=30,
        verbosity=0,
    )
end
