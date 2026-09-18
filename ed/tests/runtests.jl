using Test
using LinearAlgebra
using SparseArrays

include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")
include("../shared/optical_response.jl")
include("../shared/structure_factor.jl")

include("hopping_derivatives_test.jl")
include("optical_response_test.jl")
include("large_sector_optimizations_test.jl")
include("case4_optical_integration_test.jl")
