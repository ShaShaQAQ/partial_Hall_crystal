using Test
using LinearAlgebra
using SparseArrays

include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")
include("../shared/optical_response.jl")

include("hopping_derivatives_test.jl")
include("optical_response_test.jl")
