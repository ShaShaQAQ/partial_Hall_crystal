module CylinderModelCore

using LinearAlgebra

include("../shared/lattice.jl")
include("../shared/hoppings.jl")

export CylinderModelParams,
    NB_DISPS,
    reciprocal_vectors,
    hopping_matrices_from_get_Hk

struct CylinderModelParams
    t1::Float64
    t3::Float64
    V1::Float64
    V2::Float64
    V3::Float64
end

CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.0, V3=0.0) =
    CylinderModelParams(Float64(t1), Float64(t3), Float64(V1), Float64(V2), Float64(V3))

const NB_DISPS = Dict(
    1 => [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, 1), (1, -1)],
    2 => [(1, 1), (-1, -1), (-1, 2), (1, -2), (-2, 1), (2, -1)],
    3 => [(2, 0), (-2, 0), (0, 2), (0, -2), (-2, 2), (2, -2)],
)

function reciprocal_vectors(a1::NTuple{2,Float64}, a2::NTuple{2,Float64})
    area = a1[1] * a2[2] - a1[2] * a2[1]
    b1 = (2π * a2[2] / area, -2π * a2[1] / area)
    b2 = (-2π * a1[2] / area, 2π * a1[1] / area)
    return b1, b2
end

function hopping_matrices_from_get_Hk(
    t1::Float64,
    t3::Float64;
    nx_aux::Int=16,
    ny_aux::Int=16,
    max_dx::Int=4,
    max_dy::Int=4,
    tol::Float64=1e-10,
)
    a1 = [1.0, 0.0]
    a2 = [0.5, sqrt(3) / 2]
    b1, b2 = reciprocal_vectors((a1[1], a1[2]), (a2[1], a2[2]))
    b1v = collect(b1)
    b2v = collect(b2)
    sublat = [[0.0, 0.0], copy(a2)]
    out = Dict{NTuple{2,Int},Matrix{ComplexF64}}()
    Nk = nx_aux * ny_aux
    for n1 in -max_dx:max_dx, n2 in -max_dy:max_dy
        R = n1 .* a1 .+ n2 .* (2 .* a2)
        tmat = zeros(ComplexF64, 2, 2)
        for m1 in 0:(nx_aux - 1), m2 in 0:(ny_aux - 1)
            k = (m1 / nx_aux) .* b1v .+ (m2 / (2 * ny_aux)) .* b2v
            Hk = get_Hk(k, t1, t3)
            for alpha in 1:2, beta in 1:2
                delta = sublat[alpha] .- sublat[beta]
                tmat[alpha, beta] += exp(-1im * dot(k, R .+ delta)) * Hk[alpha, beta]
            end
        end
        tmat ./= Nk
        for alpha in 1:2, beta in 1:2
            abs(real(tmat[alpha, beta])) < tol && (tmat[alpha, beta] -= real(tmat[alpha, beta]))
            abs(imag(tmat[alpha, beta])) < tol && (tmat[alpha, beta] -= 1im * imag(tmat[alpha, beta]))
        end
        maximum(abs, tmat) > tol && (out[(n1, n2)] = tmat)
    end
    return out
end

end
