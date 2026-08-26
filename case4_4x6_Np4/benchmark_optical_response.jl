if !isdefined(Main, :GenLat)
    include("../shared/lattice.jl")
    include("../shared/hoppings.jl")
    include("../shared/basis.jl")
    include("../shared/ksector.jl")
    include("../shared/hamiltonian.jl")
end
if !isdefined(Main, :ResponseLanczosKernel)
    include("../shared/optical_response.jl")
end

using LinearAlgebra
using Printf


function hermiticity_error(operator)
    return norm(operator - operator') / max(norm(operator), 1.0)
end


function materialize_sector_operator(
        sec::KSector,
        lat::GenLat,
        hops::Vector{Tuple{Int,Int,ComplexF64}},
        V1::Float64,
        V2::Float64,
        V3::Float64)
    n = length(sec.reps)
    operator = zeros(ComplexF64, n, n)
    input = zeros(ComplexF64, n)
    diagonal = precompute_diag_H(sec, lat, V1, V2, V3)
    for column in 1:n
        input[column] = 1
        Hv!(view(operator, :, column), input, sec, lat, hops,
            V1, V2, V3, diagonal)
        input[column] = 0
    end
    return operator
end


function evaluate_curve(G, Kexp, area, omegas, eta)
    parts = [conductivity_parts(G, Kexp, area, omega, eta)
             for omega in omegas]
    return (total=ComplexF64[p.total for p in parts],
            regular=ComplexF64[p.regular for p in parts],
            drude=ComplexF64[p.drude for p in parts],
            drude_weight=parts[1].drude_weight)
end


function run_optical_benchmark(
        ; omegas=collect(0.0:0.001:10.0),
        eta::Float64=0.065,
        requested_m::Vector{Int}=[50, 100, 200, 400, 800],
        write_outputs::Bool=true,
        make_plots::Bool=true,
        output_dir::String=joinpath(@__DIR__, "optical_response_output"))
    all(m -> m > 0, requested_m) ||
        throw(ArgumentError("requested_m values must be positive"))

    lat = RectLat4x6()
    basis = gen_basis(lat.Ns, 4)
    sec = build_ksector(basis, lat, 2)
    hops, jx_hops, kxx_hops = build_hops_x_derivatives(lat, 1.0, 0.2)
    H = materialize_sector_operator(sec, lat, hops, 1.0, 0.0, 0.0)
    Jx = materialize_sector_operator(sec, lat, jx_hops, 0.0, 0.0, 0.0)
    Kxx = materialize_sector_operator(sec, lat, kxx_hops, 0.0, 0.0, 0.0)

    herm = (H=hermiticity_error(H),
            Jx=hermiticity_error(Jx),
            Kxx=hermiticity_error(Kxx))
    maximum(values(herm)) < 1e-10 ||
        error("operator Hermiticity check failed: $herm")

    eig = eigen(Hermitian(H))
    Eg = eig.values[1]
    g = eig.vectors[:, 1]
    residual = norm(H * g - Eg * g)
    f = projected_response_source(Jx, g)
    overlap = dot(g, f)
    Kexp_complex = dot(g, Kxx * g)
    abs(imag(Kexp_complex)) < 1e-10 ||
        error("Kxx expectation is not real")
    Kexp = real(Kexp_complex)
    area = system_area(lat)

    exact_data = exact_spectral_data(eig.values, eig.vectors, Eg, f)
    exact_G(z) = exact_resolvent(exact_data, z)
    exact_curve = evaluate_curve(exact_G, Kexp, area, omegas, eta)

    requested = sort(unique(min.(requested_m, length(f))))
    full_kernel = response_lanczos(H, Eg, f; mmax=maximum(requested))
    m_values = sort(unique(min.(requested, length(full_kernel.alpha))))
    kernels = [truncate_kernel(full_kernel, m) for m in m_values]
    curves = map(kernels) do kernel
        G(z) = lanczos_resolvent(kernel, z)
        evaluate_curve(G, Kexp, area, omegas, eta)
    end
    errors = [scaled_max_error(curve.regular, exact_curve.regular)
              for curve in curves]

    result = (sector=2, Eg=Eg, ground_residual=residual,
              hermiticity=herm, source_overlap=overlap,
              source_norm2=real(dot(f, f)), Kexp=Kexp, area=area,
              exact_data=exact_data, exact_curve=exact_curve,
              m_values=m_values, kernels=kernels, curves=curves,
              errors=errors, omegas=Float64.(omegas), eta=eta)

    write_outputs &&
        write_benchmark_outputs(result, output_dir; make_plots=make_plots)
    return result
end


function write_benchmark_outputs(result, output_dir; make_plots=true)
    error("output writer is implemented in Task 5")
end
