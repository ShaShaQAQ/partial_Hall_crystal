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

const OPTICAL_OPERATOR_CONVENTION =
    "H_A(k)=H(k+A_x*xhat); Jx=dH/dA_x; Kxx=d2H/dA_x2"


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

    result = (sector=2,
              operator_convention=OPTICAL_OPERATOR_CONVENTION,
              Eg=Eg, ground_residual=residual,
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
    @eval using JLD2
    mkpath(output_dir)
    final_index = length(result.m_values)
    final_curve = result.curves[final_index]
    final_kernel = result.kernels[final_index]

    Base.invokelatest(
        JLD2.jldsave,
        joinpath(output_dir, "optical_response_benchmark.jld2");
        sector=result.sector,
        operator_convention=result.operator_convention,
        Eg=result.Eg,
        ground_residual=result.ground_residual,
        hermiticity=result.hermiticity,
        source_overlap=result.source_overlap,
        source_norm2=result.source_norm2,
        Kexp=result.Kexp,
        area=result.area,
        eta=result.eta,
        omegas=result.omegas,
        m_values=result.m_values,
        errors=result.errors,
        alpha=final_kernel.alpha,
        beta=final_kernel.beta,
        exact_deltas=result.exact_data.deltas,
        exact_weights=result.exact_data.weights,
        exact_total=result.exact_curve.total,
        exact_regular=result.exact_curve.regular,
        lanczos_total=final_curve.total,
        lanczos_regular=final_curve.regular,
        exact_drude_weight=result.exact_curve.drude_weight,
        lanczos_drude_weight=final_curve.drude_weight,
    )

    open(joinpath(output_dir, "optical_response_curves.dat"), "w") do io
        println(io, "# operator_convention: $(result.operator_convention)")
        println(io, "# omega Re_sigma_exact Im_sigma_exact " *
                    "Re_sigma_L Im_sigma_L Re_reg_exact Im_reg_exact " *
                    "Re_reg_L Im_reg_L")
        for i in eachindex(result.omegas)
            values = (result.omegas[i],
                      real(result.exact_curve.total[i]),
                      imag(result.exact_curve.total[i]),
                      real(final_curve.total[i]),
                      imag(final_curve.total[i]),
                      real(result.exact_curve.regular[i]),
                      imag(result.exact_curve.regular[i]),
                      real(final_curve.regular[i]),
                      imag(final_curve.regular[i]))
            println(io, join(values, ' '))
        end
    end

    open(joinpath(output_dir, "lanczos_convergence.dat"), "w") do io
        println(io, "# operator_convention: $(result.operator_convention)")
        println(io, "# M scaled_max_error_regular")
        for (m, err) in zip(result.m_values, result.errors)
            println(io, "$m $err")
        end
    end

    if make_plots
        @eval using Plots
        p1 = Base.invokelatest(
            Plots.plot,
            result.omegas,
            real.(result.exact_curve.regular);
            label="exact",
            xlabel="omega",
            ylabel="Re sigma_xx^reg",
        )
        for (m, curve) in zip(result.m_values, result.curves)
            Base.invokelatest(
                Plots.plot!,
                p1,
                result.omegas,
                real.(curve.regular);
                label="M=$m",
            )
        end
        Base.invokelatest(
            Plots.savefig,
            p1,
            joinpath(output_dir, "optical_conductivity_regular.pdf"),
        )

        p2 = Base.invokelatest(
            Plots.plot,
            result.omegas,
            imag.(result.exact_curve.total);
            label="exact",
            xlabel="omega",
            ylabel="Im sigma_xx total",
        )
        Base.invokelatest(
            Plots.plot!,
            p2,
            result.omegas,
            imag.(final_curve.total);
            label="Lanczos M=$(result.m_values[end])",
        )
        Base.invokelatest(
            Plots.savefig,
            p2,
            joinpath(output_dir, "optical_conductivity_total.pdf"),
        )

        p3 = Base.invokelatest(
            Plots.plot,
            result.m_values,
            result.errors;
            marker=:circle,
            xscale=:log10,
            yscale=:log10,
            legend=false,
            xlabel="Lanczos M",
            ylabel="scaled max error",
        )
        Base.invokelatest(
            Plots.savefig,
            p3,
            joinpath(output_dir, "lanczos_convergence.pdf"),
        )
    end

    @printf("sector=%d  E0=%.12f  residual=%.3e\n",
            result.sector, result.Eg, result.ground_residual)
    @printf("<Kxx>=%.12g  ||QJxg||^2=%.12g  area=%.12g\n",
            result.Kexp, result.source_norm2, result.area)
    @printf("Dxx exact=%.12g  Lanczos=%.12g\n",
            result.exact_curve.drude_weight, final_curve.drude_weight)
    for (m, err) in zip(result.m_values, result.errors)
        @printf("M=%4d  scaled max error=%.3e\n", m, err)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    quick = "--quick" in ARGS
    no_plot = "--no-plot" in ARGS
    result = run_optical_benchmark(
        omegas=quick ? collect(0.0:0.05:10.0) :
                       collect(0.0:0.001:10.0),
        requested_m=quick ? [50, 100, 200, 400, 600] :
                            [50, 100, 200, 400, 800, 10_000],
        make_plots=!no_plot,
    )
    target_error = quick ? 1e-6 : 1e-8
    result.errors[end] < target_error ||
        error("final Lanczos curve did not reach target accuracy " *
              "$target_error")
end
