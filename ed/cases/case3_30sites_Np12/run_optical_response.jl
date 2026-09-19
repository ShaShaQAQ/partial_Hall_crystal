#!/usr/bin/env julia

include("../../shared/lattice.jl")
include("../../shared/hoppings.jl")
include("../../shared/basis.jl")
include("../../shared/ksector.jl")
include("../../shared/hamiltonian.jl")
include("../../shared/optical_response.jl")

using Dates
using JLD2
using LinearAlgebra
using Printf
using SparseArrays

# Sparse row multiplication uses Julia threads, so keep BLAS from competing
# for the same cores.
BLAS.set_num_threads(parse(Int, get(ENV, "PHC_BLAS_THREADS", "1")))

function parse_response_argument(name, default)
    index = findfirst(==(name), ARGS)
    index === nothing && return default
    index < length(ARGS) || error("missing value after $name")
    return ARGS[index + 1]
end

const RESPONSE_SECTOR = parse(
    Int, parse_response_argument("--sector", "5"))
const RESPONSE_MMAX = parse(
    Int, parse_response_argument("--mmax", "600"))
const RESPONSE_ETA = parse(
    Float64, parse_response_argument("--eta", "0.065"))
const RESPONSE_OMEGA_MAX = parse(
    Float64, parse_response_argument("--omega-max", "10.0"))
const RESPONSE_OMEGA_STEP = parse(
    Float64, parse_response_argument("--omega-step", "0.002"))
const RESPONSE_T1 = parse(
    Float64, parse_response_argument("--t1", "1.0"))
const RESPONSE_T3 = parse(
    Float64, parse_response_argument("--t3", "0.2"))
const RESPONSE_V1 = parse(
    Float64, parse_response_argument("--V1", "100.0"))
const RESPONSE_V2 = parse(
    Float64, parse_response_argument("--V2", "0.0"))
const RESPONSE_V3 = parse(
    Float64, parse_response_argument("--V3", "0.0"))
const RESPONSE_LATTICE = parse_response_argument(
    "--lattice", "legacy")
const RESPONSE_DATA_DIR = parse_response_argument(
    "--data-dir", joinpath(@__DIR__, "output"))
const RESPONSE_OUT_DIR = parse_response_argument(
    "--out-dir", joinpath(@__DIR__, "optical_response_output"))
const RESPONSE_PREFLIGHT =
    "--preflight" in ARGS || get(ENV, "PHC_PREFLIGHT", "0") == "1"

0 <= RESPONSE_SECTOR <= 14 || error("sector must be in 0:14")
RESPONSE_MMAX > 0 || error("mmax must be positive")
RESPONSE_OMEGA_STEP > 0 || error("omega step must be positive")
RESPONSE_LATTICE in ("legacy", "corrected") ||
    error("lattice must be legacy or corrected")

function response_lattice(convention::AbstractString=RESPONSE_LATTICE)
    convention == "legacy" && return LegacyTiltedLat30()
    convention == "corrected" && return TiltedLat30()
    throw(ArgumentError("unknown lattice convention: $convention"))
end

function response_partial_path(
        sector::Int, convention::AbstractString=RESPONSE_LATTICE)
    first_sector = convention == "corrected" ? sector : 5 * div(sector, 5)
    return joinpath(RESPONSE_DATA_DIR, "partial_$(first_sector).jld2")
end

function load_response_ground_state(sector::Int)
    path = response_partial_path(sector)
    isfile(path) || error("missing saved data: $path")
    data = load(path)
    Int(data["Np"]) == 12 || error("this driver expects Np=12")
    saved_parameters = (
        t1=Float64(data["t1"]),
        t3=Float64(data["t3"]),
        V1=Float64(data["V1"]),
        V2=Float64(data["V2"]),
        V3=Float64(data["V3"]),
    )
    requested_parameters = (
        t1=RESPONSE_T1,
        t3=RESPONSE_T3,
        V1=RESPONSE_V1,
        V2=RESPONSE_V2,
        V3=RESPONSE_V3,
    )
    saved_parameters == requested_parameters ||
        error("CLI parameters do not match saved data: $saved_parameters")
    saved_lattice = haskey(data, "lattice_convention") ?
        String(data["lattice_convention"]) : "legacy"
    saved_lattice == RESPONSE_LATTICE || error(
        "saved lattice $saved_lattice does not match requested " *
        RESPONSE_LATTICE)

    ground_states = data["gs_vecs"]
    haskey(ground_states, sector) ||
        error("sector $sector not present in $path")
    sector_energies = [Float64(energy)
                       for (momentum, energy) in data["ev_pairs"]
                       if Int(momentum) == sector]
    isempty(sector_energies) &&
        error("no energy for sector $sector in $path")
    return (state=ComplexF64.(ground_states[sector]),
            energy=minimum(sector_energies), path=path,
            lattice_convention=saved_lattice)
end

function apply_response_operator(
        vector,
        sector,
        lattice,
        hops,
        diagonal;
        skip_tol=1e-20,
        scratch=nothing)
    output = zeros(ComplexF64, length(vector))
    if scratch === nothing
        Hv!(output, vector, sector, lattice, hops, 0.0, 0.0, 0.0,
            diagonal; skip_tol=skip_tol)
    else
        Hv_threaded!(output, vector, sector, lattice, hops, diagonal,
                     scratch; skip_tol=skip_tol)
    end
    return output
end

function save_response_curve(
        path,
        sector,
        ground_energy,
        residual,
        diamagnetic_expectation,
        source_norm2,
        area,
        eta,
        frequencies,
        curve,
        kernel,
        lattice_convention;
        model_parameters,
        requested_mmax)
    operator_convention =
        "$lattice_convention TiltedLat30 Fourier mesh; " *
        "H_A(k)=H(k+A_x*xhat); " *
        "Jx=dH/dA_x; Kxx=d2H/dA_x2"
    Np = model_parameters.Np
    t1 = model_parameters.t1
    t3 = model_parameters.t3
    V1 = model_parameters.V1
    V2 = model_parameters.V2
    V3 = model_parameters.V3
    lanczos_steps = length(kernel.alpha)
    jldsave(
        path;
        sector,
        Np,
        t1,
        t3,
        V1,
        V2,
        V3,
        requested_mmax,
        lanczos_steps,
        ground_energy,
        residual,
        diamagnetic_expectation,
        source_norm2,
        area,
        eta,
        frequencies,
        alpha=kernel.alpha,
        beta=kernel.beta,
        drude_weight=curve.drude_weight,
        total=curve.total,
        regular=curve.regular,
        drude=curve.drude,
        lattice_convention,
        operator_convention,
    )

    data_path = replace(path, ".jld2" => ".dat")
    open(data_path, "w") do output
        println(output,
                "# sector=$sector Eg=$ground_energy eta=$eta area=$area " *
                "Kexp=$diamagnetic_expectation source_norm2=$source_norm2 " *
                "lattice=$lattice_convention Np=$Np t1=$t1 t3=$t3 " *
                "V1=$V1 V2=$V2 V3=$V3 requested_mmax=$requested_mmax " *
                "lanczos_steps=$lanczos_steps residual=$residual")
        println(output,
                "# omega Re_total Im_total Re_regular Im_regular " *
                "Re_drude Im_drude")
        for i in eachindex(frequencies)
            @printf(output,
                    "%.10g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                    frequencies[i], real(curve.total[i]),
                    imag(curve.total[i]), real(curve.regular[i]),
                    imag(curve.regular[i]), real(curve.drude[i]),
                    imag(curve.drude[i]))
        end
    end
end

function run_large_optical_response()
    println("="^72)
    println("30-site ED optical response, sector $RESPONSE_SECTOR")
    println("start: ", now())
    println("parameters: t1=$RESPONSE_T1 t3=$RESPONSE_T3 " *
            "V1=$RESPONSE_V1 V2=$RESPONSE_V2 V3=$RESPONSE_V3")
    println("lattice: $RESPONSE_LATTICE")
    println("response: M=$RESPONSE_MMAX eta=$RESPONSE_ETA " *
            "omega=0:$RESPONSE_OMEGA_STEP:$RESPONSE_OMEGA_MAX")
    flush(stdout)

    saved = load_response_ground_state(RESPONSE_SECTOR)
    ground_state = saved.state
    ground_energy = saved.energy
    lattice = response_lattice()
    println("[1] lattice: Ns=$(lattice.Ns) Nuc=$(lattice.Nuc) " *
            "source=$(saved.path)")
    println("[2] generating C($(lattice.Ns),12) basis")
    basis = gen_basis(lattice.Ns, 12)
    println("[3] building momentum sector $RESPONSE_SECTOR")
    translation_maps = translation_site_maps(lattice)
    sector_time = @elapsed sector = build_ksector(
        basis, lattice, RESPONSE_SECTOR; trans_maps=translation_maps)
    length(ground_state) == length(sector.reps) ||
        error("saved vector dimension does not match reconstructed sector")
    @printf("    dimension=%d sector_time=%.2f s\n",
            length(ground_state), sector_time)

    hops = build_hops(lattice, RESPONSE_T1, RESPONSE_T3, 0.0)
    scratch = threaded_Hv_workspace(length(ground_state))
    diagonal_time = @elapsed diagonal = precompute_diag_H(
        sector, lattice, RESPONSE_V1, RESPONSE_V2, RESPONSE_V3)
    h_ground = similar(ground_state)
    h_time = @elapsed Hv_threaded!(
        h_ground, ground_state, sector, lattice, hops, diagonal, scratch;
        skip_tol=0.0)
    residual = norm(h_ground .- ground_energy .* ground_state)
    @printf("[4] saved-state residual=%.6e diag=%.2f s Hv=%.2f s\n",
            residual, diagonal_time, h_time)
    residual < 1e-7 ||
        error("saved ground state is incompatible with reconstructed H")

    println("[5] applying analytic Jx and Kxx")
    _, current_hops, diamagnetic_hops = build_hops_x_derivatives(
        lattice, RESPONSE_T1, RESPONSE_T3)
    zero_diagonal = zeros(Float64, length(ground_state))
    current_time = @elapsed source = apply_response_operator(
        ground_state, sector, lattice, current_hops, zero_diagonal;
        skip_tol=0.0, scratch=scratch)
    source .-= ground_state .* dot(ground_state, source)
    source_norm2 = real(dot(source, source))
    source_norm2 > 1e-14 || error("projected current source is zero")
    diamagnetic_time = @elapsed diamagnetic_vector = apply_response_operator(
        ground_state, sector, lattice, diamagnetic_hops, zero_diagonal;
        skip_tol=0.0, scratch=scratch)
    diamagnetic_expectation = real(dot(
        ground_state, diamagnetic_vector))
    @printf("    ||QJxg||^2=%.12g <Kxx>=%.12g Jx=%.2f s Kxx=%.2f s\n",
            source_norm2, diamagnetic_expectation,
            current_time, diamagnetic_time)

    if RESPONSE_PREFLIGHT
        println("preflight complete; sparse H and response Lanczos skipped")
        return nothing
    end

    scratch = nothing
    h_ground = nothing
    diamagnetic_vector = nothing
    zero_diagonal = nothing
    diagonal = nothing
    empty!(basis)
    GC.gc()

    println("[6] building sparse Hamiltonian")
    sparse_time = @elapsed sparse_hamiltonian = build_sparse_H(
        sector, lattice, hops, RESPONSE_V1, RESPONSE_V2, RESPONSE_V3)
    empty!(sector.orbit_data)
    empty!(sector.fock2rep)
    GC.gc()
    csr_time = @elapsed response_hamiltonian = threaded_csr(
        sparse_hamiltonian)
    sparse_hamiltonian = nothing
    GC.gc()
    @printf("    size=%s nnz=%d sparse=%.2f s CSR=%.2f s\n",
            string(size(response_hamiltonian)), nnz(response_hamiltonian),
            sparse_time, csr_time)

    println("[7] response Lanczos")
    kernel = response_lanczos(
        response_hamiltonian, ground_energy, source;
        mmax=min(RESPONSE_MMAX, length(source)),
        ground_state=ground_state,
        progress=(step, args...) ->
            (step % 50 == 0 && println("    response step=$step")))
    println("    steps=$(length(kernel.alpha)) " *
            "breakdown=$(kernel.breakdown)")

    area = system_area(lattice)
    frequencies = collect(
        0.0:RESPONSE_OMEGA_STEP:RESPONSE_OMEGA_MAX)
    resolvent(z) = lanczos_resolvent(kernel, z)
    components = [conductivity_parts(
                      resolvent, diamagnetic_expectation, area,
                      frequency, RESPONSE_ETA)
                  for frequency in frequencies]
    curve = (
        total=ComplexF64[part.total for part in components],
        regular=ComplexF64[part.regular for part in components],
        drude=ComplexF64[part.drude for part in components],
        drude_weight=components[1].drude_weight,
    )

    mkpath(RESPONSE_OUT_DIR)
    output_path = joinpath(
        RESPONSE_OUT_DIR,
        "sector_$(RESPONSE_SECTOR)_optical_response.jld2")
    save_response_curve(
        output_path, RESPONSE_SECTOR, ground_energy, residual,
        diamagnetic_expectation, source_norm2, area, RESPONSE_ETA,
        frequencies, curve, kernel, RESPONSE_LATTICE;
        model_parameters=(
            Np=12,
            t1=RESPONSE_T1,
            t3=RESPONSE_T3,
            V1=RESPONSE_V1,
            V2=RESPONSE_V2,
            V3=RESPONSE_V3,
        ),
        requested_mmax=RESPONSE_MMAX)
    println("[8] saved: $output_path")
    println("done: ", now())
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_large_optical_response()
end
