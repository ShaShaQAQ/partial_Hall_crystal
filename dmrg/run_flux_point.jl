#!/usr/bin/env julia

include("DMRGFullModel.jl")
using .DMRGFullModel
using Dates
using LinearAlgebra
using Printf
using ITensors

function parse_cli(args)
    opts = Dict{String,String}()
    for arg in args
        startswith(arg, "--") || continue
        keyval = split(arg[3:end], "="; limit=2)
        length(keyval) == 2 || error("Arguments must use --key=value, got $arg")
        opts[keyval[1]] = keyval[2]
    end
    return opts
end

parse_int(opts, key, default) = parse(Int, get(opts, key, string(default)))
parse_float(opts, key, default) = parse(Float64, get(opts, key, string(default)))
parse_int_list(s) = [parse(Int, x) for x in split(s, ",") if !isempty(x)]
parse_bool(opts, key, default) = lowercase(get(opts, key, string(default))) in ("1", "true", "yes", "on")

function configure_threading(; threaded_blocksparse::Bool, blas_threads::Int, disable_strided_threads::Bool)
    blas_threads > 0 || throw(ArgumentError("blas_threads must be positive"))
    BLAS.set_num_threads(blas_threads)
    disable_strided_threads && ITensors.Strided.disable_threads()
    if threaded_blocksparse
        ITensors.enable_threaded_blocksparse()
    else
        ITensors.disable_threaded_blocksparse()
    end
    @printf(
        "Thread config: julia_threads=%d BLAS=%d threaded_blocksparse=%s strided_threads=%d\n",
        Threads.nthreads(),
        BLAS.get_num_threads(),
        string(ITensors.using_threaded_blocksparse()),
        ITensors.Strided.get_num_threads(),
    )
end

function first_sweep_at_maxdim(maxdim::Vector{Int}, target_maxdim::Int)
    target_maxdim > 0 || throw(ArgumentError("target_maxdim must be positive"))
    idx = findfirst(>=(target_maxdim), maxdim)
    return idx === nothing ? length(maxdim) : idx
end

function write_point_summary(
    path::AbstractString,
    result::DMRGRunResult,
    lat::CylinderLat,
    params::CylinderModelParams;
    Lx::Int,
    Ly::Int,
    Np::Int,
    step::Int,
    steps::Int,
    phi::Float64,
    nsweeps::Int,
    max_sweeps::Int,
    maxdim::Vector{Int},
    cutoff::Float64,
    min_sweeps::Int,
    target_sweep::Int,
    stable_sweeps::Int,
    energy_tol::Float64,
    density_tol::Float64,
    truncerr_tol::Float64,
    start_time::DateTime,
    end_time::DateTime,
)
    open(path, "w") do io
        @printf(io, "energy %.16g\n", result.energy)
        @printf(io, "phi %.16g\n", phi)
        @printf(io, "step %d\nsteps %d\n", step, steps)
        @printf(io, "Lx %d\nLy %d\nNs %d\nNuc %d\nNp %d\n", Lx, Ly, lat.Ns, lat.Nuc, Np)
        @printf(io, "t1 %.16g\nt3 %.16g\nV1 %.16g\nV2 %.16g\nV3 %.16g\n", params.t1, params.t3, params.V1, params.V2, params.V3)
        @printf(io, "nsweeps %d\nmax_sweeps %d\nmin_sweeps %d\n", nsweeps, max_sweeps, min_sweeps)
        @printf(io, "target_sweep %d\nstable_sweeps %d\n", target_sweep, stable_sweeps)
        @printf(io, "maxdim %s\ncutoff %.16g\n", join(maxdim, ","), cutoff)
        @printf(io, "energy_tol %.16g\ndensity_tol %.16g\ntruncerr_tol %.16g\n", energy_tol, density_tol, truncerr_tol)
        if !isempty(result.convergence)
            last = result.convergence[end]
            @printf(io, "actual_sweeps %d\n", last.sweep)
            @printf(io, "converged %s\nconvergence_reason %s\n", string(last.converged), last.reason)
            @printf(io, "final_delta_energy %.16g\nfinal_max_density_delta %.16g\nfinal_maxerr %.16g\n", last.delta_energy, last.max_density_delta, last.maxerr)
        else
            @printf(io, "actual_sweeps %d\nconverged false\nconvergence_reason fixed_sweeps\n", nsweeps)
        end
        @printf(io, "run_started %s\nrun_finished %s\n", Dates.format(start_time, dateformat"yyyy-mm-ddTHH:MM:SS"), Dates.format(end_time, dateformat"yyyy-mm-ddTHH:MM:SS"))
        @printf(io, "julia_version %s\njulia_num_threads %d\n", string(VERSION), Threads.nthreads())
        @printf(io, "blas_threads %d\nthreaded_blocksparse %s\nstrided_threads %d\n", BLAS.get_num_threads(), string(ITensors.using_threaded_blocksparse()), ITensors.Strided.get_num_threads())
        @printf(io, "pbs_jobid %s\nhostname %s\n", get(ENV, "PBS_JOBID", "none"), get(ENV, "HOSTNAME", "unknown"))
    end
end

function main(args)
    opts = parse_cli(args)
    Lx = parse_int(opts, "Lx", 15)
    Ly = parse_int(opts, "Ly", 6)
    steps = parse_int(opts, "steps", 12)
    step = parse_int(opts, "step", 0)
    0 <= step <= steps || throw(ArgumentError("step must be between 0 and steps"))
    phi = 2π * step / steps
    lat = cylinder_lattice(Lx; Ly)
    Np = parse_int(opts, "Np", filling_np(lat))
    params = CylinderModelParams(
        t1=parse_float(opts, "t1", 1.0),
        t3=parse_float(opts, "t3", 0.2),
        V1=parse_float(opts, "V1", 1.0),
        V2=parse_float(opts, "V2", 0.0),
        V3=parse_float(opts, "V3", 0.0),
    )
    nsweeps = parse_int(opts, "nsweeps", 30)
    max_sweeps = parse_int(opts, "max_sweeps", nsweeps)
    maxdim = parse_int_list(get(opts, "maxdim", "200,400,800,1200,1600,2000"))
    cutoff = parse_float(opts, "cutoff", 1e-9)
    min_sweeps = parse_int(opts, "min_sweeps", 6)
    target_maxdim = parse_int(opts, "target_maxdim", maximum(maxdim))
    target_sweep = parse_int(opts, "target_sweep", first_sweep_at_maxdim(maxdim, target_maxdim))
    stable_sweeps = parse_int(opts, "stable_sweeps", 3)
    energy_tol = parse_float(opts, "energy_tol", 1e-6)
    density_tol = parse_float(opts, "density_tol", 2e-5)
    truncerr_tol = parse_float(opts, "truncerr_tol", 5e-8)
    threaded_blocksparse = parse_bool(opts, "threaded_blocksparse", true)
    blas_threads = parse_int(opts, "blas_threads", 1)
    disable_strided_threads = parse_bool(opts, "disable_strided_threads", threaded_blocksparse)
    output = get(opts, "output", @sprintf("dmrg/flux_parallel/phi_%03d", step))
    mkpath(output)

    configure_threading(;
        threaded_blocksparse,
        blas_threads,
        disable_strided_threads,
    )

    start_time = now()
    @printf("Flux point step=%d/%d phi=%.16g\n", step, steps, phi)
    result = run_dmrg(
        lat,
        params,
        Np;
        phi_y=phi,
        nsweeps,
        maxdim,
        cutoff,
        adaptive=true,
        max_sweeps,
        min_sweeps,
        target_sweep,
        stable_sweeps,
        energy_tol,
        density_tol,
        truncerr_tol,
    )
    end_time = now()
    @printf("DMRG energy %.16g\n", result.energy)
    write_density(joinpath(output, "density.dat"), lat, result.density)
    write_vector(joinpath(output, "density_x.dat"), collect(0:(lat.Lx - 1)), y_averaged_density(lat, result.density); header="# x y_averaged_density")
    write_convergence(joinpath(output, "convergence.dat"), result.convergence)
    write_point_summary(
        joinpath(output, "summary.dat"),
        result,
        lat,
        params;
        Lx,
        Ly,
        Np,
        step,
        steps,
        phi,
        nsweeps,
        max_sweeps,
        maxdim,
        cutoff,
        min_sweeps,
        target_sweep,
        stable_sweeps,
        energy_tol,
        density_tol,
        truncerr_tol,
        start_time,
        end_time,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
