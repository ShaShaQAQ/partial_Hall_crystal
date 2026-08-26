#!/usr/bin/env julia

include("DMRGFullModel.jl")
using .DMRGFullModel
using Dates
using LinearAlgebra
using Printf
using Random
using Serialization
using ITensors
using ITensorMPS

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

function shell_output(cmd::Cmd)
    try
        return readchomp(cmd)
    catch err
        return "unavailable: $(typeof(err))"
    end
end

function write_summary(
    path::AbstractString,
    result::DMRGRunResult,
    lat::CylinderLat,
    params::CylinderModelParams;
    Lx::Int,
    Ly::Int,
    Np::Int,
    phi::Float64,
    nsweeps::Int,
    max_sweeps::Int,
    maxdim::Vector{Int},
    cutoff::Float64,
    adaptive::Bool,
    min_sweeps::Int,
    target_sweep::Int,
    stable_sweeps::Int,
    energy_tol::Float64,
    density_tol::Float64,
    truncerr_tol::Float64,
    seed::Int,
    init_linkdim::Int,
    init_source::AbstractString,
    checkpoint_in::AbstractString,
    checkpoint_out::AbstractString,
    checkpoint_saved::Bool,
    start_time::DateTime,
    end_time::DateTime,
)
    git_status = replace(shell_output(`git status --short`), '\n' => " | ")
    open(path, "w") do io
        @printf(io, "energy %.16g\n", result.energy)
        @printf(io, "Lx %d\nLy %d\nNs %d\nNuc %d\nNp %d\n", Lx, Ly, lat.Ns, lat.Nuc, Np)
        @printf(io, "t1 %.16g\nt3 %.16g\nV1 %.16g\nV2 %.16g\nV3 %.16g\nphi %.16g\n", params.t1, params.t3, params.V1, params.V2, params.V3, phi)
        @printf(io, "adaptive %s\n", string(adaptive))
        @printf(io, "nsweeps %d\nmax_sweeps %d\nmin_sweeps %d\n", nsweeps, max_sweeps, min_sweeps)
        @printf(io, "target_sweep %d\nstable_sweeps %d\n", target_sweep, stable_sweeps)
        @printf(io, "maxdim %s\ncutoff %.16g\n", join(maxdim, ","), cutoff)
        @printf(io, "energy_tol %.16g\ndensity_tol %.16g\ntruncerr_tol %.16g\n", energy_tol, density_tol, truncerr_tol)
        @printf(io, "seed %d\ninit_linkdim %d\ninit_source %s\n", seed, init_linkdim, init_source)
        @printf(io, "checkpoint_in %s\ncheckpoint_out %s\ncheckpoint_saved %s\n", checkpoint_in, checkpoint_out, string(checkpoint_saved))
        if !isempty(result.convergence)
            last = result.convergence[end]
            @printf(io, "actual_sweeps %d\n", last.sweep)
            @printf(io, "converged %s\nconvergence_reason %s\n", string(last.converged), last.reason)
            @printf(io, "final_delta_energy %.16g\nfinal_max_density_delta %.16g\nfinal_maxerr %.16g\n", last.delta_energy, last.max_density_delta, last.maxerr)
        else
            @printf(io, "actual_sweeps %d\n", nsweeps)
            @printf(io, "converged false\nconvergence_reason fixed_sweeps\n")
        end
        @printf(io, "run_started %s\nrun_finished %s\n", Dates.format(start_time, dateformat"yyyy-mm-ddTHH:MM:SS"), Dates.format(end_time, dateformat"yyyy-mm-ddTHH:MM:SS"))
        @printf(io, "julia_version %s\njulia_num_threads %d\n", string(VERSION), Threads.nthreads())
        @printf(io, "blas_threads %d\n", BLAS.get_num_threads())
        @printf(io, "threaded_blocksparse %s\n", string(ITensors.using_threaded_blocksparse()))
        @printf(io, "strided_threads %d\n", ITensors.Strided.get_num_threads())
        @printf(io, "pbs_jobid %s\n", get(ENV, "PBS_JOBID", "none"))
        @printf(io, "hostname %s\n", get(ENV, "HOSTNAME", shell_output(`hostname`)))
        @printf(io, "git_commit %s\n", shell_output(`git rev-parse HEAD`))
        @printf(io, "git_status_short %s\n", isempty(git_status) ? "clean" : git_status)
    end
end

function checkpoint_property(obj, name::Symbol, default=nothing)
    return hasproperty(obj, name) ? getproperty(obj, name) : default
end

function model_params_tuple(params::CylinderModelParams)
    return (t1=params.t1, t3=params.t3, V1=params.V1, V2=params.V2, V3=params.V3)
end

function validate_checkpoint(obj, path::AbstractString, lat::CylinderLat, params::CylinderModelParams; Lx::Int, Ly::Int, Np::Int, phi::Float64)
    checkpoint_property(obj, :format) == "ground_dmrg_checkpoint_v1" || error("unsupported checkpoint format in $path")
    checkpoint_property(obj, :Lx) == Lx || error("checkpoint Lx mismatch in $path")
    checkpoint_property(obj, :Ly) == Ly || error("checkpoint Ly mismatch in $path")
    checkpoint_property(obj, :Ns) == lat.Ns || error("checkpoint Ns mismatch in $path")
    checkpoint_property(obj, :Np) == Np || error("checkpoint Np mismatch in $path")
    abs(checkpoint_property(obj, :phi, NaN) - phi) <= 1e-12 || error("checkpoint phi mismatch in $path")
    checkpoint_property(obj, :params) == model_params_tuple(params) || error("checkpoint model parameters mismatch in $path")
    return obj
end

function load_ground_checkpoint(path::AbstractString, lat::CylinderLat, params::CylinderModelParams; Lx::Int, Ly::Int, Np::Int, phi::Float64)
    isfile(path) || error("missing checkpoint: $path")
    obj = deserialize(path)
    validate_checkpoint(obj, path, lat, params; Lx, Ly, Np, phi)
    return obj
end

function save_ground_checkpoint(path::AbstractString, result::DMRGRunResult, lat::CylinderLat, params::CylinderModelParams; Lx::Int, Ly::Int, Np::Int, phi::Float64)
    isempty(path) && return false
    mkpath(dirname(path))
    tmp = path * ".tmp"
    data = (
        format = "ground_dmrg_checkpoint_v1",
        saved_at = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
        Lx = Lx,
        Ly = Ly,
        Ns = lat.Ns,
        Nuc = lat.Nuc,
        Np = Np,
        phi = phi,
        energy = result.energy,
        density = copy(result.density),
        sites = result.sites,
        psi = result.psi,
        params = model_params_tuple(params),
    )
    serialize(tmp, data)
    mv(tmp, path; force=true)
    return true
end

function read_density_seed_state(path::AbstractString, lat::CylinderLat, Np::Int)
    isfile(path) || error("missing density seed file: $path")
    0 <= Np <= lat.Ns || throw(ArgumentError("Np must be between 0 and Ns"))
    density_by_site = fill(NaN, lat.Ns)
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        fields = split(stripped)
        length(fields) >= 6 || error("malformed density row in $path: $line")
        site = parse(Int, fields[1])
        1 <= site <= lat.Ns || error("density site index out of range in $path: $site")
        density_by_site[site] = parse(Float64, fields[6])
    end
    any(isnan, density_by_site) && error("density seed file $path does not contain all $(lat.Ns) sites")
    ranked = sort(collect(1:lat.Ns); by=i -> (-density_by_site[i], i))
    occupied = Set(ranked[1:Np])
    return [i in occupied ? "Occ" : "Emp" for i in 1:lat.Ns]
end

function density_seed_mps(path::AbstractString, lat::CylinderLat, Np::Int; seed::Int, init_linkdim::Int)
    state = read_density_seed_state(path, lat, Np)
    sites = siteinds("Fermion", lat.Ns; conserve_qns=true)
    rng = MersenneTwister(seed)
    psi0 = random_mps(rng, ComplexF64, sites, state; linkdims=init_linkdim)
    return sites, psi0
end

function first_sweep_at_maxdim(maxdim::Vector{Int}, target_maxdim::Int)
    target_maxdim > 0 || throw(ArgumentError("target_maxdim must be positive"))
    idx = findfirst(>=(target_maxdim), maxdim)
    return idx === nothing ? length(maxdim) : idx
end

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

function main(args)
    opts = parse_cli(args)
    Lx = parse_int(opts, "Lx", 4)
    Ly = parse_int(opts, "Ly", 6)
    lat = cylinder_lattice(Lx; Ly)
    Np = parse_int(opts, "Np", filling_np(lat))
    params = CylinderModelParams(
        t1=parse_float(opts, "t1", 1.0),
        t3=parse_float(opts, "t3", 0.2),
        V1=parse_float(opts, "V1", 1.0),
        V2=parse_float(opts, "V2", 0.0),
        V3=parse_float(opts, "V3", 0.0),
    )
    phi = parse_float(opts, "phi", 0.0)
    nsweeps = parse_int(opts, "nsweeps", 8)
    max_sweeps = parse_int(opts, "max_sweeps", nsweeps)
    maxdim = parse_int_list(get(opts, "maxdim", "20,50,100,200"))
    cutoff = parse_float(opts, "cutoff", 1e-8)
    adaptive = parse_bool(opts, "adaptive", false)
    min_sweeps = parse_int(opts, "min_sweeps", min(4, max_sweeps))
    target_maxdim = parse_int(opts, "target_maxdim", maximum(maxdim))
    target_sweep = parse_int(opts, "target_sweep", first_sweep_at_maxdim(maxdim, target_maxdim))
    stable_sweeps = parse_int(opts, "stable_sweeps", 1)
    energy_tol = parse_float(opts, "energy_tol", 1e-7)
    density_tol = parse_float(opts, "density_tol", 1e-5)
    truncerr_tol = parse_float(opts, "truncerr_tol", Inf)
    seed = parse_int(opts, "seed", 1234)
    init_linkdim = parse_int(opts, "init_linkdim", 8)
    init_density = get(opts, "init_density", "")
    checkpoint_in = get(opts, "checkpoint_in", "")
    threaded_blocksparse = parse_bool(opts, "threaded_blocksparse", false)
    blas_threads = parse_int(opts, "blas_threads", 1)
    disable_strided_threads = parse_bool(opts, "disable_strided_threads", threaded_blocksparse)
    output = get(opts, "output", "dmrg/output_Lx$(Lx)")
    checkpoint_out = get(opts, "checkpoint_out", joinpath(output, "checkpoint_final.jls"))
    mkpath(output)

    configure_threading(;
        threaded_blocksparse,
        blas_threads,
        disable_strided_threads,
    )

    start_time = now()
    @printf("Cylinder DMRG: Lx=%d Ly=%d Ns=%d Nuc=%d Np=%d phi=%.12g\n", Lx, Ly, lat.Ns, lat.Nuc, Np, phi)
    @printf("Parameters: t1=%.6g t3=%.6g V1=%.6g V2=%.6g V3=%.6g\n", params.t1, params.t3, params.V1, params.V2, params.V3)
    @printf("Sweeps=%d maxdim=%s cutoff=%.3g\n", adaptive ? max_sweeps : nsweeps, join(maxdim, ","), cutoff)
    adaptive && @printf(
        "Adaptive convergence: min_sweeps=%d target_sweep=%d stable_sweeps=%d energy_tol=%.3g density_tol=%.3g truncerr_tol=%.3g\n",
        min_sweeps,
        target_sweep,
        stable_sweeps,
        energy_tol,
        density_tol,
        truncerr_tol,
    )

    sites = nothing
    psi0 = nothing
    init_source = "default_random_product"
    if !isempty(checkpoint_in)
        cp = load_ground_checkpoint(checkpoint_in, lat, params; Lx, Ly, Np, phi)
        sites = cp.sites
        psi0 = cp.psi
        init_source = "checkpoint:$checkpoint_in"
        @printf("Loaded ground-state checkpoint %s: energy=%.16g saved_at=%s\n", checkpoint_in, cp.energy, checkpoint_property(cp, :saved_at, "unknown"))
        flush(stdout)
    elseif !isempty(init_density)
        sites, psi0 = density_seed_mps(init_density, lat, Np; seed, init_linkdim)
        init_source = "density_product:$init_density"
        @printf("Built density-biased product MPS from %s with init_linkdim=%d seed=%d\n", init_density, init_linkdim, seed)
        flush(stdout)
    end

    result = run_dmrg(
        lat,
        params,
        Np;
        phi_y=phi,
        nsweeps,
        maxdim,
        cutoff,
        seed,
        sites,
        psi0,
        init_linkdim,
        adaptive,
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
    nx = collect(0:(lat.Lx - 1))
    write_vector(joinpath(output, "density_x.dat"), nx, y_averaged_density(lat, result.density); header="# x y_averaged_density")
    write_complex_matrix(joinpath(output, "green.dat"), green_function(result.psi); header="# i j Re(<c_i^dag c_j>) Im(<c_i^dag c_j>)")
    write_matrix(joinpath(output, "connected_density.dat"), connected_density_correlation(result.psi); header="# i j <n_i n_j>-<n_i><n_j>")
    write_convergence(joinpath(output, "convergence.dat"), result.convergence)
    checkpoint_saved = save_ground_checkpoint(checkpoint_out, result, lat, params; Lx, Ly, Np, phi)

    write_summary(
        joinpath(output, "summary.dat"),
        result,
        lat,
        params;
        Lx,
        Ly,
        Np,
        phi,
        nsweeps,
        max_sweeps,
        maxdim,
        cutoff,
        adaptive,
        min_sweeps,
        target_sweep,
        stable_sweeps,
        energy_tol,
        density_tol,
        truncerr_tol,
        seed,
        init_linkdim,
        init_source,
        checkpoint_in,
        checkpoint_out,
        checkpoint_saved,
        start_time,
        end_time,
    )
end

main(ARGS)
