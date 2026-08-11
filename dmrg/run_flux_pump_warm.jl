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
    flush(stdout)
end

function first_sweep_at_maxdim(maxdim::Vector{Int}, target_maxdim::Int)
    target_maxdim > 0 || throw(ArgumentError("target_maxdim must be positive"))
    idx = findfirst(>=(target_maxdim), maxdim)
    return idx === nothing ? length(maxdim) : idx
end

function point_dir(root::AbstractString, step::Int)
    return joinpath(root, @sprintf("phi_%03d", step))
end

function convergence_status(records)
    isempty(records) && return (0, false, "no_records", Inf, Inf, Inf)
    last = records[end]
    return (last.sweep, last.converged, last.reason, last.delta_energy, last.max_density_delta, last.maxerr)
end

function write_convergence_attempts(path::AbstractString, attempts)
    open(path, "w") do io
        println(io, "# total_sweep attempt local_sweep energy delta_energy max_density_delta maxlinkdim maxerr converged reason")
        total_offset = 0
        for (attempt_idx, records) in enumerate(attempts)
            for rec in records
                @printf(
                    io,
                    "%d %d %d %.16g %.16g %.16g %d %.16g %s %s\n",
                    total_offset + rec.sweep,
                    attempt_idx,
                    rec.sweep,
                    rec.energy,
                    rec.delta_energy,
                    rec.max_density_delta,
                    rec.maxlinkdim,
                    rec.maxerr,
                    string(rec.converged),
                    rec.reason,
                )
            end
            total_offset += isempty(records) ? 0 : records[end].sweep
        end
    end
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
    total_steps::Int,
    phi::Float64,
    cycles::Int,
    steps_per_2pi::Int,
    direction::String,
    total_sweeps::Int,
    attempts_count::Int,
    converged::Bool,
    reason::String,
    delta_energy::Float64,
    max_density_delta::Float64,
    maxerr::Float64,
    max_sweeps::Int,
    retry_sweeps::Int,
    max_retries::Int,
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
        @printf(io, "trajectory_step %d\ntrajectory_total_steps %d\n", step, total_steps)
        @printf(io, "cycles %d\nsteps_per_2pi %d\ndirection %s\n", cycles, steps_per_2pi, direction)
        @printf(io, "Lx %d\nLy %d\nNs %d\nNuc %d\nNp %d\n", Lx, Ly, lat.Ns, lat.Nuc, Np)
        @printf(io, "t1 %.16g\nt3 %.16g\nV1 %.16g\nV2 %.16g\nV3 %.16g\n", params.t1, params.t3, params.V1, params.V2, params.V3)
        @printf(io, "actual_sweeps %d\nattempts %d\n", total_sweeps, attempts_count)
        @printf(io, "max_sweeps %d\nretry_sweeps %d\nmax_retries %d\nmin_sweeps %d\n", max_sweeps, retry_sweeps, max_retries, min_sweeps)
        @printf(io, "target_sweep %d\nstable_sweeps %d\n", target_sweep, stable_sweeps)
        @printf(io, "maxdim %s\ncutoff %.16g\n", join(maxdim, ","), cutoff)
        @printf(io, "energy_tol %.16g\ndensity_tol %.16g\ntruncerr_tol %.16g\n", energy_tol, density_tol, truncerr_tol)
        @printf(io, "converged %s\nconvergence_reason %s\n", string(converged), reason)
        @printf(io, "final_delta_energy %.16g\nfinal_max_density_delta %.16g\nfinal_maxerr %.16g\n", delta_energy, max_density_delta, maxerr)
        @printf(io, "run_started %s\nrun_finished %s\n", Dates.format(start_time, dateformat"yyyy-mm-ddTHH:MM:SS"), Dates.format(end_time, dateformat"yyyy-mm-ddTHH:MM:SS"))
        @printf(io, "julia_version %s\njulia_num_threads %d\n", string(VERSION), Threads.nthreads())
        @printf(io, "blas_threads %d\nthreaded_blocksparse %s\nstrided_threads %d\n", BLAS.get_num_threads(), string(ITensors.using_threaded_blocksparse()), ITensors.Strided.get_num_threads())
        @printf(io, "pbs_jobid %s\nhostname %s\n", get(ENV, "PBS_JOBID", "none"), get(ENV, "HOSTNAME", "unknown"))
    end
end

function run_until_converged(
    lat,
    params,
    Np;
    phi,
    sites,
    psi0,
    seed,
    max_sweeps,
    retry_sweeps,
    max_retries,
    maxdim,
    cutoff,
    min_sweeps,
    target_sweep,
    stable_sweeps,
    energy_tol,
    density_tol,
    truncerr_tol,
)
    attempts = Vector{Vector{ConvergenceRecord}}()
    current_sites = sites
    current_psi = psi0
    result = nothing
    total_sweeps = 0
    for attempt in 0:max_retries
        sweeps = attempt == 0 ? max_sweeps : retry_sweeps
        @printf("  attempt %d/%d max_sweeps=%d\n", attempt + 1, max_retries + 1, sweeps)
        flush(stdout)
        result = run_dmrg(
            lat,
            params,
            Np;
            phi_y=phi,
            nsweeps=sweeps,
            maxdim,
            cutoff,
            outputlevel=1,
            seed,
            sites=current_sites,
            psi0=current_psi,
            adaptive=true,
            max_sweeps=sweeps,
            min_sweeps,
            target_sweep,
            stable_sweeps,
            energy_tol,
            density_tol,
            truncerr_tol,
        )
        push!(attempts, result.convergence)
        local_sweeps, converged, reason, delta_energy, max_density_delta, maxerr = convergence_status(result.convergence)
        total_sweeps += local_sweeps
        @printf(
            "  attempt %d finished: energy=%.16g local_sweeps=%d total_sweeps=%d converged=%s reason=%s |dE|=%.3e max|dn|=%.3e maxerr=%.3e\n",
            attempt + 1,
            result.energy,
            local_sweeps,
            total_sweeps,
            string(converged),
            reason,
            delta_energy,
            max_density_delta,
            maxerr,
        )
        flush(stdout)
        current_sites = result.sites
        current_psi = result.psi
        converged && return result, attempts, total_sweeps, converged, reason, delta_energy, max_density_delta, maxerr
    end
    local_sweeps, converged, reason, delta_energy, max_density_delta, maxerr = convergence_status(result.convergence)
    return result, attempts, total_sweeps, converged, reason, delta_energy, max_density_delta, maxerr
end

function main(args)
    opts = parse_cli(args)
    Lx = parse_int(opts, "Lx", 15)
    Ly = parse_int(opts, "Ly", 6)
    cycles = parse_int(opts, "cycles", 3)
    steps_per_2pi = parse_int(opts, "steps_per_2pi", 16)
    direction = lowercase(get(opts, "direction", "forward"))
    direction in ("forward", "backward") || error("direction must be forward or backward")
    total_steps = cycles * steps_per_2pi
    sign = direction == "forward" ? 1.0 : -1.0
    phi0 = direction == "forward" ? 0.0 : 2π * cycles
    output = get(opts, "output", @sprintf("dmrg/flux_Lx%d_warm_%s_spp%d", Lx, direction, steps_per_2pi))
    mkpath(joinpath(output, "logs"))

    threaded_blocksparse = parse_bool(opts, "threaded_blocksparse", true)
    blas_threads = parse_int(opts, "blas_threads", 1)
    disable_strided_threads = parse_bool(opts, "disable_strided_threads", threaded_blocksparse)
    configure_threading(;
        threaded_blocksparse,
        blas_threads,
        disable_strided_threads,
    )

    lat = cylinder_lattice(Lx; Ly)
    Np = parse_int(opts, "Np", filling_np(lat))
    params = CylinderModelParams(
        t1=parse_float(opts, "t1", 1.0),
        t3=parse_float(opts, "t3", 0.2),
        V1=parse_float(opts, "V1", 1.0),
        V2=parse_float(opts, "V2", 0.0),
        V3=parse_float(opts, "V3", 0.0),
    )
    max_sweeps = parse_int(opts, "max_sweeps", 40)
    retry_sweeps = parse_int(opts, "retry_sweeps", 40)
    max_retries = parse_int(opts, "max_retries", 2)
    maxdim = parse_int_list(get(opts, "maxdim", "200,400,800,1200,1600,2000"))
    cutoff = parse_float(opts, "cutoff", 1e-9)
    min_sweeps = parse_int(opts, "min_sweeps", 4)
    target_maxdim = parse_int(opts, "target_maxdim", maximum(maxdim))
    target_sweep = parse_int(opts, "target_sweep", first_sweep_at_maxdim(maxdim, target_maxdim))
    stable_sweeps = parse_int(opts, "stable_sweeps", 3)
    energy_tol = parse_float(opts, "energy_tol", 1e-6)
    density_tol = parse_float(opts, "density_tol", 2e-5)
    truncerr_tol = parse_float(opts, "truncerr_tol", 5e-8)
    seed = parse_int(opts, "seed", 1234)
    require_converged = parse_bool(opts, "require_converged", true)

    @printf("Warm flux pump trajectory\n")
    @printf("Lx=%d Ly=%d Np=%d cycles=%d steps_per_2pi=%d total_steps=%d direction=%s output=%s\n", Lx, Ly, Np, cycles, steps_per_2pi, total_steps, direction, output)
    @printf("DMRG max_sweeps=%d retry_sweeps=%d max_retries=%d maxdim=%s seed=%d\n", max_sweeps, retry_sweeps, max_retries, join(maxdim, ","), seed)
    flush(stdout)

    sites = nothing
    psi = nothing
    density0 = nothing
    xs = collect(0:(lat.Lx - 1))
    pump_path = joinpath(output, "pumping.dat")
    open(pump_path, "w") do pumpio
        println(pumpio, "# step phi energy cumulative_left_edge cumulative_mid cumulative_right_edge actual_sweeps attempts converged reason")
        for step in 0:total_steps
            phi = phi0 + sign * 2π * step / steps_per_2pi
            dir = point_dir(output, step)
            mkpath(dir)
            start_time = now()
            @printf("\n=== trajectory step %d/%d phi=%.16g ===\n", step, total_steps, phi)
            flush(stdout)
            result, attempts, total_sweeps, converged, reason, delta_energy, max_density_delta, maxerr = run_until_converged(
                lat,
                params,
                Np;
                phi,
                sites,
                psi0=psi,
                seed=seed + step,
                max_sweeps,
                retry_sweeps,
                max_retries,
                maxdim,
                cutoff,
                min_sweeps,
                target_sweep,
                stable_sweeps,
                energy_tol,
                density_tol,
                truncerr_tol,
            )
            end_time = now()
            sites = result.sites
            psi = result.psi
            dens = result.density
            density0 === nothing && (density0 = copy(dens))
            delta = dens .- density0
            cumulative = cumulative_charge(lat, delta)
            mid = cumulative[clamp(lat.Lx ÷ 2, 1, lat.Lx)]

            write_density(joinpath(dir, "density.dat"), lat, dens)
            write_vector(joinpath(dir, "density_x.dat"), xs, y_averaged_density(lat, dens); header="# x y_averaged_density")
            write_vector(joinpath(output, @sprintf("cumulative_phi_%03d.dat", step)), xs, cumulative; header="# x cumulative_delta_charge_from_trajectory_start")
            write_convergence_attempts(joinpath(dir, "convergence.dat"), attempts)
            write_point_summary(
                joinpath(dir, "summary.dat"),
                result,
                lat,
                params;
                Lx,
                Ly,
                Np,
                step,
                total_steps,
                phi,
                cycles,
                steps_per_2pi,
                direction,
                total_sweeps,
                attempts_count=length(attempts),
                converged,
                reason,
                delta_energy,
                max_density_delta,
                maxerr,
                max_sweeps,
                retry_sweeps,
                max_retries,
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
            @printf(
                pumpio,
                "%d %.16g %.16g %.16g %.16g %.16g %d %d %s %s\n",
                step,
                phi,
                result.energy,
                cumulative[1],
                mid,
                cumulative[end],
                total_sweeps,
                length(attempts),
                string(converged),
                reason,
            )
            flush(pumpio)
            @printf(
                "trajectory step %d done: energy=%.16g mid_charge=%.12g total_sweeps=%d attempts=%d converged=%s reason=%s\n",
                step,
                result.energy,
                mid,
                total_sweeps,
                length(attempts),
                string(converged),
                reason,
            )
            flush(stdout)
            if require_converged && !converged
                error("trajectory step $step did not converge after $(length(attempts)) attempts and $total_sweeps total sweeps")
            end
        end
    end
    @printf("Warm flux pump complete: %s\n", pump_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
