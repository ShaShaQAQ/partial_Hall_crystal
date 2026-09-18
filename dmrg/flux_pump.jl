#!/usr/bin/env julia

include("DMRGFullModel.jl")
using .DMRGFullModel
using Printf

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

function main(args)
    opts = parse_cli(args)
    Lx = parse_int(opts, "Lx", 4)
    Ly = parse_int(opts, "Ly", 6)
    lat = cylinder_lattice(Lx; Ly)
    Np = parse_int(opts, "Np", filling_np(lat))
    steps = parse_int(opts, "steps", 12)
    params = CylinderModelParams(
        t1=parse_float(opts, "t1", 1.0),
        t3=parse_float(opts, "t3", 0.2),
        V1=parse_float(opts, "V1", 1.0),
        V2=parse_float(opts, "V2", 0.0),
        V3=parse_float(opts, "V3", 0.0),
    )
    nsweeps = parse_int(opts, "nsweeps", 6)
    maxdim = parse_int_list(get(opts, "maxdim", "20,50,100"))
    cutoff = parse_float(opts, "cutoff", 1e-8)
    output = get(opts, "output", "dmrg/flux_Lx$(Lx)_steps$(steps)")
    mkpath(output)

    sites = nothing
    psi = nothing
    density0 = nothing

    open(joinpath(output, "pumping.dat"), "w") do pumpio
        println(pumpio, "# step phi energy cumulative_left_edge cumulative_mid cumulative_right_edge")
        for m in 0:steps
            phi = 2π * m / steps
            @printf("Flux step %d/%d phi=%.12g\n", m, steps, phi)
            result = run_dmrg(lat, params, Np; phi_y=phi, nsweeps, maxdim, cutoff, outputlevel=1, sites, psi0=psi)
            sites = result.sites
            psi = result.psi
            density0 === nothing && (density0 = copy(result.density))
            delta = result.density .- density0
            cumulative = cumulative_charge(lat, delta)
            write_density(joinpath(output, @sprintf("density_phi_%03d.dat", m)), lat, result.density)
            write_vector(joinpath(output, @sprintf("cumulative_phi_%03d.dat", m)), collect(0:(lat.Lx - 1)), cumulative; header="# x cumulative_delta_charge_from_left")
            mid = cumulative[clamp(lat.Lx ÷ 2, 1, lat.Lx)]
            @printf(pumpio, "%d %.16g %.16g %.16g %.16g %.16g\n", m, phi, result.energy, cumulative[1], mid, cumulative[end])
            flush(pumpio)
        end
    end
end

main(ARGS)
