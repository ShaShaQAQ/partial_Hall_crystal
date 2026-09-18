module FluxPumpMerge

include("DMRGFullModel.jl")
using .DMRGFullModel
using Printf

export FluxPumpRow, merge_flux_data, write_summary

struct FluxPumpRow
    step::Int
    phi::Float64
    energy::Float64
    left_edge::Float64
    mid::Float64
    right_edge::Float64
end

function read_density(path::AbstractString)
    vals = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            startswith(strip(line), "#") && continue
            isempty(strip(line)) && continue
            parts = split(line)
            length(parts) >= 6 || error("bad density row in $path: $line")
            push!(vals, parse(Float64, parts[6]))
        end
    end
    return vals
end

function read_summary_value(path::AbstractString, key::AbstractString)
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        startswith(line, "#") && continue
        parts = split(line)
        length(parts) >= 2 || continue
        parts[1] == key && return parse(Float64, parts[2])
    end
    error("missing $key in $path; contents:\n$(read(path, String))")
end

function write_summary(path::AbstractString, phi::Real, energy::Real)
    open(path, "w") do io
        @printf(io, "phi %.16g\n", Float64(phi))
        @printf(io, "energy %.16g\n", Float64(energy))
    end
end

function point_dir(root::AbstractString, step::Int)
    return joinpath(root, @sprintf("phi_%03d", step))
end

function merge_flux_data(root::AbstractString; Lx::Int, Ly::Int, steps::Int)
    lat = cylinder_lattice(Lx; Ly)
    density0 = read_density(joinpath(point_dir(root, 0), "density.dat"))
    length(density0) == lat.Ns || error("density length $(length(density0)) != Ns $(lat.Ns)")

    rows = FluxPumpRow[]
    open(joinpath(root, "pumping.dat"), "w") do pumpio
        println(pumpio, "# step phi energy cumulative_left_edge cumulative_mid cumulative_right_edge")
        for m in 0:steps
            dir = point_dir(root, m)
            dens = read_density(joinpath(dir, "density.dat"))
            length(dens) == lat.Ns || error("density length $(length(dens)) != Ns $(lat.Ns) for step $m")
            delta = dens .- density0
            cumulative = cumulative_charge(lat, delta)
            xs = collect(0:(lat.Lx - 1))
            write_vector(
                joinpath(root, @sprintf("cumulative_phi_%03d.dat", m)),
                xs,
                cumulative;
                header="# x cumulative_delta_charge_from_left",
            )
            summary = joinpath(dir, "summary.dat")
            phi = read_summary_value(summary, "phi")
            energy = read_summary_value(summary, "energy")
            mid = cumulative[clamp(lat.Lx ÷ 2, 1, lat.Lx)]
            row = FluxPumpRow(m, phi, energy, cumulative[1], mid, cumulative[end])
            push!(rows, row)
            @printf(
                pumpio,
                "%d %.16g %.16g %.16g %.16g %.16g\n",
                row.step,
                row.phi,
                row.energy,
                row.left_edge,
                row.mid,
                row.right_edge,
            )
        end
    end
    return rows
end

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

function main(args)
    opts = parse_cli(args)
    root = get(opts, "input", "dmrg/flux_parallel")
    Lx = parse_int(opts, "Lx", 15)
    Ly = parse_int(opts, "Ly", 6)
    steps = parse_int(opts, "steps", 12)
    rows = merge_flux_data(root; Lx, Ly, steps)
    last = rows[end]
    @printf("Merged %d flux points into %s\n", length(rows), joinpath(root, "pumping.dat"))
    @printf("Final pumped charge at mid cut: %.12g\n", last.mid)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

end
