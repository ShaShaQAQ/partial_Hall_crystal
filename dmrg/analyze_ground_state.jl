#!/usr/bin/env julia

module GroundStateDiagnostics

using LinearAlgebra
using Printf
using Serialization
using Statistics
using ITensors
using ITensorMPS

export DensityPoint,
    ConvergencePoint,
    StructureFactorPoint,
    GreenDecayPoint,
    EntanglementCut,
    DiagnosticReport,
    analyze_output,
    main

struct DensityPoint
    site::Int
    x::Int
    y::Int
    rx::Float64
    ry::Float64
    density::Float64
end

struct ConvergencePoint
    sweep::Int
    energy::Float64
    delta_energy::Float64
    max_density_delta::Float64
    maxlinkdim::Int
    maxerr::Float64
    converged::Bool
    reason::String
end

struct StructureFactorPoint
    mx::Int
    my::Int
    value::Float64
end

struct GreenDecayPoint
    dx::Int
    count::Int
    mean_abs::Float64
    min_abs::Float64
    max_abs::Float64
end

struct EntanglementCut
    bond::Int
    entropy::Float64
    schmidt_rank::Int
end

struct DensitySummary
    total_density::Float64
    mean_density::Float64
    min_density::Float64
    max_density::Float64
    cdw_contrast::Float64
    sublattice::Dict{Int,NamedTuple}
    x_profile::Vector{NamedTuple}
end

struct DiagnosticReport
    summary::Dict{String,String}
    density::DensitySummary
    convergence::Vector{ConvergencePoint}
    structure_factor::Vector{StructureFactorPoint}
    structure_factor_peaks::Vector{StructureFactorPoint}
    green_decay::Vector{GreenDecayPoint}
    green_correlation_length::Union{Nothing,Float64}
    entanglement::Vector{EntanglementCut}
    entanglement_spectrum::Vector{NamedTuple}
    entanglement_error::String
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

function parse_float_token(token::AbstractString)
    s = strip(token)
    lowercase(s) == "inf" && return Inf
    lowercase(s) == "-inf" && return -Inf
    lowercase(s) == "nan" && return NaN
    return parse(Float64, s)
end

parse_bool_token(token::AbstractString) = lowercase(strip(token)) in ("1", "true", "yes", "on")

function read_key_value_file(path::AbstractString)
    data = Dict{String,String}()
    isfile(path) || return data
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        fields = split(stripped; limit=2)
        length(fields) == 2 || continue
        data[fields[1]] = strip(fields[2])
    end
    return data
end

function read_density_file(path::AbstractString)
    isfile(path) || error("missing density file: $path")
    points = DensityPoint[]
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        fields = split(stripped)
        length(fields) >= 6 || error("malformed density row in $path: $line")
        push!(
            points,
            DensityPoint(
                parse(Int, fields[1]),
                parse(Int, fields[2]),
                parse(Int, fields[3]),
                parse_float_token(fields[4]),
                parse_float_token(fields[5]),
                parse_float_token(fields[6]),
            ),
        )
    end
    sort!(points; by=p -> p.site)
    return points
end

function read_convergence_file(path::AbstractString)
    isfile(path) || return ConvergencePoint[]
    rows = ConvergencePoint[]
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        fields = split(stripped)
        length(fields) >= 8 || error("malformed convergence row in $path: $line")
        push!(
            rows,
            ConvergencePoint(
                parse(Int, fields[1]),
                parse_float_token(fields[2]),
                parse_float_token(fields[3]),
                parse_float_token(fields[4]),
                parse(Int, fields[5]),
                parse_float_token(fields[6]),
                parse_bool_token(fields[7]),
                fields[8],
            ),
        )
    end
    sort!(rows; by=r -> r.sweep)
    return rows
end

function read_real_matrix_file(path::AbstractString, n::Int)
    isfile(path) || return zeros(Float64, 0, 0)
    matrix = zeros(Float64, n, n)
    seen = falses(n, n)
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        fields = split(stripped)
        length(fields) >= 3 || error("malformed matrix row in $path: $line")
        i = parse(Int, fields[1])
        j = parse(Int, fields[2])
        1 <= i <= n || error("row index $i out of range 1:$n in $path")
        1 <= j <= n || error("column index $j out of range 1:$n in $path")
        matrix[i, j] = parse_float_token(fields[3])
        seen[i, j] = true
    end
    any(!, seen) && @warn "matrix file does not contain every i,j entry" path
    return matrix
end

function read_complex_matrix_file(path::AbstractString, n::Int)
    isfile(path) || return zeros(ComplexF64, 0, 0)
    matrix = zeros(ComplexF64, n, n)
    seen = falses(n, n)
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        fields = split(stripped)
        length(fields) >= 4 || error("malformed complex matrix row in $path: $line")
        i = parse(Int, fields[1])
        j = parse(Int, fields[2])
        1 <= i <= n || error("row index $i out of range 1:$n in $path")
        1 <= j <= n || error("column index $j out of range 1:$n in $path")
        matrix[i, j] = complex(parse_float_token(fields[3]), parse_float_token(fields[4]))
        seen[i, j] = true
    end
    any(!, seen) && @warn "complex matrix file does not contain every i,j entry" path
    return matrix
end

function group_summary(values::AbstractVector{<:Real})
    isempty(values) && return (count=0, total=0.0, mean=NaN, min=NaN, max=NaN)
    vals = Float64.(values)
    return (
        count=length(vals),
        total=sum(vals),
        mean=mean(vals),
        min=minimum(vals),
        max=maximum(vals),
    )
end

function density_summary(points::Vector{DensityPoint})
    isempty(points) && error("density file contains no points")
    dens = [p.density for p in points]
    sub = Dict{Int,NamedTuple}()
    for key in sort(unique(mod(p.y, 2) for p in points))
        sub[key] = group_summary([p.density for p in points if mod(p.y, 2) == key])
    end
    xs = NamedTuple[]
    for x in sort(unique(p.x for p in points))
        push!(
            xs,
            merge((x=x,), group_summary([p.density for p in points if p.x == x])),
        )
    end
    return DensitySummary(
        sum(dens),
        mean(dens),
        minimum(dens),
        maximum(dens),
        maximum(dens) - minimum(dens),
        sub,
        xs,
    )
end

function dimensions_from_points(points::Vector{DensityPoint})
    isempty(points) && error("cannot infer dimensions from empty density points")
    Lx = maximum(p.x for p in points) + 1
    Ly = maximum(p.y for p in points) + 1
    return Lx, Ly
end

function structure_factor(points::Vector{DensityPoint}, corr::AbstractMatrix{<:Real})
    size(corr, 1) == 0 && return StructureFactorPoint[]
    n = length(points)
    size(corr) == (n, n) || error("connected density matrix has size $(size(corr)), expected ($n, $n)")
    Lx, Ly = dimensions_from_points(points)
    rows = StructureFactorPoint[]
    for mx in 0:(Lx - 1), my in 0:(Ly - 1)
        accum = 0.0 + 0.0im
        for i in 1:n, j in 1:n
            phase = 2π * (mx * (points[i].x - points[j].x) / Lx + my * (points[i].y - points[j].y) / Ly)
            accum += corr[i, j] * cis(phase)
        end
        push!(rows, StructureFactorPoint(mx, my, max(real(accum) / n, 0.0)))
    end
    sort!(rows; by=row -> (row.mx, row.my))
    return rows
end

function top_structure_peaks(rows::Vector{StructureFactorPoint}, n::Int)
    candidates = [row for row in rows if !(row.mx == 0 && row.my == 0)]
    sort!(candidates; by=row -> (-row.value, row.mx, row.my))
    return candidates[1:min(n, length(candidates))]
end

function green_decay_by_x(points::Vector{DensityPoint}, green::AbstractMatrix{<:Complex})
    size(green, 1) == 0 && return GreenDecayPoint[]
    n = length(points)
    size(green) == (n, n) || error("Green matrix has size $(size(green)), expected ($n, $n)")
    grouped = Dict{Int,Vector{Float64}}()
    for i in 1:n, j in 1:n
        i == j && continue
        points[i].y == points[j].y || continue
        dx = abs(points[i].x - points[j].x)
        dx == 0 && continue
        push!(get!(grouped, dx, Float64[]), abs(green[i, j]))
    end
    rows = GreenDecayPoint[]
    for dx in sort(collect(keys(grouped)))
        values = grouped[dx]
        push!(rows, GreenDecayPoint(dx, length(values), mean(values), minimum(values), maximum(values)))
    end
    return rows
end

function estimate_green_correlation_length(rows::Vector{GreenDecayPoint})
    xs = Float64[]
    ys = Float64[]
    for row in rows
        row.dx > 0 || continue
        row.mean_abs > 0 || continue
        push!(xs, Float64(row.dx))
        push!(ys, log(row.mean_abs))
    end
    length(xs) >= 2 || return nothing
    xbar = mean(xs)
    ybar = mean(ys)
    denom = sum((x - xbar)^2 for x in xs)
    denom == 0 && return nothing
    slope = sum((x - xbar) * (y - ybar) for (x, y) in zip(xs, ys)) / denom
    slope < 0 || return nothing
    return -1 / slope
end

function parse_bond_selector(selector::AbstractString, nsites::Int)
    s = lowercase(strip(selector))
    s in ("", "none", "false", "0") && return Int[]
    if s == "center"
        return [max(1, nsites ÷ 2)]
    elseif s == "all"
        return collect(1:(nsites - 1))
    elseif occursin(":", s)
        parts = split(s, ":")
        length(parts) == 2 || error("bad entanglement bond range: $selector")
        lo = parse(Int, parts[1])
        hi = parse(Int, parts[2])
        return collect(max(1, lo):min(nsites - 1, hi))
    else
        return sort(unique(parse(Int, x) for x in split(s, ",") if !isempty(x)))
    end
end

function checkpoint_property(obj, name::Symbol, default=nothing)
    return hasproperty(obj, name) ? getproperty(obj, name) : default
end

function resolve_checkpoint(input_dir::AbstractString, summary::Dict{String,String}, checkpoint::AbstractString)
    raw = strip(checkpoint)
    raw == "" && return ""
    raw == "auto" && (raw = get(summary, "checkpoint_out", ""))
    raw == "" && return ""
    candidates = String[]
    push!(candidates, raw)
    push!(candidates, abspath(raw))
    push!(candidates, joinpath(input_dir, raw))
    push!(candidates, joinpath(input_dir, basename(raw)))
    for candidate in candidates
        isfile(candidate) && return candidate
    end
    return ""
end

function left_inds_for_bond(psi, bond::Int)
    if bond == 1
        return (siteind(psi, 1),)
    end
    return (linkind(psi, bond - 1), siteind(psi, bond))
end

function checkpoint_entanglement(checkpoint::AbstractString, selector::AbstractString; spectrum_levels::Int=12)
    obj = deserialize(checkpoint)
    psi0 = checkpoint_property(obj, :psi, nothing)
    psi0 === nothing && error("checkpoint does not contain psi")
    psi = copy(psi0)
    nsites = length(psi)
    cuts = EntanglementCut[]
    spectrum_rows = NamedTuple[]
    for bond in parse_bond_selector(selector, nsites)
        1 <= bond < nsites || error("entanglement bond $bond out of range 1:$(nsites - 1)")
        orthogonalize!(psi, bond)
        decomp = svd(psi[bond], left_inds_for_bond(psi, bond)...; cutoff=0.0)
        spec = decomp.spec
        weights = collect(eigs(spec))
        ent = entropy(spec)
        push!(cuts, EntanglementCut(bond, ent, length(weights)))
        for (level, weight) in enumerate(weights[1:min(spectrum_levels, length(weights))])
            energy = weight > 0 ? -log(weight) : Inf
            push!(
                spectrum_rows,
                (
                    bond=bond,
                    level=level,
                    weight=weight,
                    entanglement_energy=energy,
                ),
            )
        end
        decomp = nothing
        GC.gc(false)
    end
    return cuts, spectrum_rows
end

function analyze_output(
    input_dir::AbstractString;
    top_peaks::Int=12,
    checkpoint::AbstractString="",
    entanglement_bonds::AbstractString="none",
    spectrum_levels::Int=12,
)
    summary = read_key_value_file(joinpath(input_dir, "summary.dat"))
    points = read_density_file(joinpath(input_dir, "density.dat"))
    dens = density_summary(points)
    convergence = read_convergence_file(joinpath(input_dir, "convergence.dat"))
    corr = read_real_matrix_file(joinpath(input_dir, "connected_density.dat"), length(points))
    sf = structure_factor(points, corr)
    peaks = top_structure_peaks(sf, top_peaks)
    green = read_complex_matrix_file(joinpath(input_dir, "green.dat"), length(points))
    green_decay = green_decay_by_x(points, green)
    green_xi = estimate_green_correlation_length(green_decay)

    ent = EntanglementCut[]
    ent_spec = NamedTuple[]
    ent_error = ""
    resolved_checkpoint = resolve_checkpoint(input_dir, summary, checkpoint)
    if !isempty(resolved_checkpoint) && lowercase(strip(entanglement_bonds)) != "none"
        try
            ent, ent_spec = checkpoint_entanglement(
                resolved_checkpoint,
                entanglement_bonds;
                spectrum_levels,
            )
        catch err
            ent_error = sprint(showerror, err)
        end
    elseif !isempty(strip(checkpoint)) && strip(checkpoint) != "auto"
        ent_error = "checkpoint not found: $checkpoint"
    end

    return DiagnosticReport(
        summary,
        dens,
        convergence,
        sf,
        peaks,
        green_decay,
        green_xi,
        ent,
        ent_spec,
        ent_error,
    )
end

function csv_escape(value)
    s = string(value)
    if occursin(",", s) || occursin("\"", s) || occursin("\n", s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_csv(path::AbstractString, header::Vector{String}, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(csv_escape.(row), ","))
        end
    end
end

function fmt(x::Real)
    if !isfinite(x)
        return string(x)
    end
    return @sprintf("%.12g", x)
end

function write_density_by_sublattice(path::AbstractString, report::DiagnosticReport)
    rows = []
    for key in sort(collect(keys(report.density.sublattice)))
        row = report.density.sublattice[key]
        push!(rows, [key, row.count, fmt(row.total), fmt(row.mean), fmt(row.min), fmt(row.max)])
    end
    write_csv(path, ["sublattice", "count", "total", "mean", "min", "max"], rows)
end

function write_density_x(path::AbstractString, report::DiagnosticReport)
    rows = [[row.x, row.count, fmt(row.total), fmt(row.mean), fmt(row.min), fmt(row.max)] for row in report.density.x_profile]
    write_csv(path, ["x", "count", "total", "mean", "min", "max"], rows)
end

function write_structure_factor(path::AbstractString, rows::Vector{StructureFactorPoint})
    write_csv(path, ["mx", "my", "S"], [[row.mx, row.my, fmt(row.value)] for row in rows])
end

function write_green_decay(path::AbstractString, rows::Vector{GreenDecayPoint})
    write_csv(path, ["dx", "count", "mean_abs", "min_abs", "max_abs"], [[row.dx, row.count, fmt(row.mean_abs), fmt(row.min_abs), fmt(row.max_abs)] for row in rows])
end

function write_entanglement(path::AbstractString, rows::Vector{EntanglementCut})
    write_csv(path, ["bond", "entropy", "schmidt_rank"], [[row.bond, fmt(row.entropy), row.schmidt_rank] for row in rows])
end

function write_entanglement_spectrum(path::AbstractString, rows::Vector{NamedTuple})
    write_csv(path, ["bond", "level", "weight", "entanglement_energy"], [[row.bond, row.level, fmt(row.weight), fmt(row.entanglement_energy)] for row in rows])
end

function write_diagnostics_summary(path::AbstractString, report::DiagnosticReport)
    open(path, "w") do io
        for key in ("energy", "Lx", "Ly", "Ns", "Nuc", "Np", "t1", "t3", "V1", "V2", "V3", "phi", "converged", "convergence_reason", "actual_sweeps")
            haskey(report.summary, key) && println(io, "$key $(report.summary[key])")
        end
        println(io, "density_total $(fmt(report.density.total_density))")
        println(io, "density_mean $(fmt(report.density.mean_density))")
        println(io, "density_min $(fmt(report.density.min_density))")
        println(io, "density_max $(fmt(report.density.max_density))")
        println(io, "density_cdw_contrast $(fmt(report.density.cdw_contrast))")
        for key in sort(collect(keys(report.density.sublattice)))
            row = report.density.sublattice[key]
            println(io, "sublattice_$(key)_mean $(fmt(row.mean))")
            println(io, "sublattice_$(key)_total $(fmt(row.total))")
        end
        if !isempty(report.convergence)
            last = report.convergence[end]
            println(io, "last_sweep $(last.sweep)")
            println(io, "last_delta_energy $(fmt(last.delta_energy))")
            println(io, "last_max_density_delta $(fmt(last.max_density_delta))")
            println(io, "last_maxerr $(fmt(last.maxerr))")
        end
        for (i, peak) in enumerate(report.structure_factor_peaks)
            println(io, "structure_peak_$(i)_mx $(peak.mx)")
            println(io, "structure_peak_$(i)_my $(peak.my)")
            println(io, "structure_peak_$(i)_S $(fmt(peak.value))")
        end
        if report.green_correlation_length !== nothing
            println(io, "green_correlation_length_x $(fmt(report.green_correlation_length))")
        end
        for cut in report.entanglement
            println(io, "entanglement_bond_$(cut.bond)_entropy $(fmt(cut.entropy))")
            println(io, "entanglement_bond_$(cut.bond)_schmidt_rank $(cut.schmidt_rank)")
        end
        !isempty(report.entanglement_error) && println(io, "entanglement_error $(report.entanglement_error)")
    end
end

html_escape(s) = replace(replace(replace(string(s), "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")

function color_hex(t::Real)
    x = clamp(Float64(t), 0.0, 1.0)
    r = round(Int, 255 * x)
    g = round(Int, 255 * (1 - abs(2x - 1) * 0.35))
    b = round(Int, 255 * (1 - x))
    return @sprintf("#%02x%02x%02x", r, g, b)
end

function draw_density_panel(parts, points::Vector{DensityPoint}, x0, y0, width, height)
    Lx, Ly = dimensions_from_points(points)
    vals = [p.density for p in points]
    vmin, vmax = minimum(vals), maximum(vals)
    span = max(vmax - vmin, eps(Float64))
    cellw = width / Lx
    cellh = height / Ly
    push!(parts, @sprintf("<text x=\"%.1f\" y=\"%.1f\" class=\"title\">density heatmap</text>", x0, y0 - 12))
    for p in points
        t = (p.density - vmin) / span
        rx = x0 + p.x * cellw
        ry = y0 + (Ly - 1 - p.y) * cellh
        push!(
            parts,
            @sprintf(
                "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\" stroke=\"#fff\" stroke-width=\"0.6\"><title>site %d x=%d y=%d n=%.6g</title></rect>",
                rx,
                ry,
                cellw,
                cellh,
                color_hex(t),
                p.site,
                p.x,
                p.y,
                p.density,
            ),
        )
    end
    push!(parts, @sprintf("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#333\"/>", x0, y0, width, height))
    push!(parts, @sprintf("<text x=\"%.1f\" y=\"%.1f\" class=\"note\">min %s, max %s</text>", x0, y0 + height + 18, fmt(vmin), fmt(vmax)))
end

function draw_polyline_panel(parts, title, rows, xfield, yfield, x0, y0, width, height; logy=false)
    isempty(rows) && return
    xs = [Float64(getproperty(row, xfield)) for row in rows]
    ys_raw = [Float64(getproperty(row, yfield)) for row in rows]
    valid = [isfinite(y) && (!logy || y > 0) for y in ys_raw]
    any(valid) || return
    ys = [logy ? log10(y) : y for y in ys_raw]
    xs = xs[valid]
    ys = ys[valid]
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    xmax == xmin && (xmax += 1)
    ymax == ymin && (ymax += 1)
    xmap(x) = x0 + (x - xmin) / (xmax - xmin) * width
    ymap(y) = y0 + height - (y - ymin) / (ymax - ymin) * height
    pts = join((@sprintf("%.2f,%.2f", xmap(x), ymap(y)) for (x, y) in zip(xs, ys)), " ")
    push!(parts, @sprintf("<text x=\"%.1f\" y=\"%.1f\" class=\"title\">%s</text>", x0, y0 - 12, html_escape(title)))
    push!(parts, @sprintf("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"white\" stroke=\"#333\"/>", x0, y0, width, height))
    push!(parts, @sprintf("<polyline points=\"%s\" fill=\"none\" stroke=\"#1f5f99\" stroke-width=\"2\"/>", pts))
    for (x, y) in zip(xs, ys)
        push!(parts, @sprintf("<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3\" fill=\"#1f5f99\"/>", xmap(x), ymap(y)))
    end
    label = logy ? "log10 scale" : "linear scale"
    push!(parts, @sprintf("<text x=\"%.1f\" y=\"%.1f\" class=\"note\">%s</text>", x0, y0 + height + 18, label))
end

function draw_peak_panel(parts, peaks::Vector{StructureFactorPoint}, x0, y0, width, height)
    isempty(peaks) && return
    n = min(8, length(peaks))
    shown = peaks[1:n]
    vmax = maximum(row.value for row in shown)
    vmax <= 0 && (vmax = 1)
    push!(parts, @sprintf("<text x=\"%.1f\" y=\"%.1f\" class=\"title\">top S(q) peaks</text>", x0, y0 - 12))
    push!(parts, @sprintf("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"white\" stroke=\"#333\"/>", x0, y0, width, height))
    barw = width / n * 0.62
    for (idx, peak) in enumerate(shown)
        cx = x0 + (idx - 0.5) * width / n
        bh = height * peak.value / vmax
        push!(parts, @sprintf("<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"#a34d42\"/>", cx - barw / 2, y0 + height - bh, barw, bh))
        push!(parts, @sprintf("<text x=\"%.2f\" y=\"%.2f\" class=\"tick\" text-anchor=\"middle\">(%d,%d)</text>", cx, y0 + height + 15, peak.mx, peak.my))
    end
end

function write_svg(path::AbstractString, points::Vector{DensityPoint}, report::DiagnosticReport)
    width = 1120
    height = 780
    parts = String[
        @sprintf("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">", width, height, width, height),
        "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111}.heading{font-size:22px;font-weight:700}.title{font-size:15px;font-weight:700}.note{font-size:12px;fill:#555}.tick{font-size:11px;fill:#333}</style>",
        "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>",
    ]
    title = get(report.summary, "converged", "unknown") == "true" ? "Converged DMRG ground-state diagnostics" : "DMRG ground-state diagnostics"
    push!(parts, @sprintf("<text x=\"42\" y=\"42\" class=\"heading\">%s</text>", title))
    push!(parts, @sprintf("<text x=\"42\" y=\"68\" class=\"note\">energy=%s, total density=%s, CDW contrast=%s</text>", get(report.summary, "energy", "unknown"), fmt(report.density.total_density), fmt(report.density.cdw_contrast)))
    draw_density_panel(parts, points, 48.0, 118.0, 440.0, 270.0)
    draw_polyline_panel(parts, "density averaged along y", report.density.x_profile, :x, :mean, 560.0, 118.0, 500.0, 185.0)
    draw_peak_panel(parts, report.structure_factor_peaks, 560.0, 405.0, 500.0, 190.0)
    if !isempty(report.convergence)
        draw_polyline_panel(parts, "density convergence", report.convergence, :sweep, :max_density_delta, 48.0, 470.0, 440.0, 190.0; logy=true)
    end
    if !isempty(report.green_decay)
        draw_polyline_panel(parts, "Green decay by dx", report.green_decay, :dx, :mean_abs, 48.0, 690.0, 440.0, 58.0; logy=true)
    end
    push!(parts, "</svg>")
    open(path, "w") do io
        println(io, join(parts, "\n"))
    end
end

function write_artifacts(input_dir::AbstractString, outdir::AbstractString, report::DiagnosticReport)
    mkpath(outdir)
    points = read_density_file(joinpath(input_dir, "density.dat"))
    write_diagnostics_summary(joinpath(outdir, "diagnostics_summary.dat"), report)
    write_density_by_sublattice(joinpath(outdir, "density_by_sublattice.csv"), report)
    write_density_x(joinpath(outdir, "density_x_summary.csv"), report)
    write_structure_factor(joinpath(outdir, "structure_factor.csv"), report.structure_factor)
    write_structure_factor(joinpath(outdir, "structure_factor_peaks.csv"), report.structure_factor_peaks)
    write_green_decay(joinpath(outdir, "green_decay_x.csv"), report.green_decay)
    write_entanglement(joinpath(outdir, "entanglement_entropy.csv"), report.entanglement)
    write_entanglement_spectrum(joinpath(outdir, "entanglement_spectrum.csv"), report.entanglement_spectrum)
    write_svg(joinpath(outdir, "ground_state_diagnostics.svg"), points, report)
end

function main(args=ARGS)
    opts = parse_cli(args)
    input = get(opts, "input", "")
    isempty(input) && error("missing --input=/path/to/dmrg/output")
    outdir = get(opts, "outdir", joinpath(input, "diagnostics"))
    top_peaks = parse(Int, get(opts, "top_peaks", "12"))
    checkpoint = get(opts, "checkpoint", "auto")
    entanglement_bonds = get(opts, "entanglement_bonds", "none")
    spectrum_levels = parse(Int, get(opts, "spectrum_levels", "12"))
    report = analyze_output(
        input;
        top_peaks,
        checkpoint,
        entanglement_bonds,
        spectrum_levels,
    )
    write_artifacts(input, outdir, report)
    println("input=$input")
    println("outdir=$outdir")
    println("energy=$(get(report.summary, "energy", "unknown"))")
    println("converged=$(get(report.summary, "converged", "unknown"))")
    println("density_total=$(fmt(report.density.total_density))")
    println("density_cdw_contrast=$(fmt(report.density.cdw_contrast))")
    if !isempty(report.structure_factor_peaks)
        peak = report.structure_factor_peaks[1]
        println("top_structure_peak=($(peak.mx),$(peak.my)) S=$(fmt(peak.value))")
    end
    if report.green_correlation_length !== nothing
        println("green_correlation_length_x=$(fmt(report.green_correlation_length))")
    end
    !isempty(report.entanglement_error) && println("entanglement_error=$(report.entanglement_error)")
    println("wrote=$(joinpath(outdir, "ground_state_diagnostics.svg"))")
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    GroundStateDiagnostics.main(ARGS)
end
