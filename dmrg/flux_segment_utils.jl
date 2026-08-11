using Printf

struct PumpingRow
    step::Int
    line::String
end

function point_dir(root::AbstractString, step::Int)
    return joinpath(root, @sprintf("phi_%03d", step))
end

function checkpoint_path(root::AbstractString, step::Int)
    return joinpath(root, "checkpoints", @sprintf("state_%03d.jls", step))
end

function default_phi_for_step(step::Int; cycles::Int, steps_per_2pi::Int, direction::AbstractString)
    cycles > 0 || throw(ArgumentError("cycles must be positive"))
    steps_per_2pi > 0 || throw(ArgumentError("steps_per_2pi must be positive"))
    direction = lowercase(direction)
    direction in ("forward", "backward") || throw(ArgumentError("direction must be forward or backward"))
    phi0 = direction == "forward" ? 0.0 : 2pi * cycles
    sign = direction == "forward" ? 1.0 : -1.0
    return phi0 + sign * 2pi * step / steps_per_2pi
end

function segment_phi_for_step(
    step::Int;
    cycles::Int,
    steps_per_2pi::Int,
    direction::AbstractString,
    start_step::Int,
    phi_start::Union{Nothing,Float64}=nothing,
    dphi::Union{Nothing,Float64}=nothing,
)
    if phi_start !== nothing || dphi !== nothing
        phi_start === nothing && throw(ArgumentError("phi_start is required when dphi is set"))
        dphi === nothing && throw(ArgumentError("dphi is required when phi_start is set"))
        return phi_start + dphi * (step - start_step)
    end
    return default_phi_for_step(step; cycles, steps_per_2pi, direction)
end

function segment_pumping_path(root::AbstractString, start_step::Int, end_step::Int)
    return joinpath(root, @sprintf("pumping_segment_%03d_%03d.dat", start_step, end_step))
end

function pumping_header()
    return "# step phi energy cumulative_left_edge cumulative_mid cumulative_right_edge actual_sweeps attempts converged reason"
end

function collect_pumping_rows(root::AbstractString)
    isdir(root) || return PumpingRow[]
    rows = Dict{Int,PumpingRow}()
    for name in sort(readdir(root))
        occursin(r"^pumping_segment_[0-9]{3}_[0-9]{3}\.dat$", name) || continue
        path = joinpath(root, name)
        for line in eachline(path)
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(stripped, "#") && continue
            fields = split(stripped)
            length(fields) >= 10 || throw(ArgumentError("malformed pumping row in $path: $line"))
            step = parse(Int, fields[1])
            rows[step] = PumpingRow(step, stripped)
        end
    end
    return [rows[step] for step in sort(collect(keys(rows)))]
end

function refresh_pumping_index(root::AbstractString)
    rows = collect_pumping_rows(root)
    open(joinpath(root, "pumping.dat"), "w") do io
        println(io, pumping_header())
        for row in rows
            println(io, row.line)
        end
    end
    return length(rows)
end
