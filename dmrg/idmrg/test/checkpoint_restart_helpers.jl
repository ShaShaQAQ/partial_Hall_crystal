using TOML
using ITensors
using ITensorMPS
using ITensorInfiniteMPS
using InfiniteCylinderDMRG

function restart_config(phi_y::Real=0.0)
    return InfiniteCylinderConfig(;
        geometry=:paper_straight,
        Ny=2,
        x_period=3,
        filling_num=1,
        filling_den=3,
        phi_y,
    )
end

restart_model() = CylinderModelParams(; t1=1.0, t3=0.2, V1=1.0, V2=0.0, V3=0.0)

function index_identity(index)
    blocks = [
        Dict("qn" => string(qn(index, block)), "dimension" => blockdim(index, block))
        for block in 1:nblocks(index)
    ]
    return Dict(
        "id" => string(id(index)),
        "tags" => string(tags(index)),
        "direction" => string(dir(index)),
        "prime_level" => plev(index),
        "dimension" => dim(index),
        "qn_blocks" => blocks,
    )
end

function state_index_table(psi::InfiniteCanonicalMPS)
    rows = Dict{String,Any}[]
    for component_name in ("AL", "C", "AR")
        component = getproperty(psi, Symbol(component_name))
        for site in 1:nsites(psi)
            for (position, index) in enumerate(inds(component[site]))
                push!(
                    rows,
                    Dict(
                        "component" => component_name,
                        "site" => site,
                        "position" => position,
                        "index" => index_identity(index),
                    ),
                )
            end
        end
    end
    return rows
end

site_identity_table(psi::InfiniteCanonicalMPS) =
    [index_identity(index) for index in siteinds(only, psi.AL)]

function schmidt_qn_table(psi::InfiniteCanonicalMPS, config::InfiniteCylinderConfig)
    data = entanglement_data(psi, config; cut_x=1)
    return [
        Dict(
            "qn" => string(level.qn),
            "raw_charge" => level.raw_charge,
            "physical_charge" => level.physical_charge,
        ) for level in data.levels
    ]
end

function write_restart_toml(path::AbstractString, data)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return path
end
