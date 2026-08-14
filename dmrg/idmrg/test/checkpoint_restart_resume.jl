include(joinpath(@__DIR__, "checkpoint_restart_helpers.jl"))

length(ARGS) == 1 || error("usage: checkpoint_restart_resume.jl OUTPUT_DIRECTORY")
output_directory = abspath(only(ARGS))
config = restart_config()
checkpoint = joinpath(output_directory, "state.h5")

psi = load_checkpoint(checkpoint, config)
sites = siteinds(only, psi.AL)
H = build_infinite_mpo(config, restart_model(), sites)
pre_maxlinkdim = maximum(link_dimensions(psi))
pre_energy = energy_data(psi, H, config).per_cell
pre_density = [row.density for row in density_data(psi, config)]
pre_schmidt_qns = schmidt_qn_table(psi, config)
write_restart_toml(
    joinpath(output_directory, "loaded_indices.toml"),
    Dict(
        "site_indices" => site_identity_table(psi),
        "state_indices" => state_index_table(psi),
        "link_dimensions" => link_dimensions(psi),
    ),
)

no_expansion = vumps_iteration(H, psi; vumps_tol=1e-4)
write_restart_toml(
    joinpath(output_directory, "iterated_indices.toml"),
    Dict(
        "site_indices" => site_identity_table(no_expansion.psi),
        "state_indices" => state_index_table(no_expansion.psi),
        "link_dimensions" => link_dimensions(no_expansion.psi),
    ),
)
continued = InfiniteCylinderDMRG._canonicalize_vumps_state(no_expansion.psi)
expanded = expand_subspace(continued, H, 8; cutoff=1e-8)

next_config = restart_config(0.1)
next_H = build_infinite_mpo(next_config, restart_model(), sites)
next_step = vumps_iteration(next_H, expanded; vumps_tol=1e-4)
next_psi = InfiniteCylinderDMRG._canonicalize_vumps_state(next_step.psi)

resumed = Dict(
    "pre_energy" => pre_energy,
    "pre_density" => pre_density,
    "pre_schmidt_qns" => pre_schmidt_qns,
    "pre_site_indices" => site_identity_table(psi),
    "pre_state_indices" => state_index_table(psi),
    "pre_maxlinkdim" => pre_maxlinkdim,
    "post_maxlinkdim" => maximum(link_dimensions(next_psi)),
    "next_flux_iteration_completed" => true,
)
write_restart_toml(joinpath(output_directory, "resumed.toml"), resumed)
