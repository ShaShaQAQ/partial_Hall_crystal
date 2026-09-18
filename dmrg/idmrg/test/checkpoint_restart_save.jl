include(joinpath(@__DIR__, "checkpoint_restart_helpers.jl"))

length(ARGS) == 1 || error("usage: checkpoint_restart_save.jl OUTPUT_DIRECTORY")
output_directory = abspath(only(ARGS))
mkpath(output_directory)

config = restart_config()
sites, _, product_state = initial_infinite_mps(config)
H = build_infinite_mpo(config, restart_model(), sites)
expanded = expand_subspace(product_state, H, 4; cutoff=1e-8)
step = vumps_iteration(H, expanded; vumps_tol=1e-4)
psi = InfiniteCylinderDMRG._canonicalize_vumps_state(step.psi)

checkpoint = joinpath(output_directory, "state.h5")
save_checkpoint(checkpoint, psi, config)
saved = Dict(
    "energy" => energy_data(psi, H, config).per_cell,
    "density" => [row.density for row in density_data(psi, config)],
    "schmidt_qns" => schmidt_qn_table(psi, config),
    "site_indices" => site_identity_table(psi),
    "state_indices" => state_index_table(psi),
    "maxlinkdim" => maximum(link_dimensions(psi)),
)
write_restart_toml(joinpath(output_directory, "saved.toml"), saved)
