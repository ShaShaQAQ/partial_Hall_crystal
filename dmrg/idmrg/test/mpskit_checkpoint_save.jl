include(joinpath(@__DIR__, "mpskit_checkpoint_helpers.jl"))

length(ARGS) == 1 || error("usage: mpskit_checkpoint_save.jl OUTPUT_DIRECTORY")
output_directory = abspath(only(ARGS))
mkpath(output_directory)

config, hamiltonian, initial_state, number_operator =
    mpskit_checkpoint_test_fixture()
result = run_mpskit_idmrg(
    hamiltonian,
    initial_state;
    maxdim_schedule=[4],
    cutoff=1e-10,
    idmrg_maxiter=20,
    vumps_maxiter=50,
    galerkin_tol=1e-7,
)
result.converged || error("checkpoint fixture did not converge")

checkpoint = joinpath(output_directory, "state.h5")
metadata = (
    candidate_id="connected_open_chain",
    generation=1,
    requested_maxdim=4,
)
save_mpskit_checkpoint(checkpoint, result.state, config, metadata)
fingerprints = mpskit_state_space_fingerprints(result.state, config)
observations = mpskit_checkpoint_observations(
    result.state,
    hamiltonian,
    number_operator,
)
write_mpskit_checkpoint_toml(
    joinpath(output_directory, "before.toml"),
    Dict(
        "configuration_signature" => configuration_signature(config),
        "physical_space_fingerprint" => fingerprints.physical,
        "virtual_space_fingerprint" => fingerprints.virtual,
        "energy_per_site" => observations.energy_per_site,
        "galerkin_residual" => observations.galerkin_residual,
        "densities" => observations.densities,
        "sector_weights" => observations.sector_weights,
        "checkpoint_sha256" => checkpoint_sha256(checkpoint),
    ),
)
