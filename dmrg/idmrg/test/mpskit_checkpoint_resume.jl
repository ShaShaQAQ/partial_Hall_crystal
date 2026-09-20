include(joinpath(@__DIR__, "mpskit_checkpoint_helpers.jl"))

length(ARGS) == 1 || error("usage: mpskit_checkpoint_resume.jl OUTPUT_DIRECTORY")
output_directory = abspath(only(ARGS))
config, hamiltonian, _, number_operator = mpskit_checkpoint_test_fixture()
checkpoint = joinpath(output_directory, "state.h5")

loaded = load_mpskit_checkpoint(checkpoint, config)
state = loaded.state
original_sha256 = checkpoint_sha256(checkpoint)
pre_fingerprints = mpskit_state_space_fingerprints(state, config)
pre = mpskit_checkpoint_observations(state, hamiltonian, number_operator)

refinement_iterations = Ref(0)
refinement_environments = MPSKit.environments(state, hamiltonian, state)
refined_state, _, returned_residual = MPSKit.find_groundstate(
    state,
    hamiltonian,
    MPSKit.VUMPS(;
        tol=0.0,
        maxiter=1,
        verbosity=0,
        alg_eigsolve=static_checkpoint_lanczos(),
        finalize=(iteration, psi, operator, environments) -> begin
            refinement_iterations[] += 1
            return psi, environments
        end,
    ),
    refinement_environments,
)
refinement_iterations[] == 1 || error("restart did not execute one VUMPS iteration")
post = mpskit_checkpoint_observations(
    refined_state,
    hamiltonian,
    number_operator,
)

resumed_checkpoint = joinpath(output_directory, "state-resumed.h5")
save_mpskit_checkpoint(
    resumed_checkpoint,
    refined_state,
    config,
    (
        candidate_id=loaded.metadata.candidate_id,
        generation=loaded.metadata.generation + 1,
        requested_maxdim=loaded.metadata.requested_maxdim,
    ),
)
write_mpskit_checkpoint_toml(
    joinpath(output_directory, "after.toml"),
    Dict(
        "configuration_signature" => configuration_signature(config),
        "physical_space_fingerprint" => pre_fingerprints.physical,
        "virtual_space_fingerprint" => pre_fingerprints.virtual,
        "pre_energy_per_site" => pre.energy_per_site,
        "pre_galerkin_residual" => pre.galerkin_residual,
        "pre_densities" => pre.densities,
        "pre_sector_weights" => pre.sector_weights,
        "post_energy_per_site" => post.energy_per_site,
        "post_galerkin_residual" => post.galerkin_residual,
        "returned_residual" => returned_residual,
        "refinement_iterations" => refinement_iterations[],
        "original_checkpoint_sha256" => original_sha256,
        "resumed_checkpoint_sha256" => checkpoint_sha256(resumed_checkpoint),
    ),
)
