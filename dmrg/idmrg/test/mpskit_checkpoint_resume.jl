include(joinpath(@__DIR__, "mpskit_checkpoint_helpers.jl"))

length(ARGS) == 1 || error("usage: mpskit_checkpoint_resume.jl OUTPUT_DIRECTORY")
output_directory = abspath(only(ARGS))
config, hamiltonian, _, number_operator = mpskit_checkpoint_restart_fixture()
checkpoint = joinpath(output_directory, "state.h5")

loaded = load_mpskit_checkpoint(checkpoint, config)
state = loaded.state
original_sha256 = checkpoint_sha256(checkpoint)
pre_fingerprints = mpskit_state_space_fingerprints(state, config)
pre = mpskit_checkpoint_observations(state, hamiltonian, number_operator)
pre_maxlinkdim = maximum(
    InfiniteCylinderDMRG._mpskit_solver_link_dimensions(state),
)
pre_maxlinkdim == 4 || error("loaded checkpoint does not have maxlinkdim=4")

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
no_expansion_maxlinkdim = maximum(
    InfiniteCylinderDMRG._mpskit_solver_link_dimensions(refined_state),
)
no_expansion_maxlinkdim == pre_maxlinkdim || error(
    "no-expansion VUMPS iteration changed the maximum link dimension",
)

expanded = run_mpskit_idmrg(
    hamiltonian,
    refined_state;
    maxdim_schedule=[8],
    cutoff=1e-10,
    idmrg_maxiter=1,
    vumps_maxiter=1,
    galerkin_tol=1e-6,
)
post_maxlinkdim = maximum(expanded.link_dimensions)
post_maxlinkdim == 8 || error(
    "controlled checkpoint expansion did not reach maxlinkdim=8",
)
post = mpskit_checkpoint_observations(
    expanded.state,
    hamiltonian,
    number_operator,
)

resumed_checkpoint = joinpath(output_directory, "state-resumed.h5")
save_mpskit_checkpoint(
    resumed_checkpoint,
    expanded.state,
    config,
    (
        candidate_id=loaded.metadata.candidate_id,
        generation=loaded.metadata.generation + 1,
        requested_maxdim=8,
    ),
)

next_flux_phi_y = 0.1
next_config, next_hamiltonian, _, _ =
    mpskit_checkpoint_restart_fixture(next_flux_phi_y)
next_flux_iterations = Ref(0)
next_environments = MPSKit.environments(
    expanded.state,
    next_hamiltonian,
    expanded.state,
)
next_flux_state, _, _ = MPSKit.find_groundstate(
    expanded.state,
    next_hamiltonian,
    MPSKit.VUMPS(;
        tol=0.0,
        maxiter=1,
        verbosity=0,
        alg_eigsolve=static_checkpoint_lanczos(),
        finalize=(iteration, psi, operator, environments) -> begin
            next_flux_iterations[] += 1
            return psi, environments
        end,
    ),
    next_environments,
)
next_flux_iterations[] == 1 || error(
    "restart did not execute one next-flux VUMPS iteration",
)
next_flux_checkpoint = joinpath(output_directory, "state-next-flux.h5")
save_mpskit_checkpoint(
    next_flux_checkpoint,
    next_flux_state,
    next_config,
    (
        candidate_id=loaded.metadata.candidate_id,
        generation=loaded.metadata.generation + 2,
        requested_maxdim=8,
        phi_y=next_flux_phi_y,
    ),
)
write_mpskit_checkpoint_toml(
    joinpath(output_directory, "after.toml"),
    Dict(
        "configuration_signature" => configuration_signature(config),
        "process_id" => string(getpid()),
        "load_completed" => true,
        "no_expansion_iteration_completed" => refinement_iterations[] == 1,
        "controlled_expansion_completed" => post_maxlinkdim == 8,
        "next_flux_iteration_completed" => next_flux_iterations[] == 1,
        "physical_space_fingerprint" => pre_fingerprints.physical,
        "virtual_space_fingerprint" => pre_fingerprints.virtual,
        "pre_maxlinkdim" => pre_maxlinkdim,
        "no_expansion_post_maxlinkdim" => no_expansion_maxlinkdim,
        "post_maxlinkdim" => post_maxlinkdim,
        "pre_energy_per_site" => pre.energy_per_site,
        "pre_galerkin_residual" => pre.galerkin_residual,
        "pre_densities" => pre.densities,
        "pre_sector_weights" => pre.sector_weights,
        "pre_entanglement_spectrum" => pre.entanglement_spectrum,
        "post_energy_per_site" => post.energy_per_site,
        "post_galerkin_residual" => post.galerkin_residual,
        "returned_residual" => returned_residual,
        "refinement_iterations" => refinement_iterations[],
        "next_flux_phi_y" => next_flux_phi_y,
        "next_flux_iterations" => next_flux_iterations[],
        "original_checkpoint_sha256" => original_sha256,
        "resumed_checkpoint_sha256" => checkpoint_sha256(resumed_checkpoint),
        "next_flux_checkpoint_sha256" =>
            checkpoint_sha256(next_flux_checkpoint),
    ),
)
