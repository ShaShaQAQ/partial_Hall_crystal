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
iteration_validation = try
    InfiniteCylinderDMRG._validate_checkpoint_state(no_expansion.psi, config)
    "valid"
catch error
    sprint(showerror, error)
end
_, right_probe, _ = ITensorInfiniteMPS.right_orthogonalize(
    no_expansion.psi.AL;
    left_tags=ts"",
    right_tags=ts"Right",
)
write_restart_toml(
    joinpath(output_directory, "iterated_indices.toml"),
    Dict(
        "site_indices" => site_identity_table(no_expansion.psi),
        "state_indices" => state_index_table(no_expansion.psi),
        "link_dimensions" => link_dimensions(no_expansion.psi),
        "checkpoint_validation" => iteration_validation,
        "right_orthogonalized_indices" => state_index_table(
            InfiniteCanonicalMPS(right_probe, no_expansion.psi.C, right_probe)
        ),
    ),
)
continued = InfiniteCylinderDMRG._canonicalize_vumps_state(no_expansion.psi)

probe_n1, probe_n2 = 2, 3
probe_left_common = commoninds(continued.AL[probe_n1], continued.C[probe_n1])
probe_right_common = commoninds(continued.AR[probe_n2], continued.C[probe_n1])
probe_NL = nullspace(continued.AL[probe_n1], probe_left_common; atol=1e-2)
probe_NR = nullspace(continued.AR[probe_n2], probe_right_common; atol=1e-2)
probe_nL = uniqueinds(probe_NL, continued.AL[probe_n1])
probe_nR = uniqueinds(probe_NR, continued.AR[probe_n2])
probe_HN = ITensorInfiniteMPS.generate_twobody_nullspace(
    continued, H, (probe_n1, probe_n2); atol=1e-2
) * probe_NL * probe_NR
probe_U, _, probe_V = svd(
    probe_HN,
    probe_nL;
    maxdim=8 - dim(only(probe_left_common)),
    cutoff=1e-8,
)
probe_NL *= dag(probe_U)
probe_NR *= dag(probe_V)
write_restart_toml(
    joinpath(output_directory, "expansion_direction_probe.toml"),
    Dict(
        "bond" => [probe_n1, probe_n2],
        "left_common" => index_identity.(collect(probe_left_common)),
        "right_common" => index_identity.(collect(probe_right_common)),
        "left_existing" => index_identity.(
            collect(uniqueinds(continued.AL[probe_n1], probe_NL))
        ),
        "left_addition" => index_identity.(
            collect(uniqueinds(dag(probe_NL), continued.AL[probe_n1]))
        ),
        "right_existing" => index_identity.(
            collect(uniqueinds(continued.AR[probe_n2], probe_NR))
        ),
        "right_addition" => index_identity.(
            collect(uniqueinds(dag(probe_NR), continued.AR[probe_n2]))
        ),
    ),
)
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
